import DFMesh
import numpy as np
import scipy
from scipy.optimize import minimize
from scipy.optimize import LinearConstraint
from itertools import groupby
from operator import itemgetter
from sortedcontainers import SortedList


# Regularization inputs
n_element_reg = 10              # Nb elem in the reg length
l = n_element_reg * DFMesh.hun  # Regularization length
w = 2.0  # Weight quadrature

# Energy release rate
Yc = [DFMesh.sigmac[el] ** 2 / (2.0 * DFMesh.E) for el in range(DFMesh.n_el)]

# Constant lambda
lamb_lip = [2.0 * Yc[el] * l / DFMesh.Gc for el in range(DFMesh.n_el)]
lamb_czm = [
    DFMesh.sigmac[el] * DFMesh.hun / (DFMesh.E * DFMesh.deltac[el])
    for el in range(DFMesh.n_el)
]
# lamb = lamb_czm
lamb = lamb_lip

if max(lamb) > 1./3.:
    raise Exception("lambda > 1/3 -> h(d) not convex!")

# Softening function
def h_lip(lamb, d):
    return (2.0 * d - d**2) / (1.0 - d + lamb * d**2) ** 2


def h_czm(lamb, d):
    return 1.0 / (1.0 - lamb) * ((1.0 / ((1.0 - lamb) * (1.0 - d) ** 2 + lamb)) - 1.0)


# h = h_czm
h = h_lip


# Functional given by Moes, Le and Stershic (2022)
def getFunctional(strain):
    """Returns a lambda function of damage to be use in the optimization problem for the whole domain.\n
    Arguments: \n
    strain -- strain computed using the displacement at the next time-step(n+1)."""

    functional_whole_domain = lambda damage: w * sum(
        [
            (0.5 * (1.0 - damage[el]) ** 2 * DFMesh.E * strain[el] ** 2
            + Yc[el] * h(lamb[el], damage[el]))
            * DFMesh.hun
            * 0.5
            for el in range(DFMesh.n_el)
        ]
    )

    return functional_whole_domain


def getFunctionalSubdomain(strain, region_optimization):
    """Returns a lambda function of damage to be use in the optimization problem for a specific region in the domain (subdomain).\n
    Arguments: \n
    strain -- strain computed using the displacement at the next time-step(n+1).\n
    region -- an array with the element indexes of the subdomain to compute damage at next time-step(n+1)."""

    functional_subdomain = lambda damage: w * sum(
        [
            (
                0.5
                * (1.0 - damage[i]) ** 2
                * DFMesh.E
                * strain[region_optimization[i]] ** 2
                + Yc[region_optimization[i]]
                * h(lamb[region_optimization[i]], damage[i])
            )
            * DFMesh.hun
            * 0.5
            for i in range(len(region_optimization))
        ]
    )

    return functional_subdomain


def computeDamagePredictor(u_next, dn):
    """Returns an array with the damage predictor for the whole domain.\n
    Arguments: \n
    u_next -- displacement at the next time-step(n+1); \n
    dn -- damage at time-step(n)"""

    functional = getFunctional(u_next)

    damage_predictor_opt = minimize(
        fun=functional,
        x0=dn,
        method="SLSQP",
        bounds=zip(dn, [1.0] * DFMesh.n_el),
        tol=1e-6,
    )
    if damage_predictor_opt.success == False:
        raise Exception("optimization damage predictor failed")

    damage_predictor = damage_predictor_opt.x

    return damage_predictor


def computeProjections(damage_prediction):
    """Returns the damage prediction upper and lower projections."""

    lower = np.zeros(DFMesh.n_el)
    upper = np.zeros(DFMesh.n_el)

    for el in range(DFMesh.n_el):
        upper_opt = minimize(
            lambda y: -damage_prediction[np.searchsorted(DFMesh.node_coord, y[0]) - 1]
            + abs(DFMesh.x[el] - y[0]) / l,
            x0=0.,
            method="SLSQP",
            bounds=[(DFMesh.x0, DFMesh.xf)],
            tol=1e-5
        )
        if upper_opt.success == False:
            raise Exception("upper projection of damage predictor failed")
        upper[el] = -upper_opt.fun
        

        lower_opt = minimize(
            lambda y: damage_prediction[np.searchsorted(DFMesh.node_coord, y[0])-1] + abs(DFMesh.x[el]-y[0])/l,
            x0=0.,
            method='SLSQP',
            bounds=[(DFMesh.x0, DFMesh.xf)],
            tol=1e-5
        )
        if lower_opt.success == False:
            raise Exception('lower projection damage predictor failed')
        lower[el] = lower_opt.fun


    return upper, lower

def get_neighbour(index):
    if index == 0:
        return [0]
    if index == DFMesh.n_el - 1:
        return [index - 1]
    return [index - 1, index + 1]


def computeProjectionsUsingFM_lip_projector_1D(damage_predictor, regularization_lenght, flank):
    
    # Let's assume initially the projection equal to the damage predictor
    projection = damage_predictor.copy()

    # Configure key for trial set according to the projection if upper (flank=max) or lower (flank=min)
    if flank == "min":
        trial_set = SortedList(key=lambda x: (-x[0], x[1]))
    else:
        trial_set = SortedList()

    # Add all elements to the trial set
    for index, projection_value in enumerate(projection):
        trial_set.add((projection_value, index))

    # Initialize the frozen set 
    frozen_set = set()

    while len(trial_set) > 0:

        projection_current_index, index = trial_set.pop()
        frozen_set.add(index)
        neighbours = get_neighbour(index)

        for index_neighbour in neighbours:
            if index_neighbour not in frozen_set:      
                update_projection_value = False

                delta_projection = (projection[index_neighbour] - projection_current_index) / DFMesh.hun * regularization_lenght 

                if delta_projection < -1.0:
                    update_projection_value = True
                    new_projection = projection_current_index - DFMesh.hun /  regularization_lenght

                elif delta_projection > 1.0:
                    update_projection_value = True
                    new_projection = projection_current_index + DFMesh.hun / regularization_lenght

                if update_projection_value == True:
                    projection[index_neighbour] = new_projection
                    trial_set.discard((projection[index_neighbour], index_neighbour))
                    trial_set.add((new_projection, index_neighbour))

    return projection



def computeDamageLipConstraint(strain, region_optimization, dn):
    """Returns an array with the damage at the next time-step for one subdomain imposing the Lipischitz constraint.\n
    Arguments:
    strain -- array with the strain for next time-step for all the elements;\n
    region_optimization -- an array with the indexes of elements in one subdomain; \n
    dn -- damage for all of the domain at time-step (n)."""

    functional = getFunctionalSubdomain(strain, region_optimization)
    # Size of the domain
    size = len(region_optimization)
    # Inputs for LinearConstraint
    A = scipy.sparse.eye(size - 1, size) - scipy.sparse.eye(size - 1, size, 1)
    b = DFMesh.hun / l
    constraints = LinearConstraint(A, -b * np.ones(size - 1), b * np.ones(size - 1))
    # Bounds
    bound_inf = [dn[region_optimization[i]] for i in range(size)]
    bound_sup = [1.0 for i in range(size)]

    dlip_opt = minimize(
        fun=functional,
        x0=bound_inf,
        method="SLSQP",
        bounds=zip(bound_inf, bound_sup),
        tol=1e-6,
        constraints=constraints,
    )
    if dlip_opt.success == False:
        raise Exception("optimization failed")

    dlip = dlip_opt.x

    return dlip


def groupSubregion(region_lip):
    regions = []
    for i, subgroup in groupby(
        enumerate(region_lip), lambda index: index[0] - index[1]
    ):
        regions.append(list(map(itemgetter(1), subgroup)))
    return regions


def computeDamageNextTimeStep(u_next, dn, use_FM=False):
    """Returns the damage at the next time-step (n+1) for all the domain.\n
    Arguments:\n
    u_next -- displacement at the next time-step (n+1);\n
    dn -- damage at time-step (n)."""

    d_next = np.zeros(DFMesh.n_el)
    region_lip = []
    small_number = 10e-5

    # Compute strain using u_next
    strain = [
        (u_next[DFMesh.connect[el][1]] - u_next[DFMesh.connect[el][0]])
        / DFMesh.ElemLength(el)
        for el in range(DFMesh.n_el)
    ]

    # Compute damage predictor (dp) -- Only imposition of D space
    dp = computeDamagePredictor(strain, dn)

    # Compute upper and lower projections of the damage predictor
    if use_FM == True:
        upper = computeProjectionsUsingFM_lip_projector_1D(dp, l, flank="max")
        lower = computeProjectionsUsingFM_lip_projector_1D(dp, l, flank="min")
    else:
        upper, lower = computeProjections(dp)
    
    # Verify if the projections are supperposed
    for el in range(DFMesh.n_el):

        # If projections are superposed
        if upper[el] - lower[el] < small_number:
            d_next[el] = dp[el]

        else:
            # Add element to region to impose the Lipschitz constraint
            region_lip.append(el)

    if region_lip:
        # Separate the consecutive elements in subregions
        regions = groupSubregion(region_lip)
        # print(region_lip)

        # Solve the optimization problem for each subregion
        for subregion in regions:
            if len(subregion) > 1:
                dlip = computeDamageLipConstraint(strain, subregion, dn)
                i = 0
                for intpoint in subregion:
                    d_next[intpoint] = dlip[i]
                    i = +1

    return d_next



    