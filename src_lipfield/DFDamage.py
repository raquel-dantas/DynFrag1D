import DFMesh
import numpy as np
import scipy
from scipy.optimize import minimize
from scipy.optimize import LinearConstraint
from itertools import groupby
from operator import itemgetter


# Regularization inputs
n_element_reg = 10  # Nb elem in the reg length
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
            (
                0.5 * (1.0 - damage[el]) ** 2 * DFMesh.E * strain[el] ** 2
                + Yc[el] * h(lamb[el], damage[el])
            )
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
            x0=0.5 * DFMesh.L,
            method="SLSQP",
            bounds=[(DFMesh.x0, DFMesh.xf)],
            tol=1e-6,
        )
        if upper_opt.success == False:
            raise Exception("upper projection of damage predictor failed")
        upper[el] = -upper_opt.fun

        lower_opt = minimize(
            lambda y: damage_prediction[np.searchsorted(DFMesh.node_coord, y[0]) - 1]
            + abs(DFMesh.x[el] - y[0]) / l,
            x0=0.5 * DFMesh.L,
            method="SLSQP",
            bounds=[(DFMesh.x0, DFMesh.xf)],
            tol=1e-6,
        )
        if lower_opt.success == False:
            raise Exception("lower projection of damage predictor failed")
        lower[el] = lower_opt.fun

    return upper, lower


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


def computeDamageNextTimeStep(u_next, dn):
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

        # Solve the optimization problem for each subregion
        for subregion in regions:
            if len(subregion) > 1:
                dlip = computeDamageLipConstraint(strain, subregion, dn)
                i = 0
                for intpoint in subregion:
                    d_next[intpoint] = dlip[i]
                    i = +1

    return d_next
