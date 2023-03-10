import numpy as np
import scipy
from scipy.optimize import minimize
from scipy.optimize import LinearConstraint
from itertools import groupby
from operator import itemgetter
from sortedcontainers import SortedList
import copy

import DFMesh
import DFFem


if DFMesh.use_lipfield == True:

    # Regularization inputs
    # n_elements_regularization = 10
    # regularization_length = n_elements_regularization * DFMesh.h_uniform
    regularization_length = 2.21 * 10**-6
    weight_quadrature = 2.0
    tolerance_opt = 10e-5


    # Energy release rate (Yc)
    Yc = [
        DFMesh.stress_critical[el] ** 2 / (2.0 * DFMesh.young_modulus)
        for el in range(DFMesh.n_elements)
    ]


    # Constant lambda
    lamb_lip = [
        2.0 * Yc[el] * regularization_length / DFMesh.fracture_energy
        for el in range(DFMesh.n_elements)
    ]
    lamb_czm_mimic = [
        DFMesh.stress_critical[el]
        * DFMesh.h_uniform
        / (DFMesh.young_modulus * DFMesh.crack_critical[el])
        for el in range(DFMesh.n_elements)
    ]
    lamb_const = lamb_lip
    if max(lamb_const) > 1.0 / 3.0:
        raise Exception("lambda > 1/3 -> h(d) not convex!")


    # Softening function of damage h(d)
    def h_lip(lamb, d):
        return (2.0 * d - d**2) / (1.0 - d + lamb * d**2) ** 2


    def h_czm(lamb, d):
        return 1.0 / (1.0 - lamb) * ((1.0 / ((1.0 - lamb) * (1.0 - d) ** 2 + lamb)) - 1.0)


    h = h_lip


    # Functional given by Moes, Le and Stershic (2022)
    def getFunctionalWholeDomain(strain):
        """Returns a lambda function of damage to be use in the optimization problem for the whole domain.\n
        Arguments: \n
        strain -- strain computed using the displacement at the next time-step(n+1)."""

        functional_whole_domain = lambda damage: weight_quadrature * sum(
            [
                (
                    0.5 * (1.0 - damage[el]) ** 2 * DFMesh.young_modulus * strain[el] ** 2
                    + Yc[el] * h(lamb_const[el], damage[el])
                )
                * DFMesh.getElemLength(el)
                * 0.5
                for el in range(DFMesh.n_elements)
            ]
        )

        return functional_whole_domain


    def getFunctionalSubdomain(strain, region_lip_optimization):
        """Returns a lambda function of damage to be use in the optimization problem for a specific region in the domain (subdomain).\n
        Arguments: \n
        strain -- strain computed using the displacement at the next time-step(n+1).\n
        region_lip_optimization -- an array with the element indexes of the subdomain to compute damage at next time-step(n+1)."""

        functional_subdomain = lambda damage: weight_quadrature * sum(
            [
                (
                    0.5
                    * (1.0 - damage[i]) ** 2
                    * DFMesh.young_modulus
                    * strain[region_lip_optimization[i]] ** 2
                    + Yc[region_lip_optimization[i]]
                    * h(lamb_const[region_lip_optimization[i]], damage[i])
                )
                * DFMesh.getElemLength(region_lip_optimization[i])
                * 0.5
                for i in range(len(region_lip_optimization))
            ]
        )

        return functional_subdomain


    def computeDamagePredictor_useSLSQP(strain, dn):
        """Returns an array with the damage predictor for the whole domain.\n
        Arguments: \n
        strain -- strain computed with the displacement at the next time-step(n+1); \n
        dn -- damage at time-step(n)"""

        functional = getFunctionalWholeDomain(strain)

        damage_predictor_opt = minimize(
            fun=functional,
            x0=dn,
            method="SLSQP",
            bounds=zip(dn, [1.0] * DFMesh.n_elements),
            tol=tolerance_opt,
        )
        if damage_predictor_opt.success == False:
            raise Exception("optimization damage predictor failed")

        damage_predictor = damage_predictor_opt.x

        return damage_predictor


    def getGradFunctional(strain, damage):

        first_derivative_g = [-2.0 * (1.0 - damage[i]) for i in range(len(damage))]
        first_derivative_h = [
            (2.0 + 2.0 * (-3.0 + damage[i]) * damage[i] ** 2 * lamb_const[i])
            / (1.0 + damage[i] * (-1.0 + damage[i] * lamb_const[i])) ** 3
            for i in range(len(damage))
        ]
        grad = np.array(
            [
                0.5 * first_derivative_g[i] * DFMesh.young_modulus * strain[i] ** 2
                + Yc[i] * first_derivative_h[i]
                for i in range(len(damage))
            ]
        )

        return grad


    def getHessFunctional(strain, damage):

        second_derivative_g = 2.0
        second_derivative_h = [
            (
                6.0
                + 6.0
                * damage[i]
                * lamb_const[i]
                * (-4.0 + damage[i])
                * damage[i] ** 2
                * lamb_const[i]
            )
            / (1.0 + damage[i] * (-1.0 + damage[i] * lamb_const[i])) ** 4
            for i in range(len(damage))
        ]
        hess = np.array(
            [
                0.5 * second_derivative_g * DFMesh.young_modulus * strain[i] ** 2
                + Yc[i] * second_derivative_h[i]
                for i in range(len(damage))
            ]
        )

        return hess


    def computeDamagePredictor_useNewton(strain, damage_previous_step):

        d_previous = damage_previous_step

        grad = getGradFunctional(strain, d_previous)
        hess = getHessFunctional(strain, d_previous)

        dp = d_previous - grad / hess
        dp = np.clip(dp, damage_previous_step, 1.)

        n_interrupt = 0
        n_max_interations = 500
        error = tolerance_opt
        norm_residual = np.linalg.norm(dp - d_previous) / np.sqrt(DFMesh.n_elements)

        while n_interrupt < n_max_interations and norm_residual > error:

            n_interrupt += 1
            d_previous = dp

            grad = getGradFunctional(strain, d_previous)
            hess = getHessFunctional(strain, d_previous)

            dp = d_previous - grad / hess
            dp = np.clip(dp, damage_previous_step, 1.)

            norm_residual = np.linalg.norm(dp - d_previous) / np.sqrt(DFMesh.n_elements)

        if norm_residual > error:
            print(norm_residual)
            raise Exception("optimization damage predictor failed ")

        return dp


    def computeProjections_useSLSQP(damage_prediction):
        """Returns the damage_prediction upper and lower projections by solving an optimization problem with scipy.optimize.minimize method=SLSQP."""

        lower = np.zeros(DFMesh.n_elements)
        upper = np.zeros(DFMesh.n_elements)

        for el in range(DFMesh.n_elements):
            upper_opt = minimize(
                lambda y: -damage_prediction[np.searchsorted(DFMesh.node_coord, y[0]) - 1]
                + abs(DFMesh.intpoint_coord[el] - y[0]) / regularization_length,
                x0=0.0,
                method="SLSQP",
                bounds=[(DFMesh.x0, DFMesh.xf)],
                tol=tolerance_opt,
            )
            if upper_opt.success == False:
                raise Exception("upper projection of damage predictor failed")
            upper[el] = -upper_opt.fun

            lower_opt = minimize(
                lambda y: damage_prediction[np.searchsorted(DFMesh.node_coord, y[0]) - 1]
                + abs(DFMesh.intpoint_coord[el] - y[0]) / regularization_length,
                x0=0.0,
                method="SLSQP",
                bounds=[(DFMesh.x0, DFMesh.xf)],
                tol=tolerance_opt,
            )
            if lower_opt.success == False:
                raise Exception("lower projection damage predictor failed")
            lower[el] = lower_opt.fun

        return upper, lower


    def get_neighbour(index):
        if index == 0:
            return [0]
        if index == DFMesh.n_elements - 1:
            return [index - 1]
        return [index - 1, index + 1]


    # Use FM to compute projections based in function provide by GEM code
    def computeProjections_useFM(damage_predictor, flank):
        """Returns the damage_predictor projection by usign FastMarching.\n
        For the uuper prediction set flank='max'\n
        For the lower prediction set flank='min'."""

        # Let's assume initially the projection equal to the damage predictor
        projection = damage_predictor.copy()
        slope_limit = 1.0 / regularization_length
        # Configure key for trial set according to the projection if upper (flank=max) or lower (flank=min)
        if flank == "min":
            trial_set = SortedList(key=lambda x: (-x[0], x[1]))
        if flank == "max":
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

                    delta_projection = (
                        projection[index_neighbour] - projection_current_index
                    ) / DFMesh.h_uniform

                    if delta_projection < -slope_limit:
                        update_projection_value = True
                        new_projection = (
                            projection_current_index - DFMesh.h_uniform / regularization_length
                        )

                    elif delta_projection > slope_limit:
                        update_projection_value = True
                        new_projection = (
                            projection_current_index + DFMesh.h_uniform / regularization_length
                        )

                    if update_projection_value == True:
                        trial_set.discard((projection[index_neighbour], index_neighbour))
                        trial_set.add((new_projection, index_neighbour))
                        projection[index_neighbour] = new_projection

        return projection



    def computeDamageLipConstraint(strain, region_optimization, dn, upper, lower):
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
        b = DFMesh.h_uniform / regularization_length
        constraints = LinearConstraint(A, -b * np.ones(size - 1), b * np.ones(size - 1))
        # Bounds
        bound_inf = [dn[region_optimization[i]] for i in range(size)]
        # bound_inf = [lower[region_optimization[i]] for i in range(size)]
        bound_sup = [1.0 for i in range(size)]
        # bound_sup = [upper[region_optimization[i]] for i in range(size)]

        dlip_opt = minimize(
            fun=functional,
            x0=bound_inf,
            method="SLSQP",
            bounds=zip(bound_inf, bound_sup),
            tol=tolerance_opt,
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


    def computeDamageNextStep_useProjection(
        u_next, dn, predictor_method, projection_method
    ):
        """Returns the damage at the next time-step (n+1) for all the domain.\n
        Arguments:\n
        u_next -- displacement at the next time-step (n+1);\n
        dn -- damage at time-step (n);\n
        predictor_method -- choose the method 'SLSQP' or 'Newton' to solve damage predictor; \n
        projection_method -- choose the method 'SLSQP' or 'FM' to solve the projections."""

        d_next = np.zeros(DFMesh.n_elements)
        region_lip = []
        small_number = 10e-5

        # Compute strain using u_next
        strain = [
            (u_next[DFMesh.connect[el][1]] - u_next[DFMesh.connect[el][0]])
            / DFMesh.getElemLength(el)
            for el in range(DFMesh.n_elements)
        ]

        # Compute damage predictor (dp) -- Only imposition of D space
        if predictor_method == "SLSQP":
            dp = computeDamagePredictor_useSLSQP(strain, dn)
        if predictor_method == "Newton":
            dp = computeDamagePredictor_useNewton(strain, dn)

        # dp_copy = copy.deepcopy(dp)
        # Compute upper and lower projections of the damage predictor
        if projection_method == "SLSQP":
            upper, lower = computeProjections_useSLSQP(dp)
        if projection_method == "FM":
            upper = computeProjections_useFM(dp, flank="max")
            lower = computeProjections_useFM(dp, flank="min")
        # assert np.linalg.norm(dp_copy - dp) < small_number

        # Verify if the projections are supperposed
        for el in range(DFMesh.n_elements):

            d_next[el] = dp[el]

            # If projections are superposed
            if abs(upper[el] - lower[el]) > small_number:
                region_lip.append(el)

        if region_lip:
            # Separate the consecutive elements in subregions
            # region_lip_copy = copy.deepcopy(region_lip)
            regions = groupSubregion(region_lip)
            # assert np.linalg.norm(np.array(region_lip_copy) - np.array(region_lip)) < small_number

            # Solve the optimization problem for each subregion
            for subregion in regions:
                if len(subregion) > 1:
                    dlip = computeDamageLipConstraint(strain, subregion, dn, upper, lower)
                    d_next[subregion] = dlip

                    # for el in subregion:
                    #     tolerance = 10e-8
                    #     assert d_next[el] - dn[el] > - tolerance
                    #     assert dp[el] - dn[el] > - tolerance
                    #     assert upper[el] - dp[el] > - tolerance
                    #     assert lower[el] - dp[el] < tolerance

        # test derivative
        # if __debug__:
        #     for el in range(1, DFMesh.n_elements):
        #         assert abs( d_next[el] - d_next[el - 1]) < abs( dn[el] - dn[el - 1])

        return d_next


    def computeDamageNextStep_noProjection(u_next, dprevious_step):

        # Compute strain using u_next
        strain = [
            (u_next[DFMesh.connect[el][1]] - u_next[DFMesh.connect[el][0]])
            / DFMesh.getElemLength(el)
            for el in range(DFMesh.n_elements)
        ]

        functional = getFunctionalWholeDomain(strain)
        # Size of the domain
        size = len(dprevious_step)

        # Inputs for LinearConstraint
        A = scipy.sparse.eye(size - 1, size) - scipy.sparse.eye(size - 1, size, 1)
        b = DFMesh.h_uniform / regularization_length
        constraints = LinearConstraint(A, -b * np.ones(size - 1), b * np.ones(size - 1))

        dlip_opt = minimize(
            fun=functional,
            x0=dprevious_step,
            method="SLSQP",
            bounds=zip(dprevious_step, [1.0] * size),
            tol=tolerance_opt,
            constraints=constraints,
        )
        if dlip_opt.success == False:
            raise Exception("optimization failed")

        d_next = dlip_opt.x

        return d_next


    def internalForce(u, d):
        """Returns the internal force vector (kd(d)u)\n
        Arguments:\n
        u -- displacemnt vector for all dofs; \n
        d -- for all elements"""

        n_dofs = u.shape[0]
        fint = np.zeros(n_dofs)

        for el in range(DFMesh.n_elements):
            g = (1.0 - d[el]) ** 2
            # u_loc returns a vector contained u for a local dof
            u_loc = np.array(
                [u[DFFem.getGlobalIndex(el, 0)], u[DFFem.getGlobalIndex(el, 1)]]
            )
            fint_loc = g * np.matmul(DFFem.k_elem, u_loc) / DFMesh.getElemLength(el)
            # Contribution of each dof in the internal force vector
            for i_loc in range(2):
                i_gl = DFFem.getGlobalIndex(el, i_loc)
                fint[i_gl] += fint_loc[i_loc]

        return fint
