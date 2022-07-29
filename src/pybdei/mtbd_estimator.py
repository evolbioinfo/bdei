import os
from itertools import chain
from multiprocessing.pool import ThreadPool

import numpy as np
import pandas as pd
from ete3 import Tree
from pybdei import initial_rate_guess
from scipy.integrate import odeint
from scipy.optimize import minimize, fsolve

from treesimulator.mtbd_models import BirthDeathExposedInfectiousModel

N_U_STEPS = int(1e7)

STATE_K = 'state'
TI = 'ti'

MIN_VALUE = np.log(np.finfo(np.float64).eps + 1)
MAX_VALUE = np.log(np.finfo(np.float64).max)


def state_frequencies(MU, LA, PSI):
    """
    Calculates equilibrium state frequencies for given rate values.

    :param MU: an array of state transition rates
    :param LA: an array of transmission rates
    :param PSI:  an array of removal rates
    :return: an array of equilibrium state frequencies [pi_0, ..., pi_m]
    """
    m = len(PSI)

    def func(PI):
        SIGMA = PI.dot(LA.sum(axis=1) - PSI)
        res = [PI.sum() - 1]
        for k in range(m - 1):
            pi_k = PI[k]
            res.append(pi_k * SIGMA + pi_k * (PSI[k] + MU[k, :].sum()) - PI.dot(MU[:, k] + LA[:, k]))
        return res

    PI = fsolve(func, np.ones(m) / m)
    if np.any(PI < 0) or np.any(PI > 1):
        return np.ones(m) / m
    return PI


def find_time_index(v, tt, start, stop):
    """
    Searches for an index i in time array tt, such that tt[i] >= v > tt[i + 1], using a binary search.

    :param v: a time value for which the index is searched for
    :param tt: a time array [t_n, ..., t_0] such that t_{i + 1} > t_i.
    :param start: start index for the search (inclusive)
    :param stop: stop index for the search (exclusive)
    :return: an index i, such that tt[i] >= v > tt[i + 1].
    """
    i = start + ((stop - start) // 2)
    if i == start or i == stop - 1:
        return i
    if tt[i] >= v:
        if tt[i + 1] < v:
            return i
        return find_time_index(v, tt, i + 1, stop)
    if tt[i - 1] >= v:
        return i - 1
    return find_time_index(v, tt, start, i)


def compute_U(T, MU, LA, PSI, RHO, SIGMA, nsteps=N_U_STEPS):
    """
    Calculates a function get_U which for a given time t: 0 <= t <= T, would return
    an array of unobserved probabilities [U_1(t), ..., U_m(t)].

    U_k(t) are calculated by
    (1) solving their ODEs numerically for an array tt of nsteps times equally spaced between t=T and t=0,
    producing an array of solutions sol of length nstep (for corresponding times in tt)s.
    (2) creating a linear approximation which for a given time t (2a) find an index i such that tt[i] >= t > tt[i+1];
    (2b) returns sol[i + 1] + (sol[i] - sol[i + 1]) * (tt[i] - t) / (tt[i] - tt[i + 1]).


    :param T: time at end of the sampling period
    :param MU: an array of state transition rates
    :param LA: an array of transmission rates
    :param PSI: an array of removal rates
    :param RHO: an array of sampling probabilities
    :param SIGMA: an array of rate sums: MU.sum(axis=1) + LA.sum(axis=1) + PSI
    :return: a function that for a given time t returns the array of corresponding unsampled probabilities:
        t ->  [U_1(t), ..., U_m(t)].
    """
    tt = np.linspace(T, 0, nsteps)
    y0 = np.ones(LA.shape[0], np.float64)
    PSI_NOT_RHO = PSI * (1 - RHO)

    def pdf_U(U, t):
        dU = (SIGMA - LA.dot(U)) * U - MU.dot(U) - PSI_NOT_RHO
        return dU

    sol = odeint(pdf_U, y0, tt)

    def get_U(t):
        t = max(0, min(t, T))
        tt_len = len(tt)
        i = find_time_index(t, tt, 0, tt_len)
        sol_prev = sol[i, :]
        if i == (tt_len - 1) or t == tt[i]:
            return sol_prev
        sol_next = sol[i + 1, :]
        if t == tt[i + 1]:
            return sol_next
        return sol_next + (sol_prev - sol_next) * (tt[i] - t) / (tt[i] - tt[i + 1])

    return get_U


def find_index_within_bounds(sol, start, stop, upper_bound, lower_bound=0):
    """
    Find an index i in a given array sol (of shape n x m),
    such that all the values of sol[i, :] are withing the given bounds
    and at least one of them is above the lower bound.

    Note that such index might not be unique. The search starts with the middle of the array
    and if needed proceeds with the binary search in the lower half of the array.

    :param sol: an array where to search
    :param start: start position in the array (inclusive)
    :param stop: stop position in the array (exclusive)
    :param upper_bound: upper bound
    :param lower_bound: lower bound
    :return: the index i that satisfies the above conditions
    """
    i = start + ((stop - start) // 2)
    if i == start or i == stop - 1:
        return i
    value = sol[i, :]
    if np.all(value >= lower_bound) and np.all(value <= upper_bound) and np.any(value > lower_bound):
        return i
    return find_index_within_bounds(sol, start, i, upper_bound)


def get_P(ti, l, t0, get_U, MU, LA, SIGMA):
    """
    Calculates P_{kl}^{(i)}(t0) for k in 1:m, where the initial condition is specified at time ti >= t0 (time of node i):
    P_{kl}^{(i)}(ti) = 0 for all k=l;
    P_{ll}^{(i)}(ti) = 1.

    :param ti: time for the initial condition (at node i)
    :param l: state of node i (the only state for which the initial condition is non-zero)
    :param t0: time to calculate the values at (t0 <= ti)
    :param get_U: a function to calculate an array of unsampled probabilities for a given time: t -> [U_1, .., U_m]
    :param MU: an array of state transition rates
    :param LA: an array of transmission rates
    :param SIGMA:  an array of rate sums: MU.sum(axis=1) + LA.sum(axis=1) + PSI, where PSI is the array of removal rates
    :return: a tuple containing an array of (potentially rescaled) branch evolution probabilities at time t0:
        [CP^{(i)}_{0l}(t0), .., CP^{(i)}_{ml}(t0)] and a log of the scaling factor: logC
    """

    def pdf_Pany_l(P, t):
        U = get_U(t)
        return (SIGMA - LA.dot(U)) * P - (MU + U * LA).dot(P)

    y0 = np.zeros(LA.shape[0], np.float64)
    y0[l] = 1

    nsteps = 100
    tt = np.linspace(ti, t0, nsteps)
    sol = odeint(pdf_Pany_l, y0, tt)

    # If there was an underflow during P_{kl}^{(i)}(t) calculations, we find a time tt[i] before the problem happened
    # and use its values sol[i, :] as new initial values for a rescaled ODE calculation,
    # which we solve for CP_{kl}^{(i)}(t). The new initial values become:
    # CP_{kl}^{(i)}(tt[i]) = C sol[i, k],
    # where C = 1 / min_positive(sol[i, :]).
    cs = [1]
    while np.any(sol[-1, :] < 0) or np.all(sol[-1, :] == 0) or np.any(sol[-1, :] > np.prod(cs)):
        i = find_index_within_bounds(sol, 0, len(tt), np.prod(cs))
        if i == 0:
            nsteps *= 10
            tt = np.linspace(ti, t0, nsteps)
            sol = odeint(pdf_Pany_l, y0, tt)
            continue

        vs = sol[i, :]

        c = 1 / min(sol[i, sol[i, :] > 0])
        cs.append(c)
        y0 = vs * c
        ti = tt[i]
        nsteps = 100
        tt = np.linspace(ti, t0, nsteps)
        sol = odeint(pdf_Pany_l, y0, tt)

    return np.maximum(sol[-1, :], 0), np.log(cs).sum()


def rescale_log_array(loglikelihood_array):
    """
    Rescales the likelihood array if it gets too small/large, by multiplying it by a factor of 10.
    :param loglikelihood_array: numpy array containing the loglikelihood to be rescaled
    :return: float, factor of e by which the likelihood array has been multiplied.
    """

    max_limit = MAX_VALUE
    min_limit = MIN_VALUE

    non_zero_loglh_array = loglikelihood_array[loglikelihood_array > -np.inf]
    min_lh_value = np.min(non_zero_loglh_array)
    max_lh_value = np.max(non_zero_loglh_array)

    factors = 0

    if max_lh_value > max_limit:
        factors = max_limit - max_lh_value - 1
    elif min_lh_value < min_limit:
        factors = min(min_limit - min_lh_value + 1, max_limit - max_lh_value - 1)
    loglikelihood_array += factors
    return factors


def loglikelihood_known_states(forest, T, MU, LA, PSI, RHO, u=0, threads=1):
    """
    Calculates loglikelihood for a given forest of trees,
    whose nodes are annotated with their state ids in 1:m (via feature STATE_K),
    and treir times (via feature TI), under given MTBD model parameters.

    :param forest: a list of ete3.Tree trees
    :param T: time at end of the sampling period
    :param MU: an array of state transition rates
    :param LA: an array of transmission rates
    :param PSI: an array of removal rates
    :param RHO: an array of sampling probabilities
    :param u: number of hidden trees, where no tip got sampled
    :return: the value of loglikelihood
    """
    PI = state_frequencies(MU, LA, PSI)
    PSI_RHO = PSI * RHO
    SIGMA = MU.sum(axis=1) + LA.sum(axis=1) + PSI
    get_U = compute_U(T, MU=MU, LA=LA, PSI=PSI, RHO=RHO, SIGMA=SIGMA)


    if threads < 1:
        threads = max(os.cpu_count(), 1)

    def _work(n):
        ti = getattr(n, TI)
        P, factors = get_P(t0=ti - n.dist, l=getattr(n, STATE_K), ti=ti, get_U=get_U, MU=MU, LA=LA, SIGMA=SIGMA)
        # if np.all(P == 0):
        #     plot_P(getattr(n, STATE_K), get_U, ti, ti - n.dist, MU, LA, PSI)
        return n, (P, factors)

    if threads > 1:
        with ThreadPool(processes=threads) as pool:
            acr_results = \
                pool.map(func=_work, iterable=chain(*(tree.traverse() for tree in forest)))
    else:
        acr_results = [_work(n) for n in chain(*(tree.traverse() for tree in forest))]

    n2p = {n: _ for (n, _) in acr_results}
    # for tree in forest:
    #     for n in tree.traverse('preorder'):
    #         ti = getattr(n, TI)
    #         n2p[n] = get_P(t0=ti - n.dist, l=getattr(n, STATE_K), ti=ti, get_U=get_U, MU=MU, LA=LA, SIGMA=SIGMA)
    #
    #         if np.all(n2p[n][0] == 0):
    #             plot_P(getattr(n, STATE_K), get_U, ti, ti - n.dist, MU, LA, PSI)
    #             return -np.inf

    res = 0 if not u else u * np.log(PI.dot(get_U(0)))
    for tree in forest:
        for n in tree.traverse('preorder'):
            k = getattr(n, STATE_K)
            if n.is_root() and n.ti:
                P, factors = n2p[n]
                res += np.log(PI.dot(P)) - factors
            if n.is_leaf():
                res += np.log(PSI_RHO[k])
            else:
                lc, rc = n.children
                (P_lc, factors_lc), (P_rc, factors_rc) = n2p[lc], n2p[rc]
                max_factors = max(factors_lc, factors_lc)
                diff_lc, diff_rc = max_factors - factors_lc, max_factors - factors_rc
                P_lc *= np.exp(diff_lc)
                P_rc *= np.exp(diff_rc)
                LA_k = LA[k, :]
                res += np.log(P_lc[k] * LA_k.dot(P_rc) + P_rc[k] * LA_k.dot(P_lc)) - 2 * max_factors
            if np.isnan(res):
                break

    if np.isnan(res):
        res = -np.inf
    return res #- len(forest) * np.log(1 - prob_evolve_unobserved)


def rescale_forest(forest):
    """
    Rescales forest branches and node times so that the average branch length in the rescaled forest is one.

    :param forest: a list of ete3.Tree trees
    :return: the scaling factor, which is equal to 1 / average_branch_length
    """
    br_lengths = []
    for tree in forest:
        for _ in tree.traverse():
            br_lengths.append(_.dist)
    avg_br_len = np.mean(br_lengths)
    scaling_factor = 1 / avg_br_len
    for tree in forest:
        for n in tree.traverse():
            n.dist *= scaling_factor
            n.add_feature(TI, getattr(n, TI) * scaling_factor)
    return scaling_factor


def optimize_likelihood_params(forest, model, T, optimise, u=0):
    """
    Optimizes the likelihood parameters for a given forest and a given MTBD model.


    :param forest: a list of ete3.Tree trees, annotated with node states and times via features STATE_K and TI.
    :param T: time at end of the sampling period
    :param model: MTBD model containing starting parameter values
    :param optimise: MTBD model whose rates indicate which parameters need to optimized:
        positive rates correspond to optimized parameters
    :param u: number of hidden trees, where no tip got sampled
    :return: the values of optimised parameters and the corresponding loglikelihood: (MU, LA, PSI, RHO, best_log_lh)
    """
    scaling_factor = 1
    # scaling_factor = rescale_forest(forest)
    MU, LA, PSI, RHO = model.transition_rates, model.transmission_rates, model.removal_rates, model.ps
    MU, LA, PSI = MU * scaling_factor, LA * scaling_factor, PSI * scaling_factor
    print('rescaled everything by {}'.format(scaling_factor))

    bounds = []
    opt_transition_rates = (optimise.transition_rates > 0)
    opt_transmission_rates = (optimise.transmission_rates > 0)
    opt_removal_rates = (optimise.removal_rates > 0)
    opt_ps = (optimise.ps > 0)
    n_mu = opt_transition_rates.sum()
    n_la = opt_transmission_rates.sum()
    n_psi = opt_removal_rates.sum()
    n_p = opt_ps.sum()

    bounds.extend([[1e-3, 10e3]] * (n_mu + n_la + n_psi))
    bounds.extend([[1e-3, 1 - 1e-3]] * n_p)
    bounds = np.array(bounds, np.float64)

    def get_real_params_from_optimised(ps):
        ps = np.maximum(np.minimum(ps, bounds[:, 1]), bounds[:, 0])
        MU[opt_transition_rates] = ps[:n_mu]
        LA[opt_transmission_rates] = ps[n_mu: n_mu + n_la]
        PSI[opt_removal_rates] = ps[n_mu + n_la: n_mu + n_la + n_psi]
        RHO[opt_ps] = ps[n_mu + n_la + n_psi:]

    def get_optimised_params_from_real():
        return np.append(
            np.append(
                np.append(MU[opt_transition_rates],
                          LA[opt_transmission_rates]),
                PSI[opt_removal_rates]),
            RHO[opt_ps])

    def get_v(ps):
        if np.any(pd.isnull(ps)):
            return np.nan
        get_real_params_from_optimised(ps)
        res = loglikelihood_known_states(forest, T, MU, LA, PSI, RHO, u)
        print("mu=", MU[0, 1] / scaling_factor, "la=", LA[1, 0] / scaling_factor, "psi=", PSI[1] / scaling_factor, "p=", RHO[1], "\t-->\t", res)
        return -res

    x0 = get_optimised_params_from_real()
    best_log_lh = -get_v(x0)

    for i in range(10):
        if i == 0:
            vs = x0
        else:
            vs = np.random.uniform(bounds[:, 0], bounds[:, 1])

        fres = minimize(get_v, x0=vs, method='TNC', bounds=bounds)
        if fres.success and not np.any(np.isnan(fres.x)) and -fres.fun >= best_log_lh:
            x0 = fres.x
            best_log_lh = -fres.fun
            break
        print('Attempt {} of trying to optimise the parameters: {}.'.format(i, -fres.fun))
    get_real_params_from_optimised(x0)
    MU, LA, PSI = MU / scaling_factor, LA / scaling_factor, PSI / scaling_factor
    return MU, LA, PSI, RHO, best_log_lh


if __name__ == '__main__':
    tree = Tree('/home/azhukova/projects/bdei/simulations/medium/trees/tree.1.nwk')
    T = 0
    for n in tree.traverse('preorder'):
        ti = (0 if n.is_root() else getattr(n.up, TI)) + n.dist
        n.add_feature(TI, ti)
        T = max(T, ti)
        n.add_feature(STATE_K, 1)

    forest = [tree]
    p = 0.6225767763228239
    real_mu = 0.2523725112488919
    real_la = 0.907081384137969
    real_psi = 0.2692907505391973

    rate = initial_rate_guess(forest).pop()
    model = BirthDeathExposedInfectiousModel(mu=rate, la=rate, psi=rate, p=p)
    optimise = BirthDeathExposedInfectiousModel(mu=1, la=1, psi=1, p=0)
    real_model = BirthDeathExposedInfectiousModel(mu=real_mu, la=real_la, psi=real_psi, p=p)

    print('Real values and likelihood are:\n', "mu=", real_mu, "la=", real_la, "psi=", real_psi, "p=", p, "\t-->\t",
          loglikelihood_known_states(forest, T, real_model.transition_rates, real_model.transmission_rates,
                                     real_model.removal_rates, real_model.ps))
    # # model = real_model
    MU, LA, PSI, RHO, lk = optimize_likelihood_params(forest, model, T, optimise)
    print("mu=", MU[0, 1], "la=", LA[1, 0], "psi=", PSI[1], "p=", RHO[1], "lk=", lk)
