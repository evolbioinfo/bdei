import os
from itertools import chain
from multiprocessing.pool import ThreadPool

import numpy as np

RTOL = 100 * np.finfo(np.float64).eps
import pandas as pd
from scipy.integrate import odeint
from scipy.optimize import minimize, fsolve

N_U_STEPS = int(1e6)

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


def compute_U(T, MU, LA, PSI, RHO, SIGMA=None, nsteps=N_U_STEPS):
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
    if SIGMA is None:
        SIGMA = MU.sum(axis=1) + LA.sum(axis=1) + PSI

    tt = np.linspace(T, 0, nsteps)
    y0 = np.ones(LA.shape[0], np.float64)
    PSI_NOT_RHO = PSI * (1 - RHO)

    def pdf_U(U, t):
        dU = (SIGMA - LA.dot(U)) * U - MU.dot(U) - PSI_NOT_RHO
        return dU

    sol = odeint(pdf_U, y0, tt, rtol=RTOL)
    sol = np.maximum(sol, 0)

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
    y0 = np.zeros(LA.shape[0], np.float64)
    y0[l] = 1

    if t0 == ti:
        return y0, 0

    def pdf_Pany_l(P, t):
        U = get_U(t)
        return (SIGMA - LA.dot(U)) * P - (MU + U * LA).dot(P)

    nsteps = 2
    tt = np.linspace(ti, t0, nsteps)
    sol = odeint(pdf_Pany_l, y0, tt, rtol=RTOL)

    while np.any(sol[-1, :] < 0) or np.all(sol[-1, :] == 0) or np.any(sol[-1, :] > y0[l]):
        y0[l] *= 2
        sol = odeint(pdf_Pany_l, y0, tt, rtol=RTOL)
    return np.minimum(np.maximum(sol[-1, :] / y0[l], 0), 1)


def get_rescale_factors(loglikelihood_array):
    """
    Checks if the input (log)array is too small/large, and return a factor of e to multiply it by.

    :param loglikelihood_array: numpy array containing the loglikelihood to be rescaled
    :return: float, factor of e by which the likelihood array should be multiplied.
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
    return factors


def loglikelihood_known_states(forest, T, MU, LA, PSI, RHO, u=-1, threads=1):
    """
    Calculates loglikelihood for a given forest of trees,
    whose nodes are annotated with their state ids in 1:m (via feature STATE_K),
    and their times (via feature TI), under given MTBD model parameters.

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
    LOG_PSI_RHO = np.log(PSI) + np.log(RHO)
    SIGMA = MU.sum(axis=1) + LA.sum(axis=1) + PSI
    get_U = compute_U(T, MU=MU, LA=LA, PSI=PSI, RHO=RHO, SIGMA=SIGMA)

    if threads < 1:
        threads = max(os.cpu_count(), 1)

    def _work(n):
        ti = getattr(n, TI)
        P = get_P(t0=ti - n.dist, l=getattr(n, STATE_K), ti=ti, get_U=get_U, MU=MU, LA=LA, SIGMA=SIGMA)
        return n, P

    if threads > 1:
        with ThreadPool(processes=threads) as pool:
            acr_results = \
                pool.map(func=_work, iterable=chain(*(tree.traverse() for tree in forest)))
    else:
        acr_results = [_work(n) for n in chain(*(tree.traverse() for tree in forest))]

    n2p = dict(acr_results)

    hidden_probability = PI.dot(get_U(0))
    if u < 0:
        u = len(forest) * hidden_probability / (1 - hidden_probability)

    res = 0 if not u else u * np.log(hidden_probability)
    for tree in forest:
        for n in tree.traverse('preorder'):
            k = getattr(n, STATE_K)
            if n.is_root() and n.dist:
                P = n2p[n]
                res += np.log(PI.dot(P))
            if n.is_leaf():
                res += LOG_PSI_RHO[k]
            else:
                LA_k = LA[k, :]
                lc, rc = n.children
                P_lc, P_rc = n2p[lc], n2p[rc]
                # res += np.log(P_lc[k] * LA_k.dot(P_rc) + P_rc[k] * LA_k.dot(P_lc))
                left_donor = np.log(P_lc[k]) + np.log(LA_k.dot(P_rc))
                right_donor = np.log(P_rc[k]) + np.log(LA_k.dot(P_lc))
                max_factors = max(get_rescale_factors(left_donor), get_rescale_factors(right_donor))
                res += np.log(np.exp(left_donor + max_factors) + np.exp(right_donor + max_factors)) - max_factors

            if np.isnan(res):
                break

    if np.isnan(res):
        res = -np.inf
    return res


def optimize_likelihood_params(forest, model, T, optimise, u=-1, bounds=None):
    """
    Optimizes the likelihood parameters for a given forest and a given MTBD model.


    :param forest: a list of ete3.Tree trees, annotated with node states and times via features STATE_K and TI.
    :param T: time at end of the sampling period
    :param model: MTBD model containing starting parameter values
    :param optimize: MTBD model whose rates indicate which parameters need to optimized:
        positive rates correspond to optimized parameters
    :param u: number of hidden trees, where no tip got sampled
    :return: the values of optimized parameters and the corresponding loglikelihood: (MU, LA, PSI, RHO, best_log_lh)
    """
    MU, LA, PSI, RHO = model.transition_rates, model.transmission_rates, model.removal_rates, model.ps

    opt_transition_rates = (optimise.transition_rates > 0)
    opt_transmission_rates = (optimise.transmission_rates > 0)
    opt_removal_rates = (optimise.removal_rates > 0)
    opt_ps = (optimise.ps > 0)
    n_mu = opt_transition_rates.sum()
    n_la = opt_transmission_rates.sum()
    n_psi = opt_removal_rates.sum()
    n_p = opt_ps.sum()

    if bounds is None:
        bounds = []
        bounds.extend([[1e-3, 1e2]] * (n_mu + n_la + n_psi))
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
        # print(ps, "\t-->\t", res)
        return -res

    x0 = get_optimised_params_from_real()
    best_log_lh = -get_v(x0)

    for i in range(5):
        if i == 0:
            vs = x0
        else:
            vs = np.random.uniform(bounds[:, 0], bounds[:, 1])

        fres = minimize(get_v, x0=vs, method='L-BFGS-B', bounds=bounds)
        if not np.any(np.isnan(fres.x)):
            if -fres.fun >= best_log_lh:
                x0 = fres.x
                best_log_lh = -fres.fun
            print('Attempt {} of trying to optimise the parameters: {} -> {}.'
                  .format(i, fres.x, -fres.fun))
    get_real_params_from_optimised(x0)
    return MU, LA, PSI, RHO, best_log_lh
