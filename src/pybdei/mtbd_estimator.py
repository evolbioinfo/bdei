
import numpy as np
import pandas as pd
from ete3 import Tree
from matplotlib import pyplot as plt
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


def binary(v, tt, start, stop):
    i = start + ((stop - start) // 2)
    if i == start or i == stop - 1:
        return i
    if tt[i] >= v:
        if tt[i + 1] < v:
            return i
        return binary(v, tt, i + 1, stop)
    if tt[i - 1] >= v:
        return i - 1
    return binary(v, tt, start, i)


def compute_U(T, MU, LA, PSI, RHO, SIGMA, nsteps=N_U_STEPS):
    tt = np.linspace(T, 0, nsteps)
    y0 = np.ones(LA.shape[0], np.float64)
    PSI_NOT_RHO = PSI * (1 - RHO)

    def pdf_U(U, t):
        dU = (SIGMA - LA.dot(U)) * U \
               - MU.dot(U) \
               - PSI_NOT_RHO
        return dU

    sol = odeint(pdf_U, y0, tt)

    def get_U(t):
        t = max(0, min(t, T))
        tt_len = len(tt)
        i = binary(t, tt, 0, tt_len)
        sol_prev = sol[i, :]
        if i == (tt_len - 1) or t == tt[i]:
            return sol_prev
        sol_next = sol[i + 1, :]
        if t == tt[i + 1]:
            return sol_next
        return sol_next + (sol_prev - sol_next) * (tt[i] - t) / (tt[i] - tt[i + 1])

    return get_U


def get_logP_plus_1(t, l, t0, get_U, MU, LA, SIGMA):
    y0 = np.zeros(LA.shape[0], np.float64)
    y0[l] = np.log(2)

    def pdf_logPany_l(logP_plus_1, t):
        U = get_U(t)
        A = SIGMA - LA.dot(U)
        B = (MU + U * LA)
        Q = np.exp(logP_plus_1)
        P = Q - 1
        return A - (B.dot(P) + A) / (P + 1)

    sol = odeint(pdf_logPany_l, y0, [t0, t])
    return sol[1, :]


def plot_P_simple(k, get_U, ti, t0, MU, LA, PSI):
    t = np.linspace(ti, t0, 1001)

    def pdf_Pany_l(P, t):
        U = get_U(t)
        return (MU.sum(axis=1) + LA.dot(1 - U) + PSI) * P - (MU + U * LA).dot(P)

    y0 = np.zeros(len(PSI), np.float64)
    y0[k] = 1
    sol = odeint(pdf_Pany_l, y0, t)

    plt.plot(t, sol[:, 0], 'b', label='logP_E(t)')
    plt.plot(t, sol[:, 1], 'g', label='logP_I(t)')
    plt.legend(loc='best')
    plt.xlabel('t')
    plt.grid()
    plt.show()


def get_P(ti, l, t0, get_U, MU, LA, SIGMA):

    def pdf_Pany_l(P, t):
        U = get_U(t)
        return (SIGMA - LA.dot(U)) * P - (MU + U * LA).dot(P)

    y0 = np.zeros(LA.shape[0], np.float64)
    y0[l] = 1

    nsteps = 100
    tt = np.linspace(ti, t0, nsteps)
    sol = odeint(pdf_Pany_l, y0, tt)

    cs = [1]
    while np.any(sol[-1, :] < 0) or np.all(sol[-1, :] == 0) or np.any(sol[-1, :] > np.prod(cs)):
        i = find_positive(sol, 0, len(tt), np.prod(cs))
        if i == 0:
            nsteps *= 10
            tt = np.linspace(ti, t0, nsteps)
            sol = odeint(pdf_Pany_l, y0, tt)
            continue

        vs = sol[i, :]
        # print(ti, '->', tt[i])
        # print(y0, '->', vs)

        c = 1 / min(sol[i, sol[i, :] > 0])
        cs.append(c)
        y0 = vs * c
        ti = tt[i]
        nsteps = 100
        tt = np.linspace(ti, t0, nsteps)
        sol = odeint(pdf_Pany_l, vs, tt)

    return np.maximum(sol[-1, :], 0), np.log(cs).sum()


def rescale_log_array(loglikelihood_array):
    """
    Rescales the likelihood array if it gets too small/large, by multiplying it by a factor of 10.
    :param loglikelihood_array: numpy array containing the loglikelihood to be rescaled
    :return: float, factor of 10 by which the likelihood array has been multiplied.
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


def loglikelihood_known_states_log(forest, T, MU, LA, PSI, RHO, u=0):
    """
    Each tree node is annotated with ti and state.

    :param forest:
    :param model:
    :return:
    """
    PI = state_frequencies(MU, LA, PSI)
    PSI_RHO = PSI * RHO
    SIGMA = MU.sum(axis=1) + LA.sum(axis=1) + PSI
    m = len(PI)
    get_U = compute_U(T, MU=MU, LA=LA, PSI=PSI, RHO=RHO, SIGMA=SIGMA)
    plot_U(get_U, T)

    res = [0 if not u else u * np.log(PI.dot(get_U(0)))]

    for tree in forest:
        for n in tree.traverse('preorder'):
            k = getattr(n, STATE_K)
            ti = getattr(n, TI)
            t0 = getattr(n.up, TI) if not n.is_root() else 0
            logP_plus_1 = get_logP_plus_1(t0, k, ti, get_U, MU, LA, SIGMA)
            factors = 0 #rescale_log_array(logP_plus_1)
            P = np.maximum(np.exp(logP_plus_1) - np.exp(factors), 0)

            if np.all(P == 0):
                plot_P(k, get_U, ti, t0, MU, LA, PSI)
                print(np.exp(logP_plus_1) - np.exp(factors))



    for tree in forest:
        for n in tree.traverse('preorder'):
            k = getattr(n, STATE_K)
            ti = getattr(n, TI)
            if n.is_root() and n.ti:
                logP_plus_1 = get_logP_plus_1(0, k, ti, get_U, MU, LA, SIGMA)
                factors = rescale_log_array(logP_plus_1)
                P = np.maximum(np.exp(logP_plus_1) - np.exp(factors), 0)
                res.append(np.log(PI.dot(P)) - factors)
            if n.is_leaf():
                res.append(np.log(PSI_RHO[k]))
            else:
                lc, rc = n.children
                ti_lc, ti_rc = getattr(lc, TI), getattr(rc, TI)
                s_lc, s_rc = getattr(lc, STATE_K), getattr(rc, STATE_K)
                logP_plus_1_lc, logP_plus_1_rc = get_logP_plus_1(ti, s_lc, ti_lc, get_U, MU, LA, SIGMA), \
                                                 get_logP_plus_1(ti, s_rc, ti_rc, get_U, MU, LA, SIGMA)
                logP_plus_1_both = np.append(logP_plus_1_lc, logP_plus_1_rc)
                factors = rescale_log_array(logP_plus_1_both)
                P_both = np.maximum(np.exp(logP_plus_1_both) - np.exp(factors), 0)
                P_lc, P_rc = P_both[:m], P_both[m:]
                LA_k = LA[k, :]
                res.append(np.log(P_lc[k] * LA_k.dot(P_rc) + P_rc[k] * LA_k.dot(P_lc)) - 2 * factors)
            if np.isnan(res[-1]) or res[-1] == -np.inf:
                print(LA[k, :], n.dist, getattr(n, TI), getattr(n, STATE_K), res[-1])
                if n.children:
                    print(P_lc, LA_k.dot(P_rc), P_rc[k], LA_k.dot(P_lc))
                    print(lc.dist, rc.dist, getattr(lc, TI), getattr(lc, STATE_K), getattr(rc, TI), getattr(rc, STATE_K))
                break
    print(res)
    if np.any(np.isnan(res)):
        print(res)
        res = -np.inf
    return sum(sorted(res))


def loglikelihood_known_states(forest, T, MU, LA, PSI, RHO, u=0):
    """
    Each tree node is annotated with ti and state.

    :param forest:
    :param model:
    :return:
    """
    PI = state_frequencies(MU, LA, PSI)
    PSI_RHO = PSI * RHO
    SIGMA = MU.sum(axis=1) + LA.sum(axis=1) + PSI
    get_U = compute_U(T, MU=MU, LA=LA, PSI=PSI, RHO=RHO, SIGMA=SIGMA)

    res = 0 if not u else u * np.log(PI.dot(get_U(0)))

    n2p = {}
    for tree in forest:
        for n in tree.traverse('preorder'):
            k = getattr(n, STATE_K)
            ti = getattr(n, TI)
            t0 = ti - n.dist
            P, factors = get_P(t0=t0, l=k, ti=ti, get_U=get_U, MU=MU, LA=LA, SIGMA=SIGMA)

            n2p[n] = P, factors

            if np.all(P == 0):
                plot_P_simple(k, get_U, ti, t0, MU, LA, PSI)

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
                if (P_lc[k] * LA_k.dot(P_rc) <= 0) or (P_rc[k] * LA_k.dot(P_lc) <= 0) or (P_lc[k] * LA_k.dot(P_rc) + P_rc[k] * LA_k.dot(P_lc)) <=0:
                    print(P_lc, P_rc, diff_lc, diff_rc)
                res += np.log(P_lc[k] * LA_k.dot(P_rc) + P_rc[k] * LA_k.dot(P_lc)) - 2 * max_factors
            if np.isnan(res):
                break

    if np.isnan(res):
        res = -np.inf
    return res #- len(forest) * np.log(1 - prob_evolve_unobserved)


def rescale_forest(forest):
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


def rescale_model(model, scaling_factor):
    model.transition_rates = scaling_factor * model.transition_rates
    model.transmission_rates = scaling_factor * model.transmission_rates
    model.removal_rates = scaling_factor * model.removal_rates


def optimize_likelihood_params(forest, model, T, optimise, u=0):
    """
    Optimizes the likelihood parameters for a given forest.
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

    bounds.extend([[0.01, 10]] * (n_mu + n_la + n_psi))
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

        fres = minimize(get_v, x0=vs, method='L-BFGS-B', bounds=bounds)
        if fres.success and not np.any(np.isnan(fres.x)) and -fres.fun >= best_log_lh:
            x0 = fres.x
            best_log_lh = -fres.fun
            break
        print('Attempt {} of trying to optimise the parameters: {}.'.format(i, -fres.fun))
    get_real_params_from_optimised(x0)
    MU, LA, PSI = MU / scaling_factor, LA / scaling_factor, PSI / scaling_factor
    return MU, LA, PSI, RHO, best_log_lh


def plot_U(get_U, T):
    t = np.linspace(0, T, 1001)
    esol_E = [get_U(_)[0] for _ in t]
    esol_I = [get_U(_)[1] for _ in t]
    plt.plot(t, esol_E, 'b', label='U_E(t)')
    plt.plot(t, esol_I, 'g', label='U_I(t)')
    plt.legend(loc='best')
    plt.xlabel('t')
    plt.grid()
    plt.show()


def find_positive(sol, start, stop, upper_bound):
    i = start + ((stop - start) // 2)
    if i == start or i == stop - 1:
        return i
    value = sol[i, :]
    if np.all(value >= 0) and np.all(value <= upper_bound) and np.any(value > 0):
        return i
    return find_positive(sol, start, i, upper_bound)


def plot_P(k, get_U, ti, t0, MU, LA, PSI):
    t = np.linspace(ti, t0, 1001)

    def pdf_Pany_l(P, t):
        U = get_U(t)
        return (MU.sum(axis=1) + LA.dot(1 - U) + PSI) * P - (MU + U * LA).dot(P)

    y0 = np.zeros(len(PSI), np.float64)
    y0[k] = 1
    sol = odeint(pdf_Pany_l, y0, t)

    cs = [1]
    while np.all(sol[-1, :] <= 0):
        i = find_positive(sol, 0, len(t))
        vs = sol[i, :]
        c = 1 / min(sol[i, sol[i, :] > 0])
        cs.append(c)
        vs *= c
        print(ti, '->', t[i])
        t = np.linspace(t[i], t0, 1001)
        sol = odeint(pdf_Pany_l, vs, t)
    print(np.prod(cs))


    plt.plot(t, sol[:, 0], 'b', label='logP_E(t)')
    plt.plot(t, sol[:, 1], 'g', label='logP_I(t)')
    plt.legend(loc='best')
    plt.xlabel('t')
    plt.grid()
    plt.show()

def plot_P_simple(k, get_U, ti, t0, MU, LA, PSI):
    t = np.linspace(ti, t0, 1001)

    def pdf_Pany_l(P, t):
        U = get_U(t)
        return (MU.sum(axis=1) + LA.dot(1 - U) + PSI) * P - (MU + U * LA).dot(P)

    y0 = np.zeros(len(PSI), np.float64)
    y0[k] = 1
    sol = odeint(pdf_Pany_l, y0, t)

    plt.plot(t, sol[:, 0], 'b', label='logP_E(t)')
    plt.plot(t, sol[:, 1], 'g', label='logP_I(t)')
    plt.legend(loc='best')
    plt.xlabel('t')
    plt.grid()
    plt.show()


def plot_P_logP():
    model = BirthDeathExposedInfectiousModel(mu=0.7299648936648243, la=1.421061063434435, psi=0.6069550954814479,
                                                  p=0.44732297990706804)

    T = 20
    t0 = 5
    t = np.linspace(t0, 0, 1001)
    MU, LA, PSI, RHO = model.transition_rates, model.transmission_rates, model.removal_rates, model.ps
    SIGMA = MU.sum(axis=1) + LA.sum(axis=1) + PSI
    get_U = compute_U(T, MU=MU, LA=LA, PSI=PSI, RHO=RHO, SIGMA=SIGMA)

    def pdf_Pany_l(P, t):
        U = get_U(t)
        dP = (MU.sum(axis=1) + LA.dot(1 - U) + PSI) * P \
             - (MU + U * LA).dot(P)
        return dP

    def pdf_logPany_l(logP_plus_1, t):
        U = get_U(t)
        A = (MU.sum(axis=1) + LA.dot(1 - U) + PSI)
        B = (MU + U * LA)
        Q = np.exp(logP_plus_1)
        P = Q - 1
        dlogP_plus_1 = A - (B.dot(P) + A) / (P + 1)
        return dlogP_plus_1

    y0 = np.zeros(len(model.states), np.float64)
    y0[1] = 1
    sol = odeint(pdf_Pany_l, y0, t)

    logy0 = np.zeros(len(model.states), np.float64)
    logy0[1] = np.log(2)
    logsol = odeint(pdf_logPany_l, logy0, t)
    logsol = np.exp(logsol) - 1

    plt.plot(t, sol[:, 0], 'b', label='P_E(t)')
    plt.plot(t, sol[:, 1], 'g', label='P_I(t)')
    plt.plot(t, logsol[:, 0], 'bv', label='logP_E(t)')
    plt.plot(t, logsol[:, 1], 'gv', label='logP_I(t)')
    plt.legend(loc='best')
    plt.xlabel('t')
    plt.grid()
    plt.show()


if __name__ == '__main__':

    # plot_P()
    # exit()
    tree = Tree('/home/azhukova/projects/bdei/simulations/medium/trees/tree.1.nwk')
    T = 0
    for n in tree.traverse('preorder'):
        ti = (0 if n.is_root() else getattr(n.up, TI)) + n.dist
        n.add_feature(TI, ti)
        T = max(T, ti)
        n.add_feature(STATE_K, 1)

    forest = [tree]
    p = 0.6225767763228239
    rate = initial_rate_guess(forest).pop()
    model = BirthDeathExposedInfectiousModel(mu=rate, la=rate, psi=rate / 2, p=p)
    optimise = BirthDeathExposedInfectiousModel(mu=1, la=1, psi=1, p=0)
    real_model = BirthDeathExposedInfectiousModel(mu=0.2523725112488919, la=0.907081384137969, psi=0.2692907505391973,
                                                  p=p)
    print('Real likelihood is', loglikelihood_known_states_log(forest, T, real_model.transition_rates, real_model.transmission_rates, real_model.removal_rates, real_model.ps))
    # prob_model = BirthDeathExposedInfectiousModel(mu=27.30227398105135, la=27.30227398105135, psi=27.30227398105135,
    #                                               p=p)
    # print('Prob likelihood is', loglikelihood_known_states(forest, T, prob_model.transition_rates, prob_model.transmission_rates, prob_model.removal_rates, prob_model.ps))
    # exit()
    # # model = real_model
    MU, LA, PSI, RHO, lk = optimize_likelihood_params(forest, model, T, optimise)
    print("mu=", MU[0, 1], "la=", LA[1, 0], "psi=", PSI[1], "p=", RHO[1], "lk=", lk)
