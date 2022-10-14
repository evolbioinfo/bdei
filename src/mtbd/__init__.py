import numpy as np


def initial_rate_guess(forest, mu=None, la=None, psi=None):

    fixed_rates = []
    for rate in (mu, la, psi):
        if rate and rate > 0:
            fixed_rates.append(rate)
    if fixed_rates:
        avg_rate = np.average(fixed_rates)
        min_rate = avg_rate
        max_rate = avg_rate
    else:
        internal_is, external_is = [], []
        internal_es, external_es = [], []
        for tree in forest:
            for n in tree.traverse('preorder'):
                if not n.is_leaf():
                    i, e = sorted(n.children, key=lambda c: c.dist)
                    (internal_is if not i.is_leaf() else external_is).append(i)
                    (internal_es if not i.is_leaf() else external_es).append(e)
        i_s_times = np.array([_.dist for _ in external_is])
        i_tr_times = np.array([_.dist for _ in internal_is])
        e_s_times = np.array([_.dist for _ in external_es])
        e_tr_times = np.array([_.dist for _ in internal_es])
        i_s_time = np.median(i_s_times)
        e_s_time = np.median(e_s_times)
        # if it is a corner case when we only have tips, let's use sampling times
        i_tr_time = np.median(i_tr_times) if len(i_tr_times) else i_s_time
        e_tr_time = np.median(e_tr_times) if len(e_tr_times) else e_s_time
        mu_times = []
        if e_s_time > i_s_time:
            mu_times.append(e_s_time - i_s_time)
        if e_tr_time > i_tr_time:
            mu_times.append(e_tr_time - i_tr_time)
        mu_time = np.average(mu_times) if mu_times else min(i_s_time, i_tr_time) * 0.01
        avg_rate = np.average([1 / i_s_time, 1 / i_tr_time, 1 / mu_time])
        min_rate = min([1 / i_s_time, 1 / i_tr_time, 1 / mu_time])
        max_rate = max([1 / i_s_time, 1 / i_tr_time, 1 / mu_time])

    return sorted({avg_rate, min_rate, max_rate})