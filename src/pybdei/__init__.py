import os
from collections import namedtuple

import _pybdei
import numpy as np
from ete3 import Tree

SAMPLING_PERIOD_LENGTH = 'T'

ERRORS = 0
WARNINGS = 1
INFO = 2
DEBUG = 3

PYBDEI_VERSION = 0.9

PS = (0.1, 0.4, 0.7)

MAX, MIN, MEDIAN, MEAN = 'max', 'min', 'median', 'mean'


BDEI_result = namedtuple('BDEI_result', ['mu', 'la', 'psi', 'p', 'mu_CI', 'la_CI', 'psi_CI', 'p_CI',
                                         'R_naught', 'incubation_period', 'infectious_time'])
BDEI_time = namedtuple('BDEI_time', ['CPU_time', 'iterations'])


def get_temp_file_name(prefix):
    return '{}.{}.temp'.format(prefix, np.random.random())


def parse_tree(nwk):
    """Parses an input tree in newick format and checks if it is binary and rooted"""
    tree = None
    for f in (0, 1, 2, 3, 5):
        try:
            tree = Tree(nwk, format=f)
            break
        except:
            pass
    if not tree:
        raise ValueError("Please make sure to use correct newick format and to specify branch lengths.")
    if len(tree.children) > 2:
        raise ValueError("The input tree must be rooted and binary, please root your tree.")
    for _ in tree.traverse():
        n_children = len(_.children)
        if not _.is_leaf() and n_children != 2:
            raise ValueError("The input tree must be binary, while your tree contains internal nodes with {} children, "
                             "please fix it.".format(n_children))
        if not _.is_root() and _.dist == 0:
            raise ValueError("The input tree must not contain non-root zero branches, while your tree does, "
                             "please fix it.")
    return tree


def parse_forest(nwk, T=0., u_policy=MEAN):
    forest = []
    with open(nwk, 'r') as f:
        forest = [parse_tree(_.strip('\n') + ';') for _ in f.read().split(';')[:-1]]
    if not forest:
        raise ValueError('Could not find any trees in your input file (check the newick format): {}'.format(nwk))
    Ts = [max(get_T(tree), float(getattr(tree, SAMPLING_PERIOD_LENGTH, 0))) for tree in forest]
    max_T = max(*Ts, T)
    for (tree, tree_T) in zip(forest, Ts):
        tree.add_feature(SAMPLING_PERIOD_LENGTH, tree_T if T <= 0 else max_T)
    u_T = max_T
    if T <= 0:
        if u_policy == MEDIAN:
            u_T = np.median(Ts)
        if u_policy == MEAN:
            u_T = np.mean(Ts)
        if u_policy == MIN:
            u_T = np.min(Ts)

    return forest, u_T


def save_forest(forest, nwk):
    with open(nwk, 'w+') as f:
        for tree in forest:
            _ = tree.write(format=5, format_root_node=True, features=[SAMPLING_PERIOD_LENGTH])\
                .replace('&&NHX:{}='.format(SAMPLING_PERIOD_LENGTH), '')
            f.write('{}\n'.format(_))


def get_T(tree):
    """
    Calculates the time between the tree start and the last sampled tip.
    """
    TIME = '_time'
    for n in tree.traverse('preorder'):
        if n.is_root():
            p_time = 0
        else:
            p_time = getattr(n.up, TIME)
            # delattr(n.up, TIME)
        n.add_feature(TIME, p_time + n.dist)
    time = 0
    for t in tree:
        time = max(time, getattr(t, TIME))
        # delattr(t, TIME)
    return time



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


def infer(nwk, start=None, upper_bounds=None, pi_E=-1,
          mu=-1, la=-1, psi=-1, p=-1, T=0., u=-1, CI_repetitions=0, threads=1, log_level=INFO, u_policy=MEAN, **kwargs):
    """Infer BDEI parameters from a phylogenetic tree."""

    forest, u_T = parse_forest(nwk, T=T, u_policy=u_policy)

    if isinstance(start, BDEI_result):
        start = np.array([start.mu, start.la, start.psi, start.p])
    if start is None:
        if not p or p <= 0 or p >= 1:
            rate = initial_rate_guess(forest, mu, la, psi).pop()
            starts = [[rate, rate, rate, pp] for pp in PS]
        else:
            starts = [[rate, rate, rate, p] for rate in initial_rate_guess(forest, mu, la, psi)]
    else:
        starts = [start]
    for s in starts:
        s[-1] = min(0.99, max(0.001, s[-1]))

    if upper_bounds is None:
        upper_bounds = np.array([np.inf, np.inf, np.inf, 1])
    elif isinstance(upper_bounds, BDEI_result):
        upper_bounds = np.array([upper_bounds.mu, upper_bounds.la, upper_bounds.psi, upper_bounds.p])
    else:
        upper_bounds = np.array(upper_bounds)
    if any(upper_bounds <= 0):
        raise ValueError('Upper bound must be positive.')
    upper_bounds[-1] = min(upper_bounds[-1], 1)

    if pi_E is not None and pi_E >= 0:
        if pi_E > 1:
            raise ValueError('Frequency pi_E must be between 0 and 1.')
    else:
        pi_E = -1

    something_is_fixed = False
    all_is_fixed = True
    if mu >= 0:
        upper_bounds[0] = mu
        for s in starts:
            s[0] = mu
        something_is_fixed = True
    else:
        all_is_fixed = False
    if la >= 0:
        upper_bounds[1] = la
        for s in starts:
            s[1] = la
        something_is_fixed = True
    else:
        all_is_fixed = False
    if psi >= 0:
        upper_bounds[2] = psi
        for s in starts:
            s[2] = psi
        something_is_fixed = True
    else:
        all_is_fixed = False
    if 0 < p <= 1:
        upper_bounds[3] = p
        for s in starts:
            s[3] = p
        something_is_fixed = True
    else:
        all_is_fixed = False

    if not something_is_fixed:
        raise ValueError('To be identifiable, the BDEI model needs one of its parameters to be fixed. '
                         'Please do so, using one of the following arguments: mu, la, psi, p.')
    if all_is_fixed:
        raise ValueError('At least one of the following arguments: mu, la, psi, p, should be left to be optimised.')
    starts = [np.minimum(s, upper_bounds * 0.999) for s in starts]
    nstarts = len(starts)
    starts = np.reshape(starts, (4 * nstarts,))

    temp_nwk = get_temp_file_name(nwk)
    save_forest(forest, temp_nwk)
    res = _pybdei.infer(f=temp_nwk, start=starts, ub=upper_bounds, pie=pi_E,
                        mu=mu, la=la, psi=psi, p=p, u=u, ut=u_T, nt=threads, nbiter=CI_repetitions,
                        debug=log_level, nstarts=nstarts)
    try:
        os.remove(temp_nwk)
    except OSError:
        pass

    return BDEI_result(mu=res[0], la=res[1], psi=res[2], p=res[3],
                       mu_CI=(res[4], res[5]) if CI_repetitions > 0 else None,
                       la_CI=(res[6], res[7]) if CI_repetitions > 0 else None,
                       psi_CI=(res[8], res[9]) if CI_repetitions > 0 else None,
                       p_CI=(res[10], res[11]) if CI_repetitions > 0 else None,
                       R_naught=res[1] / res[2], incubation_period=1 / res[0], infectious_time=1 / res[2]), \
           BDEI_time(CPU_time=res[13], iterations=res[14])


def get_loglikelihood(nwk, mu=-1, la=-1, psi=-1, p=-1, pi_E=-1, T=0., u=-1, threads=1, log_level=INFO, params=None,
                      u_policy=MEAN, **kwargs):
    """Calculate loglikelihood for given BDEI parameters from a phylogenetic tree."""
    if params is not None:
        if isinstance(params, BDEI_result):
            mu = params.mu
            la = params.la
            psi = params.psi
            p = params.p
        else:
            [mu, la, psi, p] = params
    for par in [mu, la, psi, p]:
        if par is None or par < 0:
            raise ValueError('All the parameters (mu, la, psi, p) must be specified, '
                             'either via dedicated arguments or via the params argument')

    forest, u_T = parse_forest(nwk, T=T, u_policy=u_policy)

    temp_nwk = get_temp_file_name(nwk)
    save_forest(forest, temp_nwk)
    res = _pybdei.likelihood(f=temp_nwk, mu=mu, la=la, psi=psi, p=p, pie=pi_E, u=u, ut=u_T, nt=threads, debug=log_level)
    try:
        os.remove(temp_nwk)
    except OSError:
        pass

    return res

