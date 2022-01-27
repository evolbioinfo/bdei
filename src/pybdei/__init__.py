import os
import shutil
from collections import namedtuple

import _pybdei
import numpy
import numpy as np
from ete3 import Tree


BDEI_result = namedtuple('BDEI_result', ['mu', 'la', 'psi', 'p', 'mu_CI', 'la_CI', 'psi_CI', 'p_CI',
                                         'R_naught', 'incubation_period', 'infectious_time'])
BDEI_time = namedtuple('BDEI_time', ['CPU_time', 'iterations'])


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


def parse_forest(nwk):
    with open(nwk, 'r') as f:
        return [parse_tree(_.strip('\n') + ';') for _ in f.read().split(';')[:-1]]


def save_tree(tree, nwk):
    with open(nwk, 'w+') as f:
        f.write('{}\n'.format(tree.write(format=5, format_root_node=True)))


def save_forest(forest, nwk):
    with open(nwk, 'w+') as f:
        for tree in forest:
            _ = tree.write(format=5, format_root_node=True)
            f.write('{}\n'.format(_))


def initial_guess(forest):
    internal_es, internal_is, external_es, external_is = [], [], [], []
    for tree in forest:
        for n in tree.traverse('preorder'):
            if not n.is_leaf():
                i, e = sorted(n.children, key=lambda c: c.dist)
                (internal_is if not i.is_leaf() else external_is).append(i)
                (internal_es if not e.is_leaf() else external_es).append(e)
    i_s_times = numpy.array([_.dist for _ in external_is])
    e_s_times = numpy.array([_.dist for _ in external_es])
    i_tr_times = numpy.array([_.dist for _ in internal_is])
    e_tr_times = numpy.array([_.dist for _ in internal_es])

    i_s_time = numpy.median(i_s_times)
    e_s_time = numpy.median(e_s_times)

    i_tr_time = numpy.median(i_tr_times)
    e_tr_time = numpy.median(e_tr_times)

    m_time_s = e_s_time - i_s_time
    m_time_tr = e_tr_time - i_tr_time
    m_time = numpy.mean([m_time_s, m_time_tr])
    tr_time_n = numpy.mean(list(i_tr_times) + list(e_tr_times - m_time))
    s_time_n = numpy.mean(list(i_s_times) + list(e_s_times - m_time))

    return np.array([1 / m_time, 1 / tr_time_n, 1 / s_time_n, 0.5])


def infer(nwk, start=None, upper_bounds=None, mu=-1, la=-1, psi=-1, p=-1, T=0.0, u=0, CI_repetitions=0, threads=1, **kwargs):
    """Infer BDEI parameters from a phylogenetic tree."""

    forest = parse_forest(nwk)
    # Run bdei from module library
    if start is None:
        start = initial_guess(forest)
    elif isinstance(start, BDEI_result):
        start = np.array([start.mu, start.la, start.psi, start.p])

    if upper_bounds is None:
        upper_bounds = np.array([np.inf, np.inf, np.inf, 1])
    elif isinstance(upper_bounds, BDEI_result):
        upper_bounds = np.array([upper_bounds.mu, upper_bounds.la, upper_bounds.psi, upper_bounds.p])

    something_is_fixed = False
    if mu >= 0:
        upper_bounds[0] = mu
        start[0] = mu
        something_is_fixed = True
    if la >= 0:
        upper_bounds[1] = la
        start[1] = la
        something_is_fixed = True
    if psi >= 0:
        upper_bounds[2] = psi
        start[2] = psi
        something_is_fixed = True
    if 0 < p <= 1:
        upper_bounds[3] = p
        start[3] = p
        something_is_fixed = True

    if not something_is_fixed:
        raise ValueError('To be identifiable, the BDEI model needs one of its parameters to be fixed. '
                         'Please do so, using one of the following arguments: mu, la, psi, p.')

    start = np.minimum(start, upper_bounds)
    try:
        res = _pybdei.infer(f=nwk, start=start, ub=upper_bounds, mu=mu, la=la, psi=psi, p=p, T=T, u=u,
                            nt=threads, nbiter=CI_repetitions)
    except:
        temp_nwk = nwk + '.temp'
        save_forest(forest, temp_nwk)
        res = _pybdei.infer(f=temp_nwk, start=start, ub=upper_bounds, mu=mu, la=la, psi=psi, p=p, T=T, u=u,
                            nt=threads, nbiter=CI_repetitions)
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


def get_loglikelihood(nwk, mu=None, la=None, psi=None, p=None, T=0.0, u=0, params=None, **kwargs):
    """Calculate loglikelihood for given BDEI parameters from a phylogenetic tree."""

    forest = parse_forest(nwk)
    temp_nwk = nwk + '.temp'
    save_forest(forest, temp_nwk)
    if params is not None:
        if isinstance(params, BDEI_result):
            mu = params.mu
            la = params.la
            psi = params.psi
            p = params.p
        else:
            [mu, la, psi, p] = params
    for par in [mu, la, psi, p]:
        if par is None:
            raise ValueError('All the parameters (mu, la, psi, p) must be specified, '
                             'either via dedicated arguments or via the params argument')

    res = _pybdei.likelihood(f=temp_nwk, mu=mu, la=la, psi=psi, p=p, T=T, u=u)
    shutil.rmtree(temp_nwk, ignore_errors=True)
    return res
