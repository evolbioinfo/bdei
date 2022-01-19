import logging
import os

from bd_models import BirthDeathExposedModel
from tree_generator import generate_forest, simulate_tree_gillespie, TIME_TILL_NOW
import numpy as np

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Simulates a bunch of trees for given model parameters.")
    parser.add_argument('--tips', required=True, type=int, help="desired total number of simulated leaves")
    parser.add_argument('--T', required=False, default=np.inf, type=float, help="total simulation time")
    parser.add_argument('--mu', required=True, type=float)
    parser.add_argument('--la', required=True, type=float)
    parser.add_argument('--psi', required=True, type=float)
    parser.add_argument('--p', required=True, type=float, help='sampling probability')
    parser.add_argument('--log', required=False, default=None, type=str, help="log file")
    parser.add_argument('--nwk', required=False, default=None, type=str, help="forest file")
    params = parser.parse_args()
    logging.getLogger().handlers = []
    logging.basicConfig(level=logging.INFO, format='%(asctime)s: %(message)s', datefmt="%Y-%m-%d %H:%M:%S")

    logging.info('Parameters are:\n\tmu={}\n\tlambda={}\n\tpsi={}\n\tp={}'.format(params.mu, params.la, params.psi, params.p))
    logging.info('Total time T={}'.format(params.T))

    model = BirthDeathExposedModel(p=params.p)
    model.params2rates([params.mu, params.la, params.psi])

    max_tips = 199 if params.tips < 200 else (500 if params.tips <= 500 else 10_000)
    min_tips = 50 if params.tips < 200 else (200 if params.tips <= 500 else 5_000)
    T = params.T

    if T < np.inf:
        T = T * 0.75
        while True:
            forest = generate_forest(model, max_time=T, min_tips=min_tips,
                                     keep_nones=True, root_state=model.states[1])
            total_trees = len(forest)
            forest = [tree for tree in forest if tree is not None]
            fl = len(forest)
            u = total_trees - fl
            total_tips = sum(len(list(t.iter_leaves())) for t in forest)
            logging.info('Generated a forest of {} visible and {} hidden trees with {} sampled tips over time {} '
                         '(max time is {}).'.format(fl, u, total_tips, T, params.T))
            if total_tips > max_tips:
                T *= 0.9
            else:
                break
    else:
        while True:
            max_sampled = int(min_tips + np.random.random() * (max_tips - min_tips))
            tree = simulate_tree_gillespie(model, max_time=np.inf, max_sampled=max_sampled)
            total_tips = len(tree) if tree else 0
            logging.info('Generated a tree with {} sampled tips over time {}.'
                         .format(total_tips, getattr(tree, TIME_TILL_NOW) if tree else np.inf))
            if total_tips >= min_tips:
                forest = [tree]
                fl = 1
                u = 0
                T = getattr(tree, TIME_TILL_NOW) + tree.dist
                break
    params.T = T
    params.tips = total_tips

    if params.nwk:
        os.makedirs(os.path.dirname(os.path.abspath(params.nwk)), exist_ok=True)
        with open(params.nwk, 'w+') as f:
            for tree in forest:
                nwk = tree.write(format=5, format_root_node=True)
                f.write('{}\n'.format(nwk))
    if params.log:
        os.makedirs(os.path.dirname(os.path.abspath(params.log)), exist_ok=True)
        with open(params.log, 'w+') as f:
            f.write('tips\tf\tu\tT\tmu\tla\tpsi\tp\n')
            f.write('{tips}\t{f}\t{u}\t{T}\t{mu}\t{la}\t{psi}\t{p}\n'.format(f=fl, u=u, **vars(params)))
