import logging

import numpy as np
import pandas as pd
from treesimulator import save_log, save_forest
from treesimulator.generator import simulate_tree_gillespie
from treesimulator.mtbd_models import BirthDeathExposedInfectiousModel


def random_float(min_value=0, max_value=1):
    """
    Generate a random float in ]min_value, max_value]
    :param max_value: max value
    :param min_value: min value
    :return: the generated float
    """
    return min_value + (1 - np.random.random(size=1)[0]) * (max_value - min_value)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Generates random BDPN tree.")
    parser.add_argument('--in_log', type=str, help="input log file with parameters")
    parser.add_argument('--out_log', type=str, help="output log file with parameters")
    parser.add_argument('--nwk', type=str, help="output nwk file with the forest")
    parser.add_argument('--min_tips', default=5000, type=int, help="Min number of tips")
    parser.add_argument('--max_tips', default=10000, type=int, help="Max number of tips")
    parser.add_argument('--min_T', default=0.5, type=float, help="Min proportion of the total time")
    parser.add_argument('--max_T', default=0.7, type=float, help="Max proportion of the total time")

    params = parser.parse_args()

    df = pd.read_csv(params.in_log, sep='\t')
    T, mu, la, psi, p = df.loc[next(iter(df.index)), ["T", "mu", "la", "psi", "p"]]
    min_T = T * params.min_T
    max_T = T * params.max_T
    model = BirthDeathExposedInfectiousModel(p=p, la=la, psi=psi, mu=mu)

    total_n_tips = np.inf
    while total_n_tips > params.max_tips:
        total_n_tips = 0
        u = 0
        forest = []
        while total_n_tips < params.min_tips:
            max_time = random_float(min_T, max_T) if min_T < max_T else min_T
            tree = simulate_tree_gillespie(model, max_time=max_time)
            if tree:
                total_n_tips += len(tree)
                forest.append(tree)
                for n in tree.traverse():
                    n.del_feature('T')
                tree.add_feature('T', max_time)
                print('Simulated a tree with {} sampled tips, total: {}'.format(len(tree) if tree else 0, total_n_tips))
            else:
                u += 1
        print('Simulated a forest with {} sampled tips'.format(total_n_tips))

    with open(params.nwk, 'w+') as f:
        for tree in forest:
            nwk = tree.write(format=5, format_root_node=True, features=['T'])
            f.write('{}\n'.format(nwk))
    df.loc[next(iter(df.index)), ["tips", "f", "u", "T"]] = total_n_tips, len(forest), u, '{:.2f}-{:.2f}'.format(min_T, max_T)
    df.to_csv(params.out_log, sep='\t', index=False)


