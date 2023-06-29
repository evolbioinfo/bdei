import pandas as pd
from ete3 import Tree
import numpy as np


def random_float(min_value=0, max_value=1):
    """
    Generate a random float in ]min_value, max_value]
    :param max_value: max value
    :param min_value: min value
    :return: the generated float
    """
    return min_value + (1 - np.random.random(size=1)[0]) * (max_value - min_value)


def get_unbanned_descendants(n, T):
    time = getattr(n, 'time')
    if time >= T:
        if not getattr(n, 'banned', False):
            yield n
    else:
        for c in n.children:
            for _ in get_unbanned_descendants(c, T):
                yield _


def get_subtree_at_T(n, T):
    nodes = list(get_unbanned_descendants(n, T))
    if nodes:
        n = np.random.choice(nodes, size=1, replace=False)[0]
        parent = n
        while parent:
            parent.add_feature('banned', True)
            parent = parent.up
        n.detach()
        n.dist = getattr(n, 'time') - T
        return n
    return None


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Splits input trees into sub-epidemics.")
    parser.add_argument('--in_nwk', required=True, type=str, help="Nwk tree to cut")
    parser.add_argument('--in_log', required=True, type=str, help="Nwk tree parameters")
    parser.add_argument('--out_nwk', type=str, help="Nwk containing a forest of cut subtrees")
    parser.add_argument('--out_log', type=str, help="Output parameters")
    parser.add_argument('--min_T', default=0.8, type=float, help="Min proportion of the total time")
    parser.add_argument('--max_T', default=0.7, type=float, help="Max proportion of the total time")
    parser.add_argument('--min_tips', default=5000, type=int, help="Min number of tips")
    params = parser.parse_args()

    df = pd.read_csv(params.in_log, sep='\t')
    tree = Tree(params.in_nwk)
    for n in tree.traverse('preorder'):
        parent_time = 0 if n.is_root() else getattr(n.up, 'time')
        n.add_feature('time', parent_time + n.dist)
    height = max(getattr(t, 'time') for t in tree)
    min_T = height * params.min_T
    max_T = height * params.max_T

    total_tips = np.round(random_float(params.min_tips, len(tree)), 0)

    df.loc[0, 'T'] = '{:.2f}-{:.2f}'.format(min_T, max_T)
    f = 0
    tips = 0
    u = 0
    with open(params.out_nwk, 'w+') as file:
        while tips < total_tips:
            T = random_float(min_T, max_T)
            n = get_subtree_at_T(tree, T)
            if n:
                n.add_feature('T', T)
                file.write(n.write(format_root_node=True, format=5, features=['T']) + '\n')
                f += 1
                tips += len(n)
            else:
                u += 1
    df.loc[0, 'f'] = f
    df.loc[0, 'tips'] = tips
    df.loc[0, 'u'] = u
    df.to_csv(params.out_log, sep='\t', index=False)
