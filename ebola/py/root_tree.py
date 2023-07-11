from collections import Counter

import numpy as np
import pandas as pd
from ete3 import Tree
from pastml.tree import read_forest, name_tree, remove_certain_leaves

from collections import defaultdict


def collapse_zero_branches(tree, threshold):
    n_collapsed = 0
    for n in list(tree.traverse('postorder')):
        children = list(n.children)
        for child in children:
            if child.dist < threshold and not child.is_leaf():
                    n.remove_child(child)
                    for grandchild in child.children:
                        n.add_child(grandchild, dist=grandchild.dist + child.dist)
                    n_collapsed += 1
    print('Collapsed {} internal branches'.format(n_collapsed))


def resolve(tree, mp_df):
    mp_df.index = mp_df.index.map(str)
    node2country = mp_df.idxmax(axis=1)

    name_tree(tree)

    polytomy_counter = Counter()
    for n in tree.traverse('postorder'):
        n_my = len(n.children)
        if n_my > 2:
            polytomy_counter[n_my] += 1
            this_country = node2country.loc[n.name]
            other_country2children = defaultdict(list)
            for _ in n.children:
                country = node2country.loc[_.name]
                if country != this_country:
                    other_country2children[country].append(_)
            for country, children in other_country2children.items():
                if len(children) > 1:
                    parent = n.add_child(dist=0)
                    for c in children:
                        n.remove_child(c)
                        parent.add_child(c, dist=c.dist)
            while len(n.children) > 2:
                child1, child2 = np.random.choice(n.children, 2, replace=False)
                n.remove_child(child1)
                n.remove_child(child2)
                parent = n.add_child(dist=0)
                parent.add_child(child1, dist=child1.dist)
                parent.add_child(child2, dist=child2.dist)
    print('Resolved {} polytomies in a tree of {} tips: {}'
          .format(sum(polytomy_counter.values()), len(tree),
                  ', '.join('{} of {}'.format(v, k)
                            for (k, v) in sorted(polytomy_counter.items(), key=lambda _: -_[0]))))


if '__main__' == __name__:
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument('--in_nwk', default='/home/azhukova/projects/bdei_main/ebola/data/tree.nwk', type=str)
    parser.add_argument('--in_nex', default='/home/azhukova/projects/bdei_main/ebola/data/pastml/named.tree_timetree.nwk', type=str)
    parser.add_argument('--mp', default='/home/azhukova/projects/bdei_main/ebola/data/pastml/marginal_probabilities.character_country.model_F81.tab', type=str)
    parser.add_argument('--aln_len', default=18996, type=int)
    parser.add_argument('--out_nwk', default='/home/azhukova/projects/bdei_main/ebola/data/rooted.tree.10.nwk', type=str)
    params = parser.parse_args()


    tree = Tree(params.in_nwk)
    timetree = read_forest(params.in_nex)[0]

    one_mutation = 1 / params.aln_len
    tree_len = sum(_.dist for _ in tree.traverse())


    # Fix the tree naming and topology to be compatible with the timetree
    tt_tips = {_.name for _ in timetree}
    tree = remove_certain_leaves(tree, lambda _: _.name not in tt_tips)

    tree.unroot()
    tree.set_outgroup(next(_ for _ in timetree.children[1]).name)
    tree.unroot()
    og = tree.get_common_ancestor(*(_.name for _ in timetree.children[0]))
    tree.set_outgroup(og)

    tips2name = {}
    for n in timetree.traverse('postorder'):
        if n.is_leaf():
            n.add_feature('tips', {n.name})
        else:
            n.add_feature('tips', set.union(*(getattr(_, 'tips') for _ in n.children)))
            tips2name[':'.join(sorted(getattr(n, 'tips')))] = n.name
    for n in tree.traverse('postorder'):
        if n.is_leaf():
            n.add_feature('tips', {n.name})
        else:
            n.add_feature('tips', set.union(*(getattr(_, 'tips') for _ in n.children)))
            tips = ':'.join(sorted(getattr(n, 'tips')))
            if tips not in tips2name:
                if n.dist < one_mutation / 2:
                    continue
            n.name = tips2name[tips]

    collapse_zero_branches(tree, one_mutation / 2)

    resolve(tree, pd.read_csv(params.mp, index_col=0, sep='\t'))
    d_len = 0
    for n in tree.traverse():
        addition = np.random.random(1)[0] * one_mutation
        n.dist += addition
        d_len += addition
    smoothing_factor = tree_len / (tree_len + d_len)
    for n in tree.traverse():
        n.dist *= smoothing_factor

    tree.write(outfile=params.out_nwk, format=5)

