from collections import Counter

import numpy as np
from ete3 import Tree
from pastml.tree import read_forest


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


if '__main__' == __name__:
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument('--in_nwk', required=True, type=str)
    parser.add_argument('--in_nex', required=True, type=str)
    parser.add_argument('--aln_len', required=True, type=int)
    parser.add_argument('--out_nwk', required=True, type=str)
    params = parser.parse_args()

    tree = Tree(params.in_nwk)
    one_mutation = 1 / params.aln_len
    tree.unroot()
    collapse_zero_branches(tree, one_mutation / 2)

    timetree = read_forest(params.in_nex)[0]
    tree.set_outgroup(next(_ for _ in timetree.children[1]).name)
    tree.unroot()
    og = tree.get_common_ancestor(*(_.name for _ in timetree.children[0]))
    tree.set_outgroup(og)

    polytomy_counter = Counter()
    for n in tree.traverse('postorder'):
        n_my = len(n.children)
        if n_my > 2:
            polytomy_counter[n_my] += 1
            while len(n.children) > 2:
                child1, child2 = np.random.choice(n.children, 2, replace=False)
                n.remove_child(child1)
                n.remove_child(child2)
                # print('Old child distances: {:.2f} and {:.2f}'.format(child1.dist * params.aln_len, child2.dist * params.aln_len))
                dist = np.random.random(1)[0] * min(child2.dist, child1.dist, one_mutation)
                parent = n.add_child(dist=dist)
                child1 = parent.add_child(child1, dist=child1.dist - dist)
                child2 = parent.add_child(child2, dist=child2.dist - dist)
                # print('\t new child and parent distances: {:.2f}, {:.2f}, and {:.2f}'
                #       .format(child2.dist * params.aln_len, child1.dist * params.aln_len, dist * params.aln_len))
        elif n.dist < one_mutation / 2:
            n.dist = np.random.random(1)[0] * one_mutation

    print('Resolved {} polytomies in a tree of {} tips: {}'
          .format(sum(polytomy_counter.values()), len(tree),
                  ', '.join('{} of {}'.format(v, k)
                            for (k, v) in sorted(polytomy_counter.items(), key=lambda _: -_[0]))))
    tree.write(outfile=params.out_nwk, format=5)

