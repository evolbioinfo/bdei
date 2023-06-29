import pandas as pd
from ete3 import Tree
from pastml import numeric2datetime
from pastml.tree import DATE, read_forest, annotate_dates

if '__main__' == __name__:
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument('--in_tree', required=True, type=str)
    parser.add_argument('--tab', required=True, type=str)
    parser.add_argument('--out_forest', required=True, type=str)
    parser.add_argument('--date', required=True, type=float)
    params = parser.parse_args()

    tree = read_forest(params.in_tree)[0]
    annotate_dates([tree])
    df = pd.read_csv(params.tab, sep='\t', index_col=0)
    df.index = df.index.map(str)
    max_date = max(getattr(_, DATE) for _ in tree)
    print('Root date is {}, last sampled time date is {}.'
          .format(numeric2datetime(getattr(tree, DATE)), numeric2datetime(max_date)))
    print('Full tree contains {} tips.'.format(len(tree)))

    roots = []
    todo = [(tree, getattr(tree, DATE))]
    while todo:
        node, date = todo.pop()
        if params.date > date:
            todo.extend([(child, date + child.dist) for child in node.children])
            continue
        if getattr(node, 'country') == 'SLE':
                # print('Preselected a {}-tip SLE (prob. {}) cluster'.format(len(node),  df.loc[node.name, 'SLE']))
                node.detach()
                node.dist = date - params.date
                roots.append(node)
                # annotate the total time till the end of sampling
                node.add_feature('T', max_date - date + node.dist)

    nwks = []
    for tree in roots:
        for node in tree.traverse('postorder'):
            if getattr(node, 'country') != 'SLE' or (node.is_leaf() and not node.is_root()):
                continue
            if node.is_leaf():
                if node.is_root():
                    nwks.append(node)
                continue
            non_SLE_children = [c for c in node.children if getattr(c, 'country') != 'SLE']
            for child in non_SLE_children:
                # print('Removing a {}-rooted cluster of size {}'.format(getattr(child, 'country'), len(child)))
                node.remove_child(child)
            # if our node became a leaf now,
            # it means it corresponded to an internal node with two non-SLE children in the initial tree
            if node.is_leaf():
                node.add_feature('country', 'non-SLE')
            elif len(node.children) == 1:
                child = node.children[0]
                if node.is_root():
                    child.dist += node.dist
                    child.detach()
                    nwks.append(child)
                    child.add_feature('T', getattr(node, 'T'))
                else:
                    parent = node.up
                    parent.add_child(child, dist=node.dist + child.dist)
                    parent.remove_child(node)
            elif node.is_root():
                nwks.append(node)

    with open(params.out_forest, 'w+') as f:
        f.write('\n'.join(_.write(format=3, format_root_node=True, features=['T']) for _ in nwks) + '\n')

    print('SLE forest contains {} tips'.format(sum(len(tree) for tree in nwks)))
