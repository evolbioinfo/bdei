import pandas as pd
from pastml import numeric2datetime
from pastml.tree import DATE, read_forest, annotate_dates

if '__main__' == __name__:
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument('--in_tree', default='/home/azhukova/projects/bdei/ebola/data/pastml/6/named.tree_timetree.6.nwk', type=str)
    parser.add_argument('--tab', default='/home/azhukova/projects/bdei/ebola/data/pastml/6/marginal_probabilities.character_country.model_F81.tab', type=str)
    parser.add_argument('--out_root_tree', default='/home/azhukova/projects/bdei/ebola/data/SLE/SLE.6.root.nwk', type=str)
    parser.add_argument('--date', default=2014.6, type=float, help="Aug 8 2014")
    params = parser.parse_args()

    tree = read_forest(params.in_tree)[0]
    annotate_dates([tree])
    df = pd.read_csv(params.tab, sep='\t', index_col=0)
    df.index = df.index.map(str)
    print('Root year is {}.'.format(numeric2datetime(getattr(tree, DATE))))
    print('Full tree contains {} tips.'.format(len(tree)))

    roots = []
    todo = [(tree, getattr(tree, DATE))]
    while todo:
        node, date = todo.pop()
        if date < params.date:
            todo.extend([(child, date + child.dist) for child in node.children])
        elif getattr(node, 'country') == 'SLE':
                # print('Preselected a {}-tip SLE (prob. {}) cluster'.format(len(node),  df.loc[node.name, 'SLE']))

                parent = node.up
                node.detach()
                node.dist = date - params.date
                roots.append(node)

                assert len(parent.children) == 1
                bro = parent.children[0]
                grandpa = parent.up
                grandpa.add_child(bro, dist=bro.dist + parent.dist)
                grandpa.remove_child(parent)

    # Extract SLE root tree
    todo = [(tree, getattr(tree, DATE))]
    nwk = None
    while todo:
        node, date = todo.pop()
        if getattr(node, 'country') == 'SLE' and len(node) > 100 and node.dist >= 2 / 18996 \
                and df.loc[node.name, 'SLE'] > 0.9:
            node.detach()
            node.dist = 0
            tree = node
            print(df.loc[node.name, 'SLE'], node.name, len(node))
            for node in tree.traverse('postorder'):
                if getattr(node, 'country') != 'SLE' or node.is_leaf() and not node.is_root():
                    continue
                non_SLE_children = [c for c in node.children if getattr(c, 'country') != 'SLE']
                for child in non_SLE_children:
                    # print('Removing a {}-rooted cluster of size {}'.format(getattr(child, 'country'), len(child)))
                    node.remove_child(child)
                if node.is_leaf():
                    node.add_feature('country', 'non-SLE')
                elif len(node.children) == 1:
                    child = node.children[0]
                    if node.is_root():
                        child.dist += node.dist
                        child.detach()
                        nwk = child
                    else:
                        parent = node.up
                        parent.add_child(child, dist=node.dist + child.dist)
                        parent.remove_child(node)
                elif node.is_root():
                    nwk = node
            break
        else:
            todo.extend([(child, date + child.dist) for child in node.children])

    with open(params.out_root_tree, 'w+') as f:
        f.write(nwk.write(format=3, format_root_node=True))
    print('SLE root tree contains {} tips'.format(len(nwk)))
