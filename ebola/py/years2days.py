from pastml.tree import read_forest

if '__main__' == __name__:
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument('--in_nwk', required=True, type=str)
    parser.add_argument('--out_nwk', required=True, type=str)
    params = parser.parse_args()

    trees = read_forest(params.in_nwk)
    for tree in trees:
        for n in tree.traverse():
            n.dist *= 365
        tree.add_feature('T', float(getattr(tree, 'T')) * 365)
    nwks = [tree.write(format=5, format_root_node=True, features=['T']) for tree in trees]
    with open(params.out_nwk, 'w+') as f:
        f.write('\n'.join(nwks) + '\n')

