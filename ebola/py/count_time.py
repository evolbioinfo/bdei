from pastml.tree import read_forest

if '__main__' == __name__:
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument('--in_nwk', required=True, type=str)
    parser.add_argument('--out_log', required=True, type=str)
    params = parser.parse_args()

    T = 0
    for tree in read_forest(params.in_nwk):
        for n in tree.traverse('preorder'):
            cur_T = (getattr(n.up, 'T') if not n.is_root() else 0) + n.dist
            n.add_feature('T', cur_T)
            T = max(T, cur_T)

    with open(params.out_log, 'w+') as f:
        f.write('{}'.format(T))

