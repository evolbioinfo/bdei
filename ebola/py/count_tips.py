from pastml.tree import read_forest

if '__main__' == __name__:
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument('--in_nwk', required=True, type=str)
    parser.add_argument('--out_log', required=True, type=str)
    params = parser.parse_args()

    n = sum(len(_) for _ in read_forest(params.in_nwk))
    with open(params.out_log, 'w+') as f:
        f.write('{}'.format(n))

