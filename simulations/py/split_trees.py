import pandas as pd
from ete3 import Tree

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Splits PhyloDeep trees and outputs their parameters.")
    parser.add_argument('--in_nwk', required=True, type=str, help="PhyloDeep nwk to split")
    parser.add_argument('--in_log', required=True, type=str, help="PhyloDeep log to process")
    parser.add_argument('--out_nwks', nargs='+', type=str, help="Split nwks")
    parser.add_argument('--out_logs', nargs='+', type=str, help="Logs")
    params = parser.parse_args()

    df = pd.read_csv(params.in_log, sep=',')
    with open(params.in_nwk, 'r') as f:
        nwks = f.read().split(';')[:-1]
    for nwk, (_, row), out_nwk, out_log in zip(nwks, df.iterrows(), params.out_nwks, params.out_logs):
        root = Tree(nwk.strip('\n') + ';', format=3)
        for j, t in enumerate(root.iter_leaves()):
            t.name = str(j)
        T = 0
        for n in root.traverse('preorder'):
            n.add_feature('T', n.dist + (getattr(n.up, 'T') if not n.is_root() else 0))
            if n.is_leaf() and getattr(n, 'T') > T:
                T = getattr(n, 'T')
        with open(out_nwk, 'w+') as f:
            f.write(root.write(format_root_node=True, format=5) + '\n')
        n_tips = len(root)
        R0, inc_time, inf_time, n_tips_declared, p = \
            row['R_naught'], row['incubation_period'], row['infectious_time'], row['tree_size'], row['sampling_proba']
        assert n_tips_declared == n_tips
        with open(out_log, 'w+') as f:
            f.write('tips\tf\tu\tT\tmu\tla\tpsi\tp\n')
            mu = 1 / inc_time
            psi = 1 / inf_time
            la = R0 * psi
            f.write('{tips}\t{f}\t{u}\t{T}\t{mu}\t{la}\t{psi}\t{p}\n'
                    .format(tips=n_tips, T=T, mu=mu, la=la, psi=psi, p=p, f=1, u=0))
        print('Tips={tips}\tT={T}\tmu={mu}\tla={la}\tpsi={psi}\tp={p}\n'
              .format(tips=n_tips, T=T, mu=mu, la=la, psi=psi, p=p))