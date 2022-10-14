import numpy as np
from ete3 import Tree

from mtbd.mtbd_estimator import TI, STATE_K, optimize_likelihood_params, loglikelihood_known_states
from treesimulator.mtbd_models import BirthDeathModel

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description="Estimates BD parameters.")
    parser.add_argument('--nwk', required=False, default='/home/azhukova/projects/bdei/bd/small/trees/tree.0.nwk', type=str, help="Input tree")
    parser.add_argument('--p', required=False, default=0.305336903275534, type=float, help="Input p")
    parser.add_argument('--est', type=str, default=None, help="Estimate")
    params = parser.parse_args()

    tree = Tree(params.nwk)
    T = 0
    int_len, ext_len = [], []
    for n in tree.traverse('preorder'):
        ti = (0 if n.is_root() else getattr(n.up, TI)) + n.dist
        n.add_feature(TI, ti)
        T = max(T, ti)
        n.add_feature(STATE_K, 0)
        (ext_len if n.is_leaf() else int_len).append(n.dist)

    forest = [tree]

    median_int_len = np.median(int_len)
    model = BirthDeathModel(la=1 / median_int_len, psi=1 / np.median(ext_len), p=params.p)
    optimise = BirthDeathModel(la=1, psi=1, p=0)

    real_la, real_psi = 0.3969150932340265,	0.294334586386472
    real_model = BirthDeathModel(la=real_la, psi=real_psi, p=params.p)

    print('Real values and likelihood are:\n', "la=", real_la, "psi=", real_psi, "p=", params.p, "\t-->\t",
          loglikelihood_known_states(forest, T, real_model.transition_rates, real_model.transmission_rates,
                                     real_model.removal_rates, real_model.ps))

    # # model = real_model
    _, LA, PSI, RHO, lk = optimize_likelihood_params(forest, model, T, optimise)
    la = LA[0, 0]
    psi = PSI[0]
    p = RHO[0]
    print("la=", la, "psi=", psi, "p=", p, "lk=", lk)
    if params.est:
        with open(params.est, 'w+') as f:
            f.write('la\tpsi\tp\n')
            f.write('{la}\t{psi}\t{p}\n'.format(la=la, psi=psi, p=p))
