import numpy as np
from ete3 import Tree

from mtbd import initial_rate_guess
from mtbd.mtbd_estimator import TI, STATE_K, optimize_likelihood_params, loglikelihood_known_states
from treesimulator.mtbd_models import BirthDeathExposedInfectiousModel


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description="Estimates BD parameters.")
    parser.add_argument('--nwk', required=False, default='/home/azhukova/projects/bdei_main/simulations/tiny/trees/tree.0.nwk', type=str, help="Input tree")
    parser.add_argument('--p', required=False, default=0.44732297990706804, type=float, help="Input p")
    parser.add_argument('--est', type=str, default=None, help="Estimate")
    params = parser.parse_args()

    tree = Tree(params.nwk)
    T = 0
    for n in tree.traverse('preorder'):
        ti = (0 if n.is_root() else getattr(n.up, TI)) + n.dist
        n.add_feature(TI, ti)
        T = max(T, ti)
        n.add_feature(STATE_K, 1)

    forest = [tree]
    p = 0.44732297990706804
    real_mu = 0.7299648936648243
    real_la = 1.421061063434435
    real_psi = 0.6069550954814479

    rate = initial_rate_guess(forest)[1]
    model = BirthDeathExposedInfectiousModel(mu=rate, la=rate, psi=rate, p=p)
    optimise = BirthDeathExposedInfectiousModel(mu=1, la=1, psi=1, p=0)
    real_model = BirthDeathExposedInfectiousModel(mu=real_mu, la=real_la, psi=real_psi, p=p)

    print('Real values and likelihood are:\n', "mu=", real_mu, "la=", real_la, "psi=", real_psi, "p=", p, "\t-->\t",
          loglikelihood_known_states(forest, T, real_model.transition_rates, real_model.transmission_rates,
                                     real_model.removal_rates, real_model.ps))

    bounds = np.array([[0.02, 5], [0.1, 5], [0.1, 1]])

    # # model = real_model

    MU, LA, PSI, RHO, lk = optimize_likelihood_params(forest, model, T, optimise, bounds=bounds)
    print("mu=", MU[0, 1], "la=", LA[1, 0], "psi=", PSI[1], "p=", RHO[1], "lk=", lk)
