from ete3 import Tree

from mtbd import initial_rate_guess
from mtbd.mtbd_estimator import TI, STATE_K, optimize_likelihood_params, loglikelihood_known_states
from treesimulator.mtbd_models import BirthDeathExposedInfectiousModel

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description="Estimates BD parameters.")
    parser.add_argument('--nwk', required=False, default='/home/azhukova/projects/bdei/simulations/medium/trees/tree.1.nwk', type=str, help="Input tree")
    parser.add_argument('--p', required=False, default=0.305336903275534, type=float, help="Input p")
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
    p = 0.6225767763228239
    real_mu = 0.2523725112488919
    real_la = 0.907081384137969
    real_psi = 0.2692907505391973

    rate = initial_rate_guess(forest).pop()
    model = BirthDeathExposedInfectiousModel(mu=rate, la=rate, psi=rate, p=p)
    optimise = BirthDeathExposedInfectiousModel(mu=1, la=1, psi=1, p=0)
    real_model = BirthDeathExposedInfectiousModel(mu=real_mu, la=real_la, psi=real_psi, p=p)

    print('Real values and likelihood are:\n', "mu=", real_mu, "la=", real_la, "psi=", real_psi, "p=", p, "\t-->\t",
          loglikelihood_known_states(forest, T, real_model.transition_rates, real_model.transmission_rates,
                                     real_model.removal_rates, real_model.ps))
    # # model = real_model
    MU, LA, PSI, RHO, lk = optimize_likelihood_params(forest, model, T, optimise)
    print("mu=", MU[0, 1], "la=", LA[1, 0], "psi=", PSI[1], "p=", RHO[1], "lk=", lk)
