import numpy as np
from ete3 import Tree

from treesimulator.mtbd_models import BirthDeathExposedInfectiousModel

from mtbd import initial_rate_guess
from mtbd.mtbd_estimator import \
    loglikelihood_known_states, TI, STATE_K, compute_U, state_frequencies, optimize_likelihood_params
from pybdei import get_loglikelihood, infer
from pybdei.u_calculator import get_u


if __name__ == '__main__':

    # p = 0.6225767763228239
    # real_mu = 0.2523725112488919
    # real_la = 0.907081384137969
    # real_psi = 0.2692907505391973
    # T = 40
    # model = BirthDeathExposedInfectiousModel(mu=real_mu, la=real_la, psi=real_psi, p=p)
    #
    # MU, LA, PSI, RHO = model.transition_rates, model.transmission_rates, model.removal_rates, model.ps
    # PI = state_frequencies(MU, LA, PSI)
    # PSI_RHO = PSI * RHO
    # SIGMA = MU.sum(axis=1) + LA.sum(axis=1) + PSI
    # get_U = compute_U(T, MU=MU, LA=LA, PSI=PSI, RHO=RHO, SIGMA=SIGMA)
    #
    # plot_P(1, get_U, ti=10, t0=5, MU=MU, LA=LA, PSI=PSI)
    
    NWK = '/home/azhukova/projects/bdei_main/simulations/tiny/trees/tree.0.nwk'
    mu, la, psi, p = 0.7299648936648243, 1.421061063434435, 0.6069550954814479, 0.44732297990706804
    tree = Tree(NWK)
    T = 0
    for n in tree.traverse('preorder'):
        ti = (0 if n.is_root() else getattr(n.up, TI)) + n.dist
        n.add_feature(TI, ti)
        T = max(T, ti)
        n.add_feature(STATE_K, 1)

    print(T)
    forest = [tree]
    model = BirthDeathExposedInfectiousModel(mu=mu, la=la, psi=psi, p=p)
    MU, LA, PSI, RHO = model.transition_rates, model.transmission_rates, model.removal_rates, model.ps
    print(model.state_frequencies, state_frequencies(MU, LA, PSI))

    # get_U = compute_U(T, MU, LA, PSI, RHO)
    # for t in np.linspace(0.01, T - 0.01, 10):
    #     Us = get_U(t)
    #     U = model.state_frequencies.dot(Us)
    #     print(U / (1 - U), get_u(mu=mu, la=la, psi=psi, p=p, T=T-t, f=1, log_level=0))

    print(get_loglikelihood(NWK, mu, la, psi, p),
          loglikelihood_known_states(forest, T, MU, LA, PSI, RHO))
    bounds = np.array([[0.02, 5], [0.1, 5], [0.1, 1], [0.001, 0.999]])
    bdei_res, _ = infer(NWK, upper_bounds=bounds[:, 1], p=p, log_level=0)
    model = BirthDeathExposedInfectiousModel(mu=bdei_res.mu, la=bdei_res.la, psi=bdei_res.psi, p=bdei_res.p)
    MU, LA, PSI, RHO = model.transition_rates, model.transmission_rates, model.removal_rates, model.ps
    print("mu=", bdei_res.mu, "la=", bdei_res.la, "psi=", bdei_res.psi, "p=", bdei_res.p, "lk=",
          loglikelihood_known_states(forest, T, MU, LA, PSI, RHO),
          get_loglikelihood(NWK, bdei_res.mu, bdei_res.la, bdei_res.psi, bdei_res.p))

    rate = initial_rate_guess(forest)[1]
    model = BirthDeathExposedInfectiousModel(mu=rate, la=rate, psi=rate, p=p)
    optimise = BirthDeathExposedInfectiousModel(mu=1, la=1, psi=1, p=0)

    MU, LA, PSI, RHO, lk = optimize_likelihood_params(forest, model, T, optimise, bounds=bounds[:-1])
    print("mu=", MU[0, 1], "la=", LA[1, 0], "psi=", PSI[1], "p=", RHO[1], "lk=", lk,
          get_loglikelihood(NWK, MU[0, 1], LA[1, 0], PSI[1], RHO[1]))

