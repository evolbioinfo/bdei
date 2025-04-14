import numpy as np
from scipy.integrate import odeint
from treesimulator.mtbd_models import BirthDeathExposedInfectiousModel

from pybdei import parse_forest, PYBDEI_VERSION, ERRORS, WARNINGS, INFO, DEBUG, BDEI_result, get_T, \
    SAMPLING_PERIOD_LENGTH, MEDIAN, MEAN, MIN, MAX

RTOL = 100 * np.finfo(np.float64).eps
N_U_STEPS = int(1e7)


def compute_U(T, MU, LA, PSI, RHO, SIGMA=None, nsteps=N_U_STEPS):
    """
    Calculates a function get_U which for a given time t: 0 <= t <= T, would return
    an array of unobserved probabilities [U_1(t), ..., U_m(t)].

    U_k(t) are calculated by
    (1) solving their ODEs numerically for an array tt of nsteps times equally spaced between t=T and t=0,
    producing an array of solutions sol of length nstep (for corresponding times in tt)s.
    (2) creating a linear approximation which for a given time t (2a) find an index i such that tt[i] >= t > tt[i+1];
    (2b) returns sol[i + 1] + (sol[i] - sol[i + 1]) * (tt[i] - t) / (tt[i] - tt[i + 1]).


    :param T: time at end of the sampling period
    :param MU: an array of state transition rates
    :param LA: an array of transmission rates
    :param PSI: an array of removal rates
    :param RHO: an array of sampling probabilities
    :param SIGMA: an array of rate sums: MU.sum(axis=1) + LA.sum(axis=1) + PSI
    :return: a function that for a given time t returns the array of corresponding unsampled probabilities:
        t ->  [U_1(t), ..., U_m(t)].
    """

    if SIGMA is None:
        SIGMA = MU.sum(axis=1) + LA.sum(axis=1) + PSI

    tt = np.linspace(T, 0, nsteps)
    y0 = np.ones(LA.shape[0], np.float64)
    PSI_NOT_RHO = PSI * (1 - RHO)

    def pdf_U(U, t):
        dU = (SIGMA - LA.dot(U)) * U - MU.dot(U) - PSI_NOT_RHO
        return dU

    sol = odeint(pdf_U, y0, tt, rtol=RTOL)
    sol = np.maximum(sol, 0)
    return sol


def get_u(mu=-1, la=-1, psi=-1, p=-1, pi_E=-1, T=0.0, f=0, nwk=None, log_level=INFO, params=None, u_policy=MEAN,
          **kwargs):
    """Calculates u for given BDEI parameters and a given forest."""
    if params is not None:
        if isinstance(params, BDEI_result):
            mu = params.mu
            la = params.la
            psi = params.psi
            p = params.p
        else:
            [mu, la, psi, p] = params
    for par in [mu, la, psi, p]:
        if par is None or par < 0:
            raise ValueError('All the parameters (mu, la, psi, p) must be specified, '
                             'either via dedicated arguments or via the params argument')

    if pi_E is None:
        pi_E = -1

    if not nwk:
        if f is None or f <= 0 or T is None or T <= 0:
            raise ValueError('Either the forest file (via nwk argument), '
                             'or both the number of observed trees (f) and the total time (T) must be specified.')

    if nwk:
        forest, T = parse_forest(nwk, T, u_policy=u_policy)
        f = len(forest)
        if log_level >= INFO:
            print('The input forest contains {} trees.'.format(f))
    if log_level >= INFO:
        print('T is set to {}.'.format(T))

    model = BirthDeathExposedInfectiousModel(mu=mu, la=la, psi=psi, p=p)
    Us = compute_U(T, model.transition_rates, model.transmission_rates, model.removal_rates, model.ps)[-1]
    state_frequencies = np.array([pi_E, 1 - pi_E]) if 0 <= pi_E <= 1 else model.state_frequencies
    U = state_frequencies.dot(Us)
    u = f * U / (1 - U)

    if log_level >= INFO:
        print("The number of hidden trees u is estimated {}.".format(u))

    return u


def main():
    """
    Entry point, calling :py:func:`pydbei.u_calculator.get_u` with command-line arguments.
    :return: void
    """
    import argparse

    parser = argparse.ArgumentParser(
        description="Calculator of the number of hidden trees in the forest under given BDEI model parameters.",
        prog='bdei_u')

    tree_group = parser.add_argument_group('tree-related arguments')
    tree_group.add_argument('--nwk', help="Observed tree(s) in newick format "
                                          "to calculate the total time T and the number of observed trees. "
                                          "If not given, --f and --T must be specified.",
                            type=str, required=False, default=None)
    tree_group.add_argument('--f', help="Number of observed trees. If not given, --nwk must be specified.",
                            type=int, required=False, default=0)

    parameter_group = parser.add_argument_group('parameter-related arguments')
    parameter_group.add_argument('--mu', required=True, type=float, default=None,
                                 help="Value to fix BDEI becoming-infectious rate mu. "
                                      "If not given, will be estimated.")
    parameter_group.add_argument('--la', required=True, type=float, default=None,
                                 help="Value to fix BDEI transmission rate lambda. "
                                      "If not given, will be estimated.")
    parameter_group.add_argument('--psi', required=True, type=float, default=None,
                                 help="Value to fix BDEI removal rate psi. "
                                      "If not given, will be estimated.")
    parameter_group.add_argument('-p', '--p', required=True, type=float, default=None,
                                 help="Value to fix BDEI sampling probability. "
                                      "If not given, will be estimated.")
    parameter_group.add_argument('--pi_E', required=False, type=float, default=-1,
                                 help="Frequency of E at time 0, "
                                      "should be between 0 and 1. "
                                      "If not given, will be estimated from the model parameters.")
    parameter_group.add_argument('--T', default=0, type=float,
                                 help="Total time between the hidden tree roots and the end of the sampling period. "
                                      "If the observed forest is not specified as --nwk, "
                                      "this time must be positive and will be taken as is."
                                      "If the observed forest is specified as --nwk and a positive --T value is given, "
                                      "the total time will be set to the maximum "
                                      "between this value and the maximal time between the start "
                                      "and the last sampled tip of all the trees. "
                                      "If a zero or negative --T value is given, "
                                      "the observed forest must be specified as --nwk: "
                                      "The time will be calculated based on the observed tree-specific times "
                                      "according to the measure specified in --u_policy. "
                                      "Observed tree-specific times are estimated as the time between the root "
                                      "and the last sampled tip for each tree. "
                                      "In the latter case, one can additionally annotate each tree root "
                                      "with a feature '{sp}' (e.g. '(a:2,b:3):1[&&NHX:{sp}=5];' "
                                      "is a tree with two tips, a and b, and the tree-specific time annotated to 5): "
                                      "then the tree-specific time will be set to the maximum "
                                      "between the annotated value and the time between the root "
                                      "and the last sampled tip of this tree.".format(sp=SAMPLING_PERIOD_LENGTH))
    parameter_group.add_argument('--u_policy', default=MEAN, choices=[MIN, MEDIAN, MEAN, MAX],
                                 help="How to estimate the time for unobserved trees "
                                      "in case of tree-specific observed tree times. "
                                      "By default, the mean of tree-specific observed times is taken.")
    parser.add_argument('--log_level',
                        help="level of logging information "
                             "(the lower, the less information will be printed to the output). "
                             "Possible levels are: {} (errors only), {} (errors+warnings), {} (errors+warnings+info), "
                             "{} (errors+warnings+info+debug).".format(ERRORS, WARNINGS, INFO, DEBUG), type=int,
                        default=INFO)
    parser.add_argument('-v', '--version', action='version',
                        version='%(prog)s {version}'.format(version=PYBDEI_VERSION))

    params = parser.parse_args()
    u = get_u(**vars(params))
    print(u)


if '__main__' == __name__:
    main()
