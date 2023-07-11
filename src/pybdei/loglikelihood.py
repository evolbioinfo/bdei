from pybdei import get_loglikelihood, ERRORS, WARNINGS, INFO, DEBUG, PYBDEI_VERSION, SAMPLING_PERIOD_LENGTH, MAX, MIN, \
    MEDIAN, MEAN


def main():
    """
    Entry point, calling :py:func:`pydbei.infer` with command-line arguments.
    :return: void
    """
    import argparse

    parser = argparse.ArgumentParser(description="BDEI model likelihood calculation on phylogenetic trees.", prog='bdei_loglikelihood')

    tree_group = parser.add_argument_group('tree-related arguments')
    tree_group.add_argument('--nwk', help="input tree(s) in newick format (must be rooted).",
                            type=str, required=True)
    tree_group.add_argument('-u', '--u', help="Number of unobserved trees. "
                                              "Specify -1 if you want u to be estimated.",
                            type=int, default=-1)

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
                                  help="Total time between the tree roots and the end of the sampling period. "
                                       "If a positive value is given, the total time will be set to the maximum "
                                       "between this value and the maximal time between the start "
                                       "and the last sampled tip of all the trees. "
                                       "If a zero or negative value is given, the time will be tree-specific "
                                       "and estimated as the time between the root "
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

    parser.add_argument('-t', '--threads', help="number of threads for parallelization.", type=int, default=1)
    parser.add_argument('--log_level',
                        help="level of logging information "
                             "(the lower, the less information will be printed to the output). "
                             "Possible levels are: {} (errors only), {} (errors+warnings), {} (errors+warnings+info), "
                             "{} (errors+warnings+info+debug).".format(ERRORS, WARNINGS, INFO, DEBUG), type=int,
                        default=INFO)
    parser.add_argument('-v', '--version', action='version', version='%(prog)s {version}'.format(version=PYBDEI_VERSION))

    params = parser.parse_args()
    if params.pi_E is None:
        params.pi_E = -1
    res = get_loglikelihood(**vars(params))
    print(res)


if '__main__' == __name__:
    main()
