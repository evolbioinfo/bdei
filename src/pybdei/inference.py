from pybdei import infer, ERRORS, DEBUG, WARNINGS, INFO, PYBDEI_VERSION


def main():
    """
    Entry point, calling :py:func:`pydbei.infer` with command-line arguments.
    :return: void
    """
    import argparse

    parser = argparse.ArgumentParser(description="BDEI model parameter inference from phylogenetic trees.", prog='bdei_infer')

    tree_group = parser.add_argument_group('tree-related arguments')
    tree_group.add_argument('--nwk', help="Input tree(s) in newick format (must be rooted).",
                            type=str, required=True)
    tree_group.add_argument('-u', '--u', help="Number of unobserved trees. "
                                              "Can be non-zero only when all trees started at the same time "
                                              "(i.e. T is given).",
                            type=int, default=0)

    parameter_group = parser.add_argument_group('parameter-related arguments')
    parameter_group.add_argument('--mu', default=-1, type=float,
                                  help="Value to fix BDEI becoming-infectious rate mu. "
                                       "If not given, will be estimated.")
    parameter_group.add_argument('--la', default=-1, type=float,
                                  help="Value to fix BDEI transmission rate lambda. "
                                       "If not given, will be estimated.")
    parameter_group.add_argument('--psi', default=-1, type=float,
                                  help="Value to fix BDEI removal rate psi. "
                                       "If not given, will be estimated.")
    parameter_group.add_argument('-p', '--p', default=-1, type=float,
                                  help="Value to fix BDEI sampling probability. "
                                       "If not given, will be estimated.")
    parameter_group.add_argument('--start', default=None, nargs=4, type=float,
                                  help="Starting values for parameter optimisation, "
                                       "should be 4 values in the following order: mu, lambda, psi, p. "
                                       "If not given, will be estimated.")
    parameter_group.add_argument('--upper_bounds', default=None, nargs=4, type=float,
                                  help="Upper bound on parameter values for parameter optimisation, "
                                       "should be in the following order: mu, lambda, psi, p. "
                                       "If not given, will be estimated.")
    parameter_group.add_argument('--pi_E', default=-1, type=float,
                                  help="Frequency of E at time 0, "
                                       "should be between 0 and 1. "
                                       "If not given, will be estimated from the model parameters.")
    parameter_group.add_argument('--T', default=0, type=float,
                                  help="Total time between the tree roots and the end of the epidemic "
                                       "(to be given if all trees start at the same time). "
                                       "If a positive value is given, the total time will be set to the maximum "
                                       "between this value and the maximal time between the start "
                                       "and the last sampled tip of all the trees. "
                                       "If a zero or negative value is given, the time will be tree-specific "
                                       "and estimated as the time between the root "
                                       "and the last sampled tip for each tree."
                                       "Note that if the number of unobserved trees (u) is given, "
                                       "all the trees are assumed to have started at the same time, "
                                       "hence T must be non-negative.")

    result_group = parser.add_argument_group('output-related arguments')
    result_group.add_argument('-c', '--CI_repetitions', default=0,
                              help="Number of repetitions for CI calculation "
                                   "(the higher, the more precise but also longer; a typical value is 100). "
                                   "If not specified, CIs will not be calculated.",
                              type=int)
    result_group.add_argument('--log', default=None, type=str,
                              help="Path to the output file where to write the estimates. "
                                   "If not given, the estimates will only be printed in the stdout")
    result_group.add_argument('--time_log', default=None, type=str,
                              help="Path to the output file where to write the time. "
                                   "If not given, the time will only be printed in the stdout")

    parser.add_argument('-t', '--threads', help="number of threads for parallelization.", type=int, default=1)
    parser.add_argument('--log_level',
                        help="level of logging information "
                             "(the lower, the less information will be printed to the output). "
                             "Possible levels are: {} (errors only), {} (errors+warnings), {} (errors+warnings+info), "
                             "{} (errors+warnings+info+debug).".format(ERRORS, WARNINGS, INFO, DEBUG), type=int,
                        default=INFO)
    parser.add_argument('-v', '--version', action='version', version='%(prog)s {version}'.format(version=PYBDEI_VERSION))

    params = parser.parse_args()
    res, time = infer(**vars(params))
    print(res)
    print(time)

    if params.log:
        with open(params.log, 'w+') as f:
            f.write('mu\tmu_CI\tla\tla_CI\tpsi\tpsi_CI\tp\tp_CI\tR_naught\tincubation_period\tinfectious_time\n'
                    .format('\tmu_CI' if params.CI_repetitions else '',
                            '\tla_CI' if params.CI_repetitions else '',
                            '\tpsi_CI' if params.CI_repetitions else '',
                            '\tp_CI' if params.CI_repetitions else ''))
            f.write('{mu}\t{mu_CI}\t{la}\t{la_CI}\t{psi}\t{psi_CI}\t{p}\t{p_CI}\t{R_naught}\t{incubation_period}\t{infectious_time}\n'
                    .format(**dict(zip(res._fields, res))))
    if params.time_log:
        with open(params.time_log, 'w+') as f:
            f.write('{}\n'.format('\t'.join(time._fields)))
            f.write('{}\n'.format('\t'.join('{}'.format(_) for _ in time)))


if '__main__' == __name__:
    main()
