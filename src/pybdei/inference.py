from pybdei import infer


def main():
    """
    Entry point, calling :py:func:`pydbei.infer` with command-line arguments.
    :return: void
    """
    import argparse

    parser = argparse.ArgumentParser(description="BDEI model parameter inference from phylogenetic trees.", prog='bdei_infer')

    tree_group = parser.add_argument_group('tree-related arguments')
    tree_group.add_argument('--nwk', help="input tree(s) in newick format (must be rooted).",
                            type=str, required=True)
    tree_group.add_argument('-u', '--u', help="number of unobserved trees.",
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
    parameter_group.add_argument('--p', default=0, type=float,
                                  help="Value to fix BDEI sampling probability. "
                                       "If not given, will be estimated.")

    result_group = parser.add_argument_group('output-related arguments')
    result_group.add_argument('-c', '--CI_repetitions', default=0,
                              help="Number of repetitions for CI calculation (the higher-the more precise). "
                                   "If not specified, CIs will not be calculated.",
                              type=int)
    result_group.add_argument('--log', required=False, default=None, type=str,
                              help="Path to the output file where to write the estimates. "
                                   "If not given, the estimates will only be printed in the stdout")

    params = parser.parse_args()
    res = infer(**vars(params))
    print(res)
    if params.log:
        with open(params.log, 'w+') as f:
            f.write('mu\tmu_CI\tla\tla_CI\tpsi\tpsi_CI\tp\tp_CI\n'
                    .format('\tmu_CI' if params.CI_repetitions else '',
                            '\tla_CI' if params.CI_repetitions else '',
                            '\tpsi_CI' if params.CI_repetitions else '',
                            '\tp_CI' if params.CI_repetitions else ''))
            f.write('{mu}\t{mu_CI}\t{la}\t{la_CI}\t{psi}\t{psi_CI}\t{p}\t{p_CI}\n'
                    .format(**dict(zip(res._fields, res))))


if '__main__' == __name__:
    main()
