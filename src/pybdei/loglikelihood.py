from pybdei import get_loglikelihood


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
    parameter_group.add_argument('--mu', required=True, type=float, default=None,
                                  help="Value to fix BDEI becoming-infectious rate mu. "
                                       "If not given, will be estimated.")
    parameter_group.add_argument('--la', required=True, type=float, default=None,
                                  help="Value to fix BDEI transmission rate lambda. "
                                       "If not given, will be estimated.")
    parameter_group.add_argument('--psi', required=True, type=float, default=None,
                                  help="Value to fix BDEI removal rate psi. "
                                       "If not given, will be estimated.")
    parameter_group.add_argument('--p', required=True, type=float, default=None,
                                  help="Value to fix BDEI sampling probability. "
                                       "If not given, will be estimated.")
    params = parser.parse_args()
    res = get_loglikelihood(**vars(params))
    print(res)


if '__main__' == __name__:
    main()
