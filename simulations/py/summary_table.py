import logging
import re

import pandas as pd


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Plots errors.")
    parser.add_argument('--estimated', nargs='+', type=str, help="estimated parameters")
    parser.add_argument('--estimated_beast2', type=str, help="estimated parameters by BEAST2", required=False,
                        default=None)
    parser.add_argument('--estimated_beast2_CI', type=str, help="estimated CIs by BEAST2", required=False,
                        default=None)
    parser.add_argument('--estimated_dl', type=str, help="estimated parameters by CNN_CBLV", required=False,
                        default=None)
    parser.add_argument('--estimated_dl_CI', type=str, help="estimated CIs by CNN_CBLV", required=False,
                        default=None)
    parser.add_argument('--real', nargs='+', type=str, help="real parameters")
    parser.add_argument('--tab', type=str, help="estimate table")
    params = parser.parse_args()

    logging.getLogger().handlers = []
    logging.basicConfig(level=logging.INFO, format='%(asctime)s: %(message)s', datefmt="%Y-%m-%d %H:%M:%S")

    df = pd.DataFrame(columns=['type', 'T', 'sampled_tips', 'observed_trees', 'unobserved_trees',
                               'mu', 'mu_min', 'mu_max',
                               'lambda', 'lambda_min', 'lambda_max',
                               'psi', 'psi_min', 'psi_max',
                               'R_naught', 'R_naught_min', 'R_naught_max',
                               'infectious_time', 'infectious_time_min', 'infectious_time_max',
                               'incubation_period', 'incubation_period_min', 'incubation_period_max',
                               'p'])

    if params.estimated_beast2:
        bdf = pd.read_csv(params.estimated_beast2, header=0)
        bdf.index = bdf.index.map(int)
        bdf.columns = ['R_naught', 'infectious_time', 'incubation_period']
        bdf.loc[pd.isna(bdf['R_naught']), 'R_naught'] = 3
        bdf.loc[pd.isna(bdf['infectious_time']), 'infectious_time'] = 5.5
        bdf.loc[pd.isna(bdf['incubation_period']), 'incubation_period'] = 25.1
        bdf_ci = pd.read_csv(params.estimated_beast2_CI, index_col=0, header=0)
        bdf_ci.index = bdf_ci.index.map(lambda _: int(_) - 1)
        bdf_ci = bdf_ci[['R_naught_CI_2_5', 'R_naught_CI_97_5',
                         'infectious_time_CI_2_5', 'infectious_time_CI_97_5',
                         'incubation_period_CI_2_5', 'incubation_period_CI_97_5']]
        bdf_ci.columns = ['R_naught_min', 'R_naught_max',
                          'infectious_time_min', 'infectious_time_max',
                          'incubation_period_min', 'incubation_period_max']
        bdf = bdf.join(bdf_ci, how='left')
        bdf['type'] = 'BEAST2'
        bdf.index = bdf.index.map(lambda _: '{}.{}'.format(_, 'BEAST2'))
        bdf['observed_trees'] = 1
        bdf['unobserved_trees'] = 0

        bdf['psi'] = 1 / bdf['infectious_time']
        bdf['lambda'] = bdf['R_naught'] * bdf['psi']
        bdf['mu'] = 1 / bdf['incubation_period']

        df = df.append(bdf)

    if params.estimated_dl:
        dldf = pd.read_csv(params.estimated_dl, header=0)
        dldf.index = dldf.index.map(int)
        dldf.columns = ['R_naught', 'infectious_time', 'incubation_period']
        if params.estimated_dl_CI:
            dldf_ci = pd.read_csv(params.estimated_dl_CI, index_col=0, header=0)
            dldf_ci.index = dldf_ci.index.map(int)

            dldf_ci = dldf_ci[['R_naught_HPD_2_5_1000', 'R_naught_HPD_97_5_1000',
                               'infectious_time_resc_HPD_2_5_1000', 'infectious_time_resc_HPD_97_5_1000',
                               'incubation_time_resc_HPD_2_5_1000', 'incubation_time_resc_HPD_97_5_1000']]
            dldf_ci.columns = ['R_naught_min', 'R_naught_max',
                               'infectious_time_min', 'infectious_time_max',
                               'incubation_period_min', 'incubation_period_max']
            dldf = dldf.join(dldf_ci, how='left')
        dldf['type'] = 'PhyloDeep'
        dldf.index = dldf.index.map(lambda _: '{}.{}'.format(_, 'PhyloDeep'))
        dldf['observed_trees'] = 1
        dldf['unobserved_trees'] = 0

        dldf['psi'] = 1 / dldf['infectious_time']
        dldf['lambda'] = dldf['R_naught'] * dldf['psi']
        dldf['mu'] = 1 / dldf['incubation_period']

        df = df.append(dldf)

    for real in params.real:
        i = int(re.findall(r'[0-9]+', real)[0])
        ddf = pd.read_csv(real, sep='\t')
        mu, la, psi, p, T, o_trees, h_trees \
            = ddf.loc[next(iter(ddf.index)), ['mu', 'la', 'psi', 'p', 'T', 'f', 'u']]
        df.loc['{}.real'.format(i),
               ['R_naught', 'infectious_time', 'incubation_period',
                'mu', 'lambda', 'psi', 'p', 'T', 'observed_trees', 'unobserved_trees', 'type']] \
            = [la / psi, 1 / psi, 1 / mu, mu, la, psi, p, T, o_trees, h_trees, 'real']


    def parse_CI(ci):
        return [float(_) for _ in ci[1:-1].split(', ')]


    for est in params.estimated:
        i = int(re.findall(r'[0-9]+', est)[0])
        ddf = pd.read_csv(est, sep='\t')
        est_label = 'PyBDEI{}'.format('' if 'tree' in est else ' (forest)')
        estimates = ddf.loc[next(iter(ddf.index)), :]
        df_i = '{}.{}'.format(i, est_label)
        df.loc[df_i, ['mu', 'psi', 'R_naught', 'incubation_period', 'infectious_time']] \
            = estimates[['mu', 'psi', 'R_naught', 'incubation_period', 'infectious_time']]
        df.loc['{}.{}'.format(i, est_label),
               ['mu_min', 'mu_max',
                'lambda', 'lambda_min', 'lambda_max',
                'psi_min', 'psi_max',
                'type']] \
            = [*parse_CI(estimates['mu_CI']),
               estimates['la'], *parse_CI(estimates['la_CI']),
               *parse_CI(estimates['psi_CI']), est_label]

    df.index = df.index.map(lambda _: int(_.split('.')[0]))
    df.sort_index(inplace=True)
    df[['type', 'T', 'sampled_tips', 'observed_trees', 'unobserved_trees',
        'mu', 'mu_min', 'mu_max',
        'lambda', 'lambda_min', 'lambda_max',
        'psi', 'psi_min', 'psi_max',
        'p',
        'R_naught', 'R_naught_min', 'R_naught_max',
        'infectious_time', 'infectious_time_min', 'infectious_time_max',
        'incubation_period', 'incubation_period_min', 'incubation_period_max']].to_csv(params.tab, sep='\t')
