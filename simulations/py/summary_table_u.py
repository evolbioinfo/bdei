import re

import pandas as pd

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Plots errors.")
    parser.add_argument('--estimated', nargs='+', type=str, help="estimated parameters")
    parser.add_argument('--real', nargs='+', type=str, help="real parameters")
    parser.add_argument('--tab', type=str, help="estimate table")
    params = parser.parse_args()

    df = pd.DataFrame(columns=['type', 'T', 'sampled_tips', 'observed_trees', 'unobserved_trees',
                               'mu', 'mu_min', 'mu_max',
                               'lambda', 'lambda_min', 'lambda_max',
                               'psi', 'psi_min', 'psi_max',
                               'R_naught', 'R_naught_min', 'R_naught_max',
                               'infectious_time', 'infectious_time_min', 'infectious_time_max',
                               'incubation_period', 'incubation_period_min', 'incubation_period_max',
                               'p'])
    for real in params.real:
        i = int(re.findall(r'[0-9]+', real)[0])
        ddf = pd.read_csv(real, sep=',')
        R0, it, ip, p, tips, T, h_trees \
            = ddf.loc[next(iter(ddf.index)), \
            ['R0', 'infectious time', 'incubation period', 'sampling probability', 'tips', 'time', 'hidden_trees']]
        o_trees = open(real.replace('.est', '.nwk'), 'r').read().count(';')
        psi = 1 / it
        mu = 1 / ip
        la = R0 * psi
        df.loc['{}.real'.format(i),
               ['R_naught', 'infectious_time', 'incubation_period',
                'mu', 'lambda', 'psi', 'p', 'T', 'observed_trees', 'unobserved_trees', 'sampled_tips', 'type']] \
            = [R0, it, ip, mu, la, psi, p, T, o_trees, h_trees, tips, 'real']


    def parse_CI(ci):
        return [float(_) for _ in ci[1:-1].split(', ')]


    for est in params.estimated:
        i = int(re.findall(r'[0-9]+', est)[0])
        ddf = pd.read_csv(est, sep='\t')
        est_label = 'PyBDEI'
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
