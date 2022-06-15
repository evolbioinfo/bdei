import logging
import re

import pandas as pd


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Plots errors.")
    parser.add_argument('--estimated_p', nargs='+', type=str, help="estimated parameters")
    parser.add_argument('--estimated_mu', nargs='+', type=str, help="estimated parameters")
    parser.add_argument('--estimated_la', nargs='+', type=str, help="estimated parameters")
    parser.add_argument('--estimated_psi', nargs='+', type=str, help="estimated parameters")
    parser.add_argument('--real', nargs='+', type=str, help="real parameters")
    parser.add_argument('--tab', type=str, help="estimate table")
    params = parser.parse_args()

    logging.getLogger().handlers = []
    logging.basicConfig(level=logging.INFO, format='%(asctime)s: %(message)s', datefmt="%Y-%m-%d %H:%M:%S")

    df = pd.DataFrame(columns=['type', 'T', 'sampled_tips', 'observed_trees', 'unobserved_trees',
                               'mu',
                               'lambda',
                               'psi',
                               'R_naught',
                               'infectious_time',
                               'incubation_period',
                               'p'])

    for real in params.real:
        i = int(re.findall(r'[0-9]+', real)[0])
        ddf = pd.read_csv(real, sep='\t')
        mu, la, psi, p, T, o_trees, h_trees \
            = ddf.loc[next(iter(ddf.index)), ['mu', 'la', 'psi', 'p', 'T', 'f', 'u']]
        df.loc['{}.real'.format(i),
               ['R_naught', 'infectious_time', 'incubation_period',
                'mu', 'lambda', 'psi', 'p', 'T', 'observed_trees', 'unobserved_trees', 'type']] \
            = [la / psi, 1 / psi, 1 / mu, mu, la, psi, p, T, o_trees, h_trees, 'real']

    for (est_list, fixed) in ((params.estimated_mu, 'mu'), (params.estimated_la, 'lambda'),
                         (params.estimated_psi, 'psi'), (params.estimated_p, 'p')):
        for est in est_list:
            i = int(re.findall(r'[0-9]+', est)[0])
            ddf = pd.read_csv(est, sep='\t')
            est_label = 'PyBDEI ({})'.format(fixed)
            estimates = ddf.loc[next(iter(ddf.index)), :]
            df_i = '{}.{}'.format(i, est_label)
            df.loc[df_i, ['mu', 'psi', 'p', 'R_naught', 'incubation_period', 'infectious_time']] \
                = estimates[['mu', 'psi', 'p', 'R_naught', 'incubation_period', 'infectious_time']]
            df.loc['{}.{}'.format(i, est_label), ['lambda', 'type']] = [estimates['la'], est_label]

    df.index = df.index.map(lambda _: int(_.split('.')[0]))
    df.sort_index(inplace=True)
    df[['type', 'T', 'sampled_tips', 'observed_trees', 'unobserved_trees',
        'mu',
        'lambda',
        'psi',
        'p',
        'R_naught',
        'infectious_time',
        'incubation_period']].to_csv(params.tab, sep='\t')
