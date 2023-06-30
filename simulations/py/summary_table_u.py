import re

import pandas as pd

from pybdei.u_calculator import get_u
from pybdei import ERRORS

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Plots errors.")
    parser.add_argument('--estimated', nargs='+', type=str, help="estimated parameters")
    parser.add_argument('--real', nargs='+', type=str, help="real parameters")
    parser.add_argument('--tab', type=str, help="estimate table")
    params = parser.parse_args()

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

    for est in params.estimated:
        i = int(re.findall(r'[0-9]+', est)[0])
        policy = re.findall(r'mean|median|min|max|zero', est)[0]
        ddf = pd.read_csv(est, sep='\t')
        estimates = ddf.loc[next(iter(ddf.index)), :]
        df_i = '{}.{}'.format(i, policy)
        df.loc[df_i, ['mu', 'psi', 'R_naught', 'incubation_period', 'infectious_time']] \
            = estimates[['mu', 'psi', 'R_naught', 'incubation_period', 'infectious_time']]
        mu, la, psi, p = estimates[['mu', 'la', 'psi', 'p']]
        u = get_u(mu=mu, la=la, psi=psi, p=p, nwk=est.replace('{}.est'.format(policy), 'nwk'),
                  u_policy=policy, log_level=ERRORS, T=0.001 if 'forest' in est else 0)
        df.loc['{}.{}'.format(i, policy), ['lambda', 'type', 'unobserved_trees']] \
            = [estimates['la'], policy, u]

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
