import os
import re

import pandas as pd
from ete3 import Tree


def parse_CI(ci):
    return [float(_) for _ in ci[1:-1].split(', ')]


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Combines the estimates from different forests/settings into one table.")
    parser.add_argument('--times', nargs='+', type=str, help="times")
    parser.add_argument('--forests', nargs='+', type=str, help="forests")
    parser.add_argument('--estimates', nargs='+', type=str, help="estimated parameters")
    parser.add_argument('--labels', nargs='+', type=str, help="estimated labels")
    parser.add_argument('--tab', type=str, help="estimate table")
    parser.add_argument('--time', type=str, help="time table")
    params = parser.parse_args()

    time_df = pd.DataFrame(columns=['CPU_time', 'iterations'])
    for file in params.times:
        rep = int(os.path.basename(file).split('_')[0].split('.')[1]) + 1
        ddf = pd.read_csv(file, sep='\t')
        time_df.loc[rep, :] = ddf.loc[next(iter(ddf.index)), :]
    me, m, M = time_df.mean(), time_df.min(), time_df.max()
    time_df.loc['mean', :] = me
    time_df.loc['min', :] = m
    time_df.loc['max', :] = M
    time_df.to_csv(params.time, sep='\t')

    df = pd.DataFrame(columns=['repetition', 'sampled_tips', 'observed_trees', 'hidden_trees',
                               'mu', 'mu_min', 'mu_max',
                               'lambda', 'lambda_min', 'lambda_max',
                               'psi',
                               'R0', 'R0_min', 'R0_max',
                               'infectious_time',
                               'incubation_period', 'incubation_period_min', 'incubation_period_max',
                               'p', 'p_min', 'p_max'])

    i2stats = {}
    for (i, nwk) in enumerate(params.forests, start=1):
        with open(nwk, 'r') as f:
            nwks = f.read().split(';')[:-1]
        i2stats[i] = (len(nwks), sum(len(Tree(_ + ';', format=5)) for _ in nwks))

    for (file, est_label) in zip(params.estimates, params.labels):
        rep = int(est_label.split('_')[0]) + 1
        N = int(re.findall(r'N=(\d+)', est_label)[0])
        o_trees, tips = i2stats[rep]

        ddf = pd.read_csv(file, sep='\t')
        estimates = ddf.loc[next(iter(ddf.index)), :]

        df.loc[est_label, ['mu', 'psi', 'incubation_period', 'infectious_time', 'p']] \
            = estimates[['mu', 'psi', 'incubation_period', 'infectious_time', 'p']]
        df.loc[est_label,
               ['R0', 'mu_min', 'mu_max',
                'lambda', 'lambda_min', 'lambda_max',
                'p_min', 'p_max']] \
            = [estimates['R_naught'], *parse_CI(estimates['mu_CI']),
               estimates['la'], *parse_CI(estimates['la_CI']),
               *parse_CI(estimates['p_CI'])]
        df.loc[est_label,
               ['repetition', 'sampled_tips', 'observed_trees', 'hidden_trees']] \
            = [rep, tips, o_trees, N - o_trees]

    df['R0_min'] = df['lambda_min'] / df['psi']
    df['R0_max'] = df['lambda_max'] / df['psi']
    df['incubation_period_min'] = 1 / df['mu_max']
    df['incubation_period_max'] = 1 / df['mu_min']

    df.sort_values(by=['repetition', 'hidden_trees', 'infectious_time'], inplace=True)

    for col in ['mu', 'mu_min', 'mu_max',
                'lambda', 'lambda_min', 'lambda_max',
                'psi',
                'R0', 'R0_min', 'R0_max',
                'p', 'p_min', 'p_max']:
        df[col] = df[col].apply(lambda _: '{:.2f}'.format(_))
    for col in ['infectious_time',
                'incubation_period', 'incubation_period_min', 'incubation_period_max']:
        df[col] = df[col].apply(lambda _: '{:.1f}'.format(_))
    df[['repetition', 'sampled_tips', 'observed_trees', 'hidden_trees',
        'infectious_time',
        'incubation_period', 'incubation_period_min', 'incubation_period_max',
        'R0', 'R0_min', 'R0_max',
        'p', 'p_min', 'p_max',
        'mu', 'mu_min', 'mu_max',
        'lambda', 'lambda_min', 'lambda_max',
        'psi'
        ]].to_csv(params.tab, sep='\t', index=False)

    df['incubation_period'] = ' $' + df['incubation_period'].apply(str) + '\;[' + df['incubation_period_min'].apply(
        str) + '-' + df['incubation_period_max'].apply(str) + ']$ '
    df['infectious_time'] = ' $' + df['infectious_time'].apply(str) + '$ '
    df['R0'] = ' $' + df['R0'].apply(str) + '\;[' + df['R0_min'].apply(str) + '-' + df['R0_max'].apply(str) + ']$ '
    df['p'] = ' $' + df['p'].apply(str) + '\;[' + df['p_min'].apply(str) + '-' + df['p_max'].apply(str) + ']$ \\\\'
    df['mu'] = ' $' + df['mu'].apply(str) + '\;[' + df['mu_min'].apply(str) + '-' + df['mu_max'].apply(str) + ']$ '
    df['lambda'] = ' $' + df['lambda'].apply(str) + '\;[' + df['lambda_min'].apply(str) + '-' + df['lambda_max'].apply(
        str) + ']$'
    df['psi'] = ' $' + df['psi'].apply(str) + '$'
    df['repetition'] = ' $' + df['repetition'].apply(str) + '$ '
    df['observed_trees'] = ' $' + df['observed_trees'].apply(str) + '$ '
    df['hidden_trees'] = ' $' + df['hidden_trees'].apply(str) + '$ '
    df['sampled_tips'] = ' $' + df['sampled_tips'].apply(str) + '$ '

    df[['repetition', 'sampled_tips', 'observed_trees', 'hidden_trees',
        'infectious_time', 'R0', 'incubation_period', 'p',]].to_csv(params.tab + '.latex', sep='&', index=False)
