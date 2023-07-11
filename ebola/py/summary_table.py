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
    parser.add_argument('--tab', type=str, help="estimate table")
    parser.add_argument('--time', type=str, help="time table")
    params = parser.parse_args()

    time_df = pd.DataFrame(columns=['CPU_time', 'iterations'])
    for file in params.times:
        rep = int(re.findall(r'\d+', os.path.basename(file))[0]) + 1
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
                               'psi', 'psi_min', 'psi_max',
                               'R0', 'R0_min', 'R0_max',
                               'infectious_time', 'infectious_time_min', 'infectious_time_max',
                               'incubation_period', 'incubation_period_min', 'incubation_period_max',
                               'p'])

    i2stats = {}
    for (i, nwk) in enumerate(params.forests, start=1):
        with open(nwk, 'r') as f:
            nwks = f.read().split(';')[:-1]
        i2stats[i] = (len(nwks), sum(len(Tree(_ + ';', format=5)) for _ in nwks))

    for file in params.estimates:
        basename = os.path.basename(file)
        rep = int(re.findall(r'[.]u=[\w\d]+[.](\d+)', basename)[0]) + 1
        o_trees, tips = i2stats[rep]
        u = re.findall(r'[.]u=([\w\d]+)[.]', basename)[0]

        ddf = pd.read_csv(file, sep='\t')
        estimates = ddf.loc[next(iter(ddf.index)), :]
        est_label = '{i}.u={u}.p={p}'.format(i=rep, u=u, p=estimates[['p']])

        df.loc[est_label, ['mu', 'psi', 'incubation_period', 'infectious_time', 'p']] \
            = estimates[['mu', 'psi', 'incubation_period', 'infectious_time', 'p']]
        df.loc[est_label,
               ['R0', 'mu_min', 'mu_max',
                'lambda', 'lambda_min', 'lambda_max',
                'psi_min', 'psi_max']] \
            = [estimates['R_naught'], *parse_CI(estimates['mu_CI']),
               estimates['la'], *parse_CI(estimates['la_CI']),
               *parse_CI(estimates['psi_CI'])]
        df.loc[est_label,
               ['repetition', 'sampled_tips', 'observed_trees', 'hidden_trees']] \
            = [rep, tips, o_trees, u]

    df['R0_min'] = df['lambda_min'] / df['psi']
    df['R0_max'] = df['lambda_max'] / df['psi']
    df['incubation_period_min'] = 1 / df['mu_max']
    df['incubation_period_max'] = 1 / df['mu_min']
    df['infectious_time_min'] = 1 / df['psi_max']
    df['infectious_time_max'] = 1 / df['psi_min']

    df.sort_values(by=['repetition', 'hidden_trees', 'p'], inplace=True)

    for col in ['mu', 'mu_min', 'mu_max',
                'lambda', 'lambda_min', 'lambda_max',
                'psi', 'psi_min', 'psi_max',
                'R0', 'R0_min', 'R0_max',
                'infectious_time', 'infectious_time_min', 'infectious_time_max']:
        df[col] = df[col].apply(lambda _: '{:.2f}'.format(_))
    df['p'] = df['p'].apply(lambda _: '{:.3f}'.format(_))
    for col in ['incubation_period', 'incubation_period_min', 'incubation_period_max']:
        df[col] = df[col].apply(lambda _: '{:.1f}'.format(_))
    df[['repetition', 'sampled_tips', 'observed_trees', 'hidden_trees',
        'p',
        'infectious_time', 'infectious_time_min', 'infectious_time_max',
        'incubation_period', 'incubation_period_min', 'incubation_period_max',
        'R0', 'R0_min', 'R0_max',
        'mu', 'mu_min', 'mu_max',
        'lambda', 'lambda_min', 'lambda_max',
        'psi', 'psi_min', 'psi_max'
        ]].to_csv(params.tab, sep='\t', index=False)

    df['p'] = ' $' + df['p'].apply(str) + '$'
    df['incubation_period'] = ' $' + df['incubation_period'].apply(str) + '\;[' + df['incubation_period_min'].apply(
        str) + '-' + df['incubation_period_max'].apply(str) + ']$'
    df['infectious_time'] = ' $' + df['infectious_time'].apply(str) + '\;[' + df['infectious_time_min'].apply(
        str) + '-' + df['infectious_time_max'].apply(str) + ']$  \\\\'
    df['R0'] = ' $' + df['R0'].apply(str) + '\;[' + df['R0_min'].apply(str) + '-' + df['R0_max'].apply(str) + ']$'
    df['mu'] = ' $' + df['mu'].apply(str) + '\;[' + df['mu_min'].apply(str) + '-' + df['mu_max'].apply(str) + ']$ '
    df['lambda'] = ' $' + df['lambda'].apply(str) + '\;[' + df['lambda_min'].apply(str) + '-' + df['lambda_max'].apply(
        str) + ']$'
    df['psi'] = ' $' + df['psi'].apply(str) + '\;[' + df['psi_min'].apply(str) + '-' + df['psi_max'].apply(
        str) + ']$ '
    df['repetition'] = ' $' + df['repetition'].apply(str) + '$ '
    df['observed_trees'] = ' $' + df['observed_trees'].apply(str) + '$ '
    df['hidden_trees'] = ' $' + df['hidden_trees'].apply(str) + '$ '
    df['sampled_tips'] = ' $' + df['sampled_tips'].apply(str) + '$ '

    df[['repetition', 'sampled_tips', 'observed_trees', 'hidden_trees', 'p',
        'R0', 'incubation_period', 'infectious_time']].to_csv(params.tab + '.latex', sep='&', index=False)
