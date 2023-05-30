import logging

import numpy as np
import pandas as pd

RATE_PARAMETERS = ['mu', 'lambda', 'psi']
EPIDEMIOLOGIC_PARAMETERS = ['R_naught', 'infectious_time', 'incubation_period']
PARAMETERS = RATE_PARAMETERS + ['p'] + EPIDEMIOLOGIC_PARAMETERS
par2greek = {'mu': u'\u03bc', 'lambda': u'\u03bb', 'psi': u'\u03c8', 'p': '\u03c1',
             'R_naught': u'\u0052\u2080' + '=' + u'\u03bb\u002F\u03c8',
             'infectious_time': 'infectious time 1' + u'\u002F\u03c8', 'incubation_period': 'incubation period 1' + u'\u002F\u03bc'}


def asses_CIs(types):
    for type in types:
        mask = df['type'] == type
        if 'PyBDEI' in type:
            print('\n================{}==============='.format(type))
            n_observations = sum(mask)
            for par in RATE_PARAMETERS:
                df.loc[mask, '{}_within_CI'.format(par)] \
                    = np.less_equal(df.loc[mask, '{}_min'.format(par)], real_df[par]) \
                      & np.less_equal(real_df[par], df.loc[mask, '{}_max'.format(par)])
                print('{}:\t{:.1f}% within CIs'
                      .format(par, 100 * sum(df.loc[mask, '{}_within_CI'.format(par)]) / n_observations))
                df.loc[mask, '{}_CI_relative_width'.format(par)] \
                    = 100 * (df.loc[mask, '{}_max'.format(par)] - df.loc[mask, '{}_min'.format(par)]) / real_df[par]
                print('{}:\t{:.1f}% median CI width'
                      .format(par, (df.loc[mask, '{}_CI_relative_width'.format(par)].median())))
        elif type != 'real':
            print('\n================{}==============='.format(type))
            n_observations = sum(mask)
            for par in EPIDEMIOLOGIC_PARAMETERS:
                df.loc[mask, '{}_within_CI'.format(par)] \
                    = np.less_equal(df.loc[mask, '{}_min'.format(par)], real_df[par]) \
                      & np.less_equal(real_df[par], df.loc[mask, '{}_max'.format(par)])
                print('{}:\t{:.1f}% within CIs'
                      .format(par, 100 * sum(df.loc[mask, '{}_within_CI'.format(par)]) / n_observations))
                df.loc[mask, '{}_CI_relative_width'.format(par)] \
                    = 100 * (df.loc[mask, '{}_max'.format(par)] - df.loc[mask, '{}_min'.format(par)]) / real_df[par]
                print('{}:\t{:.1f}% median CI width'
                      .format(par, (df.loc[mask, '{}_CI_relative_width'.format(par)].median())))


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Plots errors.")
    parser.add_argument('--estimates', default='../large/estimates.tab', type=str, help="estimated parameters")
    params = parser.parse_args()
    logging.getLogger().handlers = []
    logging.basicConfig(level=logging.INFO, format='%(asctime)s: %(message)s', datefmt="%Y-%m-%d %H:%M:%S")

    df = pd.read_csv(params.estimates, sep='\t', index_col=0)

    real_df = df.loc[df['type'] == 'real', :]
    df = df.loc[df['type'] != 'real', :]
    types = sorted(df['type'].unique(), key=lambda _: (_[:2], len(_)))

    for type in types:
        mask = df['type'] == type
        for par in PARAMETERS:
            df.loc[mask, '{}_error'.format(par)] = (df.loc[mask, par] - real_df[par]) / real_df[par]

    asses_CIs(types)