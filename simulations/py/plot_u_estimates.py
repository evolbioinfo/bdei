from glob import glob
import re
import numpy as np
import pandas as pd
from matplotlib.pyplot import plot, show, legend, tight_layout, savefig, subplots
import os
from collections import Counter, defaultdict

par2greek = {'mu': u'\u03bc', 'lambda': u'\u03bb', 'psi': u'\u03c8', 'p': '\u03c1',
             'R_naught': u'\u0052\u2080' + '=' + u'\u03bb\u002F\u03c8',
             'infectious_time': 'infectious time 1' + u'\u002F\u03c8',
             'incubation_period': 'incubation period 1' + u'\u002F\u03bc'}

colours = ['#e41a1c', '#377eb8', '#4daf4a', '#984ea3', '#ff7f00', '#ffff33', '#a65628']


def parse_CI(ci):
    return [float(_) for _ in ci[1:-1].split(', ')]


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Plots us.")
    parser.add_argument('--data', type=str, default='/home/azhukova/projects/bdei_main/simulations/u', help="data folder")
    parser.add_argument('--svg', type=str, default='/home/azhukova/projects/bdei_main/simulations/u_params.svg', help="plot")
    params = parser.parse_args()

    fig, (ax1, ax2, ax3) = subplots(3, 1, figsize=(12, 16))
    real_us = []
    real_fs = []
    real_tips = []
    mu_u2within_CIs = np.zeros(21)
    la_u2within_CIs = np.zeros(21)
    psi_u2within_CIs = np.zeros(21)
    for i in range(100):

        for file in glob(os.path.join(params.data, 'forest.{i}_*.est'.format(i=i))):
            u = int(re.findall(r'_(\d+).est', file)[0]) + 0.25 * (i - 50)
            df = pd.read_csv(file, sep='\t')
            mu = df.loc[0, 'mu']
            la = df.loc[0, 'la']
            psi = df.loc[0, 'psi']
            mu_min, mu_max = parse_CI(df.loc[0, 'mu_CI'])
            la_min, la_max = parse_CI(df.loc[0, 'la_CI'])
            psi_min, psi_max = parse_CI(df.loc[0, 'psi_CI'])

            ax1.plot([u, u], [mu_min, mu_max], color=colours[0], linestyle='solid', alpha=0.1)
            ax2.plot([u, u], [la_min, la_max], color=colours[1], linestyle='solid', alpha=0.1)
            ax3.plot([u, u], [psi_min, psi_max], color=colours[2], linestyle='solid', alpha=0.1)

            u = int(re.findall(r'_(\d+).est', file)[0]) // 50
            mu_u2within_CIs[u] += 1 if mu_min <= REAL_MU <= mu_max else 0
            la_u2within_CIs[u] += 1 if la_min <= REAL_LA <= la_max else 0
            psi_u2within_CIs[u] += 1 if psi_min <= REAL_PSI <= psi_max else 0

        df = pd.read_csv(os.path.join(params.data, 'forest.{i}.log'.format(i=i)))
        REAL_U = df.loc[0, 'hidden_trees']
        REAL_TIPS = df.loc[0, 'tips']
        REAL_R = df.loc[0, 'R0']
        REAL_IT = df.loc[0, 'infectious time']
        REAL_IP = df.loc[0, 'incubation period']
        REAL_MU = 1 / REAL_IP
        REAL_PSI = 1 / REAL_IT
        REAL_LA = REAL_R * REAL_PSI
        real_us.append(REAL_U)
        REAL_F = open(os.path.join(params.data, 'forest.{i}.nwk'.format(i=i)), 'r').read().count(';')
        real_fs.append(REAL_F)
        real_tips.append(REAL_TIPS)

        hidden = '_' if i > 0 else ''
        ax1.plot([REAL_U], [REAL_MU], color=colours[0], marker='+', label='{}{}'.format(hidden, par2greek['mu']))
        ax2.plot([REAL_U], [REAL_LA], color=colours[1], marker='+', label='{}{}'.format(hidden, par2greek['lambda']))
        ax3.plot([REAL_U], [REAL_PSI], color=colours[2], marker='+', label='{}{}'.format(hidden, par2greek['psi']))

    ax1.plot([0, 1000], [REAL_MU, REAL_MU], color=colours[0])
    ax2.plot([0, 1000], [REAL_LA, REAL_LA], color=colours[1])
    ax3.plot([0, 1000], [REAL_PSI, REAL_PSI], color=colours[2])

    fig.set_size_inches(6, 8)

    for ax in (ax1, ax2, ax3):
        ax.legend()

        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)

        ax.set_xticks(list(np.arange(0, 1100, step=100)))


    ax3.set_xlabel('number of hidden trees u')
    ax2.set_ylabel('parameter value')
    # show()
    tight_layout()
    savefig(params.svg, dpi=300)

    print('Us vary between {} and {}, mean is {}'.format(min(real_us), max(real_us), np.mean(real_us)))
    print('Fs vary between {} and {}, mean is {}'.format(min(real_fs), max(real_fs), np.mean(real_fs)))
    print('Tips vary between {} and {}, mean is {}'.format(min(real_tips), max(real_tips), np.mean(real_tips)))

    percents = np.array(real_us) / (np.array(real_fs) + np.array(real_us))
    print('Proportions of hidden trees vary between {} and {}, mean is {}'.format(min(percents), max(percents), np.mean(percents)))

    fig, (ax1, ax2, ax3) = subplots(3, 1, figsize=(12, 16))

    us = np.arange(0, 1050, step=50)
    ax1.plot(us, mu_u2within_CIs, color=colours[0], linestyle='solid')
    ax2.plot(us, la_u2within_CIs, color=colours[1], linestyle='solid')
    ax3.plot(us, psi_u2within_CIs, color=colours[2], linestyle='solid')

    for ax in (ax1, ax2, ax3):
        ax.set_yticks(list(np.arange(0, 105, step=5)))

    print(us)
    print(mu_u2within_CIs)
    print(la_u2within_CIs)
    print(psi_u2within_CIs)

    show()
