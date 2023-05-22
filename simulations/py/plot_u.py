from glob import glob
import re
import numpy as np
import pandas as pd
from matplotlib.pyplot import plot, show, legend, tight_layout, savefig, subplots
import os

par2greek = {'mu': u'\u03bc', 'lambda': u'\u03bb', 'psi': u'\u03c8', 'p': '\u03c1',
             'R_naught': u'\u0052\u2080' + '=' + u'\u03bb\u002F\u03c8',
             'infectious_time': 'infectious time 1' + u'\u002F\u03c8',
             'incubation_period': 'incubation period 1' + u'\u002F\u03bc'}

colours = ['#e41a1c', '#377eb8', '#4daf4a', '#984ea3', '#ff7f00', '#ffff33', '#a65628']

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Plots us.")
    parser.add_argument('--data', type=str, help="data folder")
    parser.add_argument('--svg', type=str, help="plot")
    params = parser.parse_args()

    fig, (ax1, ax2) = subplots(2, 1, figsize=(18, 8))
    real_us = []
    real_fs = []
    real_tips = []

    max_v, max_V = 0, 0
    for i in range(100):
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

        # lks = []
        # us = []
        # for f in glob(os.path.join(params.data, 'forest.{i}_*.lk'.format(i=i))):
        #     u = re.findall(r'_(\d+).lk', f)[0]
        #     lk = float(open(f, 'r').read())
        #     us.append(u)
        #     lks.append(lk)
        #
        # us = np.array(us, dtype=int)
        # neworder = np.argsort(us)
        # us = us[neworder]
        # ts = us + REAL_F
        # lks = np.array(lks, dtype=int)[neworder]
        # plot(us, lks / ts)


        us = []
        mus = []
        las = []
        psis = []
        for f in glob(os.path.join(params.data, 'forest.{i}_*.est'.format(i=i))):
            u = re.findall(r'_(\d+).est', f)[0]
            df = pd.read_csv(f, sep='\t')
            us.append(u)
            mus.append(df.loc[0, 'mu'])
            las.append(df.loc[0, 'la'])
            psis.append(df.loc[0, 'psi'])

        us = np.array(us, dtype=int)
        neworder = np.argsort(us)
        us = us[neworder]
        mus = np.array(mus, dtype=float)[neworder]
        las = np.array(las, dtype=float)[neworder]
        psis = np.array(psis, dtype=float)[neworder]

        rs = las / psis
        incts = 1 / mus
        inft = 1 / psis

        hidden = '_' if i > 0 else ''

        ax1.plot(us, mus, color=colours[0], linestyle='solid', alpha=0.1)
        ax1.plot(us, las, color=colours[1], linestyle='solid', alpha=0.1)
        ax1.plot(us, psis, color=colours[2], linestyle='solid', alpha=0.1)

        max_v = max(max_v, max(mus), max(las), max(psis))

        ax2.plot(us, rs, color=colours[3], linestyle='solid', alpha=0.1)
        ax2.plot(us, incts, color=colours[4], linestyle='solid', alpha=0.1)
        ax2.plot(us, inft, color=colours[5], linestyle='solid', alpha=0.1)

        max_V = max(max_V, max(rs), max(incts), max(inft))

        ax1.plot([REAL_U], [REAL_MU], color=colours[0], marker='d', label='{}{}'.format(hidden, par2greek['mu']))
        ax1.plot([REAL_U], [REAL_LA], color=colours[1], marker='d', label='{}{}'.format(hidden, par2greek['lambda']))
        ax1.plot([REAL_U], [REAL_PSI], color=colours[2], marker='d', label='{}{}'.format(hidden, par2greek['psi']))
        ax2.plot([REAL_U], [REAL_R], color=colours[3], marker='d', label='{}{}'.format(hidden, par2greek['R_naught']))
        ax2.plot([REAL_U], [REAL_IP], color=colours[4], marker='d', label='{}{}'.format(hidden, par2greek['incubation_period']))
        ax2.plot([REAL_U], [REAL_IT], color=colours[5], marker='d', label='{}{}'.format(hidden, par2greek['infectious_time']))

    fig.set_size_inches(6, 4)

    ax1.set_yticks(list(np.arange(0, max_v + 0.2, step=0.2)))
    ax2.set_yticks(list(np.arange(0, max_V + 0.5, step=0.5)))

    for ax in (ax1, ax2):
        ax.legend()

        ax.set_xlabel('number of hidden trees u')
        ax.set_ylabel('parameter value')

        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)


    # show()
    tight_layout()
    savefig(params.svg, dpi=300)

    print('Us vary between {} and {}, mean is {}'.format(min(real_us), max(real_us), np.mean(real_us)))
    print('Fs vary between {} and {}, mean is {}'.format(min(real_fs), max(real_fs), np.mean(real_fs)))
    print('Tips vary between {} and {}, mean is {}'.format(min(real_tips), max(real_tips), np.mean(real_tips)))
