from glob import glob
import re
import numpy as np
import pandas as pd
from matplotlib.pyplot import plot, show, legend, tight_layout, savefig, subplots, subplots_adjust
from matplotlib.offsetbox import TextArea, HPacker, AnchoredOffsetbox, VPacker
import os
from collections import Counter, defaultdict
from pastml.visualisation.colour_generator import get_enough_colours

MAX_U = 500

par2greek = {'mu': u'\u03bc', 'lambda': u'\u03bb', 'psi': u'\u03c8', 'p': '\u03c1'}

colours = ['#e41a1c', '#377eb8', '#4daf4a', '#984ea3', '#ff7f00', '#ffff33', '#a65628']

# colours = get_enough_colours(100)

def parse_CI(ci):
    return [float(_) for _ in ci[1:-1].split(', ')]


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Plots us.")
    parser.add_argument('--data', type=str, default='/home/azhukova/projects/bdei_main/simulations/u', help="data folder")
    parser.add_argument('--svg', type=str, default='/home/azhukova/projects/bdei_main/simulations/u_params.svg', help="plot")
    parser.add_argument('--fixed', type=str, default='p', help="parameter that was fixed during estimations")
    params = parser.parse_args()

    fig, (ax1, ax2, ax3) = subplots(3, 1, figsize=(12, 16))
    real_us = []
    real_fs = []
    real_tips = []
    mu_u2within_CIs = np.zeros(MAX_U // 50 + 1)
    la_u2within_CIs = np.zeros(MAX_U // 50 + 1)
    psi_u2within_CIs = np.zeros(MAX_U // 50 + 1)
    p_u2within_CIs = np.zeros(MAX_U // 50 + 1)

    mu_u2total = np.zeros(MAX_U // 50 + 1)
    la_u2total = np.zeros(MAX_U // 50 + 1)
    psi_u2total = np.zeros(MAX_U // 50 + 1)
    p_u2total = np.zeros(MAX_U // 50 + 1)

    mu_i2within_CIs = defaultdict(set)
    la_i2within_CIs = defaultdict(set)
    psi_i2within_CIs = defaultdict(set)
    p_i2within_CIs = defaultdict(set)

    for i in range(100):

        df = pd.read_csv(os.path.join(params.data, 'forest.{i}.log'.format(i=i)))
        REAL_U = df.loc[0, 'hidden_trees']
        REAL_TIPS = df.loc[0, 'tips']
        REAL_R = df.loc[0, 'R0']
        REAL_P = df.loc[0, 'sampling probability']
        REAL_IT = df.loc[0, 'infectious time']
        REAL_IP = df.loc[0, 'incubation period']
        REAL_MU = 1 / REAL_IP
        REAL_PSI = 1 / REAL_IT
        REAL_LA = REAL_R * REAL_PSI
        real_us.append(REAL_U)
        REAL_F = open(os.path.join(params.data, 'forest.{i}.nwk'.format(i=i)), 'r').read().count(';')
        real_fs.append(REAL_F)
        real_tips.append(REAL_TIPS)

        def round_to_fifty(value):
            return min(MAX_U, int(((value + 25) // 50) * 50))

        # print(REAL_U, round_to_fifty(REAL_U / 4), round_to_fifty(REAL_U * 4))
        # for u in range(round_to_fifty(REAL_U / 2), round_to_fifty(REAL_U * 1.5) + 50, 50):
        for u in range(0, MAX_U + 50, 50):
        # for file in glob(os.path.join(params.data, 'forest.{par}.{i}_*.est'.format(i=i, par=params.fixed))):
        #     u = int(re.findall(r'_(\d+).est', file)[0])
            file = os.path.join(params.data, 'forest.{par}.{i}_{u}.est'.format(i=i, par=params.fixed, u=u))
            df = pd.read_csv(file, sep='\t')
            mu = df.loc[0, 'mu']
            la = df.loc[0, 'la']
            psi = df.loc[0, 'psi']
            p = df.loc[0, 'p']
            mu_min, mu_max = parse_CI(df.loc[0, 'mu_CI'])
            la_min, la_max = parse_CI(df.loc[0, 'la_CI'])
            psi_min, psi_max = parse_CI(df.loc[0, 'psi_CI'])
            p_min, p_max = parse_CI(df.loc[0, 'p_CI'])

            xs = [u + 0.42 * (i - 50)] * 2
            ax1.plot(xs, [mu_min, mu_max], linestyle='solid', color=colours[0], alpha=0.2)
            ax2.plot(xs, [la_min, la_max], linestyle='solid', color=colours[1], alpha=0.2)
            if 'p' == params.fixed:
                ax3.plot(xs, [psi_min, psi_max], linestyle='solid', color=colours[2], alpha=0.2)
            else:
                ax3.plot(xs, [p_min, p_max], linestyle='solid', color=colours[3], alpha=0.2)

            mu_ok = mu_min <= REAL_MU <= mu_max
            mu_u2within_CIs[u // 50] += 1 if mu_ok else 0
            la_ok = la_min <= REAL_LA <= la_max
            la_u2within_CIs[u // 50] += 1 if la_ok else 0
            psi_ok = psi_min <= REAL_PSI <= psi_max
            psi_u2within_CIs[u // 50] += 1 if psi_ok else 0
            p_ok = p_min <= REAL_P <= p_max
            p_u2within_CIs[u // 50] += 1 if p_ok else 0

            mu_u2total[u // 50] += 1
            psi_u2total[u // 50] += 1
            la_u2total[u // 50] += 1
            p_u2total[u // 50] += 1

            if mu_ok:
                mu_i2within_CIs[i].add(u)
            if la_ok:
                la_i2within_CIs[i].add(u)
            if psi_ok:
                psi_i2within_CIs[i].add(u)
            if p_ok:
                p_i2within_CIs[i].add(u)
        hidden = '_' if i > 0 else ''
        ax1.plot([REAL_U], [REAL_MU], color=colours[0], marker='+')
        ax2.plot([REAL_U], [REAL_LA], color=colours[1], marker='+')
        if 'p' == params.fixed:
            ax3.plot([REAL_U], [REAL_PSI], color=colours[2], marker='+')
        else:
            ax3.plot([REAL_U], [REAL_P], color=colours[3], marker='+')

    xs = [0, MAX_U]
    ax1.plot(xs, [REAL_MU, REAL_MU], color=colours[0])
    ax1.set_title(par2greek['mu'], y=0.9)
    ax2.plot(xs, [REAL_LA, REAL_LA], color=colours[1])
    ax2.set_title(par2greek['lambda'], y=0.9)
    if 'p' == params.fixed:
        ax3.plot(xs, [REAL_PSI, REAL_PSI], color=colours[2])
        ax3.set_title(par2greek['psi'], y=0.9)
    else:
        ax3.plot(xs, [REAL_P, REAL_P], color=colours[3])
        ax3.set_title(par2greek['p'], y=0.9)




    def get_pbox(par, ax):
        boxes = []
        u2within_CIs = mu_u2within_CIs if par == 'mu' \
            else (psi_u2within_CIs if par == 'psi' else (la_u2within_CIs if par == 'la' else p_u2within_CIs))
        u2total = mu_u2total if par == 'mu' \
            else (psi_u2total if par == 'psi' else (la_u2total if par == 'la' else p_u2total))
        for u in range(0, MAX_U + 50, 50):
            total = u2total[int(u // 50)]
            v = int(np.round(100 * u2within_CIs[int(u // 50)] / total, 0)) if total else ' '
            boxes.append(
                        TextArea('{}'.format(v) if isinstance(v, int) and v >= 10 else ' {}'.format(v),
                                 textprops=dict(color='black', ha='center', va='center',
                                                fontsize='small', fontweight='normal', family='monospace')))
        return AnchoredOffsetbox(loc=3, child=HPacker(children=boxes, align="center", pad=0, sep=4.5 if MAX_U == 1000 else 19),
                                 pad=0, frameon=False, bbox_to_anchor=(0.05, 0.05),
                                 bbox_transform=ax.transAxes, borderpad=0.)


    for ax in (ax1, ax2, ax3):
        # ax.legend()

        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)

        ax.set_xticks(list(np.arange(0, MAX_U + 100, step=100)))

    ax3.set_xlabel('number of hidden trees u')
    ax2.set_ylabel('parameter value')
    # tight_layout()

    ax1.add_artist(get_pbox('mu', ax1))
    ax2.add_artist(get_pbox('la', ax2))
    if 'p' == params.fixed:
        ax3.add_artist(get_pbox('psi', ax3))
    else:
        ax3.add_artist(get_pbox('p', ax3))


    fig.set_size_inches(6, 8)
    # show()
    savefig(params.svg, dpi=300)

    print('Us vary between {} and {}, mean is {}, median is {} [{}-{}]'
          .format(min(real_us), max(real_us), np.mean(real_us), np.median(real_us), np.percentile(real_us, 2.5), np.percentile(real_us, 97.5)))
    print('Fs vary between {} and {}, mean is {}, median is {} [{}-{}]'
          .format(min(real_fs), max(real_fs), np.mean(real_fs), np.median(real_fs), np.percentile(real_fs, 2.5), np.percentile(real_fs, 97.5)))
    print('Tips vary between {} and {}, mean is {}, median is {} [{}-{}]'
          .format(min(real_tips), max(real_tips), np.mean(real_tips), np.median(real_tips), np.percentile(real_tips, 2.5), np.percentile(real_tips, 97.5)))

    percents = np.array(real_us) / (np.array(real_fs) + np.array(real_us))
    print('Proportions of hidden trees vary between {} and {}, mean is {}, median is {} [{}-{}]'
          .format(min(percents), max(percents), np.mean(percents), np.median(percents), np.percentile(percents, 2.5), np.percentile(percents, 97.5)))
