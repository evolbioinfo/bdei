import re
from collections import defaultdict

import numpy as np
import pandas as pd

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Combine time stats.")
    parser.add_argument('--time', nargs='+', type=str, help="times")
    parser.add_argument('--tab', type=str, help="estimate table")
    parser.add_argument('--cores', action='store_true')
    params = parser.parse_args()

    cores2results = defaultdict(list)
    cores2times = defaultdict(list)
    cores2iterations = defaultdict(list)
    for file in params.time:
        df = pd.read_csv(file, sep='\t')
        time, n_iter = df.loc[next(iter(df.index)), ['CPU_time', 'iterations']]
        cores = int(re.findall(r'threads=(\d+)', file)[0]) if params.cores else 1
        cores2results[cores].append(time / n_iter)
        cores2times[cores].append(time)
        cores2iterations[cores].append(n_iter)
    with open(params.tab, 'w+') as f:
        f.write('cores\ts\titerations\ttotal\n')
        for cores in sorted(cores2results.keys()):
            results = cores2results[cores]
            iterations = cores2iterations[cores]
            times = cores2times[cores]
            print("{} cores:\t{} s\t{} iterations\t{} s".format(cores, np.mean(results), np.mean(iterations), np.mean(times)))
            f.write('{}\t{}\t{}\t{}\n'.format(cores, np.mean(results), np.mean(iterations), np.mean(times)))
