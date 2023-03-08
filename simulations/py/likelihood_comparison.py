import pandas as pd
from pybdei import get_loglikelihood

REAL_TYPE = 'real'


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Calculates likelihoods.")
    parser.add_argument('--estimates', type=str, help="estimated parameter", required=True)
    parser.add_argument('--tree_pattern', type=str, help="trees", required=True)
    parser.add_argument('--likelihoods', type=str, help="likelihood table")
    parser.add_argument('--log', type=str, help="likelihood stats")
    params = parser.parse_args()

    df = pd.read_csv(params.estimates, sep='\t', index_col=0)
    NON_REAL_TYPES = [_ for _ in df['type'].unique() if _ != REAL_TYPE]
    ALL_TYPES = [REAL_TYPE] + NON_REAL_TYPES

    df.index = df['type'] + '_' + df.index.map(str)
    lk_df, index = [], []
    for i in range(100):
        if 'real_{}'.format(i) not in df.index:
            continue
        p = df.loc['real_{}'.format(i), 'p']
        type2lk = {}
        for type in ALL_TYPES:
            mu, la, psi = df.loc['{}_{}'.format(type, i), ['mu', 'lambda', 'psi']]
            type2lk[type] = get_loglikelihood(nwk=params.tree_pattern.format(i), mu=mu, la=la, psi=psi, p=p)
            index.append(i)
        best_lk = max(type2lk.values())
        for _ in ALL_TYPES:
            lk_df.append([_, type2lk[_], best_lk - type2lk[_]])
    lk_df = pd.DataFrame(data=lk_df, columns=['type', 'loglk', 'diff with max'], index=index)
    lk_df.to_csv(params.likelihoods, sep='\t')

    log_lk_r = lk_df[lk_df['type'] == REAL_TYPE]

    with open(params.log, 'w+') as f:
        for type in ALL_TYPES:
            log_lk_t = lk_df[lk_df['type'] == type]
            log = '{} likelihood is the highest in {}% of cases.'.format(type, sum(log_lk_t['diff with max'] == 0))
            print(log)
            f.write(log + '\n')

        for type in NON_REAL_TYPES:
            log_lk_t = lk_df[lk_df['type'] == type]
            log = '{} likelihood is higher than real in {}% of cases'\
                .format(type, sum(log_lk_t['loglk'] > log_lk_r['loglk']))
            print(log)
            f.write(log + '\n')

        for i in range(len(NON_REAL_TYPES) - 1):
            type1 = NON_REAL_TYPES[i]
            log_lk_t1 = lk_df[lk_df['type'] == type1]
            for j in range(i + 1, len(NON_REAL_TYPES)):
                type2 = NON_REAL_TYPES[j]
                log_lk_t2 = lk_df[lk_df['type'] == type2]
                log = '{} likelihood is higher than {} in {}% of cases'\
                    .format(type1, type2, sum(log_lk_t1['loglk'] > log_lk_t2['loglk']))
                print(log)
                print(log)
                f.write(log + '\n')
