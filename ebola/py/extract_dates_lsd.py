import pandas as pd


def date2years(d, default=None, default_min_date=1900, default_max_date=2018):
    if pd.notnull(d):
        first_jan_this_year = pd.datetime(year=d.year, month=1, day=1)
        day_of_this_year = d - first_jan_this_year
        first_jan_next_year = pd.datetime(year=d.year + 1, month=1, day=1)
        days_in_this_year = first_jan_next_year - first_jan_this_year
        date = d.year + day_of_this_year / days_in_this_year
        max_date, min_date = date, date
        if date == d.year:
            max_date = d.year + (days_in_this_year.days - 1) / days_in_this_year.days
            return None, date, max_date
        return date, min_date, max_date
    else:
        return default, default_min_date, default_max_date


if '__main__' == __name__:
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument('--data', required=True, type=str)
    parser.add_argument('--dates', required=True, type=str)
    parser.add_argument('--date_col', required=True, type=str)
    params = parser.parse_args()

    df = pd.read_csv(params.data, index_col=0)[[params.date_col]]
    try:
        df[params.date_col] = pd.to_datetime(df[params.date_col], format='%d/%m/%Y')
    except ValueError:
        try:
            df[params.date_col] = pd.to_datetime(df[params.date_col], infer_datetime_format=True)
        except ValueError:
            try:
                df[params.date_col] = pd.to_datetime(df[params.date_col], format='%Y.0')
            except ValueError:
                raise ValueError('Could not infer the date format for column "{}", please check it.'
                                 .format(params.date_col))

    m_date = date2years(df[~pd.isna(df[params.date_col])][params.date_col].min())[1]
    M_date = date2years(df[~pd.isna(df[params.date_col])][params.date_col].max())[2]
    df[['date', 'lower', 'upper']] = df[params.date_col].apply(lambda _: date2years(_, None, m_date, M_date)).apply(
        pd.Series)
    df.loc[pd.isna(df['date']), 'date'] = 'b(' + df.loc[pd.isna(df['date']), 'lower'].astype(str) \
                                          + ',' + df.loc[pd.isna(df['date']), 'upper'].astype(str) + ')'
    if params.dates:
        with open(params.dates, 'w+') as f:
            f.write('%d\n' % df.shape[0])
        df['date'].to_csv(params.dates, sep='\t', header=False, mode='a')
