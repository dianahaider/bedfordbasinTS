import numpy as np
import pandas as pd
from joblib import Parallel, delayed

def rarefy_replicate(df_wide_filtered, sampling_depth, rep, random_state=None):
    if random_state is not None:
        np.random.seed(random_state + rep)
    rarefied_rows = {}
    for sample in df_wide_filtered.index:
        counts = df_wide_filtered.loc[sample].values
        total = counts.sum()
        proportions = counts / total
        # np.random.multinomial simulates sampling with replacement.
        rarefied = np.random.multinomial(sampling_depth, proportions)
        rarefied_rows[sample] = rarefied
    rep_df = pd.DataFrame.from_dict(rarefied_rows, orient='index', columns=df_wide_filtered.columns)
    return rep_df


def repeat_rarefy_df(df, sampling_depth: int, repeat_times: int,
                      with_replacement: bool = True, random_state: int = None, n_jobs: int = -1) -> pd.DataFrame:

    df_wide = df.pivot_table(index='sampleid', columns='feature_id',
                             values='feature_frequency', fill_value=0)


    sample_totals = df_wide.sum(axis=1)
    samples_to_keep = sample_totals[sample_totals >= sampling_depth].index
    df_wide_filtered = df_wide.loc[samples_to_keep]

    print(f"Removed {df_wide.shape[0] - df_wide_filtered.shape[0]} samples from {df_wide.shape[0]} below the rarefaction depth of {sampling_depth}.")


    sample_totals = df_wide_filtered.sum(axis=1)
    if (sample_totals < sampling_depth).any():
        missing = sample_totals[sample_totals < sampling_depth]
        raise ValueError(f"Some samples have total counts below the sampling depth:\n{missing}")


    replicate_list = Parallel(n_jobs=n_jobs)(
        delayed(rarefy_replicate)(df_wide_filtered, sampling_depth, rep, random_state)
        for rep in range(repeat_times)
    )


    mean_rarefied = sum(replicate_list) / repeat_times

    mean_rarefied = np.ceil(mean_rarefied).astype(int)

    return mean_rarefied
