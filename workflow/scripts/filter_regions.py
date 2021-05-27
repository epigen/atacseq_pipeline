#!/bin/env python

#### filter regions by support and mean/variance distribution ####

#### libraries
import pybedtools
import os
import time
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

import warnings
warnings.filterwarnings("ignore")


# utility functions
def cpm(counts_matrix, feature_len=1, norm_factors=1, log=False, pseudocount=0.5, libsize_pseudocount=1, million=1e6):
    """
    Counts per million normalization
    If feature_len=1 then you get counts per million
    Setting pseudocount=0.5, libsize_pseudocount=1 ensures results equivalent to LIMMA-voom
    """
    rpk = counts_matrix.astype(float) / feature_len
    cpm = million * (rpk + pseudocount) / (rpk.sum(axis=0) * norm_factors + libsize_pseudocount)

    return np.log2(cpm) if log else cpm

def filter_by_reads(counts_df, min_group_n, cpm_df=None, min_count=10, min_total_count=15, large_min_n=10,
                    large_min_n_proportion=0.7, verbose=True):

    if min_group_n > large_min_n:
        # edgeR
        min_n = large_min_n + (min_group_n - large_min_n) * large_min_n_proportion
    else:
        min_n = min_group_n

    _million, _tolerance = 1e6, 1e-14
    median_lib_size = counts_df.sum(axis=0).median()
    cpm_cutoff = min_count / median_lib_size * _million
    keep_cpm = (cpm_df >= cpm_cutoff).sum(axis=1) >= min_n - _tolerance
    keep_min_total_count = counts_df.sum(axis=1) >= min_total_count - _tolerance
    if verbose:
        print('min_n', min_n)
        print('median_lib_size', median_lib_size)
        print('cpm_cutoff', cpm_cutoff)
        print('remove based on cpm_cutoff', (~keep_cpm).sum())
        print('additionally remove based on keep_min_total_count', (~keep_min_total_count[keep_cpm]).sum())
    return keep_cpm & keep_min_total_count

def plot_sequenced_samples(df, n_samples=30, ax=None, xlabel='Read count', ylabel='Density', title=None, xlim=None, ylim=None,
                           kde=True, hist=False, legend=False, samples_as_rows=False):
    if samples_as_rows:
        df = df.T
    for i in np.random.choice(list(range(df.shape[1])), size=n_samples, replace=False) if n_samples is not None else range(df.shape[1]):
        sample = df.iloc[:, i]
        ax = sns.distplot(sample, ax=ax, kde=kde, hist=hist, label=df.columns[i] if legend else None)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    if xlim is not None:
        ax.set_xlim(xlim)
    if ylim is not None:
        ax.set_ylim(ylim)
    if title is None:
        ax.set_title('{} samples{}{}'.format(n_samples, ', xlim: {}'.format(xlim) if xlim is not None else '', ', ylim: {}'.format(ylim) if ylim is not None else ''))
    else:
        ax.set_title(title)
    if legend:
        ax.legend(bbox_to_anchor=(1, 1))
    sns.despine()
    return ax


#### configurations

# input
data_counts=pd.read_csv(snakemake.input[0], index_col=0)
support=pd.read_csv(snakemake.input[1], index_col=0)
annot=pd.read_csv(snakemake.input[2], index_col=0)

# parameters
peak_support_threshold = snakemake.params["peak_support_threshold"]
large_min_n_proportion = snakemake.params["proportion"]
min_group = snakemake.params["min_group"]

# output
output_data=snakemake.output[0]
output_plot=snakemake.output[1]

#####  filter regions by peak support
support = support.loc[data_counts.index,:]
support_sum = support.sum(axis=1)
region_filter_by_support = (support_sum>=peak_support_threshold)
data_filtered=data_counts.loc[region_filter_by_support,:]

print('before filtering by support', data_counts.shape[0])
print('after filtering by support', data_filtered.shape[0])

##### filter regions by mean/variance distribution
data_filtered_cpm = cpm(data_filtered, feature_len=1, norm_factors=1, log=False, pseudocount=0.5, libsize_pseudocount=0)

assert data_filtered_cpm.columns.equals(data_filtered.columns)
assert data_filtered_cpm.index.equals(data_filtered.index)

# determine min_group_n
if min_group=='':
    min_group_n=data_filtered.shape[1]
else:
    min_group_n=annot.groupby(min_group)[min_group].count().min()

min_count_mask = filter_by_reads(
    data_filtered,
    min_group_n=min_group_n,
    cpm_df=data_filtered_cpm,
    large_min_n_proportion=large_min_n_proportion,
)

print('before filtering', data_filtered.shape[0])
print('after filtering', data_filtered.loc[min_count_mask, :].shape[0])

lcpm = np.log2(data_filtered_cpm)

fig, axs = plt.subplots(2, 2, figsize=(10, 10))
plt.subplots_adjust(wspace=0.5, hspace=0.5)

ax = sns.scatterplot(lcpm.mean(axis=1),
                     lcpm.std(axis=1),
                     ax=axs[0, 0], rasterized=True)
ax.set_xlabel('Mean')
ax.set_ylabel('Std')
sns.despine()
ax.set_title('before filtering\n{} peaks'.format(lcpm.shape[0]))

ax = sns.scatterplot(lcpm.loc[min_count_mask, :].mean(axis=1),
                     lcpm.loc[min_count_mask, :].std(axis=1),
                     ax=axs[0, 1], rasterized=True)
ax.set_xlabel('Mean')
ax.set_ylabel('Std')
sns.despine()
ax.set_title('after filtering\n{} peaks'.format(min_count_mask.sum()))

plot_sequenced_samples(lcpm, n_samples=None, ax=axs[1, 0], xlabel='Log$2$ CPM', ylabel='Density', title='before filtering', samples_as_rows=False)
plot_sequenced_samples(lcpm.loc[min_count_mask, :], n_samples=None, ax=axs[1, 1], xlabel='Log$2$ CPM', ylabel='Density', title='after filtering', samples_as_rows=False)

fig.suptitle('Region filtering with proportion={} and min group size={}'.format(large_min_n_proportion,min_group_n), fontsize=10)

plt.savefig(output_plot, dpi=300)
plt.show()
plt.close()

# apply filter on pre filtered counts
filtered_counts = data_filtered.loc[min_count_mask, :]
# save filtered counts
filtered_counts.to_csv(output_data)