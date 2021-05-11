#!/bin/env python

#### select and plot highly variable regions (HVR) ####

#### libraries
import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

# for robust dispersion norm calculation
from statsmodels import robust
import warnings
warnings.filterwarnings("ignore")

#### configurations

# input
data=pd.read_csv(snakemake.input[0], index_col=0)

# parameters
HVR_percentage=snakemake.params['HVR_percentage']
top_n = int((data.shape[0]/100)*HVR_percentage)

# output
output_data=snakemake.output[0]
output_plot=snakemake.output[1]

#####  determine highly variable regions (HRV) based on normalized dispersion
## inspired by “highly_variable_genes” function from scanpy 1.5.1 (Python 3.6.10) with the “flavor” argument set to “cell_ranger”
## but instead of dispersion=var/mean we use dispersion=std to keep it meaningful for negative mean values
## overlap of region selecrtion in test data was 96.6%

mean = data.mean(axis=1)
# var = data.var(axis=1)

# now actually compute the dispersion
mean[mean == 0] = 1e-12  # set entries equal to zero to small value
dispersion = data.std(axis=1) #var / mean

# all of the following quantities are "per-gene" here
df = pd.DataFrame()
df['means'] = mean
df['dispersions'] = dispersion

# bin regions by their mean in 20 bins
df['mean_bin'] = pd.cut(df['means'], np.r_[
    -np.inf,
    np.percentile(df['means'], np.arange(10, 105, 5)),
    np.inf
])
# group regions by mean_bin
disp_grouped = df.groupby('mean_bin')['dispersions']
# determine median (robust) dispersion by group
disp_median_bin = disp_grouped.median()

# the next line raises the warning: "Mean of empty slice"
with warnings.catch_warnings():
    warnings.simplefilter('ignore')
    # calculate median absolute deviation (mad; robust measure of variance) per group
    disp_mad_bin = disp_grouped.apply(robust.mad)
    # determine "normalized" dispersion via (dispersion-dispersion_group_median)/dispersion_mad
    df['dispersions_norm'] = (df['dispersions'].values
        - disp_median_bin[df['mean_bin'].values].values
        ) / disp_mad_bin[df['mean_bin'].values].values
        
dispersion_norm = df['dispersions_norm'].values.astype('float32')



HVR_idx=df.index[df['dispersions_norm'].rank(ascending=False)<=top_n]
# HVR_idx=std.index[std.rank(ascending=False)<=top_n]
data_HVR = data.loc[HVR_idx, :]
data_HVR.to_csv(output_data)

# make plots describing the HVR selection
fig, axs = plt.subplots(1, 2, figsize=(10, 5))
plt.subplots_adjust(wspace=0.5, hspace=0.5)

# plot region variability
ax = sns.scatterplot(df['dispersions_norm'].rank(ascending=False),
                     df['dispersions_norm'],
                     ax=axs[0], 
                     rasterized=True,
                    )
ax = sns.scatterplot(df.loc[HVR_idx,'dispersions_norm'].rank(ascending=False),
                     df.loc[HVR_idx,'dispersions_norm'],
                     ax=axs[0], 
                     rasterized=True,
                     color='r'
                    )
ax.set_xlabel('rank')
ax.set_ylabel('normalized dispersion')
sns.despine()
ax.set_title("Regions ranked by normalized dispersion")

# mean-dispersion plot (highlighting selected regions)
ax = sns.scatterplot(df.loc[:,'means'],
                     df.loc[:,'dispersions_norm'],
                     ax=axs[1], 
                     rasterized=True,
                    )
ax = sns.scatterplot(df.loc[HVR_idx,'means'], 
                     df.loc[HVR_idx,'dispersions_norm'],
                     ax=axs[1], 
                     rasterized=True,
                     color='r'
                    )

ax.set_xlabel('mean')
ax.set_ylabel('normalized dispersion')
sns.despine()
ax.set_title("Mean to normalized dispersion relationship")

fig.suptitle('Top {} HVRs by binned normalized dispersion'.format(top_n), fontsize=14)

plt.savefig(
            fname=output_plot,
            format="svg",
            dpi=300,
            bbox_inches="tight",
            )

plt.show()
plt.close()



