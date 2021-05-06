#!/bin/env python

#### select HVR and plot their standard deviation (std) ####

#### libraries
import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

#### configurations

# input
data=pd.read_csv(snakemake.input[0], index_col=0)

# parameters
HVR_percentage=snakemake.params['HVR_percentage']
top_n = int((data.shape[0]/100)*HVR_percentage)

# output
output_data=snakemake.output[0]
output_plot=snakemake.output[1]


#####  determine highly variable regions (HRV) by standard deviation (std)
std=data.std(axis=1)

HVR_idx=std.index[std.rank(ascending=False)<=top_n]
data_HVR = data.loc[HVR_idx, :]
data_HVR.to_csv(output_data)

# plot region variability
plt.scatter(std.rank(ascending=False), std, marker='.')
plt.vlines(x=top_n, ymin=0, ymax=max(std), color='r')
plt.title("Ranked region variability with HVR threshold at {}".format(top_n))
plt.xlabel('rank')
plt.ylabel('standard deviation') 
plt.savefig(
            fname=output_plot,
            format="svg",
            dpi=300,
            bbox_inches="tight",
            )