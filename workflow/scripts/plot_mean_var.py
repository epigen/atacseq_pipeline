#!/bin/env python

#### plot mean variance relationship to inform on variance stabilizing measures ####

#### libraries

# general
import os
import pandas as pd
import numpy as np

# visualization
import seaborn as sns
import matplotlib.pyplot as plt

# mean-variance relationship
from sklearn.metrics import r2_score
from scipy import stats


#### configurations

# load data
data=pd.read_csv(snakemake.input[0], index_col=0)

mean=data.mean(axis=1)
std=data.std(axis=1)

xaxis='mean'
yaxis='std'
title='Mean-StD Relationship'
r2=r2_score(mean, std)
spr=stats.spearmanr(mean, std)[0]

plt.scatter(mean, std, marker='.', alpha=0.5)
plt.xlabel(xaxis)
plt.ylabel(yaxis) 
plt.title("{} R2={:.2f} rho={:.2f}".format(title,r2,spr))

# save plot
plt.savefig(fname=snakemake.output[0],
            format="svg",
            dpi=300,
            bbox_inches="tight",
           )