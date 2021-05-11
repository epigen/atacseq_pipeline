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
import statsmodels.formula.api as sm
from scipy import stats


#### configurations

# load data
data=pd.read_csv(snakemake.input[0], index_col=0)

if not os.path.exists(snakemake.params['results_dir']):
    os.mkdir(snakemake.params['results_dir'])

mean=data.mean(axis=1)
std=data.std(axis=1)

xaxis='mean'
yaxis='std'
title='Mean-StD Relationship'

lm_result = sm.ols(formula="std ~ mean", data=pd.DataFrame({'std': std, 'mean': mean})).fit()
# print(lm_result.summary())
# print(lm_result.rsquared, lm_result.rsquared_adj)
p = lm_result.params
x_dummy = np.arange(min(mean), max(mean))

spr=stats.spearmanr(mean, std)

plt.scatter(mean, std, marker='.', alpha=0.5)
plt.plot(x_dummy, p.Intercept + p['mean'] * x_dummy, c='red')

plt.xlabel(xaxis)
plt.ylabel(yaxis) 
plt.title("{} R2={:.2f} rho={:.2f} p={:.2f}".format(title,lm_result.rsquared,spr[0],spr[1]))

# save plot
plt.savefig(fname=snakemake.output[0],
            format="svg",
            dpi=300,
            bbox_inches="tight",
           )