#!/bin/env python

#### generic dimensionality reduction visualization function (UMAP & PCA) for discrete and continous variables ####

#### libraries

# general
import os
import pandas as pd
import numpy as np
from itertools import compress

# visualization
import seaborn as sns
import matplotlib.pyplot as plt

# # utils for enumerating string lists
# from collections import defaultdict
# from itertools import count
# from functools import partial


#### configurations

# annot=observations x annotations
annot=pd.read_csv(snakemake.input[0], index_col=0)

# type of dimensionality reduction
dimred=snakemake.params["dimred"]

# variables=columns of annot to plot
variables=snakemake.params['variables']

#  label=plot title and file name
label = snakemake.params["label"]

# results folder
results_dir = snakemake.params["results_dir"]

# if 3 samples or less UMAP could not be performed
if (dimred=="UMAP" and annot.shape[0]<4):
    from pathlib import Path
    for variable in variables:
        Path(os.path.join(results_dir, dimred+"_"+label+"_"+variable+".svg")).touch()
    import sys
    sys.exit()

# if 2 samples or less PCA only consists of one PC
if (dimred=="PCA" and annot.shape[0]<3):
    from pathlib import Path
    for variable in variables:
        Path(os.path.join(results_dir, dimred+"_"+label+"_"+variable+".svg")).touch()
    import sys
    sys.exit()

# data=observations x features
data_df=pd.read_csv(snakemake.input[1], index_col=0)

#plot parameters
color_dict=None
alpha=1


# plot data in 2D
# sns.set(color_codes=True)
for variable in variables:
    fig, axes = plt.subplots(nrows=1, ncols=1, figsize=(5, 5))
    unique_vals = annot[variable].unique()
    
    if len(unique_vals)==1 or (sum(pd.isna(unique_vals))>0):
        print("{} only one value or contains NaNs".format(variable))
        continue

    if all([isinstance(i, (str, bool)) for i in unique_vals]):
        print('discrete variable ', variable)
        
        if color_dict==None:
#             colors = plt.cm.get_cmap("tab20").colors
            cm = plt.get_cmap('gist_rainbow')
            colors=[cm(1.*i/len(unique_vals)) for i in range(len(unique_vals))]
        else:
            colors = [color_dict[x] for x in unique_vals]

        # plot data by unique value
        for g, c in zip(unique_vals, colors):
            c = [c]  # to remove the warnings
            axes.scatter(
                data_df.loc[annot.index[annot[variable] == g].tolist(), '0',],
                data_df.loc[annot.index[annot[variable] == g].tolist(), '1',],
                label=g,
                marker=".",
                c=c,
                alpha=alpha,
            )

            # centroids by mean
            axes.scatter(
                data_df.loc[annot.index[annot[variable] == g].tolist(), '0',].mean(),
                data_df.loc[annot.index[annot[variable] == g].tolist(), '1',].mean(),
                marker="s",
                alpha=0.5,
                label=str(g) + " centroid",
                c=c,
            )
            axes.text(
                data_df.loc[annot.index[annot[variable] == g].tolist(), '0',].mean(),
                data_df.loc[annot.index[annot[variable] == g].tolist(), '1',].mean(),
                s=g,
                horizontalalignment="center",
                verticalalignment="bottom",
                alpha=0.5,
            )

        # edit and position legend
        handles, labels = axes.get_legend_handles_labels()
        legend_idx = ["centroid" in label for label in labels]
        handles = list(compress(handles, legend_idx))
        labels = list(compress(labels, legend_idx))
        axes.legend(handles, labels, loc="center left", bbox_to_anchor=(1.05, 0.5))

    elif all([isinstance(i, (int, float, np.int, np.float)) for i in unique_vals]):
        print('continous variable ', variable)


        # check if enumerated ie numerical -> probably not needed anymore
#         if annot[variable].dtype != 'float64' and annot[variable].dtype != 'int64':
#             label_to_number = defaultdict(partial(next, count(1)))
#             annot[variable]=[label_to_number[label] for label in annot[variable]]

        cmap = sns.cubehelix_palette(as_cmap=True)
        points = axes.scatter(
            data_df.loc[:, '0',],
            data_df.loc[:, '1',],
            marker=".",
            c=annot[variable],
            s=50, 
            cmap=cmap,
            alpha=alpha,
        )
        fig.colorbar(points)

    else:
        print("variable type not-detected for {}".format(variable))
        continue
    
    # show and save figure
    x_postfix=''
    y_postfix=''

    if dimred=='PCA':
        explained_variance=list(pd.read_csv(os.path.join(results_dir,"PCA_{}_{}_explained_variance.csv".format(snakemake.wildcards["split"],snakemake.wildcards["step"])), index_col=0)['0'])
        x_postfix=' ({:.1f}%)'.format(100*explained_variance[0])
        y_postfix=' ({:.1f}%)'.format(100*explained_variance[1])

    axes.set_xlabel(dimred+" 1"+x_postfix)
    axes.set_ylabel(dimred+" 2 "+y_postfix)
    plt.title(dimred+" of " + label + " "+variable)
    plt.show()

    # save plot
    fig.savefig(
        fname=os.path.join(results_dir, dimred+"_"+label+"_"+variable+".svg"),
        format="svg",
        dpi=300,
        bbox_inches="tight",
    )