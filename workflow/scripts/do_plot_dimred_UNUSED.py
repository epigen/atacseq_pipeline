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

# dimensionality reduction
import umap
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA

# utils for enumerating string lists
from collections import defaultdict
from itertools import count
from functools import partial


#### configurations

# annot=observations x annotations
annot=pd.read_csv(snakemake.input[0], index_col=0)

# data=observations x features
data=pd.read_csv(snakemake.input[1], index_col=0).T
data=data.set_index(data.columns[0])

# variables=columns of annot to plot
variables=snakemake.params['variables']

#  label=plot title and file name
label = snakemake.params["label"]

# results folder
results_dir = snakemake.params["results_dir"]

#plot parameters
color_dict=None
alpha=1

#### Dimensionality reduction via unsupervised UMAP for visualization purpose (init with PCA was horrible)
data_umap = umap.UMAP(
    n_components=2,
    random_state=42,
    metric="correlation",
).fit(data)
#     print(data_umap.embedding_.shape)

# dimensionality redcution via PCA so that 99.9% of variance is preserved
pca_obj = PCA(0.999, 
               random_state=42,
              )
data_pca=pca_obj.fit_transform(StandardScaler().fit_transform(data))

#     print(pca_obj.explained_variance_ratio_)
#     print(data_pca.shape)

for dimred in ['PCA','UMAP']:

    # make dataframes for plotting
    if dimred=='PCA':
        data_df = pd.DataFrame(data_pca, index=data.index,)
        data_df = data_df.rename_axis(("sample_name"))
    #     print(data_df.shape)
    #     print(data_df.head())
    if dimred=='UMAP':
        data_df = pd.DataFrame(data_umap.embedding_, index=data.index,)
        data_df = data_df.rename_axis(("sample_name"))
    #     print(data_df.shape)
    #     print(data_df.head())



    # plot 2D UMAP of data
    # sns.set(color_codes=True)
    for variable in variables:
        fig, axes = plt.subplots(nrows=1, ncols=1, figsize=(5, 5))  # figsize=(3, 3))
        unique_vals = annot[variable].unique()

        if len(unique_vals)<=20:
            print('discrete variable ', variable)
            if color_dict==None:
                colors = plt.cm.get_cmap("tab20").colors
            else:
                colors = [color_dict[x] for x in unique_vals]

            # plot data by unique value
            for g, c in zip(unique_vals, colors):
                c = [c]  # to remove the warnings
                axes.scatter(
                    data_df.loc[annot.index[annot[variable] == g].tolist(), 0,],
                    data_df.loc[annot.index[annot[variable] == g].tolist(), 1,],
                    label=g,
                    marker=".",
                    c=c,
                    alpha=alpha,
                )

                # centroids by mean
                axes.scatter(
                    data_df.loc[annot.index[annot[variable] == g].tolist(), 0,].mean(),
                    data_df.loc[annot.index[annot[variable] == g].tolist(), 1,].mean(),
                    marker="s",
                    alpha=0.5,
                    label=str(g) + " centroid",
                    c=c,
                )
                axes.text(
                    data_df.loc[annot.index[annot[variable] == g].tolist(), 0,].mean(),
                    data_df.loc[annot.index[annot[variable] == g].tolist(), 1,].mean(),
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

        else:
            print('continous variable ', variable)


            # check if enumerated ie numerical
            if annot[variable].dtype != 'float64' and annot[variable].dtype != 'int64':
                label_to_number = defaultdict(partial(next, count(1)))
                annot[variable]=[label_to_number[label] for label in annot[variable]]

            cmap = sns.cubehelix_palette(as_cmap=True)
            points = axes.scatter(
                data_df.loc[:, 0,],
                data_df.loc[:, 1,],
                marker=".",
                c=annot[variable],
                s=50, 
                cmap=cmap,
                alpha=alpha,
            )
            fig.colorbar(points)

        # show and save figure
        x_postfix=''
        y_postfix=''

        if dimred=='PCA':
            x_postfix=' ({:.1f}%)'.format(100*pca_obj.explained_variance_ratio_[0])
            y_postfix=' ({:.1f}%)'.format(100*pca_obj.explained_variance_ratio_[1])

        axes.set_xlabel(dimred+" 1"+x_postfix)
        axes.set_ylabel(dimred+" 2 "+y_postfix)
        plt.title(dimred+" of " + label + " "+variable)
        plt.show()

        # save plot if directory is provided
        if results_dir!=None:
            fig.savefig(
                fname=os.path.join(results_dir, dimred+"_"+label+"_"+variable+".svg"),
                format="svg",
                dpi=300,
                bbox_inches="tight",
            )