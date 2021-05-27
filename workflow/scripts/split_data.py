#!/bin/env python

#### split data by metadata variable ####

#### libraries
import os
import pandas as pd
import numpy as np

#### configurations

# input
data_counts=pd.read_csv(snakemake.input[0], index_col=0)
annotation=pd.read_csv(snakemake.input[1], index_col=0)

# parameters
project_path = snakemake.params["project_path"]
split_by = snakemake.params["split_by"]

#####  make output folders, split and save data
for split in list(annotation[split_by].unique()):
    # make split folder
    if not os.path.exists(os.path.join(project_path,split)):
        os.mkdir(os.path.join(project_path,split))
    
    # split counts
    split_counts = data_counts.loc[:, annotation[(annotation[split_by]==split)].index]
    # save split counts
    split_counts.to_csv(os.path.join(project_path,split,"{}_counts.csv".format(split)))
    
    # split annotatios
    split_annot = annotation.loc[annotation[split_by]==split, :]
    # save split annotatios
    split_annot.to_csv(os.path.join(project_path,split,"{}_annotation.csv".format(split)))