#!/bin/env python

#### generic dimensionality reduction function using PCA ####

#### libraries

# general
import os
import pandas as pd
import numpy as np
from itertools import compress

# dimensionality reduction
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA

#### configurations

# data=observations x features
data=pd.read_csv(snakemake.input[0], index_col=0).T

output_data=snakemake.output[0]
output_expl_var=snakemake.output[1]

if not os.path.exists(snakemake.params['results_dir']):
    os.mkdir(snakemake.params['results_dir'])

#### Dimensionality reduction via unsupervised PCA so that 99.9% of variance is preserved for visualization purpose
pca_obj = PCA(0.999, 
               random_state=42,
              )
data_pca=pca_obj.fit_transform(StandardScaler().fit_transform(data))

data_df = pd.DataFrame(data_pca, index=data.index,)
data_df = data_df.rename_axis(("sample_name"))

data_df.to_csv(output_data)

pd.DataFrame(pca_obj.explained_variance_ratio_).to_csv(output_expl_var)