#!/bin/env python

#### generic dimensionality reduction function using UMAP ####

#### libraries

# general
import os
import pandas as pd
import numpy as np
from itertools import compress

# dimensionality reduction
import umap

#### configurations

# data=observations x features
data=pd.read_csv(snakemake.input[0], index_col=0).T

output=snakemake.output[0]

if not os.path.exists(snakemake.params['results_dir']):
    os.mkdir(snakemake.params['results_dir'])
    
# if 3 samples or less UMAP can not be performed
if data.shape[0]<4:
    from pathlib import Path
    Path(output).touch()
    import sys
    sys.exit()

#### Dimensionality reduction via unsupervised UMAP for visualization purpose
data_umap = umap.UMAP(
    n_components=2,
    random_state=42,
    metric="correlation",
).fit(data)

data_df = pd.DataFrame(data_umap.embedding_, index=data.index,)
data_df = data_df.rename_axis(("sample_name"))

data_df.to_csv(output)