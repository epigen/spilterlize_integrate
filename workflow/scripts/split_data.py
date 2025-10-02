#!/bin/env python

#### split data by annotation ####

#### libraries
import os
import pandas as pd

#### helper functions
def drop_all_zero(data: pd.DataFrame, annotation: pd.DataFrame):
    
    # align data with annotion
    d = data.loc[:, annotation.index]
    # drop all-zero columns (samples)
    col_mask = (d != 0).any(axis=0)
    d = d.loc[:, col_mask]
    a = annotation.loc[d.columns]
    # Drop all-zero rows (features)
    row_mask = (d != 0).any(axis=1)
    d = d.loc[row_mask, :]

    return d,a

#### configurations

# input
data_path = snakemake.input["data"]
annotation_path = snakemake.input["annotation"]

# parameters
result_path = snakemake.params["result_path"]
split_by = snakemake.params["split_by"]


# load data
data = pd.read_csv(data_path, index_col=0)
annotation = pd.read_csv(annotation_path, index_col=0)

# split and save count and annotation data
for split_var in split_by:
    if split_var == "all":
        # drop all zero columns (samples) and rows (features)
        d, a = drop_all_zero(data, annotation)
        
        # save cleaned split
        d.to_csv(os.path.join(result_path,"all","counts.csv"))
        a.to_csv(os.path.join(result_path,"all","annotation.csv"))
        
    
    if split_var in annotation.columns:
        for split in annotation[split_var].unique():

            # split annotation
            annotation_split = annotation.loc[annotation[split_var]==split, :]
            
            # drop all zero columns (samples) and rows (features)
            d, a = drop_all_zero(data, annotation_split)

            # save cleaned split
            d.to_csv(os.path.join(result_path,f"{split_var}_{split}","counts.csv"))
            a.to_csv(os.path.join(result_path,f"{split_var}_{split}","annotation.csv"))
            