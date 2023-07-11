#!/bin/env python

#### integrate data using reComBat ####

#### libraries
import os
import pandas as pd
from reComBat import reComBat

#### configurations

# input
data_path = os.path.join(snakemake.input["normalized_data"])
annotation_path = os.path.join(snakemake.input["annotation"])

# ouput
result_path = os.path.join(snakemake.output["integrated_data"])

# parameters
recombat_parameters = snakemake.config["reComBat_parameters"]

# load data
data = pd.read_csv(data_path, index_col=0).T
annotation = pd.read_csv(annotation_path, index_col=0)

# split and adapt metadata
# batch as pandas series
batch = annotation[recombat_parameters["batch_column"]]
# desired effects as pandas data frame
X = annotation.loc[:,recombat_parameters["desired_categorical"]+recombat_parameters["desired_numerical"]]
X.columns = [col+"_numerical" if col in recombat_parameters["desired_numerical"] else col for col in X.columns]
# unwanted effects as pandas data frame
C = annotation.loc[:,recombat_parameters["unwanted_categorical"]+recombat_parameters["unwanted_numerical"]]
C.columns = [col+"_numerical" if col in recombat_parameters["unwanted_numerical"] else col for col in C.columns]

#### integrate data using reComBat

# create reComBat model
combat_model = reComBat(parametric=recombat_parameters["parametric"],
                        model=recombat_parameters["model"],
                        config=dict(recombat_parameters["config"]) if recombat_parameters["config"]!="" else None,
                        conv_criterion=recombat_parameters["conv_criterion"],
                        max_iter=recombat_parameters["max_iter"],
                        n_jobs=None,
                        mean_only=recombat_parameters["mean_only"],
                        optimize_params=recombat_parameters["optimize_params"],
                        reference_batch=recombat_parameters["reference_batch"] if recombat_parameters["reference_batch"]!="" else None,
                        verbose=recombat_parameters["verbose"]
                       )

# fit and transform data
integrated_data = combat_model.fit_transform(data = data,
                                             batches = batch,
                                             X = X if X.shape[1]>0 else None,
                                             C = C if C.shape[1]>0 else None
                                            )


# save integrated data
integrated_data.T.to_csv(result_path)
