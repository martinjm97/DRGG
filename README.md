# Directed Random Geometric Graph (DRGG) Model Evaluation
## Overview
We provide the data and implementations that were used to evaluate the practical applications of the DRGG model. Our example used is from a word ascociation study. We also found other data sets such as http://vlado.fmf.uni-lj.si/pub/networks/data/dic/fa/FreeAssoc.htm.

## To Use
Download PairsP.net from http://vlado.fmf.uni-lj.si/pub/networks/data/dic/fa/FreeAssoc.htm and add it to the cloned repo. One may start by running theoretical_degree_distribution.py to reproduce some of the graphs given in the paper. 

## Contents
- drgg.py a class that allows for the production of graphs following the DRGG model.
- real_world_data_metrics.py computes the different metrics (e.g. degree distribution) for data from a given data set.
- drgg_model_metrics.py computes the different metrics (e.g. degree distribution) for the DRGG simulation of a given data set.
- theoretical_degree_distribution.py the implementation of the expected behavior of DRGG that uses the theoretical results of the different properties (e.g. degree distribution) of the graph produced.
- utils.py some helper methods
