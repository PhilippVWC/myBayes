# Package myBayes #

## Overview ##
Package *myBayes* is a collection of probabilistic methods for bayesian inference on dynamical systems with a strong focus on the logistic map and nonparametric, choatic maps.
It is essentially a toolbox for parameter estimation and time series prediction. The developer does not make himself responsable for potential errors, code stabililty nor mathematical correctness of the underlying logic.

Philipp van Wickevoort Crommelin

## Contents ##
Among the most important functions available in *myBayes* are:

Function | Purpose
---------|----------
MH | Simmulated annealing for optimization in one dimension.

## Special remarks ##
Many routines are implemented in C++ in order to increase performance. 

## Dependencies ##
C++ integration is accomplished by Dirk Eddelbuettels ***Rcpp*** package (For more information see [Rcpp](http://rcpp.org/)).
