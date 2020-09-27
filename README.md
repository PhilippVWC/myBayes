# Package myBayes #

## Overview ##
*myBayes* is an R package that serves as a collection of probabilistic methods for bayesian inference on dynamical systems with a strong focus on parametric and nonparametric choatic maps.
It is essentially a toolbox for parameter estimation, time series filter and prediction. The developer does not make himself responsable for potential errors, code stabililty nor mathematical correctness of the underlying logic.

## Contents ##
Among the most important functions available in *myBayes* are:

Function | Purpose
---------|----------
MH | Simmulated annealing for optimization in one dimension.
Lik_gsm_cpp | A Gaussian Likelihood given time series data of the general symmetric map.
blurrVector | For statistical dispersion of deterministic vectors.
CeBePa | For approximate parameter transformation from Gaussian distribution to Beta distribution
collapse | Routine to remove equal and adjacent values of a vector.
myOpt_N | Optimization in more than two dimensions with multiple algorithms (sim. Annealing, BFGS, etc.)
axLabels3d | Individual manipulation of axis label annotations for 3d persp plots. 

## Special remarks ##
Many routines are implemented in C++ in order to increase performance.

## Dependencies ##
C++ integration is accomplished by Dirk Eddelbuettels **Rcpp** package (For more information see [Rcpp](http://rcpp.org/)). An installation of Rcpp is crucial for proper functionality of *myBayes*.

## Installation ##
Install package *devtools* and install from github repository:
```R
install.packages("devtools")
devtools::install_github("PhilippVWC/myBayes")
```


## Contact ##
For further remarks, please contact me under (philippcrommelin@wwu.de), Philipp van Wickevoort Crommelin.
