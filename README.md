
# RgeoProfile
[![Build Status](https://travis-ci.org/bobverity/Rgeoprofile.svg?branch=master)](https://travis-ci.org/bobverity/Rgeoprofile)
[![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/bobverity/RgeoProfile?branch=master&svg=true)](https://ci.appveyor.com/project/bobverity/RgeoProfile)

## Background
RgeoProfile is an R package for carrying out geographic profiling - a technique derived from criminology that uses the spatial locations of linked crimes to infer the home location (or locations) from which the criminal is operating. More recently, this method has been applied to problems in ecology and infectious disease, where "crimes" are now equivalent to animal observations or sites of infection, and "home locations" may be roosts or persistant sources of disease. Finding these sources is often a priority, for example in malaria this may provide a more efficient method of carrying out targeted larval source management.

The RgeoProfile package uses the Dirichlet Process Mixture (DPM) modelling framework of Verity et al. (REF), which was designed to place geographic profiling in a Bayesian framework and to deal with the issue of multiple sources. In simple terms, the DPM model assumes a large (strictly infinite) number of potential source locations, although some of these sources are assumed to be more "active" than others. The RgeoProfile algorithm then attempts to infer the number of active sources (i.e. the sources that gave rise to the observed data, and not all the other potential sources), along with the spatial locations of these sources. It does this by Bayesian Markov chain Monte Carlo (MCMC) - a statistical technique in which the parameters of interest are estimated by iteratively updating our guess, given the known values of all other parameters.

RgeoProfile is written in R and C++ (through the Rcpp package). The required data is nothing more than the longitude/latitude locations of observed "crimes". Other user inputs include a series of parameters for defining priors on things like the average distance of a crime from a source. A large number of plotting functions are available for exploring the MCMC output, and for evaluating the efficiency of the search strategy of the final profile compared with other methods.

## Installation

RgeoProfile uses the Rcpp package, which requires a couple of additional installation steps:

### Windows
- test bullet point
- test bullet point

Some code:
```r
install.packages("rmarkdown")
```
### Mac OS X


