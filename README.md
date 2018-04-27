
# RgeoProfile
### Version2.1.0
[![Build Status](https://travis-ci.org/bobverity/Rgeoprofile.svg?branch=master)](https://travis-ci.org/bobverity/Rgeoprofile)
[![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/bobverity/RgeoProfile?branch=master&svg=true)](https://ci.appveyor.com/project/bobverity/RgeoProfile)

## Background
RgeoProfile is an R package for carrying out geographic profiling - a technique derived from criminology that uses the spatial locations of linked crimes to infer the home location (or locations) from which the criminal is operating<sup>1</sup>. More recently, this method has been applied to problems in ecology<sup>2-5</sup> and infectious disease<sup>6,7</sup>, where "crimes" are now equivalent to animal observations or sites of infection, and "home locations" may be roosts or persistent sources of disease. Finding these sources is often a key priority - for example in malaria this may provide a more efficient method of carrying out targeted larval source management<sup>7</sup>.

The RgeoProfile package uses the Dirichlet Process Mixture (DPM) modelling framework of *Verity et al.*<sup>7</sup>, which was designed to place geographic profiling in a Bayesian framework and to deal with the issue of multiple sources. In simple terms, the DPM model assumes a large (strictly infinite) number of *potential* source locations, although some of these sources are assumed to be more "active" than others. The RgeoProfile algorithm then attempts to infer the number of active sources (i.e. the number of sources that gave rise to the observed data), along with the spatial locations of these sources. It does this using Markov chain Monte Carlo (MCMC) - a statistical technique that provides estimates of unknown parameters by iteratively updating our guess.

<p align="center">
<img src="R_ignore/LondonExample_figure1.png" width="700" align="middle">

<p align="center"> Example geoprofile produced by analysing the tutorial data set, included with the package. </p>
</p>

RgeoProfile is written in R, and links to C++ through the Rcpp package. The required data is nothing more than the longitude/latitude of observed "crimes". The user is also required to define a series of parameters, for example defining prior beliefs about the average distance of a crime from a source. After running the main MCMC algorithm, a large number of plotting functions are available for exploring the output and, for evaluating the efficiency of the search strategy.


## Installation

RgeoProfile relies on the Rcpp package, which requires the following OS-specific steps:

* Windows
    - Download and install the appropriate version of [Rtools](https://cran.rstudio.com/bin/windows/Rtools/) for your version of R. On installation, ensure you check the box to arrange your system PATH as recommended by Rtools
* Mac OS X
    - Download and install [XCode](http://itunes.apple.com/us/app/xcode/id497799835?mt=12)
    - Within XCode go to Preferences : Downloads and install the Command Line Tools
* Linux (Debian/Ubuntu)
    - Install the core software development utilities required for R package development as well as LaTeX by executing
    
            ```
            sudo apt-get install r-base-dev texlive-full
            ```

Next, ensure that you have devtools installed by running
```r
install.packages("devtools")
library(devtools)
```
Finally, install the RgeoProfile package directly from GitHub by running
```r
install_github("bobverity/RgeoProfile")
library(RgeoProfile)
```
Once installed, run `??RgeoProfile` to open the help for the package, which contains a number of worked examples.


## References

<sup>1</sup> Rossmo, D.K., 1999. Geographic profiling. CRC press.

<sup>2</sup> Le Comber, S.C., Nicholls, B., Rossmo, D.K. and Racey, P.A., 2006. Geographic profiling and animal foraging. Journal of Theoretical Biology, 240(2), pp.233-240.

<sup>3</sup> Martin, R.A., Rossmo, D.K. and Hammerschlag, N., 2009. Hunting patterns and geographic profiling of white shark predation. Journal of Zoology, 279(2), pp.111-118.

<sup>4</sup> Raine, N.E., Rossmo, D.K. and Le Comber, S.C., 2009. Geographic profiling applied to testing models of bumble-bee foraging. Journal of the Royal Society Interface, 6(32), pp.307-319.

<sup>5</sup> Le Comber, S.C. and Stevenson, M.D., 2012. From Jack the Ripper to epidemiology and ecology. Trends in ecology & evolution, 27(6), pp.307-308.

<sup>6</sup>Le Comber, S.C., Rossmo, D., Hassan, A.N., Fuller, D.O. and Beier, J.C., 2011. Geographic profiling as a novel spatial tool for targeting infectious disease control. International journal of health geographics, 10(1), p.35.

<sup>7</sup> Verity, R., Stevenson, M.D., Rossmo, D.K., Nichols, R.A. and Le Comber, S.C., 2014. Spatial targeting of infectious disease control: identifying multiple, unknown sources. Methods in Ecology and Evolution, 5(7), pp.647-655.
