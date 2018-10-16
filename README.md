
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
Next you need an API key for Google Maps. Earlier this year GoogleMaps changed the way they price API calls – see here (https://cloud.google.com/maps-platform/pricing/#details). Under the new system there is a lower limit on API calls below which you won’t be charged, but annoyingly you still need to sign up for an API key otherwise it produces a 403 error.

Instructions on signing up for a GoogleMaps API key can be found here (https://developers.google.com/maps/documentation/javascript/get-api-key). Then – just to be extra sure that you will never be charged - you can limit your budget to something very low (like £1) by following these instructions (https://cloud.google.com/billing/docs/how-to/budgets?hl=en). This way you can be sure you will only use the free API calls, and your account will cap before you start being charged any significant amount of money.

Once you have a key (which should look like a long string of letters and numbers) you need to register it in R as below, replacing 'xxxxxxx' with your key. Start by installing and loading the “ggmap” package. Note that this should be installed from GitHub since the CRAN version doesn’t currently include the function register_google(). After loading the ggmap library, load RgeoProfile as normal. Once your key is registered mapping functions in RgeoProfile should work normally.
```r
install_github("dkahle/ggmap")
library(ggmap)
register_google(key ="xxxxxxxx")
```
More details can be found here: (https://rdrr.io/github/fresques/ggmap/man/register_google.html

We know this is inconvenient and makes it more difficult for users to work with our package. For this reason the next version of RgeoProfile is likely to switch to open source alternatives.

Finally, install the RgeoProfile package directly from GitHub by running
```r
install_github("bobverity/RgeoProfile")
library(RgeoProfile)
```
Once installed, run `??RgeoProfile` to open the help for the package, which contains a number of worked examples.

## References

<sup>1</sup> Rossmo, D.K., 1999. Geographic profiling. CRC press.<br/>
The theory of geographic profiling, including the Criminal Geographic Targeting (CGT) algorithm used in criminology.

<sup>2</sup> Verity R, Stevenson MD, Rossmo DK, Nichols RA, Le Comber SC (2014). Spatial targeting of infectious disease control: Identifying multiple, unknown sources. Methods in Ecology and Evolution vol. 5 (7), 647-655. DOI: 10.1111/2041-210X.12190 <br/>
Introducing the Dirichlet Process Mixture (DPM) model of geographic profiling.

<sup>3</sup> Faulkner SC, Verity R, Roberts D, Roy SS, Robertson PA, Stevenson MD, Le Comber SC (2016). Using geographic profiling to compare the value of sightings vs trap data in a biological invasion. Diversity and Distributions 23 (1), 104-112. DOI: 10.1111/ddi.12498 <br/>
Using citizen science data with geographic profiling, and extending the DPM model to fit sigma, the dispersal parameter.

<sup>4</sup> Faulkner SC, Stevens MCA, Romañach SS, Lindsey PA, Le Comber SC (2018). A spatial approach to combatting wildlife crime. Conservation Biology 32 (3), 685-693. DOI: 10.1111/cobi.13027 <br/>
A study of poaching in Zimbabwe showing how geospatial information can be incorporated within the GP framework.

<sup>5</sup> Struebig MJ, Linkie M, Deere N, Martyr DJ, Millyanawati B, Faulkner SC, Le Comber SC, Mangunjaya FM, Fachruddin M, Leader-Williams N, McKay JE, St John FAV (2018). Addressing human-tiger conflict using socio-ecological information on tolerance and risk. Nature Communications 9 (1), article 3455. DOI: 10.1038/s41467-018-05983-y <br/>
A recent paper using geographic profiling to improve predictions of human-tiger conflict.

<sup>6</sup> Rossmo DK, Lutermann H, Stevenson MD, Le Comber SC (2014). Geographic profiling in Nazi Berlin: fact and fiction. Geospatial Intelligence Review (Fall 2014). <br/>
Applying geographic profiling to a case of German resistance to the Nazis which formed the basis of the classic novel, ‘Alone in Berlin’.
