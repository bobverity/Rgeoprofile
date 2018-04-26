#' RgeoProfile: Dirichlet Process Mixture (DPM) model of geographic profiling
#'
#' Carries out a DPM mixture model of geographic profiling (and associated functions) as described in Verity et al. (2014) For links to papers see http:/www.sbcs.qmul.ac.uk/staff/stevenlecomber.html or contact Steven Le Comber (s.c.lecomber@qmul.ac.uk). Created and maintained by Bob Verity (r.verity@imperial.ac.uk).  
#'
#' @name RgeoProfile
#' @docType package
#' @details
#' RgeoProfile: Dirichlet Process Mixture (DPM) model of geographic profiling
#'
#' Carries out a DPM mixture model of geographic profiling (and associated functions) as described in Verity et al. (2014) For links to papers see http:/www.sbcs.qmul.ac.uk/staff/stevenlecomber.html or contact Steven Le Comber (s.c.lecomber@qmul.ac.uk). Created and maintained by Bob Verity (r.verity@imperial.ac.uk).  
#'
#' @name RgeoProfile
#' @docType package
#' @examples
#' \donttest{
#' # full example of Rgeoprofile 2.1.0 workflow, illustrating all functions
#' # for details, see help for individual functions
#' 
#' #------------------------------------------------------------------
#' # data and settings
#' #------------------------------------------------------------------
#' # example data
#' d <- LondonExample_crimes
#' s <- LondonExample_sources
#' 
#' # convert d and s to correct format for geoParams()
#' # (note that in this case the example data are already in the correct
#' # format; these steps are only relevant if for example d and s are 
#' # imported as two-column matrices. They are included here for
#' # completeness)
#' d <- geoData(d$longitude, d$latitude)
#' s <- geoDataSource(s$longitude, s$latitude)
#' 
#' # set model and MCMC parameters
#' p = geoParams(data = d, sigma_mean = 1, sigma_squared_shape = 2, chains = 5, 
#'                 burnin = 1e3, samples = 1e4)
#' 
#' #------------------------------------------------------------------
#' # run model
#' #------------------------------------------------------------------
#' # run MCMC
#' m = geoMCMC(data = d, params = p)
#' 
#' #------------------------------------------------------------------
#' # output
#' #------------------------------------------------------------------
#' # plot prior and posterior of sigma
#' geoPlotSigma(params = p, mcmc = m)
#' 
#' # plot profile on map
#' mapGP <- geoPlotMap(params = p, data = d, source = s, surface = m$geoProfile)
#' mapGP
#' 
#' # get hitscores
#' hs <- geoReportHitscores(params = p, source = s, surface =m$geoProfile)
#' hs
#' 
#' # produce Lorenz plot
#' Gini <- geoPlotLorenz(hit_scores = hs, crimeNumbers = NULL)
#' Gini
#' 
#' # zoom
#' zoomLon = c(-0.1, -0.01)
#' zoomLat = c(51.51, 51.54)
#' mapZoom <- geoPlotMap(lonLimits = zoomLon, latLimits = zoomLat, params = p, 
#'                 data = d, source = s, surface = m$geoProfile)
#' mapZoom
#' 
#' # plot allocation
#' geoPlotAllocation(mcmc = m)
#' 
#' # plot co-allocation
#' geoPlotCoallocation(mcmc = m)
#' 
#' # produce perspective plots
#' # probabilities
#' geoPersp(surface = m$posteriorSurface, aggregate_size = 3, surface_type = "prob")
#' # ranked surface
#' geoPersp(surface = m$geoProfile, aggregate_size = 3)
#' 
#' # find centroids of data split by best grouping (placeholder for more thorough method)
#' ms <- geoModelSources(mcmc = m, data = d)
#' ms
#' # add peaks to map
#' # NB requires ggplot2
#' library(ggplot2)
#' mapSource <- mapGP + geom_point(aes(ms$longitude,ms$latitude), size=6, pch = 3, col="red")
#' mapSource
#' 
#' # plot surface in style of 'unknown pleasures', for fun
#' unknownPleasures(m$geoProfile, paper_ref = "RgeoProfile 2.1.0")
#' 
#' #------------------------------------------------------------------
#' # compare to alternative ring search strategy
#' #------------------------------------------------------------------
#' # compare to geoprofile based on ring search strategy
#' surface_ring <- geoRing(params = p, data = d, source = s, mcmc = m)
#' gp_ring <- geoProfile(surface = surface_ring)
#' # map of ring search geoprofile
#' mapRing <- geoPlotMap(params = p, data = d, source = s, surface = gp_ring, 
#'                 surfaceCols <- c("red", "white"))
#' mapRing
#' 
#' # hitscores of ring search geoprofile
#' hs_ring <- geoReportHitscores(params = p, source = s, surface = gp_ring)
#' hs_ring
#' 
#' #------------------------------------------------------------------
#' # incorporate GIS data
#' #------------------------------------------------------------------
#' # read in north London shapefile as mask and adjust surface
#' north_london_mask <- geoShapefile()
#' 
#' # restrict mask to Tower Hamlets
#' TH_mask <- north_london_mask[which(north_london_mask$NAME == "Tower Hamlets"),]
#' prob_masked <- geoMask(probSurface = m$posteriorSurface, params = p, mask = TH_mask,
#'                 operation = "outside", scaleValue = 1e-9)
#' gp_masked <- geoProfile(prob_masked$prob)
#' 
#' # plot new surface
#' mapMask <- geoPlotMap(params = p, data = d, source = s, surface = gp_masked, 
#'                 breakPercent = seq(0,25,l = 11))
#' mapMask
#' 
#' # hs of masked surface
#' hs_mask <- geoReportHitscores(params = p, source = s, surface = gp_masked)
#' hs_mask
#' 
#' #------------------------------------------------------------------
#' }

NULL