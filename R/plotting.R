
#------------------------------------------------
# ggplot theme that gets rid of all margins etc.
# credit to the Wilke lab and cowplot package for original code that this is taken from (https://github.com/wilkelab/cowplot)
# (not exported)

theme_nothing <- function(font_size = 14, font_family = ""){
  theme_void(base_size = font_size, base_family = font_family) %+replace%
    theme(
      line = element_blank(),
      rect = element_blank(),
      text = element_text(
        family = font_family, face = "plain",
        colour = "black", size = font_size,
        lineheight = 0.9, hjust = 0.5, vjust = 0.5, angle = 0,
        margin = margin(), debug = FALSE
      ),
      
      axis.line =          element_blank(),
      axis.line.x =        NULL,
      axis.line.y =        NULL,
      axis.text =          element_blank(),
      axis.text.x =        element_blank(),
      axis.text.x.top =    element_blank(),
      axis.text.y =        element_blank(),
      axis.text.y.right =  element_blank(),
      axis.ticks =         element_blank(),
      axis.ticks.length =  unit(0, "pt"),
      axis.title.x =       element_blank(),
      axis.title.x.top =   element_blank(),
      axis.title.y =       element_blank(),
      axis.title.y.right = element_blank(),
      
      legend.background =  element_blank(),
      legend.spacing =     unit(0.4, "cm"),
      legend.spacing.x =   NULL,
      legend.spacing.y =   NULL,
      legend.margin =      margin(0.2, 0.2, 0.2, 0.2, "cm"),
      legend.key =         element_blank(),
      legend.key.size =    unit(1.2, "lines"),
      legend.key.height =  NULL,
      legend.key.width =   NULL,
      legend.text =        element_text(size = rel(0.8)),
      legend.text.align =  NULL,
      legend.title =       element_text(hjust = 0),
      legend.title.align = NULL,
      legend.position =    "none",
      legend.direction =   NULL,
      legend.justification = "center",
      legend.box =         NULL,
      legend.box.margin =  margin(0, 0, 0, 0, "cm"),
      legend.box.background = element_blank(),
      legend.box.spacing = unit(0.4, "cm"),
      
      panel.background =   element_blank(),
      panel.border =       element_blank(),
      panel.grid.major =   element_blank(),
      panel.grid.minor =   element_blank(),
      panel.spacing =      unit(0, "pt"),
      panel.spacing.x =    NULL,
      panel.spacing.y =    NULL,
      panel.ontop    =     FALSE,
      
      strip.background =   element_blank(),
      strip.text =         element_blank(),
      strip.text.x =       element_blank(),
      strip.text.y =       element_blank(),
      strip.placement =    "inside",
      strip.placement.x =  NULL,
      strip.placement.y =  NULL,
      strip.switch.pad.grid = unit(0., "cm"),
      strip.switch.pad.wrap = unit(0., "cm"),
      
      plot.background =    element_blank(),
      plot.title =         element_blank(),
      plot.subtitle =      element_blank(),
      plot.caption =       element_blank(),
      plot.margin =        margin(0, 0, 0, 0),
      
      complete = TRUE
    )
}

#------------------------------------------------
#' Plot prior and posterior distributions of sigma.
#'
#' Plot prior distribution of sigma as defined by current parameter values. Can optionally overlay a kernel density plot of posterior draws of sigma.
#'
#' @param params a list of parameters as defined by geoParams().
#' @param mcmc stored output obtained by running geoMCMC(). Leave as NULL to plot prior only.
#' @param plotMax maximum x-axis range to plot. Leave as NULL to use default settings.
#'
#' @export
#' @examples
#' # John Snow cholera data
#' d <- geoData(Cholera$longitude, Cholera$latitude)
#' p <- geoParams(data = d, sigma_mean = 0.2, sigma_squared_shape = 2)
#' m <- geoMCMC(data = d, params = p)
#' geoPlotSigma(params = p, mcmc = m)
#' 
#' # simulated data
#' sim <-rDPM(50, priorMean_longitude = -0.04217491, priorMean_latitude = 
#' 51.5235505, alpha=1, sigma=1, tau=3)
#' d <- geoData(sim$longitude, sim$latitude)
#' p <- geoParams(data = d, sigma_mean = 1.0, sigma_squared_shape = 2)
#' m <- geoMCMC(data = d, params = p)
#' geoPlotSigma(params = p, mcmc = m)

geoPlotSigma <- function(params, mcmc=NULL, plotMax=NULL) {
  
  # check params
  geoParamsCheck(params, silent=TRUE)
  
  # check that plotMax is sensible
  if (!is.null(plotMax)) {
    stopifnot(is.numeric(plotMax))
    stopifnot(is.finite(plotMax))
    stopifnot(plotMax>0)
  }
  
  # extract sigma parameters
  sigma_mean <- params$model$sigma_mean
  sigma_var <- params$model$sigma_var
  alpha <- params$model$sigma_squared_shape
  beta <- params$model$sigma_squared_rate
  
  # stop if using fixed sigma model
  if (sigma_var==0) { stop('can only produce this plot under variable-sigma model (i.e. sigma_var>0)') }
  
  # extract sigma draws from mcmc object
  sigma_draws <- mcmc$sigma
  
  # default plotMax based on extent of prior distribution AND the extent of posterior draws if available
  if (is.null(plotMax)) {
    plotMax <- sigma_mean + 3*sqrt(sigma_var)
    if (!is.null(sigma_draws)) {
      plotMax <- max(plotMax, 2*max(sigma_draws, na.rm=TRUE))
    }
  }
  
  # produce prior distribution
  sigma_vec <- seq(0, plotMax, l=501)
  sigma_prior <- dRIG(sigma_vec, alpha, beta)
  
  # plot prior and overlay density of posterior draws if used
  if (is.null(sigma_draws)) {
    
    plot(sigma_vec, sigma_prior, type='l', lty=1, xlab='sigma (km)', ylab='probability density', main='')
    legend(x='topright', legend='prior', lty=1)
    
  } else {
    
    sigma_posterior <- density(sigma_draws, from=0, to=plotMax)
    y_max <- max(sigma_posterior$y, na.rm=TRUE)
    
    plot(sigma_vec, sigma_prior, type='l', lty=1, ylim=c(0,y_max), xlab='sigma (km)', ylab='probability density', main='')
    lines(sigma_posterior, lty=2)
    legend(x='topright', legend=c('prior','posterior'), lty=c(1,2))
  }
  
}

#------------------------------------------------
#' Plot posterior allocation
#'
#' Produces plot of posterior allocation from output of MCMC.
#'
#' @param mcmc stored output obtained by running geoMCMC().
#' @param colours vector of colours for each allocation. If NULL then use default colour scheme.
#' @param barBorderCol colour of borders around each bar. Set as NA to omit this border (useful when there are a large number of observations).
#' @param barBorderWidth line width of borders around each bar.
#' @param mainBorderCol colour of border around plot.
#' @param mainBorderWidth line width of border around plot.
#' @param yTicks_on whether to include ticks on the y-axis.
#' @param yTicks vector of y-axis tick positions.
#' @param xTicks_on whether to include ticks on the x-axis.
#' @param xTicks_size size of ticks on the x-axis.
#' @param xlab x-axis label.
#' @param ylab x-axis label.
#' @param mainTitle main title over plot.
#' @param names individual names of each observation, written horizontally below each bar.
#' @param names_size size of names under each bar.
#' @param orderBy whether to order segments within each bar by "group" or by "probability". If ordered by group, all segments of a particular group are laid down before moving to the next group. If ordered by probability the segments within each bar are ordered from large to small.
#'
#' @export
#' @examples
#' # London example data
#' d <- LondonExample_crimes
#' p <- geoParams(data = d, sigma_mean = 1.0, sigma_squared_shape = 2)
#' m <- geoMCMC(data = d, params = p)
#' geoPlotAllocation(m)
#'
#' # John Snow cholera data
#' d <- Cholera
#' p <- geoParams(data = d, sigma_mean = 1.0, sigma_squared_shape = 2)
#' m <- geoMCMC(data = d, params = p, lambda=0.05)
#' geoPlotAllocation(m, barBorderCol=NA)	# (should allocate all to a single source!)
#' 
#' # simulated data
#' sim <-rDPM(50, priorMean_longitude = -0.04217491, priorMean_latitude = 
#'                51.5235505, alpha=1, sigma=1, tau=3)
#' d <- geoData(sim$longitude, sim $latitude)
#' p <- geoParams(data = d, sigma_mean = 1.0, sigma_squared_shape = 2)
#' m <- geoMCMC(data = d, params = p)
#' geoPlotAllocation(m)

geoPlotAllocation <- function(mcmc, colours=NULL, barBorderCol=NA, barBorderWidth=0.25, mainBorderCol="black", mainBorderWidth=2, yTicks_on=TRUE, yTicks=seq(0,1,0.2), xTicks_on=FALSE, xTicks_size=1, xlab="", ylab="posterior allocation", mainTitle="", names=NA, names_size=1, orderBy="group") {
  
  # check that orderBy is either 'group' or 'probability'
  if (!(orderBy %in% c("group","probability"))) {
    stop("orderBy must equal 'group' or 'probability'")
  }
  
  # get allocation from mcmc object
  allocation <- mcmc$allocation
  
  # check that allocation is a data frame
  if (!is.data.frame(allocation)) {
    stop("allocation must be a data frame, with observations in rows and groups in columns")
  }
  
  # extract useful parameters
  n <- nrow(allocation)
  k <- ncol(allocation)
  
  # replace colours if default
  if (is.null(colours)) {
    rawCols <- RColorBrewer::brewer.pal(n=11,name="RdYlBu")
    myPal <- colorRampPalette(rawCols)
    colours <- myPal(k)
    colours <- colours[c(2*(1:ceiling(k/2))-1,2*(1:floor(k/2)))]
  }
  
  # if ordered by group
  if (orderBy=="group") {
    
    plot(0, type='n', xlim=c(0,n), ylim=c(0,1), xlab=xlab, ylab=ylab, xaxs="i", yaxs="i",axes=FALSE, main=mainTitle)
    barplot(t(allocation), col=colours, space=0, border=NA, axes=FALSE, add=TRUE)
    segments(0:n,c(0,0),0:n,c(1,1), col=barBorderCol, lwd=barBorderWidth)
    box(col=mainBorderCol, lwd=mainBorderWidth)
    axis(2, tick=yTicks_on, labels=yTicks_on, at=yTicks)
    axis(1, at=1:n-0.5, tick=xTicks_on, lwd.ticks=xTicks_size, labels=names, las=2, cex.axis=names_size)
  }
  
  # if ordered by probability
  if (orderBy=="probability") {
    
    plot(0, type='n', xlim=c(0,n), ylim=c(0,1), xlab=xlab, ylab=ylab, xaxs="i", yaxs="i",axes=FALSE, main=mainTitle)
    tM <- t(allocation)
    for (i in 1:n) {
      temp <- tM
      temp[,-i] <- NA
      temp_order <- order(temp[,i],decreasing=TRUE)
      barplot(temp[temp_order,], col=colours[temp_order], space=0, border=NA, axes=FALSE, add=TRUE)
    }
    segments(0:n,c(0,0),0:n,c(1,1), col=barBorderCol, lwd=barBorderWidth)
    box(col=mainBorderCol, lwd=mainBorderWidth)
    axis(2, tick=yTicks_on, labels=yTicks_on, at=yTicks)
    axis(1, at=1:n-0.5, tick=xTicks_on, lwd.ticks=xTicks_size, labels=names, las=2, cex.axis=names_size)
  }
  
}

#------------------------------------------------
# get optimal zoom level given x and y values
# (not exported)

getZoom <- function(x,y) {
  
  # calculate midpoint of range in x and y
  xmid <- min(x) + diff(range(x))/2
  ymid <- min(y) + diff(range(y))/2
  
  # calculate delta (half of longitude range) for a range of zoom levels
  z <- 20:2
  delta <- 445/(2^z)
  
  # calculate left and right longitude limits at all zoom levels
  long_angle_left <- xmid - delta
  long_angle_right <- xmid + delta
  
  # calculate top and bottom latitude limits at all zoom levels
  lat_angle <- ymid/360*2*pi
  projection_mid <- log(tan(pi/4+lat_angle/2))
  projection_top <- projection_mid + delta/360*2*pi
  projection_bot <- projection_mid - delta/360*2*pi
  lat_angle_top <- (2*atan(exp(projection_top))-pi/2)*360/(2*pi)
  lat_angle_bot <- (2*atan(exp(projection_bot))-pi/2)*360/(2*pi)
  
  # find the most zoomed-in level that captures all points
  zoomTest <- (min(x)>long_angle_left) & (max(x)<long_angle_right) & (min(y)>lat_angle_bot) & (max(y)<lat_angle_top)
  if (!any(zoomTest)) {
    stop("values are outside of plotting range of earth")
  }
  bestZoom <- z[which(zoomTest)[1]]
  
  return(bestZoom)
}

#------------------------------------------------
#' Plot a map and overlay data and/or geoprofile
#'
#' Plots geoprofile on map, with various customisable options.
#'
#' @param params parameters list in the format defined by geoParams().
#' @param data data object in the format defined by geoData().
#' @param source potential sources object in the format defined by geoDataSource().
#' @param surface a surface to overlay onto the map, typically a geoprofile obtained from the output of geoMCMC().
#' @param surfaceCols vector of two or more colours to plot surface. Defaults to viridis palette.
#' @param zoom zoom level of map. If NULL then choose optimal zoom from params.
#' @param latLimits optional vector setting min and max latitude for zoom view.
#' @param lonLimits optional vector setting min and max longitude for zoom view.
#' @param mapSource which online source to use when downloading the map. Options include Google Maps ("google"), OpenStreetMap ("osm"), Stamen Maps ("stamen") and CloudMade maps ("cloudmade").
#' @param mapType the specific type of map to plot. Options available are "terrain", "satellite", "roadmap" and "hybrid" (google maps), "terrain-background", "terrain", "watercolor" and "toner" (stamen maps) or a positive integer for cloudmade maps (see ?get_cloudmademap from the package ggmap for details).
#' @param opacity value between 0 and 1 givin the opacity of surface colours.
#' @param plotContours whether or not to add contours to the surface plot.
#' @param breakPercent vector of values between 0 and 100 describing where in the surface contours appear.
#' @param contourCol single colour to plot contour lines showing boundaries on surface.
#' @param smoothScale should plot legend show continuous (TRUE) or discrete (FALSE) colours.
#' @param crimeCex relative size of symbols showing crimes.
#' @param crimeCol colour of crime symbols.
#' @param crimeBorderCol border colour of crime symbols.
#' @param crimeBorderWidth width of border of crime symbols.
#' @param sourceCex relative size of symbols showing suspect sites.
#' @param sourceCol colour of suspect sites symbols.
#' @param gpLegend whether or not to add legend to plot.
#'
#' @export
#' @examples
#' # London example data
#' d <- LondonExample_crimes
#' s <- LondonExample_sources
#' p <- geoParams(data = d, sigma_mean = 1.0, sigma_squared_shape = 2)
#' m <- geoMCMC(data = d, params = p)
#' # produce simple map
#' geoPlotMap(params = p, data = d, source = s, surface = m$geoProfile,
#'                 breakPercent = seq(0, 50, 5), mapType = "hybrid",
#'                 crimeCol = "black", crimeCex = 2, sourceCol = "red", sourceCex = 2)
#'
#' # John Snow cholera data
#' d <- Cholera
#' s <- WaterPumps
#' p <- geoParams(data = d, sigma_mean = 1.0, sigma_squared_shape = 2)
#' m <- geoMCMC(data = d, params = p, lambda=0.05)
#' # produce simple map
#' geoPlotMap(params = p, data = d, source = s, surface = m$geoProfile,
#'                 breakPercent = seq(0, 50, 5), mapType = "hybrid",
#'                 crimeCol = "black", crimeCex = 2, sourceCol = "red", sourceCex = 2)
#' 
#' # simulated data
#' sim <-rDPM(50, priorMean_longitude = -0.04217491, priorMean_latitude = 
#' 51.5235505, alpha=1, sigma=1, tau=3)
#' d <- geoData(sim$longitude, sim $latitude)
#' s <- geoDataSource(sim$source_lon, sim$source_lat)
#' p <- geoParams(data = d, sigma_mean = 1.0, sigma_squared_shape = 2)
#' m <- geoMCMC(data = d, params = p)
#' # change colour palette, map type, opacity and range of geoprofile and omit legend
#' geoPlotMap(params = p, data = d, source = s, surface = m$geoProfile,
#'                 breakPercent = seq(0, 30, 5), mapType = "terrain", 
#'                 surfaceCols = c("blue","white"), crimeCol = "black", 
#'                 crimeBorderCol = "white",crimeCex = 2, sourceCol = "red", sourceCex = 2,
#'                 opacity = 0.7, gpLegend = FALSE)

geoPlotMap <- function(params, data=NULL, source=NULL, surface=NULL, surfaceCols=NULL, zoom=NULL, latLimits=NULL, lonLimits=NULL, mapSource="google", mapType="hybrid", opacity=0.6, plotContours=TRUE, breakPercent=seq(0,100,l=11), contourCol= "grey50", smoothScale=TRUE, crimeCex=1.5, crimeCol='red', crimeBorderCol='white', crimeBorderWidth=0.5, sourceCex=1.5, sourceCol='blue', gpLegend=TRUE) {
  
  # check that inputs make sense
  geoParamsCheck(params)
  if (!is.null(data)) { geoDataCheck(data) }
  
  # set defaults
  if (is.null(surfaceCols)) { surfaceCols <- viridis::plasma(100) }
  if (is.null(latLimits)) { latLimits <- params$output$latitude_minMax }
  if (is.null(lonLimits)) { lonLimits <- params$output$longitude_minMax }
  
  # if zoom=="auto" then set zoom level based on params
  if (is.null(zoom)) { 
    zoom <- getZoom(params$output$longitude_minMax, params$output$latitude_minMax)
    cat(paste0("using zoom=", zoom, "\n"))
  }
  
  # make zoom level appropriate to map source
  if (mapSource=="stamen") { zoom <- min(zoom,18) }
  
  # download map
  cat("downloading map\n")
  loc <- c(mean(params$output$longitude_minMax), mean(params$output$latitude_minMax))
  rawMap <- get_map(location=loc, zoom=zoom, source=mapSource, maptype=mapType)
  
  # get attributes from rawMap (bounding box)
  att <- unlist(attributes(rawMap)$bb)
  latVec <- seq(att[3], att[1], l=nrow(rawMap))
  lonVec <- seq(att[2], att[4], l=ncol(rawMap))
  df_rawMap <- data.frame(lat=rep(latVec, each=ncol(rawMap)), lon=rep(lonVec, times=nrow(rawMap)))
  
  # bind with colours from rawMap
  df_rawMap <- cbind(df_rawMap, col=as.vector(rawMap))
  
  # create ggplot object
  myMap <- ggplot(df_rawMap, aes_string(x='lon', y='lat', fill='col')) + geom_raster() + scale_fill_identity()
  myMap <- myMap + coord_cartesian(xlim=lonLimits, ylim=latLimits, expand=FALSE)
  
  # overlay geoprofile
  if (!is.null(surface)) {
    
    # create colour palette
    geoCols <- colorRampPalette(rev(surfaceCols))
    nbcol <- length(breakPercent)-1
    
    # extract plotting ranges and determine midpoints of cells
    longitude_minMax  <- params$output$longitude_minMax
    latitude_minMax  <- params$output$latitude_minMax
    longitude_cells  <- params$output$longitude_cells
    latitude_cells  <- params$output$latitude_cells
    longitude_cellSize <- diff(longitude_minMax)/longitude_cells
    latitude_cellSize <- diff(latitude_minMax)/latitude_cells
    longitude_midpoints <- longitude_minMax[1] - longitude_cellSize/2 + (1:longitude_cells)* longitude_cellSize
    latitude_midpoints <- latitude_minMax[1] - latitude_cellSize/2 + (1:latitude_cells)* latitude_cellSize
    
    # create data frame of x,y,z values and labels for contour level
    df <- expand.grid(x=longitude_midpoints, y=latitude_midpoints)
    df$z <- as.vector(t(surface))
    labs <- paste(round(breakPercent,1)[-length(breakPercent)],"-",round(breakPercent,1)[-1],"%",sep='')
    df$cut <- cut(df$z, breakPercent, labels=labs)
    df$col <- rev(geoCols(nbcol))[df$cut]
    
    # remove all entries outside of breakPercent range
    df_noNA <- df[!is.na(df$cut),]
    
    # convert current map into borderless background image
    background <- myMap + theme_nothing()
    myMap <- ggplot() + annotation_custom(grob=ggplotGrob(background), xmin=lonLimits[1], xmax=lonLimits[2], ymin=latLimits[1], ymax=latLimits[2])
    
    # add surface and colour scale
    if (smoothScale) {
      myMap <- myMap + geom_raster(aes_string(x='x', y='y', fill='z'), alpha=opacity, data=df_noNA)
      myMap <- myMap + scale_fill_gradientn(name="Hitscore\npercentage", colours=rev(surfaceCols))
    } else {
      myMap <- myMap + geom_raster(aes_string(x='x', y='y', fill='col'), alpha=opacity, data=df_noNA)
      myMap <- myMap + scale_fill_manual(name="Hitscore\npercentage", labels=labs, values=geoCols(nbcol))
    }
    if (!gpLegend) {
      myMap <- myMap + theme(legend.position="none")
    }
    
    # add plotting limits
    myMap <- myMap + coord_cartesian(xlim=lonLimits, ylim=latLimits, expand=FALSE)
    
    # add contours
    if (plotContours) {
      myMap <- myMap + stat_contour(aes_string(x='x', y='y', z='z'), colour=contourCol, breaks=breakPercent, size=0.3, alpha=opacity, data=df)
    }
  }
  
  # overlay data points
  if (!is.null(data)) {
    df_data <- data.frame(longitude=data$longitude, latitude=data$latitude)
    myMap <- myMap + geom_point(aes_string(x='longitude', y='latitude'), data=df_data, pch=21, stroke=crimeBorderWidth, cex=crimeCex, fill=crimeCol, col=crimeBorderCol)
  }
  
  # overlay source points
  if (!is.null(source)) {
    df_source <- data.frame(longitude=source$longitude, latitude=source$latitude)
    myMap <- myMap + geom_point(aes_string(x='longitude', y='latitude'), data=df_source, pch=15, cex=sourceCex, col=sourceCol, fill=NA)
  }
  
  # force correct aspect ratio
  centre_lat <- mean(params$output$latitude_minMax)
  centre_lon <- mean(params$output$longitude_minMax)
  scale_lat <- latlon_to_bearing(centre_lat, centre_lon, centre_lat + 0.1, centre_lon)$gc_dist
  scale_lon <- latlon_to_bearing(centre_lat, centre_lon, centre_lat, centre_lon + 0.1)$gc_dist
  asp <- diff(latLimits)*scale_lat / (diff(lonLimits)*scale_lon)
  myMap <- myMap +  theme(aspect.ratio=asp)
  
  # add labels
  myMap <- myMap +  labs(x="longitude", y="latitude")
  
  # plot map
  myMap
}

#------------------------------------------------
#' Perspective plot of geoprofile or raw probabilities
#'
#' Plots persp plot of geoprofile or posterior surface (coloured according to height), reducing matrix dimensions if necessary to avoid grid lines being too close together. NB Only works with square matrix
#'
#' @param surface surface to plot; either the geoprofile or posteriorSurface output by geoMCMC(). 
#' @param aggregate_size the number of cells to aggregate to smooth the surface.
#' @param surface_type type of surface; should be either "gp" for geoprofile or "prob" for posteriorSurface.
#' @param perspCol colour palette. Defaults to red/orange/yellow/white.
#' @param phiGP value of phi to pass to persp().
#' @param thetaGP value of theta to pass to persp().
#'
#' @export
#' @examples
#' # John Snow cholera data
#' d <- geoData(Cholera$longitude, Cholera$latitude)
#' s <- geoDataSource(WaterPumps$longitude, WaterPumps$latitude)
#' p <- geoParams(data = d, sigma_mean = 1.0, sigma_squared_shape = 2)
#' m <- geoMCMC(data = d, params = p, lambda=0.05)
#' # raw probabilities
#' geoPersp(m$posteriorSurface, surface_type = "prob")
#' # geoprofile
#' geoPersp(m$geoProfile, aggregate_size = 3, surface_type = "gp")
#' 
#' # simulated data
#' sim <-rDPM(50, priorMean_longitude = -0.04217491, priorMean_latitude = 
#' 51.5235505, alpha=1, sigma=1, tau=3)
#' d <- geoData(sim$longitude, sim $latitude)
#' s <- geoDataSource(sim$source_lon, sim$source_lat)
#' p <- geoParams(data = d, sigma_mean = 1.0, sigma_squared_shape = 2)
#' m <- geoMCMC(data = d, params = p)
#' # raw probabilities
#' geoPersp(m$posteriorSurface, surface_type = "prob")
#' # geoprofile
#' geoPersp(surface = m$geoProfile, aggregate_size = 3, surface_type = "gp")

geoPersp <- function(surface, aggregate_size=3, surface_type="gp", perspCol=c("red", "orange", "yellow", "white"), phiGP=30, thetaGP=-30) {
  
  # rescale and transform surface
  scale <- 1
  if(surface_type=="gp") { scale <- -1 }
  surface <- t(scale*surface)
  
  # reduce matrix or not
  matrix_size <- unique(dim(surface))[1]
  breaks <- seq(1, (matrix_size-(aggregate_size-1)), aggregate_size)
  nb <- length(breaks)
  output <- matrix(NA, nrow=nb, ncol=nb)
  for (i in 1:nb) {
    for (j in 1:nb) {
      output[i,j] <- mean(as.vector(surface[breaks[i]:(breaks[i]+(aggregate_size-1)), breaks[j]:(breaks[j]+(aggregate_size-1))]))
    }
  }
  
  # Generate the desired number of colors
  colpal <- colorRampPalette(perspCol)
  color <- colpal(100)
  
  # Compute the z-value at the facet centres
  ncz <- ncol(output)
  nrz <- nrow(output)
  zfacet <- output[-1, -1] + output[-1, -ncz] + output[-nrz, -1] + output[-nrz, -ncz]
  facetcol <- cut(zfacet, length(color))
  
  persp(output, col=color[facetcol], border="black", phi=phiGP, theta=thetaGP, lwd=0.2, box=FALSE)
}

#------------------------------------------------
#' Produce Lorenz Plot
#'
#' Produces a Lorenz plot showing the proportion of suspect sites or cimes identified as a function of area and calculates
#' the corresponding Gini coefficient using trapezoid rule.
#' Also allows an optional vector called crimeNumbers with numbers of crimes per suspect site; the length of this vector
#' should equal the number of suspect sites. If this is present, the function calculates and returns the Gini coefficient 
#' based on the number of crimes; otherwise, this is calculated based on the number of suspect sites.
#'
#' @param hit_scores object in the format defined by geoReportHitscores().
#' @param crimeNumbers optional vector with numbers of crimes per suspect site.
#' @param suspects_col colour to plot curve for sources.
#' @param crimes_col colour to plot curve for crimes if crimeNumbers is supplied.
#'
#' @export
#' @examples
#' # London example data
#' d <- LondonExample_crimes
#' s <- LondonExample_sources
#' p <- geoParams(data = d, sigma_mean = 1.0, sigma_squared_shape = 2)
#' m <- geoMCMC(data = d, params = p)
#' hs <- geoReportHitscores(params = p, source = s, surface = m$geoProfile)
#' hs 
#' # Lorenz plot on sources
#' geoPlotLorenz(hs)
#' # Lorenz plot on sources and crimes
#' # extract numbers of crimes allocated per source as a proxy
#' cn <- as.vector(table(m$bestGrouping))
#' geoPlotLorenz(hs, crimeNumbers = cn)

geoPlotLorenz <- function(hit_scores, crimeNumbers=NULL, suspects_col="red", crimes_col="blue") {
  
  # trapezoid rule
  tpzd <- function(x,y) {
    idx = 2:length(x)
    return (as.double( (x[idx] - x[idx-1]) %*% (y[idx] + y[idx-1])) / 2)
  }
  
  # bind crimeNumbers if using
  crimeNumbers_on <- !is.null(crimeNumbers)
  if (crimeNumbers_on) {
    hit_scores <- cbind(hit_scores, crimeNumbers)
  }
  
  # put hitscore data frame in increasing order
  hit_scores <- hit_scores[order(hit_scores[,3]),]
  hs <- hit_scores[,3]
  n <- nrow(hit_scores)
  
  # cumulative crime sites
  cum_suspect_sites <- (1:n)/n
  
  # add initial 0 and terminal 1 to each so plots will begin at (0,0) and end at (1,1)
  hs <- c(0,hs,100)
  cum_suspect_sites <- c(0,cum_suspect_sites,1)
  
  # calculate gini coefficient for crimes
  auc <- (tpzd(hs/100, cum_suspect_sites) - 0.5)/0.5
  G_sources <- round(auc,3)
  
  # produce plot
  if (!crimeNumbers_on) {
    
    # plot
    plot(0, type="n", xlim=c(0,100), ylim=c(0,1), xlab="hit score percentage", ylab="proportion identified")
    abline(h=seq(0,1,0.2), v=seq(0,100,20), col="lightgray", lwd=0.4)
    abline(0, 0.01, col="gray")
    lines(hs, cum_suspect_sites, col=suspects_col)
    
    # add Gini value text and legend
    Gini_text <- paste0("G (sources) = ",G_sources)
    legend(60, 0.4, Gini_text, cex=0.8, bty="n")
    legend(65, 0.25, c("sources"), col=c(suspects_col), lwd=1, cex=0.8)
    
    # return Gini value(s) silently
    invisible(c(G_sources=G_sources))
    
  } else {
    
    # cumulative crime numbers
    cum_crimes <- cumsum(hit_scores[,4])/sum(hit_scores[,4])
    cum_crimes <- c(0,cum_crimes,1)
    
    # calculate gini coefficient for crimes
    auc <- (tpzd(hs/100, cum_crimes) - 0.5)/0.5
    G_crimes <- round(auc,3)
    
    # plot
    plot(0, type="n", xlim=c(0,100), ylim=c(0,1), xlab="hit score percentage", ylab="proportion identified")
    abline(h=seq(0,1,0.2), v=seq(0,100,20), col="lightgray", lwd=0.4)
    abline(0, 0.01, col="gray")
    lines(hs, cum_suspect_sites, col=suspects_col)
    lines(hs, cum_crimes, col=crimes_col)
    
    # add Gini value text and legend
    Gini_text <- c(paste0("G (sources) = ",G_sources), paste0("G (incidents) = ",G_crimes))
    legend(60, 0.5, Gini_text, cex=0.8, bty="n")
    legend(65, 0.25, c("sources","incidents"), col=c(suspects_col,crimes_col), lwd=1, cex=0.8)
    
    # return Gini value(s) silently
    invisible(c(G_sources=G_sources, G_incidents=G_crimes))
  }
  
}

#------------------------------------------------
#' Calculate and plot probability of coallocation
#'
#' For all pairs of crimes, calculates the probability that both originate from the same source and plots a coloured half matrix representing these data. The data underlying these calculations can be accessed as the object $coAllocation produced by geoMCMC().
#'
#' @param mcmc object of the type output by geoMCMC().
#' @param cols colour palette to use. Defaults to viridis palette.
#'
#' @export
#' @examples
#' # London example data
#' d <- LondonExample_crimes
#' s <- LondonExample_sources
#' p <- geoParams(data = d, sigma_mean = 1.0, sigma_squared_shape = 2)
#' m <- geoMCMC(data = d, params = p)
#' # produce simple map
#' geoPlotMap(params = p, data = d, source = s, surface = m$geoProfile,
#'                 breakPercent = seq(0, 50, 5), mapType = "hybrid",
#'                 crimeCol = "black", crimeCex = 2, sourceCol = "red", sourceCex = 2)
#' # calculate coallocation matrix and plot
#' geoPlotCoallocation(m)

geoPlotCoallocation <- function(mcmc, cols=NULL) {
  
  # default colours
  if (is.null(cols)) {
    cols <- viridis::viridis(100)
  }
  
  # get co-allocation matrix from raw MCMC output
  comat <- mcmc$coAllocation
  
  # produce data frame for ggplot
  n <- nrow(comat)
  df <- data.frame(x=as.vector(row(comat)), y=as.vector(col(comat)), z=as.vector(t(comat)))
  
  # produce plot
  gg <- ggplot(df) + geom_tile(aes_string(x='x', y='y', fill='z')) + scale_fill_gradientn(colours=cols, name="Probability\nco-allocation") + coord_cartesian(xlim=c(0.5,n+0.5), ylim=c(0.5,n+0.5), expand=FALSE) + labs(x="observation", y="observation")
  
  gg
}

#------------------------------------------------
#' Unknown pleasures
#' 
#' A frivolous alternative to geoPlotMap(), this function takes the output of geoMCMC() and plots the resulting geoprofile in the style of the cover of Joy Division's 'Unknown pleasures' album.
#' 
#' @param input_matrix The surface to plot, usually the object $geoProfile produced by geoMCMC().
#' @param nlines The number of lines (defaults to the correct number of 80).
#' @param paper_ref A text string, for example a reference to a paper.
#' @param bgcol Background colour
#' @param fgcol Foreground colour
#' @param wt line weight
#' 
#' @export
#' @examples
#' # simulated data
#' sim <-rDPM(50, priorMean_longitude = -0.04217491, priorMean_latitude = 
#' 51.5235505, alpha=10, sigma=1, tau=3)
#' d <- geoData(sim$longitude, sim $latitude)
#' s <- geoDataSource(sim$source_lon, sim$source_lat)
#' p <- geoParams(data = d, sigma_mean = 1.0, sigma_squared_shape = 2)
#' m <- geoMCMC(data = d, params = p)
#' unknownPleasures(m$geoProfile, paper_ref = "Rgeoprofile v2.1.0")

unknownPleasures <- function(input_matrix, paper_ref = NULL, nlines = 80, bgcol = "black", fgcol = "white", wt = 2) {
  
  orig_prof <- input_matrix
  citation <- paper_ref
  # extract other params
  ncols <- ncol(orig_prof)
  # functions
  expandMatrix <- function(mat,output_long,output_lat)
  {
    # define function expanding vector
    expandVector <- function(input_vec,output_length)
    {
      my_vec <- input_vec
      desired_length <- output_length
      new_vec <- rep(NA, desired_length)
      
      vec_ID <- seq(1,length(my_vec),length=desired_length)
      
      for(i in 1:length(new_vec))
      {
        ifelse(vec_ID[i] %% 1 == 0,
               new_vec[i] <- my_vec[floor(vec_ID[i])],
               new_vec[i] <- mean((1-vec_ID[i] %% 1) * my_vec[floor(vec_ID[i])] + (vec_ID[i] %% 1) * my_vec[ceiling(vec_ID[i])])
        )
      }
      return(new_vec)
    }
    mat1 <- apply(mat,2, function(x) expandVector(x, output_long))
    mat2 <- apply(mat1,1, function(x) expandVector(x, output_lat))
    return(t(mat2))	
  }
  
  # reduce to manageable number of rows and columns!
  reduced_mat <- expandMatrix(orig_prof ,nlines, ncols)
  
  # scale so values fall between -0.5 and +0.5
  reduced_mat <- 1-reduced_mat/max(reduced_mat)-0.5
  
  # set y coordinates of lines
  yvals <- seq(0, 1,length = nlines)
  
  # plot
  oldPar <- par(bg = bgcol)
  on.exit(par(oldPar))
  
  plot(1:ncols, reduced_mat[nlines,]+yvals[nlines], type="l", ylim=c(-1,max(yvals)+0.5), axes=FALSE, xlab="", ylab="", col=fgcol)
  
  for (i in nlines:1) {
    polygon(c(min(reduced_mat[i,]), reduced_mat[i,]+yvals[i], min(reduced_mat[i,])), col=bgcol, border=bgcol)
    points(1:ncol(reduced_mat), reduced_mat[i,]+yvals[i], type="l", col=fgcol, lwd=wt)
  }
  text(0,-0.6,citation,adj=0,col=fgcol)
}
