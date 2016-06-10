
#------------------------------------------------
#' Draw from Dirichlet process mixture model
#'
#' This function provides random draws from a 2D Dirichlet process mixture model. Coordinates are defined in units of degrees latitude and longitude to facilitate spatial analysis.
#'
#' @param n number of draws
#' @param sigma standard deviation of dispersal distribution, in units of degrees
#'
#' @export
#' @examples
#' rDPM(10)

rDPM <- function(n, sigma=0.01, priorMean_longitude=-0.1277, priorMean_latitude=51.5074, priorSD=0.03, alpha=1) {
    
    # force n to be a scalar integer
    n <- floor(n[1])
    
    # draw grouping
	group <- rep(1,n)
	freqs <- 1
    if (n>1) {
        for (i in 2:n) {
            group[i] <- sample(length(freqs)+1,1,prob=c(freqs,alpha))
            if (group[i]>length(freqs))
                freqs <- c(freqs,0)
            freqs[group[i]] <- freqs[group[i]] + 1
        }
    }
	group <- sort(group)
	
	# draw means and data
	source_lon <- rnorm(length(freqs),priorMean_longitude,priorSD)
	source_lat <- rnorm(length(freqs),priorMean_latitude,priorSD)
    
	longitude <- rnorm(n,source_lon[group],sigma)
	latitude <- rnorm(n,source_lat[group],sigma)
	
	# return results
	return(list(longitude=longitude, latitude=latitude, group=group, source_lon=source_lon, source_lat=source_lat))
}

#------------------------------------------------
#' Create Rgeoprofile data object
#'
#' This function can be used to generate a dummy data object in the format used by the Rgeoprofile MCMC.
#'
#' @param longitude The longitudinal positions of the observed data
#' @param latitude The latitudinal positions of the observed data
#'
#' @export
#' @examples
#' geoData()

geoData <- function(longitude=NULL, latitude=NULL) {
    
    # generate dummy data if none read in
    if (is.null(longitude) & is.null(latitude)) {
        longitude <- c(-0.104545976281589, -0.102659272660916, -0.0967390020136406, -0.0996246226730725, -0.100775342233937, -0.101073477576196, -0.100932674617746, -0.0983001766339886, -0.0913571765598557, -0.100211479242536, -0.139508969429415, -0.14403082311245, -0.143607222414313, -0.137174795971723, -0.140884394738737, -0.142723755125487, -0.143380928147727, -0.136989691342132, -0.13837666855334, -0.138297288871952, -0.0773858357074935, -0.0818743917621333, -0.0738310357273188, -0.0744118149244568, -0.0757833597110897, -0.0762193916493531, -0.0810467015747727, -0.110052994420826, -0.106600836874167, -0.105104028808356, -0.101934241194567, -0.0683752111183375, -0.0758607240702608, -0.079153744918552, -0.087964365345432)
        latitude <- c(51.4996147329979, 51.4925230579844, 51.4947129689414, 51.4922683109254, 51.5007532206834, 51.4960640374896, 51.4996917836745, 51.4976936749008, 51.4977904998888, 51.4894186202378, 51.5002583182117, 51.5033510595094, 51.4984697991335, 51.5063306206839, 51.4961516950408, 51.4994464819411, 51.5067557678594, 51.4977275537675, 51.4988718377984, 51.4974782970503, 51.5137643501102, 51.5204498816501, 51.5213788858189, 51.5144343479237, 51.5212383455566, 51.5088225370868, 51.512547894056, 51.5144758355252, 51.5218865924773, 51.5218808497196, 51.5152330574081, 51.4836680563637, 51.4885211991595, 51.486842412489, 51.48845301363455)
    }
    
    # combine and return
    data <- list(longitude=longitude, latitude=latitude)
    return(data)
}

#------------------------------------------------
#' Create Rgeoprofile parameters object
#'
#' This function can be used to generate parameters in the format used by the Rgeoprofile MCMC. Parameter values can be specified as input arguments.
#'
#' @param data
#' @param sigma The standard deviation of the dispersal distribution
#' @param priorMean_longitude The position (longitude) of the mean of the prior distribution
#' @param priorMean_latitude The position (latitude) of the mean of the prior distribution
#' @param priorSD
#' @param alpha_shape
#' @param alpha_rate
#' @param chains
#' @param burnin
#' @param samples
#' @param burnin_printConsole
#' @param samples_printConsole
#' @param longitude_minMax
#' @param latitude_minMax
#' @param longitude_cells
#' @param latitude_cells
#'
#' @export
#' @examples
#' geoParams()

geoParams <- function(data=NULL, sigma=0.008, priorMean_longitude=-0.1277, priorMean_latitude=51.5074, priorSD=0.03, alpha_shape=0.1, alpha_rate=0.1, chains=10, burnin=500, samples=5000, burnin_printConsole=100, samples_printConsole=1000, longitude_minMax=NULL, latitude_minMax=NULL, longitude_cells=500, latitude_cells=500) {
    
	# if data argument used then get map limits from data
    if (!is.null(data)) {
        
        # check correct format of data
        geoDataCheck(data, silent=TRUE)
        
		# get midpoints and ranges
	    xmin <- min(data$longitude);
        xmax <- max(data$longitude)
	    ymin <- min(data$latitude);
        ymax <- max(data$latitude)
        xdiff <- diff(range(data$longitude))
        ydiff <- diff(range(data$latitude))
        xmid <- xmin + xdiff/2
        ymid <- ymin + ydiff/2
        
        # calculate latitude upper and lower limits corresponding to a square map
        delta <- xdiff/2
        projection_mid <- log(tan(pi/4+ymid/360*2*pi/2))
        projection_top <- projection_mid + delta/360*2*pi
        projection_bot <- projection_mid - delta/360*2*pi
        lat_angle_top <- (2*atan(exp(projection_top))-pi/2)*360/(2*pi)
        lat_angle_bot <- (2*atan(exp(projection_bot))-pi/2)*360/(2*pi)
		
		# if data within these limits then great. Otherwise try the reverse operation - calculate longitude left and right limits corresponding to a square map. In both cases ass a 10% buffer zone.
        if (ymin>=lat_angle_bot & ymax<=lat_angle_top) {
        	frame_xmin <- xmin-xdiff*0.1
        	frame_xmax <- xmax+xdiff*0.1
			frame_ymin <- (lat_angle_top+lat_angle_bot)/2-(lat_angle_top-lat_angle_bot)*0.6
			frame_ymax <- (lat_angle_top+lat_angle_bot)/2+(lat_angle_top-lat_angle_bot)*0.6
        } else {
        	lat_angle_bot <- ymin
        	lat_angle_top <- ymax
        	projection_bot <- log(tan((lat_angle_bot/360*(2*pi)+pi/2)/2))
        	projection_top <- log(tan((lat_angle_top/360*(2*pi)+pi/2)/2))
			projection_mid <- (projection_bot+projection_top)/2
			delta <- (projection_top-projection_mid)*360/(2*pi)
			frame_xmin <- xmid - delta*1.1
			frame_xmax <- xmid + delta*1.1
			frame_ymin <- lat_angle_bot - delta*0.1
			frame_ymax <- lat_angle_top + delta*0.1
        }
        
        # set output values
        if (is.null(longitude_minMax))
            longitude_minMax <- c(frame_xmin, frame_xmax)
        if (is.null(latitude_minMax))
            latitude_minMax <- c(frame_ymin, frame_ymax)
        
    } else {
        if (is.null(longitude_minMax))
            longitude_minMax <- -0.1277 + c(-0.1,0.1)
        if (is.null(latitude_minMax))
            latitude_minMax <- 51.5074 + c(-0.1,0.1)
    }
    
    # set model parameters
    model <- list(sigma=sigma, priorMean_longitude=priorMean_longitude, priorMean_latitude=priorMean_latitude, priorSD=priorSD, alpha_shape=alpha_shape, alpha_rate=alpha_rate)
    
    # set MCMC parameters
    MCMC <- list(chains=chains, burnin=burnin, samples=samples, burnin_printConsole=burnin_printConsole, samples_printConsole=samples_printConsole)
    
    # set output parameters
    xmin <- longitude_minMax[1];
    xmax <- longitude_minMax[2]
    ymin <- latitude_minMax[1];
    ymax <- latitude_minMax[2]
    xCells <- longitude_cells
    yCells <- latitude_cells
    xCellSize <- (xmax-xmin)/xCells
    yCellSize <- (ymax-ymin)/yCells
    xMids <- xmin - xCellSize/2 + (1:xCells)*xCellSize
    yMids <- ymin - yCellSize/2 + (1:yCells)*yCellSize
    
    output <- list(longitude_minMax=longitude_minMax, latitude_minMax=latitude_minMax, longitude_cells=longitude_cells, latitude_cells=latitude_cells, longitude_midpoints=xMids, latitude_midpoints=yMids)
    
    # combine and return
    params <- list(model=model, MCMC=MCMC, output=output)
    return(params)
}

#------------------------------------------------
#' Check data
#'
#' Check that all data for use in Rgeoprofile MCMC is in the correct format.
#'
#' @param data A list of data
#' @param silent Whether to report passing check to console
#'
#' @export
#' @examples
#' # tester
#' geoDataCheck(myData)

geoDataCheck <- function(data, silent=FALSE) {
    
    # check that data is a list
    if (!is.list(data))
    stop("data must be in list format")
    
    # check that contains longitude and latitude
    if (!("longitude"%in%names(data)))
    stop("data must contain element 'longitude'")
    if (!("latitude"%in%names(data)))
    stop("data must contain element 'latitude'")
    
    # check that data values are correct format and range
    if (!is.numeric(data$longitude) | !all(is.finite(data$longitude)))
    stop("data$longitude values must be numeric and finite")
    if (!is.numeric(data$latitude) | !all(is.finite(data$latitude)))
    stop("data$latitude values must be numeric and finite")
    
    # check same number of observations in logitude and latitude, and n>1
    if (length(data$longitude)!=length(data$latitude))
    stop("data$longitude and data$latitude must have the same length")
    if (length(data$longitude)<=1)
    stop("data must contain at least two observations")
    
    # if passed all checks
    if (!silent)
    cat("data file passed all checks\n")
}

#------------------------------------------------
#' Check parameters
#'
#' Check that all parameters for use in Rgeoprofile MCMC are OK.
#'
#' @param params A list of parameters
#' @param silent Whether to report passing check to console
#'
#' @export
#' @examples
#' # tester
#' geoParamsCheck(myParams)

geoParamsCheck <- function(params, silent=FALSE) {
    
    # check that params is a list
    if (!is.list(params))
        stop("params must be in list format")
        
    # check that contains 'model' and 'MCMC' as sublists
    if (!("model"%in%names(params) & "MCMC"%in%names(params)))
        stop("params must contain sublists 'model' and 'MCMC'")

    # check that 'model' and 'MCMC' are indeed lists
    if (!is.list(params$model))
        stop("params$model must be in list format")
    if (!is.list(params$MCMC))
        stop("params$MCMC must be in list format")

    # check that params$model contains all necessary parameters
    if (!("sigma"%in%names(params$model)))
        stop("params$model must contain parameter 'sigma'")
    if (!("priorMean_longitude"%in%names(params$model)))
        stop("params$model must contain parameter 'priorMean_longitude'")
    if (!("priorMean_latitude"%in%names(params$model)))
        stop("params$model must contain parameter 'priorMean_latitude'")
    if (!("priorSD"%in%names(params$model)))
        stop("params$model must contain parameter 'priorSD'")
    if (!("alpha_shape"%in%names(params$model)))
        stop("params$model must contain parameter 'alpha_shape'")
    if (!("alpha_rate"%in%names(params$model)))
        stop("params$model must contain parameter 'alpha_rate'")

    # check that params$MCMC contains all necessary parameters
    if (!("chains"%in%names(params$MCMC)))
        stop("params$MCMC must contain parameter 'chains'")
    if (!("burnin"%in%names(params$MCMC)))
        stop("params$MCMC must contain parameter 'burnin'")
    if (!("samples"%in%names(params$MCMC)))
        stop("params$MCMC must contain parameter 'samples'")
    if (!("burnin_printConsole"%in%names(params$MCMC)))
        stop("params$MCMC must contain parameter 'burnin_printConsole'")
    if (!("samples_printConsole"%in%names(params$MCMC)))
        stop("params$MCMC must contain parameter 'samples_printConsole'")
    
    # check that params$model values are correct format and range
    if (!is.numeric(params$model$sigma) | !is.finite(params$model$sigma))
        stop("params$model$sigma must be numeric and finite")
    if (params$model$sigma<=0)
        stop("params$model$sigma must be greater than 0")
    if (!is.numeric(params$model$priorMean_longitude) | !is.finite(params$model$priorMean_longitude))
        stop("params$model$priorMean_longitude must be numeric and finite")
    if (!is.numeric(params$model$priorMean_latitude) | !is.finite(params$model$priorMean_latitude))
        stop("params$model$priorMean_latitude must be numeric and finite")
    if (!is.numeric(params$model$priorSD) | !is.finite(params$model$priorSD))
        stop("params$model$priorSD must be numeric and finite")
    if (params$model$priorSD<=0)
        stop("params$model$priorSD must be greater than 0")
    if (!is.numeric(params$model$alpha_shape) | !is.finite(params$model$alpha_shape))
        stop("params$model$alpha_shape must be numeric and finite")
    if (params$model$alpha_shape<=0)
        stop("params$model$alpha_shape must be greater than 0")
    if (!is.numeric(params$model$alpha_rate) | !is.finite(params$model$alpha_rate))
        stop("params$model$alpha_rate must be numeric and finite")
    if (params$model$alpha_rate<=0)
        stop("params$model$alpha_rate must be greater than 0")

    # check that params$MCMC values are correct format and range
    if (!is.numeric(params$MCMC$chains) | !is.finite(params$MCMC$chains))
        stop("params$MCMC$chains must be numeric and finite")
    if (params$MCMC$chains<=1)
        stop("params$MCMC$chains must be 2 or more")
    if (!is.numeric(params$MCMC$burnin) | !is.finite(params$MCMC$burnin))
        stop("params$MCMC$burnin must be numeric and finite")
    if (params$MCMC$burnin<0)
        stop("params$MCMC$burnin must be greater than or equal to 0")
    if (!is.numeric(params$MCMC$samples) | !is.finite(params$MCMC$samples))
        stop("params$MCMC$samples must be numeric and finite")
    if (params$MCMC$samples<=0)
        stop("params$MCMC$samples must be greater than 0")
    if (!is.numeric(params$MCMC$burnin_printConsole) | !is.finite(params$MCMC$burnin_printConsole))
        stop("params$MCMC$burnin_printConsole must be numeric and finite")
    if (params$MCMC$burnin_printConsole<=0)
        stop("params$MCMC$burnin_printConsole must be greater than 0")
    if (!is.numeric(params$MCMC$samples_printConsole) | !is.finite(params$MCMC$samples_printConsole))
        stop("params$MCMC$samples_printConsole must be numeric and finite")
    if (params$MCMC$samples_printConsole<=0)
        stop("params$MCMC$samples_printConsole must be greater than 0")

    # if passed all checks
    if (!silent)
	    cat("params file passed all checks\n")
}

#------------------------------------------------
# Scaled Student's t distribution
# (not exported)

dts <- function(x,df,scale=1,log=FALSE) {
    output <- lgamma((df+1)/2)-lgamma(df/2)-0.5*log(pi*df*scale^2) - ((df+1)/2)*log(1 + x^2/(df*scale^2))
    if (!log)
    output <- exp(output)
    return(output)
}

#------------------------------------------------
#' MCMC under Rgeoprofile model
#'
#' This function carries out the main MCMC under the Rgeoprofile model.
#'
#' @param data input data in the format defined by geoData
#' @param params input parameters in the format defined by geoParams
#'
#' @export
#' @examples
#' geoMCMC()

geoMCMC <- function(data, params) {
    
    # check that data and parameters in correct format
    geoDataCheck(data)
    geoParamsCheck(params)
    cat("\n")
    
    # extract some values from params
    n <- length(data$longitude)
    sigma <- params$model$sigma
    priorMean_longitude <- params$model$priorMean_longitude
    priorMean_latitude <- params$model$priorMean_latitude
    priorSD <- params$model$priorSD
    samples <- params$MCMC$samples
    longitude_cells <- params$output$longitude_cells
    latitude_cells <- params$output$latitude_cells
    
    # carry out MCMC
    rawOutput <- C_geoMCMC(data, params)
    
    # extract raw draws and check that at least one draw in chosen region
    surface_raw <- matrix(unlist(rawOutput$geoSurface),latitude_cells,byrow=TRUE)
    if (all(surface_raw==0))
        stop('chosen lat/long window contains no posterior draws')
    
    # get some basic properties of the surface
    lon_min <- params$output$longitude_minMax[1]
    lon_max <- params$output$longitude_minMax[2]
    lat_min <- params$output$latitude_minMax[1]
    lat_max <- params$output$latitude_minMax[2]
    cells_lon <- longitude_cells
    cells_lat <- latitude_cells
    cellSize_lon <- (lon_max-lon_min)/cells_lon
    cellSize_lat <- (lat_max-lat_min)/cells_lat
    
    # set lambda (bandwidth) increment size based on cell size
    lambda_step <- min(cellSize_lon,cellSize_lat)/5
    
    # temporarily add guard rail to surface to avoid Fourier series bleeding round edges
    rail_lon <- ceiling(200*lambda_step/cellSize_lon)
    rail_lat <- ceiling(200*lambda_step/cellSize_lat)
    railMat_lon <- matrix(0,cells_lat,rail_lon)
    railMat_lat <- matrix(0,rail_lat,cells_lon+2*rail_lon)
    
    surface_normalised <- surface_raw/sum(surface_raw)
    surface_normalised <- cbind(railMat_lon, surface_normalised, railMat_lon)
    surface_normalised <- rbind(railMat_lat, surface_normalised, railMat_lat)
    
    f1 = fftw2d(surface_normalised)
    
    kernel_lon <- cellSize_lon * c(0:floor(ncol(surface_normalised)/2), floor((ncol(surface_normalised)-1)/2):1)
    kernel_lat <- cellSize_lat * c(0:floor(nrow(surface_normalised)/2), floor((nrow(surface_normalised)-1)/2):1)
    kernel_lon_mat <- outer(rep(1,length(kernel_lat)), kernel_lon)
    kernel_lat_mat <- outer(kernel_lat, rep(1,length(kernel_lon)))
    kernel_s_mat <- sqrt(kernel_lon_mat^2+kernel_lat_mat^2)
    
    logLike <- -Inf
    for (i in 1:100) {
        
        lambda <- lambda_step*i
        kernel <- dts(kernel_s_mat,df=3,scale=lambda)
        f2 = fftw2d(kernel)
        
        #### carry out fast Fourier transform, combine, and take inverse
        f3 = f1*f2
        f4 = Re(fftw2d(f3,inverse=T))/length(surface_normalised)
        f5 <- f4 - surface_normalised*dts(0,df=3,scale=lambda)
        f5[f5<0] <- 0
        f5 <- f5/sum(f4)
        f6 <- surface_normalised*log(f5)
        
        if (sum(f6,na.rm=T)<logLike)
        break()
        logLike <- sum(f6,na.rm=T)

    }
    f4 <- f4[,(rail_lon+1):(ncol(f4)-rail_lon)]
    f4 <- f4[(rail_lat+1):(nrow(f4)-rail_lat),]
    
    # produce prior matrix. Note that each cell of this matrix contains the probability density at that point multiplied by the size of that cell, meaning the total sum of the matrix from -infinity to +infinity would equal 1. As the matrix is limited to the region specified by the limits, in reality this matrix will sum to some value less than 1.
    lon_mids <- params$output$longitude_midpoints
    lat_mids <- params$output$latitude_midpoints
    lon_mids_mat <- outer(rep(1,cells_lat),lon_mids)
    lat_mids_mat <- outer(lat_mids,rep(1,cells_lon))
    
    priorMat <- dnorm(lon_mids_mat,priorMean_longitude,sd=priorSD)*dnorm(lat_mids_mat,priorMean_latitude,sd=priorSD)*(cellSize_lon*cellSize_lat)
    
    
    # finalise output format
    output <- list()
    
    # alpha
    alpha <- rawOutput$alpha
    output$alpha <- alpha
    
    # combine prior surface with stored posterior surface (the prior never fully goes away under a DPM model)
    output$surface_raw <- surface_raw
    output$surface <-  f4 + priorMat*mean(alpha/(alpha+n))
    
    # posterior allocation
    allocation <- matrix(unlist(rawOutput$allocation),n,byrow=T)
    allocation <- data.frame(allocation/samples)
    names(allocation) <- paste("group",1:ncol(allocation),sep="")
    output$allocation <- allocation
    
    return(output)
}

#------------------------------------------------
#' Calculate geoprofile from surface
#'
#' Converts surface to rank order geoprofile.
#' @param z matrix to convert to geoprofile
#' @export
#' @examples
#' geoProfile(MCMCoutput$geoSurface)

geoProfile <- function(z) {
    
    # check that z is in correct format
    if (!all(is.finite(z)) | !all(is.numeric(z)))
        stop("values in z must be finite and numeric")
    if (!is.matrix(z))
        stop("z must be a matrix")
    
    # create geoprofile from z
    z[order(z,decreasing=TRUE)] <- 1:length(z)
    
    return(z)
}

#------------------------------------------------
#' Plot posterior allocation
#'
#' Produces plot of posterior allocation, given output of MCMC.
#' @param MCMCoutput output generated from an MCMC.
#' @export
#' @examples
#' geoPlotAllocation(MCMCoutput)

geoPlotAllocation <- function(allocation, colours="default", barBorderCol="white", barBorderWidth=0.25, mainBorderCol="black", mainBorderWidth=2, yTicks_on=TRUE, yTicks=seq(0,1,0.2), xlab="", ylab="posterior allocation", mainTitle="", names=NA, names_size=1, xTicks_on=FALSE, xTicks_size=1, orderBy="group") {
    
    # check that orderBy is either 'group' or 'probability'
    if (!(orderBy%in%c("group","probability")))
        stop("orderBy must equal 'group' or 'probability'")
    
    # check that allocation is a data frame
    if (!is.data.frame(allocation))
        stop("allocation must be a data frame, with observations in rows and groups in columns")
    
    n <- nrow(allocation)
    k <- ncol(allocation)
    
    # replace colours if default
    if (colours=="default") {
        rawCols <- brewer.pal(n=11,name="RdYlBu")
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
getZoom <- function(x,y) {
	
	# calculate midpoint of range in x and y
	xmid <- min(x)+diff(range(x))/2
	ymid <- min(y)+diff(range(y))/2
	
	# calculate delta (half of longitude range) for a range of zoom levels
	z <- 20:2
	delta <- 445/(2^z)
	
	# calculate left and right longitude limits at all zoom levels
	long_angle_left <- xmid-delta
	long_angle_right <- xmid+delta
	
	# calculate top and bottom latitude limits at all zoom levels
	lat_angle <- ymid/360*2*pi
	projection_mid <- log(tan(pi/4+lat_angle/2))
	projection_top <- projection_mid + delta/360*2*pi
	projection_bot <- projection_mid - delta/360*2*pi
	lat_angle_top <- (2*atan(exp(projection_top))-pi/2)*360/(2*pi)
	lat_angle_bot <- (2*atan(exp(projection_bot))-pi/2)*360/(2*pi)
	
	# find the most zoomed-in level that captures all points
	zoomTest <- (min(x)>long_angle_left) & (max(x)<long_angle_right) & (min(y)>lat_angle_bot) & (max(y)<lat_angle_top)
    if (!any(zoomTest))
        stop("values are outside of plotting range of earth")
	bestZoom <- z[which(zoomTest)[1]]
    
	return(bestZoom)
}

#------------------------------------------------
#' Create geoprofile plot object
#'
#' Creates geoprofile plotting object, for use with other ggplot2 elements.
#' @param surface
#' @keywords bob
#' @export
#' @examples
#' geoSurface(surface)
geoSurface <- function(longitude_midpoints=1:ncol(surface), latitude_midpoints=1:nrow(surface), surface) {
    
    image(longitude_midpoints, latitude_midpoints, t(surface)^0.2, col=heat.colors(50))
    
}

#------------------------------------------------
#' Create quick geoprofile plot
#'
#' Creates quick geoprofile plot, choosing some parameters automatically.
#' @param surface
#' @export
#' @examples
#' geoQuickPlot(surface)

geoQuickPlot <- function(params, surface=NULL, data=NULL, zoom="auto", source="google", maptype="hybrid", breakPercent=seq(0,10,l=11), plotContours=TRUE, data_fillCol='black', data_borderCol='white') {
    
    # check that inputs make sense
    geoParamsCheck(params)
    if (!is.null(data))
	    geoDataCheck(data)
    
    # if zoom=="auto" then set zoom level based on params
    if (zoom=="auto")
        zoom <- getZoom(params$output$longitude_minMax, params$output$latitude_minMax)
    
    # make zoom level appropriate to map source
    if (source=="stamen")
    	zoom <- min(zoom,18)
    
    # make map
    rawMap <- get_map(location=c(mean(params$output$longitude_minMax), mean(params$output$latitude_minMax)), zoom=zoom, source=source, maptype=maptype)
    myMap <- ggmap(rawMap) + coord_cartesian(xlim=params$output$longitude_minMax, ylim=params$output$latitude_minMax)
    
    # overlay geoprofile
    if (!is.null(surface)) {
    	geoCols <- colorRampPalette(c("#00008F","#0000FF","#0070FF","#00DFFF","#50FFAF","#BFFF40","#FFCF00","#FF6000","#EF0000","#800000"))
    	
		df <- expand.grid(x=params$output$longitude_midpoints, y=params$output$latitude_midpoints)
		df$z <- as.vector(t(surface))
		labs <- paste(round(breakPercent,1)[-length(breakPercent)],"-",round(breakPercent,1)[-1],"%",sep='')
		df$cut <- cut(df$z, breakPercent/100*length(surface), labels=labs)
		df_noNA <- df[!is.na(df$cut),]
		
		myMap <- myMap + geom_tile(aes(x=x,y=y,fill=cut), alpha=0.6, data=df_noNA)
		myMap <- myMap + scale_fill_manual(name="Hitscore\npercentage", values=rev(geoCols(11)))
		
		# add contours
		if (plotContours) {
			myMap <- myMap + stat_contour(aes(x=x,y=y,z=z), breaks=breakPercent/100*length(surface), size=0.3, alpha=0.5, data=df)
		}
	}

    # overlay data points
    if (!is.null(data)) {
    	p <- data.frame(longitude=data$longitude, latitude=data$latitude)
		myMap <- myMap + geom_point(aes(x=longitude, y=latitude), data=p, cex=1.5, col=data_borderCol)
		myMap <- myMap + geom_point(aes(x=longitude, y=latitude), data=p, pch=20, cex=1.5, col=data_fillCol)
    }
    
    myMap
}