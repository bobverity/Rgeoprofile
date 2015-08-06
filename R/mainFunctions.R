
#*------------------------------------------------*
#' Draw from Dirichlet process mixture model
#'
#' This function provides draws from a 2D Dirichlet process mixture model.
#' @param n The number of draws
#' @keywords bob
#' @export
#' @examples
#' rDPM(10)
rDPM <- function(n, sigma=1,priorMean_x=0, priorMean_y=0, priorVar=10, alpha=1) {
    
    # draw grouping
	group = rep(1,n)
	freqs = 1
	for (i in 2:n) {
		group[i] = sample(length(freqs)+1,1,prob=c(freqs,alpha))
		if (group[i]>length(freqs))
        freqs = c(freqs,0)
		freqs[group[i]] = freqs[group[i]] + 1
	}
	group = sort(group)
	
	# draw means and data
	mu_x = rnorm(length(freqs),priorMean_x,sqrt(priorVar))
	mu_y = rnorm(length(freqs),priorMean_y,sqrt(priorVar))
	x = rnorm(n,mu_x[group],sigma)
	y = rnorm(n,mu_y[group],sigma)
	
	# return results
	return(list(x=x, y=y, group=group, mu_x=mu_x, mu_y=mu_y))
}

#*------------------------------------------------*
#' Create Rgeoprofile parameters object
#'
#' This function can be used to generate parameters in the format used by the Rgeoprofile MCMC. Particular parameter values can be specified as input arguments.
#' @param sigma The standard deviation of the dispersal distribution
#' @param priorMean_longitude The position (longitude) of the mean of the prior distribution
#' @param priorMean_latitude The position (latitude) of the mean of the prior distribution
#' @keywords bob
#' @export
#' @examples
#' geoParams()
geoParams <- function(sigma=1, priorMean_longitude=0, priorMean_latitude=0, priorVar=10, alpha_shape=0.1, alpha_rate=0.1, chains=10, burnin=500, samples=5000, burnin_printConsole=100, samples_printConsole=1000, longitude_minMax=c(-10,10), latitude_minMax=c(-10,10), longitude_cells=50, latitude_cells=50) {
    
    # set model parameters
    model = list(sigma=sigma, priorMean_longitude=priorMean_longitude, priorMean_latitude=priorMean_latitude, priorVar=priorVar, alpha_shape=alpha_shape, alpha_rate=alpha_rate)
    
    # set MCMC parameters
    MCMC = list(chains=chains, burnin=burnin, samples=samples, burnin_printConsole=burnin_printConsole, samples_printConsole=samples_printConsole)
    
    # set output parameters
    output = list(longitude_minMax=longitude_minMax, latitude_minMax=latitude_minMax, longitude_cells=longitude_cells, latitude_cells=latitude_cells)
    
    # combine and return
    params = list(model=model, MCMC=MCMC, output=output)
    return(params)
}

#*------------------------------------------------*
#' Create Rgeoprofile data object
#'
#' This function can be used to generate a dummy data object in the format used by the Rgeoprofile MCMC.
#' @param longitude The longitudinal positions of the observed data
#' @param latitude The latitudinal positions of the observed data
#' @keywords bob
#' @export
#' @examples
#' geoData()
geoData <- function(longitude=NULL, latitude=NULL) {
    
    # generate dummy data if none read in
    if (is.null(longitude) & is.null(latitude)) {
        longitude = c(1.3552866, 2.8332016, 1.5458939, 3.4439858, 0.1850686, 2.1763534, 0.6144457, 2.6072376, 1.4571848, 4.5698746, 4.0574786, 4.0210364, 4.0359571, 4.0343223, 4.0614972, 4.6414387, 2.5414649, -1.3519314, -2.2367752, -1.9417826, -2.2836110, -2.0883225, -1.7118926, -0.7698030, 5.3877648, 5.6224173, -16.3016122, -16.1339459, -14.6866445, -15.9828965, -15.7719747, -16.0328715, -14.4931409, -16.3640593, 0.8146278, -1.4039898, 3.5071933, 5.2664047, 4.8191548, 5.2940034, 5.6423132, 5.8105845, 6.0036396, 13.8852696, 12.8413242, -1.0211645, -1.9139261, 0.5087536, 17.7389148, 1.5007123)
        latitude = c(11.5820227, 11.4602990, 10.0136256, 12.4405047, 11.7671286, 8.4943908, 11.2519649, 10.5894781, 12.3373432, -3.4205490, -3.0006835, -3.4391243, -1.6182644, -2.6239842, -3.0284404, -2.2300918, 12.4177397, -4.1720895, -3.2941502, -3.7335080, -2.9844355, -1.9521769, -2.6494257, -2.0376150, -8.2084885, -8.9080376, -0.4101773, -1.8466111, -2.6371649, -9.3363200, -8.8889480, -7.7863900, -10.9563316, -10.7300902, 19.7024929, 18.9879279, -0.8767754, -0.7894751, -1.5367136, -0.1503125, -1.5041463, -1.5992631, -3.1876078, 1.5872134, 2.2382692, -0.2236694, 2.1864671, 10.6193734, 6.1884289, -12.5835766)
    }
    
    # combine and return
    data = list(longitude=longitude, latitude=latitude)
    return(data)
}

#*------------------------------------------------*
#' Check parameters
#'
#' Check that all parameters for use in Rgeoprofile MCMC are OK.
#' @param params A list of parameters
#' @keywords bob
#' @export
#' @examples
#' # tester
#' geoParamsCheck(myParams)
geoParamsCheck <- function(params) {
    
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
    if (!("priorVar"%in%names(params$model)))
        stop("params$model must contain parameter 'priorVar'")
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
    if (!is.numeric(params$model$priorVar) | !is.finite(params$model$priorVar))
        stop("params$model$priorVar must be numeric and finite")
    if (params$model$priorVar<=0)
        stop("params$model$priorVar must be greater than 0")
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
    cat("params file passed all checks\n")
}

#*------------------------------------------------*
#' Check data
#'
#' Check that all data for use in Rgeoprofile MCMC is in the correct format.
#' @param data A list of data
#' @keywords bob
#' @export
#' @examples
#' # tester
#' geoDataCheck(myData)
geoDataCheck <- function(data) {
    
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
    cat("data file passed all checks\n")
}

#*------------------------------------------------*
#' MCMC under Rgeoprofile model
#'
#' This function carries out the main MCMC under the Rgeoprofile model.
#' @param input A list of input parameters.
#' @keywords bob
#' @export
#' @examples
#' geoMCMC()
geoMCMC <- function(data, params) {
    
    # check that data and parameters in correct format
    geoDataCheck(data)
    geoParamsCheck(params)
    cat("\n")
    
    # extract some values from params
    n = length(data$longitude)
    sigma = params$model$sigma
    priorMean_longitude = params$model$priorMean_longitude
    priorMean_latitude = params$model$priorMean_latitude
    priorVar = params$model$priorVar
    samples = params$MCMC$samples
    longitude_cells = params$output$longitude_cells
    latitude_cells = params$output$latitude_cells
    
    # carry out MCMC
    rawOutput = .Call('RgeoProfile_C_geoMCMC', PACKAGE = 'RgeoProfile', data, params)
    
    # produce prior matrix
    xmin = params$output$longitude_minMax[1]
    xmax = params$output$longitude_minMax[2]
    ymin = params$output$latitude_minMax[1]
    ymax = params$output$latitude_minMax[2]
    xCells = longitude_cells
    yCells = latitude_cells
    xCellSize = (xmax-xmin)/xCells
    yCellSize = (ymax-ymin)/yCells
    xMids = xmin - xCellSize/2 + (1:xCells)*xCellSize
    yMids = ymin - yCellSize/2 + (1:yCells)*yCellSize
    xMids_mat = outer(rep(1,yCells),xMids)
    yMids_mat = outer(yMids,rep(1,xCells))
    priorMat = dnorm(xMids_mat,priorMean_longitude,sd=sqrt(priorVar))*dnorm(yMids_mat,priorMean_latitude,sd=sqrt(priorVar))
    
    
    # finalise output format
    output = list()
    
    alpha = rawOutput$alpha
    output$alpha = alpha
    
    allocation = matrix(unlist(rawOutput$allocation),n,byrow=T)
    allocation = data.frame(allocation/samples)
    names(allocation) = paste("group",1:ncol(allocation),sep="")
    output$allocation = allocation
    
    output$longitude_midpoints = xMids
    output$latitude_midpoints = yMids
    
    surface = priorMat*mean(alpha/(alpha+n)) + matrix(unlist(rawOutput$geoSurface),latitude_cells,byrow=TRUE)/(xCellSize*yCellSize)/samples
    surface = surface/sum(surface)
    output$surface = surface
    
    return(output)
}

#*------------------------------------------------*
#' Apply smoothing to surface
#'
#' Applies Gaussian smoothing to surface.
#' @param z matrix on which to apply smoothing.
#' @param bandwidth standard deviation of Gaussian smoothing kernel
#' @keywords bob
#' @export
#' @examples
#' geoSmooth(MCMCoutput$geoSurface)
geoSmooth <- function(x, y, z, bandwidth) {
    
    # check that input arguments are finite and numeric
    if (!all(is.finite(x)) | !all(is.numeric(x)))
        stop("values in x must be finite and numeric")
    if (!all(is.finite(y)) | !all(is.numeric(y)))
        stop("values in y must be finite and numeric")
    if (!all(is.finite(z)) | !all(is.numeric(z)))
        stop("values in z must be finite and numeric")
    if (!is.finite(bandwidth) | !is.numeric(bandwidth))
        stop("bandwidth must be finite and numeric")
    
    # check that input arguments make sense
    if (!is.vector(x) | !is.vector(y))
        stop("x and y must be vectors")
    if (!is.matrix(z))
        stop("z must be a matrix")
    if (length(x)!=ncol(z) | length(y)!=nrow(z))
        stop("length of vectors x and y must correspond to columns and rows of matrix z, respectively")
    
    # padd x and y to avoid bleeding over edges (Fourier method assumes periodic domain)
    diff_x = x[2]-x[1]
    x_pad = ceiling(3*bandwidth/diff_x)
    x_left = seq(x[1]-x_pad*diff_x, x[1]-diff_x, diff_x)
    x_right = seq(x[length(x)]+diff_x, x[length(x)]+x_pad*diff_x, diff_x)
    x = c(x_left,x,x_right)
    
    z_leftRight = matrix(0,length(y),x_pad)
    z = cbind(z_leftRight,z,z_leftRight)
    
    diff_y = y[2]-y[1]
    y_pad = ceiling(3*bandwidth/diff_y)
    y_left = seq(y[1]-y_pad*diff_y, y[1]-diff_y, diff_y)
    y_right = seq(y[length(y)]+diff_y, y[length(y)]+y_pad*diff_y, diff_y)
    y = c(y_left,y,y_right)
    
    z_topBottom = matrix(0,y_pad,length(x))
    z = rbind(z_topBottom,z,z_topBottom)
    
    # define distance matrices
    dist_x = x-x[1]
    dist_x_rev = -(x-x[length(x)])
    dist_x[dist_x_rev<dist_x] = dist_x_rev[dist_x_rev<dist_x]
    dist_x_mat = t(replicate(length(y),dist_x))
    dist_y = y-y[1]
    dist_y_rev = -(y-y[length(y)])
    dist_y[dist_y_rev<dist_y] = dist_y_rev[dist_y_rev<dist_y]
    dist_y_mat = replicate(length(x),dist_y)
    
    # define Gaussian kernel
    k = dnorm(dist_x_mat,sd=bandwidth)*dnorm(dist_y_mat,sd=bandwidth)
    
    # apply Fourier method
    f1 = fftw2d(z)
    f2 = fftw2d(k)
    f3 = f1*f2
    f4 = Re(fftw2d(f3,inverse=T))
    f4[f4<0] = 0
    
    # trim off padding
    f4 = f4[(y_pad+1):(length(y)-y_pad),(x_pad+1):(length(x)-x_pad)]
    
    return(f4)
}

#*------------------------------------------------*
#' Calculate geoprofile from surface
#'
#' Converts surface to rank order geoprofile.
#' @param z matrix to convert to geoprofile
#' @keywords bob
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
    z[order(z,decreasing=TRUE)] = 1:length(z)
    
    return(z)
}

#*------------------------------------------------*
#' Plot posterior allocation
#'
#' Produces plot of posterior allocation, given output of MCMC.
#' @param MCMCoutput output generated from an MCMC.
#' @keywords bob
#' @export
#' @examples
#' geoPlotAllocation(MCMCoutput)
geoPlotAllocation <- function(allocation, colours, barBorderCol="white", barBorderWidth=0.25, mainBorderCol="black", mainBorderWidth=2, yTicks_on=TRUE, yTicks=seq(0,1,0.2), xlab="", ylab="posterior allocation", mainTitle="", names=NA, names_size=1, xTicks_on=FALSE, xTicks_size=1, orderBy="group") {
    
    # check that orderBy is either 'group' or 'probability'
    if (!(orderBy%in%c("group","probability")))
        stop("orderBy must equal 'group' or 'probability'")
    
    # check that allocation is a data frame
    if (!is.data.frame(allocation))
        stop("allocation must be a data frame, with observations in rows and groups in columns")
    
    n = nrow(allocation)
    k = ncol(allocation)
    
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
        tM = t(MCMCoutput$allocation)
        for (i in 1:n) {
            temp = tM
            temp[,-i] = NA
            temp_order = order(temp[,i],decreasing=TRUE)
            barplot(temp[temp_order,], col=colours[temp_order], space=0, border=NA, axes=FALSE, add=TRUE)
        }
        segments(0:n,c(0,0),0:n,c(1,1), col=barBorderCol, lwd=barBorderWidth)
        box(col=mainBorderCol, lwd=mainBorderWidth)
        axis(2, tick=yTicks_on, labels=yTicks_on, at=yTicks)
        axis(1, at=1:n-0.5, tick=xTicks_on, lwd.ticks=xTicks_size, labels=names, las=2, cex.axis=names_size)
    }
    
}

#*------------------------------------------------*
#' Plot geoprofile
#'
#' Produces plot of posterior surface or geoprofile
#' @param surface
#' @keywords bob
#' @export
#' @examples
#' geoPlotSurface(surface)
geoPlotSurface <- function(longitude_midpoints=1:ncol(surface), latitude_midpoints=1:nrow(surface), surface) {
    
    image(longitude_midpoints,latitude_midpoints,t(surface)^0.2,col=heat.colors(50))
    
}
