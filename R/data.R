
#------------------------------------------------
#' @title Simulate data
#'
#' @description Simulate data from the same presence-absence model used in the
#'   inference step.
#'
#' @param sentinel_lon vector giving longitudes of sentinel sites
#' @param sentinel_lat vector giving latitudes of sentinel sites
#' @param sentinel_radius observation radius of the sentinel site (km)
#' @param K the number of sources
#' @param source_lon_min minimum limit on source longitudes
#' @param source_lon_max maximum limit on source longitudes
#' @param source_lat_min minimum limit on source latitudes
#' @param source_lat_max maximum limit on source latitudes
#' @param sigma_model set as "single" to use the same dispersal distance for all
#'   sources, or "separate" to use an independently drawn dispersal distance for
#'   each source
#' @param sigma_mean the prior mean of the parameter sigma (km)
#' @param sigma_var the prior variance of the parameter sigma (km). Set to zero 
#'   to use a fixed distance
#' @param expected_popsize the expected total number of observations (observed
#'   and unobserved) in the study area
#'
#' @export
#' @examples
#' # TODO

sim_data <- function(sentinel_lon,
                     sentinel_lat,
                     sentinel_radius = 0.1,
                     K = 3,
                     source_lon_min = -0.2,
                     source_lon_max = 0.0,
                     source_lat_min = 51.45,
                     source_lat_max = 51.55,
                     sigma_model = "single",
                     sigma_mean = 1.0,
                     sigma_var = 0.1,
                     expected_popsize = 100) {
  
  # check inputs
  assert_numeric(sentinel_lon)
  assert_numeric(sentinel_lat)
  assert_single_pos(sentinel_radius, zero_allowed = FALSE)
  assert_single_numeric(source_lon_min)
  assert_single_numeric(source_lon_max)
  assert_single_numeric(source_lat_min)
  assert_single_numeric(source_lat_max)
  assert_same_length(sentinel_lon, sentinel_lat)
  assert_single_pos_int(K, zero_allowed = FALSE)
  assert_single_pos(sigma_mean, zero_allowed = FALSE)
  assert_single_pos(sigma_var, zero_allowed = TRUE)
  assert_single_string(sigma_model)
  assert_in(sigma_model, c("single", "independent"))
  assert_single_pos(expected_popsize, zero_allowed = FALSE)
  
  # draw source locations from prior
  source_lon <- runif(K, source_lon_min, source_lon_max)
  source_lat <- runif(K, source_lat_min, source_lat_max)
  
  # draw total number of points
  N <- rpois(1, expected_popsize)
  
  # draw true allocation of all points to sources
  group <- sort(sample(K, N, replace = TRUE))
  source_N <- tabulate(group)
  
  # draw sigma
  varlog <- log(sigma_var/sigma_mean^2 + 1)
  meanlog <- log(sigma_mean) - varlog/2
  switch(sigma_model,
         "single" = {
           sigma <- rep(rlnorm(1, meanlog, sqrt(varlog)), K)
         },
         "independent" = {
           sigma <- rlnorm(K, meanlog, sqrt(varlog))
         })
  
  # draw points around sources
  df_all <- NULL
  for (k in 1:K) {
    if (source_N[k]>0) {
      rand_k <- rnorm_sphere(source_N[k], source_lon[k], source_lat[k], sigma[k])
      df_all <- rbind(df_all, as.data.frame(rand_k))
    }
  }
  
  # assign points as observed or unobserved based on distance to sentinel sites
  gc_dist <- mapply(function(x, y) {
                      lonlat_to_bearing(x, y, df_all$longitude, df_all$latitude)$gc_dist
                    }, x = sentinel_lon, y = sentinel_lat)
  counts <- colSums(gc_dist < sentinel_radius)
  df_observed <- data.frame(longitude = sentinel_lon,
                            latitude = sentinel_lat,
                            counts = counts)
  
  # add record of whether data point is observed or unobserved to df_all. Value
  # gives the number of sentinel sites that observe the data
  df_all$observed <- rowSums(gc_dist < sentinel_radius)
  
  # return simulated data and true parameter values
  ret <- list()
  ret$sentinel_radius <- sentinel_radius
  ret$source <- data.frame(longitude = source_lon,
                              latitude = source_lat)
  ret$source_N <- source_N
  ret$group <- group
  ret$sigma <- sigma
  ret$data_all <- df_all
  ret$data_observed <- df_observed
  
  # make custom class
  class(ret) <- "rgeoprofile_simdata"
  
  return(ret)
}

#------------------------------------------------
# simulate bi-allelic data
#' @noRd
sim_data_biallelic <- function(n, L, K, true_group, true_p = true_p, true_m, locus_names, samp_names, e1, e2, prop_missing, pop_col_on) {
  
  # simulate raw data
  dat <- NULL
  for (k in 1:K) {
    if (any(true_group == k)) {
      
      # draw raw numbers of REF allele for individuals in this deme
      true_m_k <- true_m[true_group == k]
      true_p_k <- mapply(function(x) {x[k,1]}, true_p)
      dat_raw_k <- t(mapply(rbinom, n = L, size = true_m_k, MoreArgs = list(prob = true_p_k)))
      
      # convert to matrix of {0.0, 0.5, 1.0}
      dat_k <- matrix(0.5, length(true_m_k), L)
      dat_k[dat_raw_k == matrix(true_m_k, length(true_m_k), L)] <- 1
      dat_k[dat_raw_k == 0] <- 0
      
      # append data
      dat <- rbind(dat, dat_k)
    }
  }
  colnames(dat) <- locus_names
  
  # add errors and missing data
  dat_uncorrupted <- NULL
  if (e1>0 || e2>0 || prop_missing>0) {
    dat_uncorrupted <- dat
    
    # error1 - homo missclassified as het
    if (e1 > 0) {
      homo1 <- dat[dat_uncorrupted == 1]
      homo1[runif(length(homo1)) < e1] <- 0.5
      dat[dat_uncorrupted == 1] <- homo1
      
      homo2 <- dat[dat_uncorrupted == 0]
      homo2[runif(length(homo2)) < e1] <- 0.5
      dat[dat_uncorrupted == 0] <- homo2
    }
    
    # error2 - het missclassified as homo
    if (e2 > 0) {
      het <- dat[dat_uncorrupted == 0.5]
      rand1 <- (runif(length(het)) < e2)
      if (any(rand1)) {
        het[rand1] <- sample(c(0,1), sum(rand1), replace = TRUE)
      }
      dat[dat_uncorrupted == 0.5] <- het
    }
    
    # missing data
    if (prop_missing > 0) {
      prop_missing_round <- round(prop_missing*n*L)
      dat[sample.int(n*L, prop_missing_round )] <- -9
    }
  }
  
  # convert dat and dat_uncorrupted to dataframe
  df <- data.frame(sample_ID = samp_names, stringsAsFactors = FALSE)
  rownames(df) <- NULL
  if (pop_col_on) {
    df$pop <- true_group
  }
  df_uncorrupted <- NULL
  if (!is.null(dat_uncorrupted)) {
    df_uncorrupted <- cbind(df, dat_uncorrupted)
  }
  df <- cbind(df, dat)
  
  # return list
  return(list(df = df, df_uncorrupted = df_uncorrupted))
}

#------------------------------------------------
# simulate multi-allelic data
#' @noRd
sim_data_multiallelic <- function(n, L, K, true_group, true_p = true_p, true_m, locus_names, samp_names, e1, e2, prop_missing, pop_col_on) {
  
  # simulate raw data
  df <- NULL
  for (l in 1:L) {
    true_p_group <- lapply(true_group, function(i) {true_p[[l]][i,]})
    haplotypes <- mapply(function(x,y) {
      sort(unique(sample.int(length(x), y, replace = TRUE, prob = x)))
    }, true_p_group, y = true_m)
    df_l <- data.frame(sample_ID = rep(samp_names, times = sapply(haplotypes,length)), locus = l, haplotype = unlist(haplotypes), stringsAsFactors = FALSE)
    df <- rbind(df, df_l)
  }
  df <- df[order(df$sample),]
  row.names(df) <- NULL
  
  # add missing data
  df_uncorrupted <- NULL
  if (prop_missing > 0) {
    df_uncorrupted <- df
    
    prop_missing_round <- round(prop_missing*n*L)
    missing_index <- expand.grid(unique(df$sample_ID), 1:L, -9)[sample.int(n*L, prop_missing_round, replace = FALSE),]
    names(missing_index) <- c("sample_ID", "locus", "haplotype")
    df <- subset(df, !( paste(df$sample_ID, df$locus, sep=".") %in% paste(missing_index$sample_ID, missing_index$locus, sep=".") ))
    df <- rbind(df, missing_index)
    df <- df[order(df$locus),]
    df <- df[order(df$sample),]
    row.names(df) <- NULL
  }
  
  # return list
  return(list(df = df, df_uncorrupted = df_uncorrupted))
}

#------------------------------------------------
#' @title Simulate genetic data subject to constraints
#'
#' @description TODO - text
#'
#' @details TODO
#'
#' @param ... TODO
#' @param data_format TODO
#' @param no_invariant_loci TODO
#' @param no_missing_samples TODO
#' @param no_missing_loci TODO
#' @param max_attempts TODO
#'
#' @export
#' @examples
#' # TODO

sim_data_safe <- function(..., data_format = "biallelic", no_invariant_loci = TRUE, no_missing_samples = TRUE, no_missing_loci = TRUE, max_attempts = 1e3) {

  # attempt to simulate satisfactory data a finite number of times
  for (i in 1:max_attempts) {

    # simulate data
    sim1 <- sim_data(..., data_format = data_format)
    data <- sim1$data
    n <- sim1$n
    L <- sim1$L

    # bi-allelic data
    if (data_format=="biallelic") {

      # check no invariant loci
      if (no_invariant_loci & any(colSums(data==1 | data==0)==nrow(data)) ) {
        next
      }

      # check no missing samples
      if ( no_missing_samples & any(colSums(data==-1)==nrow(data)) ) {
        next
      }

      # check no missing loci
      if ( no_missing_loci & any(rowSums(data==-1)==ncol(data)) ) {
        next
      }
    }

    # multi-allelic data
    if (data_format=="multiallelic") {

      # check no invariant loci
      n_haplotypes <- mapply(function(i) {
        s <- data$haplotype[data$locus==i]
        length(unique(s[s>0]))
        }, 1:L)
      if (no_invariant_loci & any(n_haplotypes==1)) {
        next
      }

      # check no missing samples
      n_nonmissing_samples <- mapply(function(i){
        sum(data$sample==i & data$haplotype>0)
        }, 1:n)
      if (no_missing_samples & any(n_nonmissing_samples==0)) {
        next
      }

      # check no missing loci
      n_nonmissing_loci <- mapply(function(i){
        sum(data$locus==i & data$haplotype>0)
        }, 1:L)
      if ( no_missing_loci & any(n_nonmissing_loci==0) ) {
        next
      }
    }

    # if made it to here then data passed all checks
    return(sim1)
  }

  stop(paste("Unable to produce data set satisfying constraints within", max_attempts, "random draws"))
}
