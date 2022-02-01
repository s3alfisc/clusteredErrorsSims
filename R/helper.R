# Code slightly adapted from code written by Gordon Burtch, published with MIT
# license. The original code can be found here https://github.com/gburtch/simulating_cluster_SEs

# Content: Helper functions for clustered regressions simulations.

### FUNCTION TO SIMULATE A DATA FRAME COMPRISED OF CLUSTERED DATA

# This function simulates clustered data.
# n = observations.
# n_cluster = number of clusters they are associated with.
# rho = intra-cluster correlation.
# param is a set of two values that reflect our regression beta coefficients.

# This function is not perfect, by any means.
# The shared cluster contribution to variance in x and e leads to a correlation between x and e,
# which can bias estimates?
gen_cluster <-
  function(param = c(.1, .5),
           n = 1000,
           n_cluster = 50,
           rho = .5,
           balanced_cluster = TRUE,
           did = FALSE,
           treatment_state = NULL) {
    # Required package: mvtnorm

    while (n %% n_cluster != 0){
      n <- n+1
    }

    # Here we are simulating two vectors, from a bivariate normal distribution.
    # These two variables will serve as our variable, X, and the error term.
    # The two vectors (variables) will have means of 0, and will not covary, i.e., we don't have endogeneity here!

    # NOTE about how rho works in this setup:

    # The simulation will essentially make the variance (sd) of x and e will both depend on cluster membership.
    # However, notice that x and e are simulated to be uncorrelated with one another - we do not have OVB problem.

    # The individual level contribution to the error will have a variance of 1-rho.
    # So, if rho = 1, the individual level draw essentially contributes nothing to X.
    # As you will see below, cluster contribution will drive all the variation in the observed variable in that situation.
    Sigma_i <- matrix(c(1, 0, 0, 1 - rho), ncol = 2)
    values_i <- rmvnorm(n = n, sigma = Sigma_i)
    x_i <- values_i[, 1]
    error_i <- values_i[, 2]

    # Now we are going to assign the above datapoints to clusters, evenly.
    if(balanced_cluster == TRUE){
      cluster_name <- rep(1:n_cluster, each = n / n_cluster)
    } else {
      state_pop <- get_state_freq()
      if(length(unique(state_pop$state)) %% n_cluster != 0){
        stop("If using 'balanced_cluster = FALSE', n_cluster must be divisible by 50.")
      }
      cluster_name <- sample(x = 1:length(state_pop$state),
                             size = n,
                             replace = TRUE,
                             prob = state_pop$population_share)
    }

    # For X and error, simulated above, we assign each observation to a cluster.
    # The assigned cluster will then shift each observation's pair X and error by a pair of values drawn from another bivariate normal distribution.
    # The cluster-specific contribution to error will have a variance of rho, again just because?
    # We are now simulating distributional draws that will serve as the cluster-specific contributions to X and error.
    Sigma_cl <- matrix(c(1, 0, 0, rho), ncol = 2)
    values_cl <- rmvnorm(n = n_cluster, sigma = Sigma_cl)
    x_cl <- rep(values_cl[, 1], each = n / n_cluster)
    error_cl <- rep(values_cl[, 2], each = n / n_cluster)

    # Finally, we simulate y as a function of our betas, and the X and error we just simulated.

    if(did == FALSE){
      # Now we stick together the individual level values and their cluster-specific contributions.
      x <- x_i + x_cl
      error <- error_i + error_cl
      y <- param[1] + param[2] * x + error
    } else if(did == TRUE){
      # dummy variable for treatment_state, "no state level"
      x <- ifelse(cluster_name %in% treatment_states, 1, 0)
      y <- param[1] + param[2] * x + error
    }

    df_tmp <-
      data.frame(x, y, cluster = cluster_name, x_i, x_cl, error, error_i, error_cl)
  return(df_tmp)
}

### SINGLE ITERATION OF THE SIM: CREATE A SAMPLE, ESTIMATE MODEL, STORE RESULTS

# This function nests the data generation function.
# For each simulation iteration, simulate clustered data with the provided params.
# Then estimate the linear regression model, store the beta, and the standard errors + confidence intervals.
# Calculate the inference naively, with regular clustering, or fast and wild cluster bootstrap (good for few clusters).
cluster_sim <- function(param = c(.1, .5),
                        n = 1000,
                        n_cluster = 50,
                        rho = .5,
                        inference = "regular",
                        FE = FALSE,
                        balanced_cluster = TRUE,
                        boot = "fwildclusterboot",
                        did = FALSE) {
  # Required packages: multiwayvcov, fwildclusterboot

  # Note we have to do some weirdness here because of the way the boottest function is coded.
  # We need df_it to be a global variable, because boottest looks for it implicitly in global memory.
  df_it <<-
    gen_cluster(
      param = param,
      n = n ,
      n_cluster = n_cluster,
      rho = rho,
      balanced_cluster = balanced_cluster,
      did = did
    )

  x_err_cor <- cor(df_it$x,df_it$error)

  if (FE == TRUE) {
    fit <- lm(data = df_it, y ~ x + factor(cluster))
  } else {
    fit <- lm(data = df_it, y ~ x)
  }

  b1 <- coef(fit)[2]

  if (inference == "regular") {
    Sigma <- stats::vcov(fit)
    b1_ci95 <- stats::confint(fit)[2, ]
  } else if (inference == "clustered") {
    Sigma <- suppressMessages(clubSandwich::vcovCR(fit,
                                  cluster = df_it$cluster,
                                  type = "CR1S"))
    conf_int <- suppressMessages(clubSandwich::conf_int(obj = fit,
                                       vcov = Sigma ,
                                       test = "z",
                                       level = 0.95,
                                       coefs = "x"))
    lower <- conf_int$CI_L
    upper <- conf_int$CI_U
    b1_ci95 <- c(lower, upper)
  } else if (inference == "fast_n_wild") {
    # Per Roodman et al. 2019, employ Webb weights if cluster count <= 12.
    if (n_cluster <= 12) {
      weights = "webb"
    } else {
      weights = "rademacher"
    }
    if(boot == "fwildclusterboot"){
      suppressWarnings(boot <-
                         fwildclusterboot::boottest(
                           fit,
                           B = 9999,
                           param = "x",
                           clustid = "cluster",
                           type = weights
                         ))
    } else if (boot == "wildboottestjlr"){
      suppressWarnings(boot <-
                         wildboottestjlr::boottest(
                           fit,
                           B = 9999,
                           param = "x",
                           clustid = "cluster",
                           type = weights
                         ))
    }

    lower <- boot$conf_int[1]
    upper <- boot$conf_int[2]
    b1_ci95 <- c(lower, upper)
  } else if(inference == "Satterthwaite"){
    Sigma <- suppressMessages(clubSandwich::vcovCR(fit,
                                  cluster = df_it$cluster,
                                  type = "CR2"))
    conf_int <- suppressMessages(clubSandwich::conf_int(obj = fit,
                                       vcov = Sigma ,
                                       test = "Satterthwaite",
                                       level = 0.95,
                                       coefs = "x"))
    lower <- conf_int$CI_L
    upper <- conf_int$CI_U
    b1_ci95 <- c(lower, upper)
  } else if(inference == "saddlepoint"){
    Sigma <- suppressMessages(clubSandwich::vcovCR(fit,
                                  cluster = df_it$cluster,
                                  type = "CR2"))
    conf_int <- suppressMessages(clubSandwich::conf_int(obj = fit,
                                       vcov = Sigma ,
                                       test = "saddlepoint",
                                       level = 0.95,
                                       coefs = "x"))
    lower <- conf_int$CI_L
    upper <- conf_int$CI_U
    b1_ci95 <- c(lower, upper)
  }

  return(c(b1, b1_ci95, x_err_cor))
}

# Lastly, a function to iterate / repeat the simulation.
# A data frame is returned containing all the betas, se's, CIs, etc.
# We are also creating, for each sim, a binary indicator of whether the CIs contain the true parameter for b1.
run_cluster_sim <-
  function(n_sims = 1000,
           param = c(.1, .5),
           n = 1000,
           n_cluster = 50,
           rho = .5,
           inference = "regular",
           FE = FALSE,
           balanced_cluster = TRUE,
           workers = 8,
           did = FALSE,
           treatment_states = NULL) {

    #'@import future
    # Required packages: dplyr

    cat(paste0("inference: ", inference), "\n")
    #if(inference != "fast_n_wild"){
      future::plan(multisession, workers = workers)
      progressr::handlers("progress")

      progressr::with_progress({
        p <- progressor(along = 1:n_sims)
        df_res <-
          future.apply::future_lapply(X = 1:n_sims,
                                      future.seed = sample(1:10000000, 1),
                                      FUN = function(x){
                                         p(sprintf("x=%g", x))
                                         cluster_sim(param = param,
                                                     n = n,
                                                     rho = rho,
                                                     n_cluster = n_cluster,
                                                     inference = inference,
                                                     FE = FE)

                                      })

      })
    #} else if(inference == "fast_n_wild") {
      # future.apply runs into problems if tryCatch happens to catch an error
      # df_res <-
      #           replicate(n_sims,
      #                     cluster_sim(param = param,
      #                                 n = n,
      #                                 rho = rho,
      #                                 n_cluster = n_cluster,
      #                                 inference = inference,
      #                                 FE = FE)
      #                     )
    #}

    df_res <- as.data.frame(Reduce(rbind, df_res))
    #df_res <- as.data.frame(t(df_res))
    names(df_res) <- c('b1', 'ci95_lower', 'ci95_upper','x_err_cor')
    df_res <- df_res %>%
      mutate(id = 1:n(),
             param_caught = ci95_lower <= param[2] &
               ci95_upper >= param[2])
    return(df_res)
  }

# source https://en.wikipedia.org/wiki/List_of_U.S._states_and_territories_by_population

get_state_freq <- function(){

  # source https://en.wikipedia.org/wiki/List_of_U.S._states_and_territories_by_population

  state_population <-
    c("California",	39538223,
      "Texas",	29145505,
      "Florida",	21538187,
      "New York",	20201249,
      "Pennsylvania",	13002700,
      "Illinois",	12812508,
      "Ohio",	11799448,
      "Georgia",	10711908,
      "North Carolina",	10439388,
      "Michigan",	10077331,
      "New Jersey",	9288994,
      "Virginia",	8631393,
      "Washington",	7705281,
      "Arizona",	7151502,
      "Massachusetts",	7029917,
      "Tennessee",	6910840,
      "Indiana",	6785528,
      "Maryland",	6177224,
      "Missouri",	6154913,
      "Wisconsin",	5893718,
      "Colorado",	5773714,
      "Minnesota",	5706494,
      "South Carolina",	5118425,
      "Alabama",	5024279,
      "Louisiana",	4657757,
      "Kentucky",	4505836,
      "Oregon",	4237256,
      "Oklahoma",	3959353,
      "Connecticut",	3605944,
      "Puerto Rico",	3285874,
      "Utah",	3271616,
      "Iowa",	3190369,
      "Nevada",	3104614,
      "Arkansas",	3011524,
      "Mississippi",	2961279,
      "Kansas",	2937880,
      "New Mexico",	2117522,
      "Nebraska",	1961504,
      "Idaho",	1839106,
      "West Virginia",	1793716,
      "Hawaii",	1455271,
      "New Hampshire",	1377529,
      "Maine",	1362359,
      "Rhode Island",	1097379,
      "Montana",	1084225,
      "Delaware",	989948,
      "South Dakota",	886667,
      "North Dakota",	779094,
      "Alaska",	733391,
      "DC",	689545,
      "Vermont", 643077,
      "Wyoming", 576851
    )

  dt = data.table(state = state_population[seq(1, length(state_population), 2)],
                  population = as.numeric(state_population[seq(2, length(state_population), 2)]))
  dt <- dt[!(state %in% c("Puerto Rico", "DC"))]
  total_population <- sum(dt$population)
  dt[, population_share := population / total_population]

  dt[, list(state, population_share)]
}


