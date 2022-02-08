wildly_different_sim <- function(n = 2000, workers = 8, n_sims = 50000, did = FALSE, treatment_states = NULL){

  # ----------------------------------------------------------------------
  # Wildly Different Cluster Sizes (US states relative sizes)
  #'@import   arm mvtnorm ggplot2 dplyr data.table RColorBrewer ggthemes fishmethods fwildclusterboot wildboottestjlr patchwork clubSandwich future.apply progressr
  #'@export

  res_wildly_different <-

    lapply(seq(0, 0.9, 0.1), function(x){

      cat(paste("rho = ", x), "\n")
      # Step 1: how does "regular" inference do for c(50, 100) with wildly different sizes?
      # mean_CI_coverage_reg <- rep(0,2)
      # sd_CI_coverage_reg <- rep(0,2)
      # i <- 1
      # for (j in c(50, 100)) {
      #   sim_clustered_fix_tmp <-
      #     run_cluster_sim(n_sims = n_sims,
      #                     n = n,
      #                     rho = x,
      #                     n_cluster = j,
      #                     inference = "regular",
      #                     workers = workers,
      #                     balanced_cluster = FALSE)
      #   mean_CI_coverage_reg[i] <- mean(sim_clustered_fix_tmp$param_caught)
      #   sd_CI_coverage_reg[i] <- sd(sim_clustered_fix_tmp$param_caught)
      #   i <- i + 1
      # }

      # regular_SE_by_clustN <- data.frame(coverage = mean_CI_coverage_reg,
      #                                    #sd_CI_coverage_reg,
      #                                    clusters = c(50, 100),
      #                                    type = "regular",
      #                                    rho = x)



      # Step 2: how does "cluster robust/sandwich" inference do for c(50, 100) with wildly different sizes?

      mean_CI_coverage_reg <- rep(0,2)
      sd_CI_coverage_reg <- rep(0,2)
      i <- 1
      for (j in c(50, 100)) {
        sim_clustered_fix_tmp <-
          run_cluster_sim(n_sims = n_sims,
                          n = n,
                          rho = x,
                          n_cluster = j,
                          inference = "clustered",
                          workers = workers,
                          balanced_cluster = FALSE)
        mean_CI_coverage_reg[i] <- mean(sim_clustered_fix_tmp$param_caught)
        sd_CI_coverage_reg[i] <- sd(sim_clustered_fix_tmp$param_caught)
        i <- i + 1
      }

      cluster_SE_by_clustN <- data.frame(coverage = mean_CI_coverage_reg,
                                         #sd_CI_coverage_reg,
                                         clusters = c(50, 100),
                                         type = "clustered",
                                         rho = x)



      # Step 3: how does "wild cluster robust" inference do for c(50, 100) with wildly different sizes?

      mean_CI_coverage_fw <- rep(0,2)
      sd_CI_coverage_fw <- rep(0,2)
      i <- 1
      for (j in c(50, 100)) {
        tryCatch(
          {
            sim_clustered_fix_tmp <-
              run_cluster_sim(n_sims = n_sims,
                              n = n,
                              rho = x,
                              n_cluster = j,
                              inference = "fast_n_wild",
                              workers = workers,
                              boot = "fwildclusterboot")
            mean_CI_coverage_fw[i] <- mean(sim_clustered_fix_tmp$param_caught)
            sd_CI_coverage_fw[i] <- sd(sim_clustered_fix_tmp$param_caught)
            i <- i + 1
          },
          error=function(cond) {
            message(paste("\nFast'n'Wild Failed.\n"))
            message("Here's the original error message:")
            message(cond)
            # Choose a return value in case of error
            return(NA)
          })

      }

      cluster_SE_by_clustN_fw <- data.frame(coverage = mean_CI_coverage_fw,
                                            #sd_CI_coverage_fw,
                                            clusters = c(50, 100),
                                            type = "fast_n_wild",
                                            rho = x)



      # Step 4: how do satterthwaite corrections do for c(50, 100) clusters?

      mean_CI_coverage <- rep(0,2)
      sd_CI_coverage <- rep(0,2)
      i <- 1
      for (j in c(50, 100)) {
        sim_clustered_fix_tmp <-
          run_cluster_sim(n_sims = n_sims,
                          n = n,
                          rho = x,
                          n_cluster = j,
                          inference = "Satterthwaite",
                          workers = workers)
        mean_CI_coverage[i] <- mean(sim_clustered_fix_tmp$param_caught)
        sd_CI_coverage[i] <- sd(sim_clustered_fix_tmp$param_caught)
        i <- i + 1
      }

      cluster_SE_by_clustN_satterthwaite <- data.frame(coverage = mean_CI_coverage,
                                                       #sd_CI_coverage,
                                                       clusters = c(50, 100),
                                                       type = "Satterthwaite",
                                                       rho = x)


      res <- rbind(#regular_SE_by_clustN,
                   cluster_SE_by_clustN,
                   cluster_SE_by_clustN_fw,
                   cluster_SE_by_clustN_satterthwaite)

      res

    })

  res_wildly_different <-
    data.table::rbindlist(res_wildly_different)

  res_wildly_different[, clusters2:= paste("G =", clusters)]
  res_wildly_different[, clusters2 := factor(clusters2, levels = c("G = 50", "G = 100"))]

  #res_wildly_different <- readRDS("res_wildly_different.rds")
  saveRDS(res_wildly_different, "res_wildly_different_10000.rds")

  #res_wildly_different <- readRDS("res_wildly_different.rds")
  res_wildly_different[type == "fast_n_wild"]$type <- "Wild Cluster Bootstrap"
  res_wildly_different[type == "clustered"]$type <- "CRVE"

  plot <-
  ggplot(data = res_wildly_different[type != "regular"], aes(x = rho, y = coverage * 100, group = type, color = type, fill = type)) +
    facet_wrap(~ clusters2) +
    geom_hline(yintercept=95,color="red",size=1,linetype="dashed")+
    geom_line(size = 1) +
    geom_point() +
    xlab("intra-cluster correlation rho (10.000 Sims Each)") +
    ylab("Avg 95% CI Coverage of True Beta") +
    theme_classic() +
    theme(
      text = element_text(size=18,family = "Economica"),
      axis.title = element_text(size = 16, margin = margin(
        t = 10,
        b = 0,
        r = 20,
        l = 0)
      ),
      axis.text = element_text(size = 14),
      panel.spacing.x = unit(1, "lines")
    ) +
    #ylim(-5,97) +
    scale_x_continuous(breaks= c(0, 0.9))+
    scale_y_continuous(breaks=seq(90, 100 ,by= 1))+
    scale_alpha_continuous(guide=FALSE)



  ggsave("wildly_different.png", width = 9, height = 4)

  plot


}
