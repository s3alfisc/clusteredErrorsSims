treatment_effect_sim <- function(n = 2000,
                                 n_sims = 1000,
                                 workers = 8
                                 ){

  #'@export

  # devtools::load_all()

  res <-
  lapply(list("increasing", "decreasing", "balanced"), function(x){

    if(x == "increasing"){
      ordered_states <- 50:1
    } else if(x == "decreasing"){
      ordered_states <- 1:50
    } else if(x == "balanced"){
      ordered_states <- 1:50
    }

    # res <-
    # lapply(0.5, function(x){


      # Step 2: how does "cluster robust/sandwich" inference do for 2:50 with wildly different sizes?

      mean_CI_coverage_reg <- rep(0,49)
      sd_CI_coverage_reg <- rep(0,49)
      #i <- 1
      for (j in 2:50) {
        treatment_states <- ordered_states[1:j]
        sim_clustered_fix_tmp <-
          run_cluster_sim(n_sims = n_sims,
                          n = n,
                          rho = 0.5,
                          n_cluster = 50,
                          inference = "clustered",
                          workers = workers,
                          balanced_cluster = ifelse(x == "balanced", TRUE, FALSE),
                          did = TRUE,
                          treatment_states = treatment_states)

        mean_CI_coverage_reg[j-1] <- mean(sim_clustered_fix_tmp$param_caught)
        sd_CI_coverage_reg[j-1] <- sd(sim_clustered_fix_tmp$param_caught)
        #i <- i + 1
      }

      cluster_SE_by_clustN <- data.frame(coverage = mean_CI_coverage_reg,
                                         #sd_CI_coverage_reg,
                                         clusters = 50,
                                         sort = x,
                                         type = "clustered",
                                         rho = 0.5)



      # Step 3: how does "wild cluster robust" inference do for 2:50 with wildly different sizes?

      mean_CI_coverage_fw <- rep(0,49)
      sd_CI_coverage_fw <- rep(0,49)
      #i <- 1
      for (j in 2:50) {
        tryCatch(
          {
            treatment_states <- ordered_states[1:j]
            sim_clustered_fix_tmp <-
              run_cluster_sim(n_sims = n_sims,
                              n = n,
                              rho = 0.5,
                              n_cluster = 50,
                              inference = "fast_n_wild",
                              workers = workers,
                              balanced_cluster = ifelse(x == "balanced", TRUE, FALSE),
                              did = TRUE,
                              treatment_states = treatment_states)

            mean_CI_coverage_fw[j-1] <- mean(sim_clustered_fix_tmp$param_caught)
            sd_CI_coverage_fw[j-1] <- sd(sim_clustered_fix_tmp$param_caught)
            #i <- i + 1
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
                                            clusters = 50,
                                            sort = x,
                                            type = "fast_n_wild",
                                            rho = 0.5)



      # Step 4: how do satterthwaite corrections do for 2:50 clusters?

      mean_CI_coverage <- rep(0,49)
      sd_CI_coverage <- rep(0,49)
      #i <- 1
      for (j in 2:50) {
        treatment_states <- ordered_states[1:j]
        sim_clustered_fix_tmp <-
          run_cluster_sim(n_sims = n_sims,
                          n = n,
                          rho = 0.5,
                          n_cluster = 50,
                          inference = "Satterthwaite",
                          workers = workers,
                          balanced_cluster = ifelse(x == "balanced", TRUE, FALSE),
                          did = TRUE,
                          treatment_states = treatment_states)

        mean_CI_coverage[j-1] <- mean(sim_clustered_fix_tmp$param_caught)
        sd_CI_coverage[j-1] <- sd(sim_clustered_fix_tmp$param_caught)
        #i <- i + 1
      }

      cluster_SE_by_clustN_satterthwaite <- data.frame(coverage = mean_CI_coverage,
                                                       clusters = 50,
                                                       sort = x,
                                                       type = "Satterthwaite",
                                                       rho = 0.5)


      res <- rbind(cluster_SE_by_clustN,
                   cluster_SE_by_clustN_fw,
                   cluster_SE_by_clustN_satterthwaite)

      res

  })

  res2 <- data.table::rbindlist(res)
  res2$order <- rep(2:50, 9) / 49

  saveRDS(res2, "treatment_simulation")

  ggplot(data = res2, aes(x = order, y = coverage * 100, group = sort, color = sort)) +
    geom_line() +
    facet_wrap(~ type, nrow = 3) +
    geom_hline(yintercept=95,color="red",size=1,linetype="dashed")+
    geom_line(size = 2) +
    xlab("Share of Treated Clusters (1000 Sims Each)") +
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
      panel.spacing.x = unit(2, "lines")
    ) +
    #ylim(-5,97) +
    scale_x_continuous(breaks= seq(0.1, 1, 0.1))+
    scale_y_continuous(breaks=seq(90, 100 ,by= 1))+
    scale_alpha_continuous(guide=FALSE)

    ggsave("treatment_simulations.png", width = 8, height = 8)

}

