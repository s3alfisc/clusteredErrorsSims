### LET'S NOW EXPLORE MORE SYSTEMATICALLY HOW THIS VARIES OVER CLUSTER COUNTS.

sim_balanced_clusters <- function(n = 1000, n_sims = 1000, rho = 0.7, workers = 8){

  #'@export


  # Steps:
  # 1) Robust inference via sandwich estimators
  # 2) Cluster robust inference via sandwich estimators
  # 3) Cluster robust inference via the wild cluster bootstrap
  # 4) Cluster robust inference via sandwich estimators with Satterthwaite correction
  # 5) Cluster robust inference via sandwich estimators with saddlepoint correction
  # Finally, Step 6 collects and plots all results

  # Step 1: how does "regular" inference do for 2:50 clusters?
  mean_CI_coverage_reg <- rep(0,49)
  sd_CI_coverage_reg <- rep(0,49)
  for (j in 2:50) {
    sim_clustered_fix_tmp <-
      run_cluster_sim(n = n,
                      n_sims = n_sims,
                      rho = rho,
                      n_cluster = j,
                      inference = "regular",
                      workers = workers)
    mean_CI_coverage_reg[j-1] <- mean(sim_clustered_fix_tmp$param_caught)
    sd_CI_coverage_reg[j-1] <- sd(sim_clustered_fix_tmp$param_caught)
  }

  regular_SE_by_clustN <- data.frame(coverage = mean_CI_coverage_reg,
                                     sd_CI_coverage_reg,
                                     clusters = seq(2:50))

  # Step 2: how does "standard clustered" inference do for 2:50 clusters?

  mean_CI_coverage <- rep(0,49)
  sd_CI_coverage <- rep(0,49)
  for (j in 2:50) {
    sim_clustered_fix_tmp <-
      run_cluster_sim(n = n,
                      n_sims = n_sims,
                      rho = rho,
                      n_cluster = j,
                      inference = "clustered",
                      workers = workers)
    mean_CI_coverage[j-1] <- mean(sim_clustered_fix_tmp$param_caught)
    sd_CI_coverage[j-1] <- sd(sim_clustered_fix_tmp$param_caught)
  }

  cluster_SE_by_clustN <- data.frame(coverage = mean_CI_coverage, sd_CI_coverage, clusters = seq(2:50))

  # Step 3: how does wild cluster bootstrapped inference do for 2:50 clusters?
  # with B = 9999 iterations & webb weights for G <= 12, else rademacher weights

  mean_CI_coverage_fw <- rep(0,49)
  sd_CI_coverage_fw <- rep(0,49)
  for (j in 2:50) {
    tryCatch(
      {
        sim_clustered_fix_tmp <-
          run_cluster_sim(n = n,
                          n_sims = n_sims,
                          rho = rho,
                          n_cluster = j,
                          inference = "fast_n_wild")
        mean_CI_coverage_fw[j-1] <- mean(sim_clustered_fix_tmp$param_caught)
        sd_CI_coverage_fw[j-1] <- sd(sim_clustered_fix_tmp$param_caught)
      },
      error=function(cond) {
        message(paste("\nFast'n'Wild Failed.\n"))
        message("Here's the original error message:")
        message(cond)
        # Choose a return value in case of error
        return(NA)
      })

  }

  cluster_SE_by_clustN_fw <- data.frame(coverage = mean_CI_coverage_fw, sd_CI_coverage_fw, clusters = seq(2:50))


  # Step 4: how do satterthwaite corrections do for 2:50 clusters?

  mean_CI_coverage <- rep(0,49)
  sd_CI_coverage <- rep(0,49)
  for (j in 2:50) {
    sim_clustered_fix_tmp <-
      run_cluster_sim(n = n,
                      n_sims = n_sims,
                      rho = rho,
                      n_cluster = j,
                      inference = "Satterthwaite",
                      workers = workers)
    mean_CI_coverage[j-1] <- mean(sim_clustered_fix_tmp$param_caught)
    sd_CI_coverage[j-1] <- sd(sim_clustered_fix_tmp$param_caught)
  }

  cluster_SE_by_clustN_satterthwaite <- data.frame(coverage = mean_CI_coverage, sd_CI_coverage, clusters = seq(2:50))


  # Step 5: how do saddlepoint corrections do for 2:50 clusters?

  #mean_CI_coverage <- rep(0,49)
  #sd_CI_coverage <- rep(0,49)

  #for (j in 2:50) {
  #  sim_clustered_fix_tmp <-
  #    run_cluster_sim(n = n,
  #                    rho = rho,
  #                    n_cluster = j,
  #                    inference = "saddlepoint",
  #                    workers = workers)
  #  mean_CI_coverage[j-1] <- mean(sim_clustered_fix_tmp$param_caught)
  #  sd_CI_coverage[j-1] <- sd(sim_clustered_fix_tmp$param_caught)
  #}

  #cluster_SE_by_clustN_saddlepoint <- data.frame(coverage = mean_CI_coverage, sd_CI_coverage, clusters = seq(2:50))

  # Step 6: Collect all results

  perfPalette_fw <- colorRampPalette(brewer.pal(9, "Greens"))(length(unique(cluster_SE_by_clustN_fw$coverage)))
  unique_coverages_fw <- cluster_SE_by_clustN_fw %>%
    group_by(coverage) %>%
    summarize(coverage = first(coverage)) %>%
    arrange(desc(coverage))
  unique_coverages_fw$heat <- perfPalette_fw

  cluster_SE_by_clustN_fw2 <- cluster_SE_by_clustN_fw %>%
    merge(unique_coverages_fw, by="coverage")

  # Now we can plot the SE performance side-by-side.
  plot <-
  ggplot(data = cluster_SE_by_clustN_fw2, aes(y=coverage*100,x=(clusters+1),alpha=0.7)) +
    geom_hline(yintercept=95,color="red",size=1,linetype="dashed")+
    #geom_vline(xintercept=30,color="black",size=1)+
    geom_line(aes(color="green"),size=2) +
    geom_line(data=cluster_SE_by_clustN, aes(color="blue"),size=2)+
    geom_line(data=regular_SE_by_clustN, aes(color="purple"),size=2)+
    #geom_line(data = cluster_SE_by_clustN_saddlepoint, aes(color = "violet"), size = 2) +
    geom_line(data = cluster_SE_by_clustN_satterthwaite, aes(color = "orange"), size = 2) +
    xlab("# of Clusters (1000 Sims Each)") +
    ylab("Avg 95% CI Coverage of True Beta") +
    theme_classic() +
    theme(
      text = element_text(size=18,family = "Economica"),
      axis.title = element_text(size = 16, margin = margin(
        t = 10,
        b = 0,
        r = 20,
        l = 0
      )),
      axis.text = element_text(size = 14)
    ) +
    ylim(-5,97) +
    scale_x_continuous(breaks=c(2,10,20,30,40,50))+
    scale_y_continuous(breaks=seq(5,95,by=10))+
    scale_alpha_continuous(guide=FALSE)+
    #ggtitle("") +
    scale_color_manual(name="SE Type",values=c("blue","green", "purple", "orange"),labels=c("Clustered","Wild Clustered", "Satterthwaite","Homoskedastic iid"))
  NULL

  ggsave(filename = "CI Coverage for Several Cluster Robust Inference Methods.png", width = 8, height = 6)

  plot

}

