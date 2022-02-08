# res <- treatment_effect_sim(n = 2000, n_sims = 50000, workers = 8)
#
# set.seed(1234)

treatment_effect_sim <- function(n,
                                 n_sims,
                                 param = c(.1, .5),
                                 n_cluster = 50,
                                 rho = .5,
                                 #inference = "regular",
                                 #balanced_cluster = TRUE,
                                 workers = 8,
                                 did = TRUE,
                                 #treatment_states = NULL,
                                 boot = "wildboottestjlr"
                                 ){

  # n_sims = 10000
  # param = c(.1, .5)
  # n = 2000
  # n_cluster = 50
  # rho = .5
  # workers = 8
  # did = TRUE
  # boot = "wildboottestjlr"
  # x = "decreasing"

  #'@export

  # devtools::load_all()

  res <-
  lapply(list("increasing", "decreasing", "balanced"), function(x){


    if(x == "increasing"){#
      ordered_states <- c(2:5, 7, 10, 15, 20, 25, 30, 35, 40, 43, 46:50)
      len_points <- length(ordered_states)
    } else {
      ordered_states <- c(1:5, 7, 10, 15, 20, 25, 30, 35, 40, 43, 46:49)
      len_points <- length(ordered_states)
    }

      cat(x, "\n")
      mean_CI_coverage_reg <- matrix(NA, len_points, 3)
      sd_CI_coverage_reg <- matrix(NA, len_points, 3)

      for (j in 16:18) {

        cat(paste0(j, "treatment states:"), "\n")

        if(x == "increasing"){
          treatment_states <- max(ordered_states):ordered_states[length(ordered_states) - j + 1]
        } else {
          treatment_states <- min(ordered_states):ordered_states[j]
        }

        df_res <-
          run_cluster_sim_wildly(n_sims = n_sims,
                          param = param,
                          n = n,
                          n_cluster = n_cluster,
                          rho = rho,
                          balanced_cluster = ifelse(x == "balanced", TRUE, FALSE),
                          workers = workers,
                          did = did,
                          treatment_states = treatment_states,
                          boot = "wildboottestjlr")

        mean_CI_coverage_reg[j,] <- df_res[, mean(param_caught), by = "type"][, V1]
        sd_CI_coverage_reg[j, ] <- df_res[, sd(param_caught), by = "type"][, V1]
      }

      #saveRDS(mean_CI_coverage_reg, paste0("mean_treatment_", n_sims, "_", x, ".rds"))
      #saveRDS(sd_CI_coverage_reg, paste0("sd_treatment_", n_sims, "_", x, ".rds"))

  })

  res <- list(readRDS("mean_treatment_10000_increasing.rds"),
              readRDS("mean_treatment_10000_decreasing.rds"),
              readRDS("mean_treatment_10000_balanced.rds"))

  res2 <- Reduce(rbind, res)
  res2 <- data.table::as.data.table(res2)
  res2$sort <- c(rep("imbalanced clusters, increasing", 18), rep("imbalanced clusters, decreasing", 18), rep("balanced clusters", 18))
  res2$order <- rep(c(1:5, 7, 10, 15, 20, 25, 30, 35, 40, 43, 46:49) / 50, 3)
  setnames(res2, c("V1", "V2", "V3"), c("CRVE", "fast_n_wild", "Satterthwaite"))

  saveRDS(res2, "treatment_simulation")
  cat("Save succesful!", "\n")

  crve <- res2[, list(CRVE, sort, order)][, type := "CRVE"]
  setnames(crve, "CRVE", "coverage")
  fast_n_wild <- res2[, list(fast_n_wild, sort, order)][, type := "Wild Cluster Bootstrap"]
  setnames(fast_n_wild, "fast_n_wild", "coverage")
  satterthwaite <- res2[, list(Satterthwaite, sort, order)][, type := "Satterthwaite"]
  setnames(satterthwaite, "Satterthwaite", "coverage")

  all <- rbind(crve, fast_n_wild, satterthwaite)

  ggplot(data = all, aes(x = order, y = coverage * 100, group = sort, color = sort)) +
    geom_line(alpha = 5) +
    geom_point() +
    facet_wrap(~ type, nrow = 3) +
    geom_hline(yintercept=95,color="red",size=1,linetype="dashed")+
    geom_line(size = 1) +
    xlab("Share of Treated Clusters (10000 Sims Each)") +
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
    scale_x_continuous(breaks= round(c(1:5, 7, 10, 15, 20, 25, 30, 35, 40, 43, 46:49) / 50, 2), guide = guide_axis(angle = 90))+
    scale_y_continuous(breaks=seq(15, 95 ,by= 40))+
    scale_alpha_continuous(guide=FALSE)

    ggsave("treatment_simulations.png", width = 12, height = 9)

}

