library(doParallel)

Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1")

run_window_lmer <- function(q, epoch, varname, use_next_trial = TRUE, time_window = c(-0.5, 0.5)) {
  library(lme4)
  library(tidyverse)
  library(car)
  library(emmeans)
  source("utils.R")

  pc <- 7
  df <- read_csv(paste0("output/PC_", pc, "_", epoch, "_trial_burst_counts.csv")) %>%
    filter(block > 0 & block < 7) %>%
    filter(quartile == q)

  print(paste0('PC',pc,' Q',q,': ',epoch,': ',varname, ' next trial=', use_next_trial))
  print(time_window)

  window_cols <- grep("^time_", names(df), value = TRUE)
  time_vals <- as.numeric(gsub("time_", "", window_cols))
  time_mask <- time_vals >= time_window[1] & time_vals <= time_window[2]
  cols_to_avg <- window_cols[time_mask]

  if (length(cols_to_avg) == 0) stop("No time points in the specified range.")

  if (use_next_trial) {
    df <- df %>% add_next_trial_value(varname)
    value_col <- paste0("next_", varname)
    df <- df %>% filter(!is.na(.data[[value_col]]) & !is.na(.data[[varname]]))

    if (grepl("err", varname)) {
      df <- df %>% filter(.data[[value_col]] < 60) %>% filter(.data[[varname]] < 60)
    }
  } else {
    df <- df %>% add_last_trial_value(varname)
    value_col <- varname
    df <- df %>% filter(!is.na(.data[[value_col]]) & !is.na(.data[[paste0("last_", varname)]]))

    if (grepl("err", varname)) {
      df <- df %>% filter(.data[[value_col]] < 60) %>% filter(.data[[paste0("last_", varname)]] < 60)
    }
  }

  df$mean_burst_count <- rowMeans(df[cols_to_avg], na.rm = TRUE)

  formula <- if (use_next_trial) {
    as.formula(paste0(value_col, " ~ group * mean_burst_count + trial + block + ", varname, " + (1 + group*mean_burst_count + ", varname, " | subject)"))
  } else {
    as.formula(paste0(value_col, " ~ group * mean_burst_count + trial + block + ", paste0('last_',varname), " + (1 + group*mean_burst_count + ", paste0('last_',varname), " | subject)"))
  }

  model <- lmer(formula, data = df)
  print(summary(model))
  res <- Anova(model,type=3)
  print(res)
  pw <- emtrends(model, pairwise~group, var = "mean_burst_count", infer = TRUE)
  print(pw)

  out <- list(Anova = res, emtrends = pw)  # <? safe
  return(out)
}


run_trial_permutation <- function(varname, use_next_trial = TRUE, min_cluster_size = 15) {
  library(permutes)
  library(tidyverse)

  source("utils.R")
  options(contrasts = c("contr.sum", "contr.poly"))

  for (epoch in c("vis", "mot")) {
    all_perm_results <- list()
    all_win_results <- data.frame()

    for (pc in 7:7){#10) {
      df <- read_csv(paste0("output/PC_", pc, "_", epoch, "_trial_burst_counts.csv")) %>%
        filter(block > 0 & block < 7)

      window_cols <- grep("^time_", names(df), value = TRUE)
      window_cols <- window_cols[3:(length(window_cols) - 2)]
      times <- as.numeric(sub("time_", "", window_cols))


      for (q in 1:4) {
        df_q <- df %>% filter(quartile == q)

        if (use_next_trial) {
          df_q <- df_q %>% add_next_trial_value(varname)
          value_col <- paste0("next_", varname)

          # Drop trials with NA in next_trial or current_trial variable
          df_q <- df_q %>% filter(!is.na(.data[[value_col]]) & !is.na(.data[[varname]]))

          # Drop trials with extreme next-trial error or current values (if applicable)
          if (grepl("err", varname)) {
            df_q <- df_q %>% filter(.data[[value_col]] < 60) %>% filter(.data[[varname]] < 60)
          }
        } else {
          df_q <- df_q %>% add_last_trial_value(varname)
          value_col <- varname

          # Drop trials with NA in last_trial or current_trial variable
          df_q <- df_q %>% filter(!is.na(.data[[value_col]]) & !is.na(.data[[paste0("last_", varname)]]))

          if (grepl("err", varname)) {
            df_q <- df_q %>% filter(.data[[value_col]] < 60) %>% filter(.data[[paste0("last_", varname)]] < 60)
          }
        }

        select_vars <- c("subject", "group", "block", "trial", value_col, all_of(window_cols))
        if (use_next_trial) {
          select_vars <- c(select_vars, varname)  # include current-trial covariate
        } else {
          select_vars <- c(select_vars, paste0("last_", varname))
        }

        df_long <- df_q %>%
          select(all_of(select_vars)) %>%
          pivot_longer(cols = all_of(window_cols), names_to = "time", values_to = "burst_count") %>%
          mutate(time = as.numeric(gsub("time_", "", time)))

        if (use_next_trial) {
          formula <- as.formula(paste0(value_col, " ~ group * burst_count + trial + block + ", varname, " + (1 + group*burst_count + ", varname, " | subject)"))
        } else {
          formula <- as.formula(paste0(value_col, " ~ group * burst_count + trial + block + ", paste0("last_", varname), " + (1 + group*burst_count + ", paste0("last_", varname), " | subject)"))
        }

        cl <- makeCluster(60, outfile = '')
        registerDoParallel(cl)

        perms <- clusterperm.lmer(
          formula,
          data = df_long,
          series.var = ~time,
          nperm = 10000,
          parallel = TRUE,
          progress = '',
          type = 'anova',
        )

        stopCluster(cl)
        perms$time<-as.numeric(perms$time)

        interaction_df<-perms[perms$Factor=='group:burst_count' & perms$p.cluster_mass<0.05 & !is.na(perms$p.cluster_mass),]
        cluster_ids<-unique(interaction_df$cluster)
        for(cluster_id in cluster_ids) {
          cluster_df<-interaction_df[interaction_df$cluster==cluster_id,]
          if(nrow(cluster_df)>=min_cluster_size) {
            start<-times[min(cluster_df$time)]
            end<-times[max(cluster_df$time)]

            win_results<-run_window_lmer(q, epoch, varname, use_next_trial = use_next_trial, time_window = c(start, end))
            # --- ANOVA (Type III) ---
            aov_tab <- win_results$Anova
            int_row <- grep("^group:.*mean_burst_count$", rownames(aov_tab))

            interaction_chisq <- aov_tab$Chisq[int_row]
            interaction_dof   <- aov_tab$Df[int_row]
            interaction_p     <- aov_tab[["Pr(>Chisq)"]][int_row]

            # --- emtrends (emmGrid -> data.frame) ---
            em_df   <- as.data.frame(win_results$emtrends$emtrends)   # by-group slopes
            ctr_df  <- as.data.frame(win_results$emtrends$contrasts)  # group contrast of slopes

            if (interaction_p < 0.05) {
              explicit_onesamp_z <- em_df$z.ratio[em_df$group == "Explicit"]
              explicit_onesamp_p <- em_df$p.value[em_df$group == "Explicit"]

              implicit_onesamp_z <- em_df$z.ratio[em_df$group == "Implicit"]
              implicit_onesamp_p <- em_df$p.value[em_df$group == "Implicit"]

              # typically first row is Explicit - Implicit
              exp_imp_z <- ctr_df$z.ratio[1]
              exp_imp_p <- ctr_df$p.value[1]

              sig_cluster <- data.frame(
                pc = pc,
                quartile = q,
                factor = "group:mean_burst_count",
                start = start,
                end = end,
                chi_sq = interaction_chisq,
                df = interaction_dof,
                p_val = interaction_p,
                explicit_onesamp_z = explicit_onesamp_z,
                explicit_onesamp_p = explicit_onesamp_p,
                implicit_onesamp_z = implicit_onesamp_z,
                implicit_onesamp_p = implicit_onesamp_p,
                explicit_implicit_z = exp_imp_z,
                explicit_implicit_p = exp_imp_p
              )

              all_win_results <- rbind(all_win_results, sig_cluster)
            }
          }
        }

        main_df<-perms[perms$Factor=='burst_count' & perms$p.cluster_mass<0.05 & !is.na(perms$p.cluster_mass),]
        cluster_ids<-unique(main_df$cluster)
        for(cluster_id in cluster_ids) {
          cluster_df<-main_df[main_df$cluster==cluster_id,]
          if(nrow(cluster_df)>=min_cluster_size) {
            start<-times[min(cluster_df$time)]
            end<-times[max(cluster_df$time)]
            win_results<-run_window_lmer(q, epoch, varname, use_next_trial = use_next_trial, time_window = c(start, end))

            # --- ANOVA (Type III) ---
            aov_tab <- win_results$Anova
            int_row <- grep("^mean_burst_count$", rownames(aov_tab))

            main_chisq <- aov_tab$Chisq[int_row]
            main_dof   <- aov_tab$Df[int_row]
            main_p     <- aov_tab[["Pr(>Chisq)"]][int_row]

            # --- emtrends (emmGrid -> data.frame) ---
            if (main_p < 0.05) {
              sig_cluster <- data.frame(
                pc = pc,
                quartile = q,
                factor = "mean_burst_count",
                start = start,
                end = end,
                chi_sq = main_chisq,
                df = main_dof,
                p_val = main_p,
                explicit_onesamp_z = NaN,
                explicit_onesamp_p = NaN,
                implicit_onesamp_z = NaN,
                implicit_onesamp_p = NaN,
                explicit_implicit_z = NaN,
                explicit_implicit_p = NaN
              )

              all_win_results <- rbind(all_win_results, sig_cluster)
            }
          }
        }

        perms$pc <- pc
        perms$quartile <- q
        all_perm_results[[paste0("PC", pc, "_Q", q)]] <- perms
      }
    }

    filename_suffix <- if (use_next_trial) paste0(varname, "_next_trial") else varname
    results <- bind_rows(all_perm_results)
    write_csv(results, paste0("output/PC_", epoch, "_", filename_suffix, "_perm.csv"))

    write_csv(all_win_results, paste0("output/PC_", epoch, "_", filename_suffix, "_win.csv"))
  }
}


run_burst_rate_trend_analysis <- function(varname, use_next_trial = TRUE) {
  library(lme4)
  library(lmerTest)
  library(tidyverse)
  library(scales)
  library(emmeans)
  source("utils.R")

  for (epoch in c("vis", "mot")) {
    all_trends <- list()

    for (pc in 7:7){#10) {
      df <- read_csv(paste0("output/PC_", pc, "_", epoch, "_trial_burst_counts.csv")) %>%
        filter(block > 0 & block < 7)

      window_cols <- grep("^time_", names(df), value = TRUE)
      window_cols <- window_cols[3:(length(window_cols) - 2)]

      for (q in 1:4) {
        df_q <- df %>% filter(quartile == q)

        if (use_next_trial) {
          df_q <- df_q %>% add_next_trial_value(varname)
          value_col <- paste0("next_", varname)

          # Drop trials with NA in next_trial or current_trial variable
          df_q <- df_q %>% filter(!is.na(.data[[value_col]]) & !is.na(.data[[varname]]))

          # Drop trials with extreme next-trial error or current values (if applicable)
          if (grepl("err", varname)) {
            df_q <- df_q %>% filter(.data[[value_col]] < 60) %>% filter(.data[[varname]] < 60)
          }
        } else {
          df_q <- df_q %>% add_last_trial_value(varname)
          value_col <- varname

          # Drop trials with NA in last_trial or current_trial variable
          df_q <- df_q %>% filter(!is.na(.data[[value_col]]) & !is.na(.data[[paste0("last_", varname)]]))

          if (grepl("err", varname)) {
            df_q <- df_q %>% filter(.data[[value_col]] < 60) %>% filter(.data[[paste0("last_", varname)]] < 60)
          }
        }

        trends <- list()

        for (col in window_cols) {
          df_q$burst_count <- df_q[[col]]
          time_val <- as.numeric(gsub("time_", "", col))

          if (use_next_trial) {
            formula <- as.formula(paste0(value_col, " ~ group * burst_count + trial + block + ", varname, " + (1 + group*burst_count + ", varname, " | subject)"))
          } else {
            formula <- as.formula(paste0(value_col, " ~ group * burst_count + trial + block + ", paste0("last_", varname), " + (1 + group*burst_count + ", paste0("last_", varname), " | subject)"))
          }
          model <- lmer(
            formula = formula,
            data = df_q
          )

          tr_group <- emtrends(model, ~group, var = "burst_count", infer = TRUE)
          tr_df <- summary(tr_group) %>%
            as_tibble() %>%
            mutate(
              time = time_val,
              quartile = q,
              pc = pc
            )

          trends[[col]] <- tr_df
        }

        trends_df <- bind_rows(trends) %>%
          rename(estimate = burst_count.trend)

        all_trends[[paste0("PC", pc, "_Q", q)]] <- trends_df
      }
    }

    results <- bind_rows(all_trends)
    suffix <- if (use_next_trial) paste0("burst_rate_", varname, "_next_trial") else paste0("burst_rate_", varname)
    write_csv(results, paste0("output/PC_", epoch, "_", suffix, ".csv"))
  }
}




vars <- c('aim_vis_abs_err', 'reach_vis_abs_err', 'reach_dur', 'reach_rt')

for (var in vars) {
  run_trial_permutation(var, use_next_trial = TRUE)
  run_burst_rate_trend_analysis(var, use_next_trial = TRUE)
  run_trial_permutation(var, use_next_trial = FALSE)
  run_burst_rate_trend_analysis(var, use_next_trial = FALSE)
}


