add_next_trial_value <- function(df, colname) {
  df <- df %>%
    arrange(subject, block, trial) %>%
    mutate(next_trial = trial + 1)
  
  next_df <- df %>%
    select(subject, block, trial, !!sym(colname)) %>%
    rename(
      next_trial = trial,
      next_value = !!sym(colname)
    )
  
  df_joined <- left_join(df, next_df, by = c("subject", "block", "next_trial"))
  
  # Rename the joined column to reflect the input
  new_colname <- paste0("next_", colname)
  df_joined <- df_joined %>%
    rename(!!new_colname := next_value)
  
  return(df_joined)
}

add_last_trial_value <- function(df, colname) {
  df <- df %>%
    arrange(subject, block, trial) %>%
    mutate(last_trial = trial - 1)

  last_df <- df %>%
    select(subject, block, trial, !!sym(colname)) %>%
    rename(
      last_trial = trial,
      last_value = !!sym(colname)
    )

  df_joined <- left_join(df, last_df, by = c("subject", "block", "last_trial"))

  # Rename the joined column to reflect the input
  new_colname <- paste0("last_", colname)
  df_joined <- df_joined %>%
    rename(!!new_colname := last_value)

  return(df_joined)
}

add_lagged_trial_value <- function(df, colname, lag = 1) {
  df <- df %>%
    arrange(subject, block, trial) %>%
    mutate(lagged_trial = trial + lag)
  
  lagged_df <- df %>%
    select(subject, block, trial, !!sym(colname)) %>%
    rename(
      lagged_trial = trial,
      lagged_value = !!sym(colname)
    )
  
  df_joined <- left_join(df, lagged_df, by = c("subject", "block", "lagged_trial"))
  
  new_colname <- paste0("lag", lag, "_", colname)
  df_joined <- df_joined %>%
    rename(!!new_colname := lagged_value) %>%
    select(-lagged_trial)
  
  return(df_joined)
}
