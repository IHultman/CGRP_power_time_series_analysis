library(ggplot2)

fused_lasso_results_dir = paste0(
  "/home",
  "/ikhultman",
  "/Desktop",
  "/power_time_series_analysis",
  "/CGRP_power_time_series_analysis",
  "/fused_lasso_analysis",
  "/fused_lasso_results240404");

day_dirs = c("/Day1_results", "/Day2_results", "/Day3_results");
n_days = length(day_dirs);

exp_phase_dirs = c("/baseline", "/post_injection");
n_exp_phases = length(exp_phase_dirs);

fitted_lambda_filenames = data.frame(
  mean=c("/mean_fitted_values_lambda_1se.csv", "/mean_fitted_values_lambda_min.csv"),
  sd=c("/sd_fitted_values_lambda_1se.csv", "/sd_fitted_values_lambda_min.csv") );

n_lambdas = nrow(fitted_lambda_filenames);

for (dx in 1:n_days) {
  next_day_dir = paste0(fused_lasso_results_dir, day_dirs[dx]);

  for (px in 1:n_exp_phases) {
    next_day_phase_dir = paste0(next_day_dir, exp_phase_dirs[px]);

    for (lx in 1:n_lambdas) {
      next_mean_filename = paste0(
        next_day_phase_dir,
        fitted_lambda_filenames$mean[lx]);

      next_sd_filename = paste0(
        next_day_phase_dir,
        fitted_lambda_filenames$sd[lx]);

      next_mean_fitted = read.csv(next_mean_filename, header=TRUE);
      next_sd_fitted = read.csv(next_sd_filename, header=TRUE);
    }
  }
}
