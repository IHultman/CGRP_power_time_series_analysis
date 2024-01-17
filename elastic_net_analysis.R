library(dplyr)
library(gglasso)

data_filename = paste0(
  "/home",
  "/ikhultman",
  "/Desktop",
  "/power_time_series_analysis",
  "/CGRP_power_time_series_analysis",
  "/elastic_net_data",
  "/mean_log_power_group_lasso.csv");

light_dark_filename = paste0(
  "/home",
  "/ikhultman",
  "/Desktop",
  "/power_time_series_analysis",
  "/CGRP_power_time_series_analysis",
  "/elastic_net_data",
  "/light_dark_phase_table.csv");

logpower_data = read.csv(data_filename, header=TRUE);
light_dark_data = read.csv(light_dark_filename, header=TRUE);
ld_cols_to_keep = c("mouse", "phase", "logodds_time_in_dark");
light_dark_data = light_dark_data[,ld_cols_to_keep];

exp_phases_of_interest = c("Baseline", "PI2");
phase_row_ixs = is.element(logpower_data$phase, exp_phases_of_interest);

freq_bands_of_interest = c(1, 3);
freq_row_ixs = is.element(logpower_data$freq_band, freq_bands_of_interest);

logpower_data = logpower_data[phase_row_ixs & freq_row_ixs,];

rownames(logpower_data) = NULL;

logpower_data$cgrp = 0;
cgrp_exp_phases = c("Post_injection", "PI2");
cgrp_phase_row_ixs = is.element(logpower_data$phase, cgrp_exp_phases);
cgrp_ixs = (logpower_data$day == 2) & cgrp_phase_row_ixs;
logpower_data$cgrp[cgrp_ixs] = 1;

logpower_phase_means_data = (
  logpower_data %>%
  group_by(mouse, freq_band, Hz, region, phase, cgrp) %>%
  summarise(logpower=mean(mean_log_power, na.rm=TRUE), .groups="drop") %>%
  as.data.frame() );

logpower_phase_means_data$phase[
  logpower_phase_means_data$phase == "Baseline"] = "baseline";

logpower_phase_means_data$phase[
  (logpower_phase_means_data$phase == "PI2") &
  (logpower_phase_means_data$cgrp == 0)] = "vehicle";

logpower_phase_means_data$phase[
  (logpower_phase_means_data$phase == "PI2") &
  (logpower_phase_means_data$cgrp == 1)] = "treatment";

order_ixs = order(
  logpower_phase_means_data$mouse,
  logpower_phase_means_data$phase,
  logpower_phase_means_data$region,
  logpower_phase_means_data$freq_band,
  logpower_phase_means_data$Hz);

logpower_phase_means_data = logpower_phase_means_data[order_ixs,];
rownames(logpower_phase_means_data) = NULL;

logpower_phase_means_data$label = paste(
  logpower_phase_means_data$region,
  logpower_phase_means_data$freq_band,
  logpower_phase_means_data$Hz,
  sep="_");

logpower_design_mat = (
  logpower_phase_means_data %>%
  group_by(mouse, phase) %>%
  summarise(x=matrix(logpower, nrow=1), .groups="drop") %>%
  as.data.frame() );

group_count_check = (
  logpower_phase_means_data %>%
  group_by(mouse, phase) %>%
  summarise(n=n(), .groups="drop") %>%
  as.data.frame() );

stopifnot(all(group_count_check$n == group_count_check$n[1]) );

n_data_pts = group_count_check$n[1];

design_mat_check = t(matrix(logpower_phase_means_data$logpower, n_data_pts) );
stopifnot(
  mean((logpower_design_mat$x - design_mat_check)^2, na.rm=TRUE) < 1e-16);

group_id_check = (
  logpower_phase_means_data %>%
  group_by(region, freq_band) %>%
  mutate(gid=cur_group_id() ) %>%
  ungroup() %>%
  select(gid) %>%
  as.matrix() %>%
  matrix(n_data_pts) );

group_ids = group_id_check[,1];

stopifnot(
  all(
    matrix(rep(group_ids, ncol(group_id_check) ), n_data_pts) ==
    group_id_check) );

logpower_colnames = (
  logpower_phase_means_data %>%
  group_by(mouse, phase) %>%
  summarise(x=matrix(label, nrow=1), .groups="drop") %>%
  as.data.frame() );

stopifnot(
  all(
    logpower_colnames$x ==
    matrix(
      rep(logpower_colnames[1,"x"], each=nrow(logpower_colnames) ),
      nrow(logpower_colnames) )));

logpower_colnames = as.vector(logpower_colnames[1,"x"]);

logpower_design_mat = cbind(logpower_design_mat[,1:2], logpower_design_mat$x);
colnames(logpower_design_mat) = c(
  colnames(logpower_design_mat)[1:2],
  logpower_colnames);

logpower_design_mat = merge(logpower_design_mat, light_dark_data);
logpower_design_mat = logpower_design_mat[,-(1:2)];

X_mat = as.matrix(logpower_design_mat[,logpower_colnames]);
y_vec = logpower_design_mat$logodds_time_in_dark;

na_vec_ixs = which(is.na(X_mat) );
na_col_ixs = 1 + floor((na_vec_ixs - 1) / nrow(X_mat) );

col_means = colMeans(X_mat, na.rm=TRUE);
na_mean_replace = col_means[na_col_ixs];
X_mat[na_vec_ixs] = na_mean_replace;

mod_gglasso = cv.gglasso(X_mat, y_vec, group_ids, nfolds=10);



#mouse_ids = unique(logpower_data$mouse);
#n_mouse = length(mouse_ids);
#region_ids = unique(logpower_data$region);
#n_region = length(region_ids);
#n_phases = length(exp_phases_of_interest);
#match_vals = c();

#for (mx in 1:n_mouse) {
#  next_mouse = mouse_ids[mx];
#  freqs_mx = unique(logpower_data$freq_band[logpower_data$mouse == next_mouse]);
#  n_freqs_mx = length(freqs_mx);

#  for (fx in 1:n_freqs_mx) {
#    next_freq_band = freqs_mx[fx];
#    hz_mx_fx = unique(
#      logpower_data$Hz[
#        (logpower_data$mouse == next_mouse) &
#        (logpower_data$freq_band == next_freq_band)]);

#    n_hz_mx_fx = length(hz_mx_fx);

#    for (hx in 1:n_hz_mx_fx) {
#      next_hz = hz_mx_fx[hx];

#      for (rx in 1:n_region) {
#        next_region = region_ids[rx];

#        for (px in 1:n_phases) {
#          next_phase = exp_phases_of_interest[px];

#          for (cx in 0:1) {
#            mean_check1 = mean(
#              logpower_data$mean_log_power[
#                (logpower_data$mouse == next_mouse) &
#                (logpower_data$freq_band == next_freq_band) &
#                (logpower_data$Hz == next_hz) &
#                (logpower_data$region == next_region) &
#                (logpower_data$phase == next_phase) &
#                (logpower_data$cgrp == cx)]);

#            mean_check2 = logpower_means_data$y[
#              (logpower_means_data$mouse == next_mouse) &
#              (logpower_means_data$freq_band == next_freq_band) &
#              (logpower_means_data$Hz == next_hz) &
#              (logpower_means_data$region == next_region) &
#              (logpower_means_data$phase == next_phase) &
#              (logpower_means_data$cgrp == cx)];

#            if (is.na(mean_check1) ) {
#              match_vals = c(match_vals, length(mean_check2) == 0);
#            } else {
#              match_vals = c(match_vals, mean_check1 == mean_check2);
#            }
#          }
#        }
#      }
#    }
#  }
#}


