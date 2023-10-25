library(lme4)

data_dir = paste0(
  "/home",
  "/ikhultman",
  "/Desktop",
  "/power_time_series_analysis",
  "/CGRP_power_time_series_analysis",
  "/time_series_dataframes/");

data_files = paste0(data_dir, list.files(data_dir) );
n_files = length(data_files);

expected_vars = c(
  "y",
  "mouse",
  "freq_band",
  "phase",
  "time",
  "day");

for (fx in 1:n_files) {
  region_data = read.csv(data_files[fx], header=TRUE);

  stopifnot(
    (length(
      setdiff(
        expected_vars,
        names(region_data) )) == 0) &&
    (length(
      setdiff(
        names(region_data),
        expected_vars) ) == 0) );

  freq_band_of_interest = 1;
  region_data = region_data[region_data$freq_band == freq_band_of_interest,];

  na_ixs = is.na(region_data$y);
  region_data = region_data[!na_ixs,];

  region_data$mouse = factor(region_data$mouse);
  region_data$phase = factor(region_data$phase);
  region_data$day = factor(region_data$day);
  region_data$time = scale(region_data$time);

  mdl_lmer = lmer(
    y ~ phase:time + day:time + (1 | mouse) + (0 + time | mouse),
    data=region_data);
}

