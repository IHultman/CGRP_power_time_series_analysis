library(emmeans)
library(lme4)
library(lmerTest)
library(multcomp)
library(nlme)
library(splines)


data_filename = paste0(
  "/home",
  "/ikhultman",
  "/Desktop",
  "/power_time_series_analysis",
  "/CGRP_power_time_series_analysis",
  "/combined_time_series_data",
  "/combined_brain_regions_3_phases_ts_data.csv");

logpow_data = read.csv(data_filename, header=TRUE);

expected_vars = c(
  "y",
  "mouse",
  "freq_band",
  "phase",
  "time",
  "day",
  "cgrp",
  "region");

stopifnot(
  (length(
    setdiff(
      expected_vars,
      names(logpow_data) )) == 0) &&
  (length(
    setdiff(
      names(logpow_data),
      expected_vars) ) == 0) );

logpow_data$mouse = factor(logpow_data$mouse);
logpow_data$freq_band = factor(logpow_data$freq_band);
logpow_data$phase = factor(logpow_data$phase);
logpow_data$day = factor(logpow_data$day);
logpow_data$cgrp = factor(logpow_data$cgrp);
logpow_data$region = factor(logpow_data$region);

time_bspline = data.frame(bs(logpow_data$time) );
colnames(time_bspline) = c("time_bs1", "time_bs2", "time_bs3");
logpow_data = cbind(logpow_data, time_bspline);

#logpow_data$mouse_day = factor(
#  paste(logpow_data$mouse, logpow_data$day, sep="_day") );

mdl_lmer_formula = as.formula(
  paste0(
    "y ~ ",
    "time_bs1 + time_bs2 + time_bs3 +",
    "(freq_band * phase * cgrp) + ",
    "(1 | mouse) + ",
    "(1 | mouse:day)") );

mdl_lmer_formula = as.formula(
  paste0(
    "y ~ ",
    "time_bs1 + time_bs2 + time_bs3 +",
    "(region * freq_band * phase * cgrp) + ",
    "(1 | mouse) + ",
    "(1 | mouse:day) + ",
    "(0 + time | mouse) + ",
    "(0 + time | mouse:day)") );

mdl_lmer = lmer(mdl_lmer_formula, data=logpow_data);

#mdl_lme = lme(
#  fixed=y ~ time + (freq_band * phase * cgrp),
#  random=list(mouse=~ mouse, mouse_day=~ mouse_day),
#  data=logpow_data);

mdl_emmeans = emmeans(mdl_lmer, ~ region:freq_band:phase:cgrp);
mdl_emmeans = emmeans(mdl_lmer, ~ phase:cgrp | freq_band);

library(dplyr)

logpow_means_data = (
  logpow_data %>%
  group_by(mouse, region, freq_band, day, phase) %>%
  summarise(y=mean(y) ) %>%
  as.data.frame() );

logpow_means_data = logpow_means_data[order(logpow_means_data$mouse),];
rownames(logpow_means_data) = NULL;

logpow_means_data$cgrp = 0;
cgrp_ixs = (
  (logpow_means_data$day == 2) &
  ((logpow_means_data$phase == "Post_injection") |
   (logpow_means_data$phase == "PI2") ));

logpow_means_data$cgrp[cgrp_ixs] = 1;

logpow_means_data$mouse = factor(logpow_means_data$mouse);
logpow_means_data$freq_band = factor(logpow_means_data$freq_band);
logpow_means_data$phase = factor(logpow_means_data$phase);
logpow_means_data$day = factor(logpow_means_data$day);
logpow_means_data$cgrp = factor(logpow_means_data$cgrp);
logpow_means_data$region = factor(logpow_means_data$region);

logpow_means_mdl_formula = as.formula(
  paste0(
    "y ~ ",
    "(region * freq_band * phase * cgrp) + ",
    "(1 | mouse) + ",
    "(1 | (mouse:day) )") );

logpow_means_mdl_re_design = lFormula(logpow_means_mdl_formula, logpow_means_data[1:36,]);

logpow_means_mdl_lmer = lmer(logpow_means_mdl_formula, logpow_means_data);

logpow_means_mdl_emmeans = emmeans(
  logpow_means_mdl_lmer,
  ~ phase:cgrp | region:freq_band,
  lmerTest.limit=5000,
  adjust="fdr");



n_region = length(unique(logpow_means_data$region) );
n_freq_band = length(unique(logpow_means_data$freq_band) );
n_days = length(unique(logpow_means_data$day) );
n_phase = length(unique(logpow_means_data$phase) );
n_cgrp = length(unique(logpow_means_data$cgrp) );

n_X_rows = n_region * n_freq_band * n_days * n_phase;
logpow_means_design_mat = Matrix(model.matrix(logpow_means_mdl_lmer)[1:n_X_rows,]);

X_region_mat = kronecker(diag(n_region), matrix(1, n_freq_band * n_days * n_phase, 1) );
X_freq_mat = kronecker(
  matrix(1, n_region, 1),
  kronecker(diag(n_freq_band), matrix(1, n_days * n_phase, 1) ));

X_day_mat = kronecker(
  matrix(1, n_region * n_freq_band),
  kronecker(diag(n_days), matrix(1, n_phase, 1) ));

X_phase_mat = kronecker(matrix(1, n_region * n_freq_band * n_days, 1), diag(n_phase) );
X_cgrp_mat = X_day_mat[,2] * X_phase_mat[,2];
X_cgrp_mat = cbind(X_cgrp_mat, -(X_cgrp_mat - 1) );

X_mat = cbind(1, X_region_mat, X_freq_mat, X_phase_mat, X_cgrp_mat);
colnames(X_mat) = NULL;

