library(genlasso)
#library(ggplot2)
library(purrr)
library(stringr)


###################################################################################
## Change this to the path of the directory where the mean log-power time series ##
## table is saved.                                                               ##
###################################################################################
mean_logpower_ts_filename = paste0(
  "/Users",
  "/ikhultman",
  "/power_time_series_analysis",
  "/CGRP_power_time_series_analysis",
  "/fused_lasso_analysis/",
  "/fused_lasso_data",
  "/combined_mean_logpower_ts.csv");

#################################################################################
## Change this to the path of the directory where the results of this analysis ##
## should be saved.                                                            ##
#################################################################################
save_dir = paste0(
  "/Users",
  "/ikhultman",
  "/power_time_series_analysis",
  "/CGRP_power_time_series_analysis",
  "/fused_lasso_analysis",
  "/fused_lasso_results/");

if (!dir.exists(save_dir) ) {
  dir.create(save_dir);

  stopifnot(dir.exists(save_dir) );
}

mean_logpower_ts_data = read.csv(mean_logpower_ts_filename, header=TRUE);

expected_nonpow_colnames = c("mouse_id", "freq_band", "region");
nonpow_col_ixs = which(
  is.element(colnames(mean_logpower_ts_data), expected_nonpow_colnames) );

stopifnot(length(nonpow_col_ixs) == length(expected_nonpow_colnames) );

power_col_ixs = sort(setdiff(1:ncol(mean_logpower_ts_data), nonpow_col_ixs) );
thinning_factor = 8;
power_col_ixs = seq(
  power_col_ixs[1],
  power_col_ixs[length(power_col_ixs)],
  thinning_factor);

n_power_pts = length(power_col_ixs);

# Construct the D penalty matrix for fused lasso regression.
D = getDtf(n_power_pts, 0);

brain_regions = unique(mean_logpower_ts_data$region);
n_regions = length(brain_regions);

for (rx in 1:n_regions) {
  next_region = brain_regions[rx];
  next_region_ixs = mean_logpower_ts_data$region == next_region;
  next_region_data = mean_logpower_ts_data[next_region_ixs,];
  rownames(next_region_data) = NULL;

  freq_bands = unique(next_region_data$freq_band);
  n_freq_bands = length(freq_bands);

  for (fx in 1:n_freq_bands) {
    next_freq_band = freq_bands[fx];
    next_freq_ixs = next_region_data$freq_band == next_freq_band;
    next_reg_freq_data = next_region_data[next_freq_ixs,];
    rownames(next_reg_freq_data) = NULL;

    n_mice = length(unique(next_reg_freq_data$mouse_id) );

    stopifnot(nrow(next_reg_freq_data) == n_mice);

    # Construct the X design matrix for fused lasso regression.
    X = diag(n_power_pts);
    X = t(matrix(rep(X, n_mice), n_power_pts) );

    y_scaled = as.vector(scale(t(as.matrix(next_reg_freq_data[,power_col_ixs]) )));
    #y = as.vector(t(as.matrix(next_reg_freq_data[,power_col_ixs]) ));
  }
}


fused_lasso_loocv = function(Y, X, D, log_lambdas=NA, gamma=0) {
# DESCRIPTION:
#   Leave-one-out cross-validation for fused lasso.
#
# PARAMETERS:
#        Y: Frequency response vector.
#
#        X: Design matrix for fused lasso.
#
#        D: Fused lasso penalty matrix.
#
#   lambda: If lambda(s) is(are) provided, then that will be used for
#           cross-validation tests; otherwise lambdas over the
#           range from zero to the largest produced during training
#           are tested in order to find which lambda value minimizes
#           the average test MSE.
#
#    gamma: Parameter used to penalize beta parameter estimates.
#
# RETURN VALUE:
#   Smallest MSE's and corresponding lambda values for each test.

  p_dim = ncol(X);
  n_dim = nrow(X) / p_dim;

  tst_inds = matrix(1:(p_dim*n_dim), ncol=n_dim);
  trn_inds = apply(
    tst_inds,
    MARGIN=2,
    FUN=function(col) setdiff(1:(p_dim*n_dim), col) );

  cat(sprintf("Fitting %d leave-one-out cross-validation models ...\n", n_dim) );

  mods_train = lapply(
    as.list(1:n_dim),
    FUN=function(gx) {
      X_train = X[trn_inds[,gx],];
      Y_train = Y[trn_inds[,gx]];
      fusedlasso(Y_train, X_train, D, gamma=gamma)
    });

  if (any(is.na(log_lambdas) )) {
    n_lambdas = 100;
    log_lam_min = 1e-6;
    log_lam_max = log(
      ceiling(
        max(
          unlist(
            lapply(
              mods_train,
              FUN=function(mod) max(mod$lambda) )))));

    log_lambdas = seq(log_lam_min, log_lam_max, length.out=n_lambdas);
  }

  n_lambdas = length(log_lambdas);

  cat(
    sprintf(
      "Computing LOOCV test MSE's at %d different lambda values ...\n\n",
      n_lambdas) );

  loocv_results = data.frame(
    log_lambda=log_lambdas,
    mean_mse=rep(NA, n_lambdas),
    sd_mse=rep(NA, n_lambdas) );

  all_mse = matrix(NA, n_lambdas, n_dim);

  for (mx in 1:n_dim) {
    next_mod_trn = mods_train[[mx]];
    betas = coef(next_mod_trn, lambda=exp(loocv_results$lambda) )$beta;
    X_test = X[tst_inds[,mx],];
    Y_test = Y[tst_inds[,mx]];
    Y_hat = X_test %*% betas;
    all_mse[,mx] = colMeans((Y_test - Y_hat)^2);
  }

  loocv_results$mean_mse = rowMeans(all_mse);
  loocv_results$sd_mse = apply(all_mse, 1, sd);

  loocv_results
};


# Fused lasso parameters
gamma = 0.5;
lambdas_1se = rep(NA, n_reg_pairs);
lambdas_min = rep(NA, n_reg_pairs);

# Location where the plots constructed by this fused lasso analysis will be saved.
figures_dir = paste0(save_dir, "figures/");

if (!dir.exists(figures_dir) ) {
  dir.create(figures_dir);

  stopifnot(dir.exists(figures_dir) );
}

# Location where the plots showing how mean test MSE changes with increasing
# lambda values will be saved. These plots also show a one standard error region
# around the smallest mean test MSE which is used to find the largest lambda
# value that provides decent prediction results.
figures_lam_choice_dir = paste0(figures_dir, "lambda_choice/");

if (!dir.exists(figures_lam_choice_dir) ) {
  dir.create(figures_lam_choice_dir);

  stopifnot(dir.exists(figures_lam_choice_dir) );
}

bg_color = rgb(245, 245, 245, maxColorValue=255);

# Using the leave-one-out cross validation function defined above, determine the
# lambda parameter value whose corresponding mean test MSE is the minimum observed
# and the largest lambda parameter whose corresponding mean test MSE falls within
# one standard error of the smallest mean test MSE for each pair of brain regions.
for (yx in 1:n_reg_pairs) {
  cat(
    paste0(
      "\nPerforming cross-validation for frequency band selection between ",
      sprintf("brain regions %s and %s ...\n", region_pairs[yx,1], region_pairs[yx,2]) ));

  results_loocv = fused_lasso_loocv(Y[,yx], X, D, gamma=gamma);
  min_tst_mse = min(results_loocv$mean_mse);
  min_mse_ix = max(which(results_loocv$mean_mse == min_tst_mse) );
  tst_mse_1se = min_tst_mse + ((results_loocv$sd_mse[min_mse_ix] / sqrt(n_mice) ) * c(-1, 1) );
  lam_1se_ix = max(which(results_loocv$mean_mse <= tst_mse_1se[2]) );
  lambdas_1se[yx] = results_loocv$lambda[lam_1se_ix];
  lambdas_min[yx] = results_loocv$lambda[min_mse_ix];

  plot_df1 = data.frame(
    log_lambda=log(results_loocv$lambda),
    mean_mse=results_loocv$mean_mse,
    min_tst_mse=rep(min_tst_mse, nrow(results_loocv) ),
    mse_1se_lower=rep(tst_mse_1se[1], nrow(results_loocv) ),
    mse_1se_upper=rep(tst_mse_1se[2], nrow(results_loocv) ));

  rm_ixs = is.infinite(plot_df1$log_lambda) | is.na(plot_df1$log_lambda);
  plot_df1 = plot_df1[!rm_ixs,];

  plot_df2 = data.frame(
    log_lambda_min=log(results_loocv$lambda[min_mse_ix]),
    log_lambda_1se=log(results_loocv$lambda[lam_1se_ix]) );

  next_mse_results_plot = (
    ggplot() +
    theme(
      panel.background=element_rect(fill=bg_color) ) +
    geom_line(
      data=plot_df1,
      aes(
        x=log_lambda,
        y=mean_mse,
        colour="mean test MSE"),
      linewidth=1) +
    geom_vline(
      data=plot_df2,
      aes(
        xintercept=log_lambda_min,
        colour="log(lambda_min)"),
      linetype="dotted",
      linewidth=0.7,
      show.legend=FALSE) +
    geom_vline(
      data=plot_df2,
      aes(
        xintercept=log_lambda_1se,
        colour="log(lambda_1se)"),
      linetype="dashed",
      linewidth=0.7,
      show.legend=FALSE) +
    geom_ribbon(
      data=plot_df1,
      aes(
        x=log_lambda,
        y=min_tst_mse,
        ymin=mse_1se_lower,
        ymax=mse_1se_upper,
        colour="1 std. error"),
      alpha=0.2,
      show.legend=FALSE) +
    scale_colour_manual(
      values=c("gray", "red", "blue", "black"),
      guide=guide_legend(
        title="Legend",
        override.aes=list(
          fill=c("gray40", bg_color, bg_color, bg_color),
          linetype=c("solid", "dashed", "dotted", "solid"),
          linewidth=c(6, 0.7, 0.7, 1) ))) +
    labs(
      x="log(lambda)",
      y="mean test MSE") );

  save_filename = paste0(
    figures_lam_choice_dir,
    sprintf(
      "lambda_choice_%s_%s.png",
      region_pairs[yx,1],
      region_pairs[yx,2]) );

  ggsave(
    save_filename,
    next_mse_results_plot,
    units="in",
    height=7,
    width=7,
    device="png");
}

# Location where the plots comparing the observed data to the fused lasso
# fitted values for each brain region will be saved.
figures_fit_v_obs_dir = paste0(figures_dir, "fitted_vs_observed/");

if (!dir.exists(figures_fit_v_obs_dir) ) {
  dir.create(figures_fit_v_obs_dir);

  stopifnot(dir.exists(figures_fit_v_obs_dir) );
}

dir.create(paste0(figures_fit_v_obs_dir, "lambda_1se") );
stopifnot(dir.exists(paste0(figures_fit_v_obs_dir, "lambda_1se") ));

dir.create(paste0(figures_fit_v_obs_dir, "lambda_min") );
stopifnot(dir.exists(paste0(figures_fit_v_obs_dir, "lambda_min") ));

# The fitted values at each frequency for each pair of brain regions will
# be stored in the following table and saved to a CSV file.
results_freq_bands_1se = as.data.frame(matrix(NA, n_freq, n_reg_pairs) );
results_freq_bands_min = as.data.frame(matrix(NA, n_freq, n_reg_pairs) );

for (yx in 1:n_reg_pairs) {
  cat(
    sprintf(
      "Brain regions: %s vs. %s\n",
      region_pairs[yx,1],
      region_pairs[yx,2]) );

  cat(paste0("Gamma: ", gamma, "\n") );
  cat(paste0("Lambda 1 std. err.: ", lambdas_1se[yx], "\n") );
  cat(paste0("Lambda min.: ", lambdas_min[yx], "\n") );

  # Fit fused lasso model
  mod = fusedlasso(Y[,yx], X, D, gamma=gamma);

  results_freq_bands_1se[,yx] = as.numeric(coef(mod, lambda=lambdas_1se[yx])$beta);
  colnames(results_freq_bands_1se)[yx] = paste0(region_pairs[yx,], collapse="-vs-");

  results_freq_bands_min[,yx] = as.numeric(coef(mod, lambda=lambdas_min[yx])$beta);
  colnames(results_freq_bands_min)[yx] = paste0(region_pairs[yx,], collapse="-vs-");

  save_filename = sprintf(
    "fused_lasso_lambda_1se_%s_%s.png",
    region_pairs[yx,1],
    region_pairs[yx,2]);

  save_filename = paste0(
    figures_fit_v_obs_dir,
    "lambda_1se/",
    save_filename);

  cat(
    sprintf(
      paste0(
        "Saving plot comparing fitted and observed values for ",
        "lambda w/ 1 std. err. test MSE to:\n%s\n\n"),
      save_filename) );

  png(filename=save_filename);

  matplot(
    1:50,
    matrix(Y[,yx], ncol=n_mice),
    lty=2,
    col="black",
    main=paste(
      region_pairs[yx,1],
      "->",
      region_pairs[yx,2],
      sep=" "),
    ylab="Scaled Phase Offset",
    xlab="Frequency (Hz)",
    type="l",
    ylim=c(-3, 3) );

  lines(results_freq_bands_1se[[yx]], col="steelblue", lwd=3);
  legend(
    34, 3,
    legend=c("Observed Data", "Fitted Values"),
    col=c("black", "steelblue"),
    lty=c(2, 1),
    lwd=c(1, 3) );

  dev.off();

  save_filename = sprintf(
    "fused_lasso_lambda_min_%s_%s.png",
    region_pairs[yx,1],
    region_pairs[yx,2]);

  save_filename = paste0(
    figures_fit_v_obs_dir,
    "lambda_min/",
    save_filename);

  cat(
    sprintf(
      paste0(
        "Saving plot comparing fitted and observed values for ",
        "lambda w/ min test MSE to:\n%s\n\n"),
      save_filename) );

  png(filename=save_filename);

  matplot(
    1:50,
    matrix(Y[,yx], ncol=n_mice),
    lty=2,
    col="black",
    main=paste(
      region_pairs[yx,1],
      "->",
      region_pairs[yx,2],
      sep=" "),
    ylab="Scaled Phase Offset",
    xlab="Frequency (Hz)",
    type="l",
    ylim=c(-3, 3) );

  lines(results_freq_bands_min[[yx]], col="steelblue", lwd=3);
  legend(
    34, 3,
    legend=c("Observed Data", "Fitted Values"),
    col=c("black", "steelblue"),
    lty=c(2, 1),
    lwd=c(1, 3) );

  dev.off();
}

results_freq_bands_1se = cbind(frequency=1:n_freq, results_freq_bands_1se);
results_freq_bands_min = cbind(frequency=1:n_freq, results_freq_bands_min);

# File where the fused lasso fitted values at each frequency for each pair of
# brain regions corresponding to the largest lambda parameter value whose
# mean test MSE falls within one standard error of the minimum will be saved.
save_filename = paste0(
  save_dir,
  "selected_frequency_bands_lambda_1se.csv");

write.table(
  results_freq_bands_1se,
  save_filename,
  sep=", ",
  row.names=FALSE);

# File where the fused lasso fitted values at each frequency for each pair of
# brain regions corresponding to the lambda parameter value with the smallest
# mean test MSE will be saved.

save_filename = paste0(
  save_dir,
  "selected_frequency_bands_lambda_min.csv");

write.table(
  results_freq_bands_min,
  save_filename,
  sep=", ",
  row.names=FALSE);

