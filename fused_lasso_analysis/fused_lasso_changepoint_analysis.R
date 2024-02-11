library(genlasso)
library(ggplot2)
library(parallel)
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

# Location where the plots comparing the observed data to the fused lasso
# fitted values for each brain region will be saved.
figures_fit_v_obs_dir = paste0(figures_dir, "fitted_vs_observed/");

if (!dir.exists(figures_fit_v_obs_dir) ) {
  dir.create(figures_fit_v_obs_dir);

  stopifnot(dir.exists(figures_fit_v_obs_dir) );
}

lambda_1se_dir = paste0(figures_fit_v_obs_dir, "lambda_1se");

if (!dir.exists(lambda_1se_dir) ) {
  dir.create(lambda_1se_dir);
  stopifnot(dir.exists(lambda_1se_dir) );
}

lambda_min_dir = paste0(figures_fit_v_obs_dir, "lambda_min");

if (!dir.exists(lambda_min_dir) ) {
  dir.create(lambda_min_dir);
  stopifnot(dir.exists(lambda_min_dir) );
}

mean_logpower_ts_data = read.csv(mean_logpower_ts_filename, header=TRUE);

expected_nonpow_colnames = c("mouse_id", "freq_band", "region");
nonpow_col_ixs = which(
  is.element(colnames(mean_logpower_ts_data), expected_nonpow_colnames) );

stopifnot(length(nonpow_col_ixs) == length(expected_nonpow_colnames) );

nonpower_info = mean_logpower_ts_data[,nonpow_col_ixs];
n_obs = nrow(nonpower_info);

power_col_ixs = sort(setdiff(1:ncol(mean_logpower_ts_data), nonpow_col_ixs) );
logpower_mat_init = as.matrix(mean_logpower_ts_data[,power_col_ixs]);
thinning_factor = 16;

n_blocks = ceiling(length(power_col_ixs) / thinning_factor);
block_sizes = rep(thinning_factor, n_blocks)
n_remaining = length(power_col_ixs) %% thinning_factor;

if (n_remaining > 0) {
  block_sizes[n_blocks] = n_remaining;
}

stopifnot(sum(block_sizes) == length(power_col_ixs) );

block_ixs_table = rbind(
  cumsum(c(1, block_sizes[1:(n_blocks-1)]) ),
  cumsum(block_sizes) );

logpower_mat = matrix(NA, n_obs, n_blocks);

for (bx in 1:n_blocks) {
  col_ax = block_ixs_table[1,bx];
  col_bx = block_ixs_table[2,bx];
  logpower_mat[,bx] = rowMeans(logpower_mat_init[,col_ax:col_bx], na.rm=TRUE);
}

mean_logpower_ts_data = cbind(nonpower_info, logpower_mat);
power_col_ixs = (1:n_blocks) + ncol(nonpower_info);

#logpower_mat_test = logpower_mat_init;

#if (n_remaining > 0) {
#  logpower_mat_test = cbind(
#    logpower_mat_test,
#    matrix(NA, n_obs, thinning_factor - n_remaining) );
#}

#n_power_pts = ncol(logpower_mat_test) / thinning_factor;
#logpower_arr = array(t(logpower_mat_test), dim=c(thinning_factor, n_power_pts, n_obs) );
#logpower_mat_test = apply(
#  logpower_arr,
#  MARGIN=c(2, 3),
#  FUN=function(nxt_block) mean(nxt_block, na.rm=TRUE) );

#stopifnot(mean((t(logpower_mat_test) - logpower_mat)^2) == 0);

#power_col_ixs = seq(
#  power_col_ixs[1],
#  power_col_ixs[length(power_col_ixs)],
#  thinning_factor);

# Construct the D penalty matrix for fused lasso regression.
D = getDtf(n_blocks, 0);

# Plots background color.
bg_color = rgb(245, 245, 245, maxColorValue=255);

n_reg_freq_combos = nrow(unique(mean_logpower_ts_data[,c("region", "freq_band")]) );

# The fitted values at each frequency for each pair of brain regions will
# be stored in the following table and saved to a CSV file.
results_fitted_1se = as.data.frame(matrix(NA, n_blocks, n_reg_freq_combos) );
results_fitted_min = as.data.frame(matrix(NA, n_blocks, n_reg_freq_combos) );

brain_regions = unique(mean_logpower_ts_data$region);
n_regions = length(brain_regions);


fused_lasso_loocv = function(Y, X, D, log_lambdas, gamma) {
# DESCRIPTION:
#   Leave-one-out cross-validation for fused lasso.
#
# PARAMETERS:
#             Y: Response vector.
#
#             X: Design matrix for fused lasso.
#
#             D: Fused lasso penalty matrix.
#
#   log_lambdas: If log-lambda(s) is(are) provided, then that will be used for
#                cross-validation tests; otherwise log-lambdas over the range from
#                log(1e-6) to the largest produced during training are used.
#
#         gamma: Parameter used to penalize beta parameter estimates such that
#                larger gamma values favor sparser models.
#
# RETURN VALUE:
#   LOOCV test MSE means, standard deviations and corresponding log-lambda values.

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
    log_lam_min = log(1e-6);
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

  for (gx in 1:n_dim) {
    next_mod_trn = mods_train[[gx]];
    betas = coef(next_mod_trn, lambda=exp(loocv_results$log_lambda) )$beta;
    X_test = X[tst_inds[,gx],];
    Y_test = Y[tst_inds[,gx]];
    Y_hat = X_test %*% betas;
    all_mse[,gx] = colMeans((Y_test - Y_hat)^2);
  }

  loocv_results$mean_mse = rowMeans(all_mse);
  loocv_results$sd_mse = apply(all_mse, 1, sd);

  loocv_results
};


fused_lasso_loocv_par = function(Y, X, D, log_lambdas, gamma, cl) {
# DESCRIPTION:
#   Leave-one-out cross-validation for fused lasso.
#
# PARAMETERS:
#             Y: Response vector.
#
#             X: Design matrix for fused lasso.
#
#             D: Fused lasso penalty matrix.
#
#   log_lambdas: If log-lambda(s) is(are) provided, then that will be used for
#                cross-validation tests; otherwise log-lambdas over the range from
#                log(1e-6) to the largest produced during training are used.
#
#         gamma: Parameter used to penalize beta parameter estimates such that
#                larger gamma values favor sparser models.
#
#            cl: Cluster object for parallel computation.
#
# RETURN VALUE:
#   LOOCV test MSE means, standard deviations and corresponding log-lambda values.

  p_dim = ncol(X);
  n_dim = nrow(X) / p_dim;

  tst_inds = matrix(1:(p_dim*n_dim), ncol=n_dim);
  trn_inds = apply(
    tst_inds,
    MARGIN=2,
    FUN=function(col) setdiff(1:(p_dim*n_dim), col) );

  cat(sprintf("Fitting %d leave-one-out cross-validation models ...\n", n_dim) );

  clusterExport(cl, "X", "Y", "D");

  mods_train = parLapply(
    cl,
    as.list(1:n_dim),
    fun=function(gx) {
      X_train = X[trn_inds[,gx],];
      Y_train = Y[trn_inds[,gx]];
      fusedlasso(Y_train, X_train, D, gamma=gamma)
    });

  if (any(is.na(log_lambdas) )) {
    n_lambdas = 100;
    log_lam_min = log(1e-6);
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

  for (gx in 1:n_dim) {
    next_mod_trn = mods_train[[gx]];
    betas = coef(next_mod_trn, lambda=exp(loocv_results$log_lambda) )$beta;
    X_test = X[tst_inds[,gx],];
    Y_test = Y[tst_inds[,gx]];
    Y_hat = X_test %*% betas;
    all_mse[,gx] = colMeans((Y_test - Y_hat)^2);
  }

  loocv_results$mean_mse = rowMeans(all_mse);
  loocv_results$sd_mse = apply(all_mse, 1, sd);

  loocv_results
};


col_ix = 0;

for (rx in 1:n_regions) {
  next_region = brain_regions[rx];
  next_region_ixs = mean_logpower_ts_data$region == next_region;
  next_region_data = mean_logpower_ts_data[next_region_ixs,];
  rownames(next_region_data) = NULL;

  freq_bands = unique(next_region_data$freq_band);
  n_freq_bands = length(freq_bands);

  for (fx in 1:n_freq_bands) {
    col_ix = col_ix + 1;

    next_freq_band = freq_bands[fx];
    next_freq_ixs = next_region_data$freq_band == next_freq_band;
    next_reg_freq_data = next_region_data[next_freq_ixs,];
    rownames(next_reg_freq_data) = NULL;

    n_mice = length(unique(next_reg_freq_data$mouse_id) );

    stopifnot(nrow(next_reg_freq_data) == n_mice);

    # Construct the X design matrix for fused lasso regression.
    X = diag(n_blocks);
    X = t(matrix(rep(X, n_mice), n_blocks) );

    y_scaled_mat = scale(t(as.matrix(next_reg_freq_data[,power_col_ixs]) ));
    Y = as.vector(y_scaled_mat);

    gamma = 0;

    start_time = Sys.time();
    loocv_results = fused_lasso_loocv(Y, X, D, NA, gamma);
    run_time = difftime(Sys.time(), start_time, units="secs")[[1]];

    start_time = Sys.time();
    loocv_par_results = fused_lasso_loocv_par(Y, X, D, NA, gamma, cl);
    run_time = difftime(Sys.time(), start_time, units="secs")[[1]];

    min_tst_mse = min(loocv_results$mean_mse);
    min_mse_ix = max(which(loocv_results$mean_mse == min_tst_mse) );
    tst_mse_1se = min_tst_mse + ((loocv_results$sd_mse[min_mse_ix] / sqrt(n_mice) ) * c(-1, 1) );
    lam_1se_ix = max(which(loocv_results$mean_mse <= tst_mse_1se[2]) );
    log_lam_1se = loocv_results$log_lambda[lam_1se_ix];
    log_lam_min = loocv_results$log_lambda[min_mse_ix];

    plot_df1 = data.frame(
      log_lambda=loocv_results$log_lambda,
      mean_mse=loocv_results$mean_mse,
      min_tst_mse=rep(min_tst_mse, nrow(loocv_results) ),
      mse_1se_lower=rep(tst_mse_1se[1], nrow(loocv_results) ),
      mse_1se_upper=rep(tst_mse_1se[2], nrow(loocv_results) ));

    rm_ixs = is.infinite(plot_df1$log_lambda) | is.na(plot_df1$log_lambda);
    plot_df1 = plot_df1[!rm_ixs,];

    plot_df2 = data.frame(
      log_lambda_min=loocv_results$log_lambda[min_mse_ix],
      log_lambda_1se=loocv_results$log_lambda[lam_1se_ix]);

    next_mse_results_plot = (
      ggplot() +
      theme(
        panel.background=element_rect(fill=bg_color) ) +
      geom_line(
        data=plot_df1,
        aes(
          x=log_lambda,
          y=mean_mse,
          colour="Mean Test MSE"),
        linewidth=1) +
      geom_vline(
        data=plot_df2,
        aes(
          xintercept=log_lambda_min,
          colour="Log(lambda_min)"),
        linetype="dotted",
        linewidth=0.7,
        show.legend=FALSE) +
      geom_vline(
        data=plot_df2,
        aes(
          xintercept=log_lambda_1se,
          colour="Log(lambda_1se)"),
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
          colour="1 Std. Error"),
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
        x="Log(lambda)",
        y="Mean Test MSE") );

    save_filename = paste0(
      figures_lam_choice_dir,
      sprintf(
        "lambda_choice_power_%s_freq_band_%d.png",
        next_region,
        next_freq_band) );

    ggsave(
      save_filename,
      next_mse_results_plot,
      units="in",
      height=7,
      width=7,
      device="png");

    mod_full_data = fusedlasso(Y, X, D, gamma=gamma);

    results_fitted_1se[,col_ix] = as.numeric(coef(mod_full_data, lambda=exp(log_lam_1se) )$beta);
    colnames(results_fitted_1se)[col_ix] = paste0(next_region, "_freq_band_", next_freq_band);

    results_fitted_min[,col_ix] = as.numeric(coef(mod_full_data, lambda=exp(log_lam_min) )$beta);
    colnames(results_fitted_min)[col_ix] = paste0(next_region, "_freq_band_", next_freq_band);

    save_filename = sprintf(
      "fused_lasso_lambda_1se_power_%s_freq_band_%d.png",
      next_region,
      next_freq_band);

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

    plot_xs = ((0:(n_blocks - 1) ) * thinning_factor) + 1;
    mean_scaled_y = as.numeric(rowMeans(y_scaled_mat) );
    y_lim = c(
      min(-1, min(mean_scaled_y) ),
      max(1, max(mean_scaled_y) ));

    plot(
      plot_xs,
      mean_scaled_y,
      lty=2,
      col="black",
      main=paste0(
        next_region,
        " - frequency band ",
        next_freq_band),
      ylab="Scaled Log-power",
      xlab="Time (seconds)",
      type="l",
      ylim=y_lim);

    lines(plot_xs, results_fitted_1se[[col_ix]], col="steelblue", lwd=3);
    abline(v=1801, col="green4", lty=6, lwd=2);
    legend(
      0.65 * n_blocks * thinning_factor, 1,
      legend=c("Mean Observed", "Fitted", "Injection"),
      col=c("black", "steelblue", "green4"),
      lty=c(2, 1, 6),
      lwd=c(1, 3, 2) );

    dev.off();

    save_filename = sprintf(
      "fused_lasso_lambda_min_power_%s_freq_band_%d.png",
      next_region,
      next_freq_band);

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

    plot(
      plot_xs,
      mean_scaled_y,
      lty=2,
      col="black",
      main=paste0(
        next_region,
        " - frequency band ",
        next_freq_band),
      ylab="Scaled Log-power",
      xlab="Time (seconds)",
      type="l",
      ylim=y_lim);

    lines(plot_xs, results_fitted_min[[col_ix]], col="steelblue", lwd=3);
    abline(v=1801, col="green4", lty=6, lwd=2);
    legend(
      0.65 * n_blocks * thinning_factor, 1,
      legend=c("Mean Observed", "Fitted", "Injection"),
      col=c("black", "steelblue", "green4"),
      lty=c(2, 1, 6),
      lwd=c(1, 3, 2) );

    dev.off();
  }
}

save_filename = paste0(
  save_dir,
  "fitted_values_lambda_1se.csv");

write.table(
  results_fitted_1se,
  save_filename,
  sep=", ",
  row.names=FALSE);

save_filename = paste0(
  save_dir,
  "fitted_values_lambda_min.csv");

write.table(
  results_fitted_min,
  save_filename,
  sep=", ",
  row.names=FALSE);

