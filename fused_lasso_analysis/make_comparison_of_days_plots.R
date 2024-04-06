library(ggplot2)

fused_lasso_results_dir = paste0(
  "/home",
  "/ikhultman",
  "/Desktop",
  "/power_time_series_analysis",
  "/CGRP_power_time_series_analysis",
  "/fused_lasso_analysis",
  "/fused_lasso_results240404");

day_info_df = data.frame(
  label=c("Day_1", "Day_2", "Day_3"),
  dir=c("/Day1_results", "/Day2_results", "/Day3_results") );

n_days = nrow(day_info_df);

exp_phase_info_df = data.frame(
  label=c("baseline", "post_injection"),
  dir=c("/baseline", "/post_injection") );

n_exp_phases = nrow(exp_phase_info_df);

lambda_info_df = data.frame(
  label=c("lam_1se", "lam_min"),
  mean_fname=c("/mean_fitted_values_lambda_1se.csv", "/mean_fitted_values_lambda_min.csv"),
  sd_fname=c("/sd_fitted_values_lambda_1se.csv", "/sd_fitted_values_lambda_min.csv") );

n_lambdas = nrow(lambda_info_df);

thinning_factor = 15;

plot_dfs = NULL;
reg_freq_vars = NULL;
n_rf_vars = NULL;

for (dx in 1:n_days) {
  next_day_dir = paste0(fused_lasso_results_dir, day_info_df$dir[dx]);
  next_day_label = day_info_df$label[dx];

  next_day_plot_df = list();

  for (lx in 1:n_lambdas) {
    next_day_plot_df[[lambda_info_df$label[lx]]] = list();
  }

  for (px in 1:n_exp_phases) {
    next_day_phase_dir = paste0(next_day_dir, exp_phase_info_df$dir[px]);

    for (lx in 1:n_lambdas) {
      next_lam_label = lambda_info_df$label[lx];

      next_mean_filename = paste0(
        next_day_phase_dir,
        lambda_info_df$mean_fname[lx]);

      next_sd_filename = paste0(
        next_day_phase_dir,
        lambda_info_df$sd_fname[lx]);

      next_mean_fitted = read.csv(next_mean_filename, header=TRUE);
      next_sd_fitted = read.csv(next_sd_filename, header=TRUE);

      if (is.null(reg_freq_vars) ) {
        reg_freq_vars = colnames(next_mean_fitted);
        n_rf_vars = length(reg_freq_vars);
      }

      stopifnot(
        (length(setdiff(colnames(next_mean_fitted), reg_freq_vars) ) == 0) &&
        (length(setdiff(reg_freq_vars, colnames(next_mean_fitted) )) == 0) );

      stopifnot(
        (length(setdiff(colnames(next_sd_fitted), reg_freq_vars) ) == 0) &&
        (length(setdiff(reg_freq_vars, colnames(next_sd_fitted) )) == 0) );

      if (px == 1) {
        for (rf_ix in 1:n_rf_vars) {
          next_rf_var = reg_freq_vars[rf_ix];

          next_day_plot_df[[next_lam_label]][[next_rf_var]] = data.frame(
            exp(next_mean_fitted[[next_rf_var]]),
            exp(next_mean_fitted[[next_rf_var]] - (2 * next_sd_fitted[[next_rf_var]]) ),
            exp(next_mean_fitted[[next_rf_var]] + (2 * next_sd_fitted[[next_rf_var]]) ));

          colnames(next_day_plot_df[[next_lam_label]][[next_rf_var]]) = paste0(
            next_day_label,
            c("_mean", "_min", "_max") );
        }
      } else {
        for (rf_ix in 1:n_rf_vars) {
          next_rf_var = reg_freq_vars[rf_ix];

          next_day_phase_plot_df = data.frame(
            exp(next_mean_fitted[[next_rf_var]]),
            exp(next_mean_fitted[[next_rf_var]] - (2 * next_sd_fitted[[next_rf_var]]) ),
            exp(next_mean_fitted[[next_rf_var]] + (2 * next_sd_fitted[[next_rf_var]]) ));

          colnames(next_day_phase_plot_df) = paste0(next_day_label, c("_mean", "_min", "_max") );

          next_day_plot_df[[next_lam_label]][[next_rf_var]] = rbind(
            next_day_plot_df[[next_lam_label]][[next_rf_var]],
            next_day_phase_plot_df);
        }
      }
    }
  }

  if (dx == 1) {
    plot_dfs = next_day_plot_df;
  } else {
    for (lx in 1:n_lambdas) {
      next_lam_label = lambda_info_df$label[lx];

      for (rf_ix in 1:n_rf_vars) {
        next_rf_var = reg_freq_vars[rf_ix];

        plot_dfs[[next_lam_label]][[next_rf_var]] = cbind(
          plot_dfs[[next_lam_label]][[next_rf_var]],
          next_day_plot_df[[next_lam_label]][[next_rf_var]]);
      }
    }
  }
}

# Plots background color.
bg_color = rgb(245, 245, 245, maxColorValue=255);

for (lx in 1:n_lambdas) {
  next_lam_label = lambda_info_df$label[lx];

  for (rf_ix in 1:n_rf_vars) {
    next_rf_var = reg_freq_vars[rf_ix];

    next_plot_df = plot_dfs[[next_lam_label]][[next_rf_var]];
    next_plot_df$x = ((0:(nrow(next_plot_df) - 1) ) * thinning_factor) + 1

    next_comparison_plot = (
      ggplot() +
      theme(panel.background=element_rect(fill=bg_color) ));

    for (dx in 1:n_days) {
      next_day_label = day_info_df$label[dx];
      legend_line_label = paste0(
        str_replace_all(next_day_label, "_", " "),
        " mean fitted");

      legend_ribbon_label = paste0(
        str_replace_all(next_day_label, "_", " "),
        " two std. dev.");

      next_comparison_plot = (
        next_comparison_plot +
        geom_line(
          data=next_plot_df,
          aes(
            x=x,
            y=.data[[paste0(next_day_label, "_mean")]],
            colour=legend_line_label),
          linewidth=1) +
        geom_ribbon(
          data=next_plot_df,
          aes(
            x=x,
            y=.data[[paste0(next_day_label, "_mean")]],
            ymin=.data[[paste0(next_day_label, "_min")]],
            ymax=.data[[paste0(next_day_label, "_max")]],
            colour=legend_ribbon_label),
          alpha=0.2,
          show.legend=FALSE) );
    }

    next_comparison_plot = (
      next_comparison_plot +
      scale_colour_manual(
        values=c("gray", "black", "red", "red", "blue", "blue"),
        guide=guide_legend(
          title="Legend",
          override.aes=list(
            fill=c("gray40", bg_color, "red", bg_color, "blue", bg_color),
            linetype=c("solid", "dashed", "solid", "dotted", "solid", "solid"),
            linewidth=c(6, 1, 6, 1, 6, 1) ))) +
      labs(
        x="Time (s)",
        y="Fitted Power") );

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
  }
}



next_lam_label = "lam_min";
next_rf_var = "Acc_freq_band_2";

plot(plot_dfs[[next_lam_label]][[next_rf_var]]$Day_1_mean, type='l', ylim=c(28, 32) );
lines(plot_dfs[[next_lam_label]][[next_rf_var]]$Day_1_min);
lines(plot_dfs[[next_lam_label]][[next_rf_var]]$Day_1_max);
lines(plot_dfs[[next_lam_label]][[next_rf_var]]$Day_2_mean, col="steelblue");
lines(plot_dfs[[next_lam_label]][[next_rf_var]]$Day_2_min, col="steelblue");
lines(plot_dfs[[next_lam_label]][[next_rf_var]]$Day_2_max, col="steelblue");
lines(plot_dfs[[next_lam_label]][[next_rf_var]]$Day_3_mean, col="green");
lines(plot_dfs[[next_lam_label]][[next_rf_var]]$Day_3_min, col="green");
lines(plot_dfs[[next_lam_label]][[next_rf_var]]$Day_3_max, col="green");
abline(v=121);
