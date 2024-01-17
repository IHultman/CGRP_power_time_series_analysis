library(stringr)


data_dir = paste0(
  "/home",
  "/ikhultman",
  "/Desktop",
  "/power_time_series_analysis",
  "/CGRP_power_time_series_analysis",
  "/time_series_dataframes/");

data_files = list.files(data_dir);
data_files_table = data.frame(
  filename=paste0(data_dir, data_files),
  region=str_extract(data_files, ".*(?=_ts_data\\.csv)") );

n_files = nrow(data_files_table);

expected_vars = c(
  "y",
  "mouse",
  "freq_band",
  "phase",
  "time",
  "day");

for (fx in 1:n_files) {
  region_data_next = read.csv(data_files_table$filename[fx], header=TRUE);

  stopifnot(
    (length(
      setdiff(
        expected_vars,
        names(region_data_next) )) == 0) &&
    (length(
      setdiff(
        names(region_data_next),
        expected_vars) ) == 0) );

  freq_bands_of_interest = c(1, 3);
  freq_row_ixs = is.element(region_data_next$freq_band, freq_bands_of_interest);

  exp_phases_of_interest = c("Baseline", "Post_injection", "PI2");
  phase_row_ixs = is.element(region_data_next$phase, exp_phases_of_interest);

  na_ixs = is.na(region_data_next$y);
  region_data_next = region_data_next[freq_row_ixs & phase_row_ixs & !na_ixs,];

  rownames(region_data_next) = NULL;

  region_data_next$cgrp = 0;
  cgrp_ixs = (
    (region_data_next$day == 2) &
    ((region_data_next$phase == "Post_injection") |
     (region_data_next$phase == "PI2") ));

  region_data_next$cgrp[cgrp_ixs] = 1;

  #region_data_next$mouse_day = factor(paste(region_data_next$mouse, region_data_next$day, sep="_day") );

  region_data_next$time = scale(region_data_next$time);
  region_data_next$region = data_files_table$region[fx];

  #time_bspline = data.frame(bs(region_data_next$time) );
  #colnames(time_bspline) = c("time_bs1", "time_bs2", "time_bs3");
  #region_data_next = cbind(region_data_next, time_bspline);

  if (exists("region_data_all") ) {
    stopifnot(all(colnames(region_data_all) == colnames(region_data_next) ));

    region_data_all = rbind(region_data_all, region_data_next);
  } else {
    region_data_all = region_data_next;
  }
}

save_filename = paste0(
  "/home",
  "/ikhultman",
  "/Desktop",
  "/power_time_series_analysis",
  "/CGRP_power_time_series_analysis",
  "/combined_time_series_data",
  "/combined_brain_regions_3_phases_ts_data.csv");

write.csv(
  region_data_all,
  save_filename,
  quote=FALSE,
  row.names=FALSE);

