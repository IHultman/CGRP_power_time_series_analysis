data_dir = [ ...
  '/mnt', ...
  '/nfs', ...
  '/clasnetappvm', ...
  '/students', ...
  '/ikhultman', ...
  '/Desktop', ...
  '/power_time_series_analysis', ...
  '/CGRP_power_time_series_analysis', ...
  '/data/'];

data_files = split(ls(data_dir) );
data_files = data_files(1:(end-1) );
n_files = length(data_files);

exp_phases = {'Baseline', 'Post_injection', 'PI2'};
n_phases = length(exp_phases);

phase_start_ixs = [1, 1802, 3603];

assert(length(phase_start_ixs) == n_phases);

brain_regions = { ...
  'Acc', ...
  'BLA', ...
  'CeA', ...
  'MD_thal', ...
  'Po', ...
  'VPM', ...
  'LPBN', ...
  'RPBN'};

n_regions = length(brain_regions);

expected_nonpower_colnames = [ ...
  {'mouse'}, ...
  {'date'}, ...
  {'freq_band'}, ...
  {'Hz'}, ...
  brain_regions];

power_colnames_prefix = 'Var';

region_dfs = struct();

for rx = 1:n_regions
  region_dfs.(brain_regions{rx}) = [];
end

for fx = 1:n_files
  next_file = [data_dir, data_files{fx}];
  next_power_means = load(next_file);
  
  tmp_region_dfs = struct();

  for rx = 1:n_regions
    tmp_region_dfs.(brain_regions{rx}) = [];
  end

  for px = 1:n_phases
    next_phase = exp_phases{px};
    next_tabl_colnames = next_power_means.(next_phase).Properties.VariableNames;
    [~, nonpow_col_ixs] = sort( ...
      intersect( ...
        next_tabl_colnames, ...
        expected_nonpower_colnames) );

    pow_col_ixs = setdiff( ...
      1:size(next_power_means.(next_phase), 2), ...
      nonpow_col_ixs);

    n_power = length(pow_col_ixs);

    % Make sure we've identified the columns corresponding to the mean
    % log-power time series.
    assert( ...
      all( ...
      1:n_power == ...
      cellfun( ...
        @str2double, ...
        strrep( ...
          next_tabl_colnames(pow_col_ixs), ...
          power_colnames_prefix, ...
          '') )));

    for rx = 1:n_regions
      next_region = brain_regions{rx};
      reg_ixs = logical(next_power_means.(next_phase).(next_region) );
      region_power_means = next_power_means.(next_phase)(reg_ixs,:);
      [reg_grp_ixs, mouse_id, freq_band] = findgroups( ...
        region_power_means.mouse, ...
        region_power_means.freq_band);

      if isrow(reg_grp_ixs)
        reg_grp_ixs = reg_grp_ixs'; 
      end

      if isrow(mouse_id)
        mouse_id = mouse_id';
      end

      if isrow(freq_band)
        freq_band = freq_band';
      end

      group_mean_ts = splitapply( ...
        @(reg_pow_mat) mean(reg_pow_mat, 1), ...
        table2array(region_power_means(:,pow_col_ixs) ), ...
        reg_grp_ixs);

      mouse_labels = repmat(mouse_id, 1, n_power);
      freq_labels = repmat(freq_band, 1, n_power);
      phase_labels = repelem({next_phase}, size(group_mean_ts, 1), n_power);
      time_labels = repmat(1:n_power, size(group_mean_ts, 1), 1);

      assert(all(size(group_mean_ts) == size(mouse_labels) ));
      assert(all(size(group_mean_ts) == size(freq_labels) ));
      assert(all(size(group_mean_ts) == size(phase_labels) ));
      assert(all(size(group_mean_ts) == size(time_labels) ));

      power_ts_table = table( ...
        reshape(group_mean_ts', [], 1), ...
        reshape(mouse_labels', [], 1), ...
        reshape(freq_labels', [], 1), ...
        reshape(phase_labels', [], 1), ...
        reshape(time_labels', [], 1), ...
        'VariableNames', ...
        [{'y'}, {'mouse'}, {'freq_band'}, {'phase'}, {'time'}]);
    
      if isempty(tmp_region_dfs.(next_region) )
        tmp_region_dfs.(next_region) = power_ts_table;
      else
        tmp_region_dfs.(next_region) = [ ...
          tmp_region_dfs.(next_region); ...
          power_ts_table];
      end
    end
  end
  
  for rx = 1:n_regions
    tmp_region_dfs.(brain_regions{rx}).day = fx;
  end
end

