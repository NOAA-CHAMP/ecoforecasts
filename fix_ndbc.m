function station = fix_ndbc(stnm)

error('This function was a one-time fix-it: not meant to be re-reun!');

  badflds = { ...
      'Date', ...
      'Name', ...
      'Hour', ...
      'InCurr', ...
      'EndTime', ...
      'wind1_dir', ...
      'wind1_speed', ...
      'wind2_dir', ...
      'wind2_speed', ...
      'air_t', ...
      'air_t_dewp', ...
      'DewPt', ...
      'barom', ...
      'barom_surf', ...
      'sea_t', ...
      'MaxWS1', ...
      'MaxWS2', ...
      'tide', ...
      'WGDir1', ...
      'wind1_gust', ...
      'WGDir2', ...
      'wind2_gust', ...
      'MeanWD1', ...
      'MeanWS1', ...
      'MeanWD2', ...
      'MeanWS2', ...
      'WSStdD1', ...
      'WSStdD2', ...
      'licor_surf_par', ...
      'licor_shallow_par', ...
      'fluorometer_vsignal', ...
      'xmissometer_vsignal', ...
      'ct_shallow_seatemp', ...
      'ct_shallow_cond', ...
      'ct_shallow_salinity', ...
      'BattM', ...
      'FJO', ...
      'FJ0', ...
      'PAR2m', ...
      'PAR1m', ...
      'WaveHgt', ...
      'TidalHgt', ...
      'IDepth', ...
      'bic_surf_305nm', ...
      'bic_surf_330nm', ...
      'bic_surf_380nm', ...
      'bic_surf_par', ...
      'bic_surf_count', ...
      'bic_shallow_305nm', ...
      'bic_shallow_330nm', ...
      'bic_shallow_380nm', ...
      'bic_shallow_par', ...
      'bic_shallow_count', ...
      'amodis_chlor_a', ...
      'amodis_seadas_sst', ...
      'tmodis_chlor_a', ...
      'tmodis_seadas_sst', ...
      'IH6WDir1', ...
      'IH6WS1', ...
      'IH6WDir2', ...
      'IH6WS2', ...
      'IH5WDir1', ...
      'IH5WS1', ...
      'IH5WDir2', ...
      'IH5WS2', ...
      'IH4WDir1', ...
      'IH4WS1', ...
      'IH4WDir2', ...
      'IH4WS2', ...
      'IH3WDir1', ...
      'IH3WS1', ...
      'IH3WDir2', ...
      'IH3WS2', ...
      'IH2WDir1', ...
      'IH2WS1', ...
      'IH2WDir2', ...
      'IH2WS2', ...
      'IH1WDir1', ...
      'IH1WS1', ...
      'IH1WDir2', ...
      'IH1WS2', ...
      'ASens1', ...
      'UVB0', ...
      'UNKD', ...
      'UNKG1', ...
      'UNKG2', ...
      'FI2', ...
      'SST2', ...
      'SSC2', ...
      'SSS2', ...
      'BattM2', ...
            };

  station = load_all_ndbc_data([],stnm);

  for bix = 1:length(badflds)
    badfld = badflds{bix};
    if ( isfield(station,badfld) )
      station = rmfield(station,badfld);
    end;
  end;

  datapath = get_ecoforecasts_path('data');
  matfname = fullfile(datapath, [stnm '-ndbc.mat']);
  disp(['Overwriting old bloated ' matfname]);
  save(matfname, 'station');

return;
