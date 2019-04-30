function var = translate_var_name(colhdr)
%function var = translate_var_name(colhdr)
%
% Convert some header strings commonly used in 10-char column data files,
% into their equivalent ICON/G2 variable names, e.g., the column heading
% 'PAR1' (assigned by ~daps/cfg/SeakeysDataParser.xml) usually translates to
% 'licor_shallow_par' for SEAKEYS/C-MAN stations, while 'PAR-1m' (from
% ~daps/cfg/AllDataParser.xml for ICON/CREWS stations) normally becomes
% 'bic_shallow_par' for the BIC instruments deployed on ICON sticks.
%
% Last Saved Time-stamp: <Tue 2015-04-28 14:06:46 Eastern Daylight Time gramer>
%

  % Column headers from 10-char col files, Web-XLS files, etc.
  vars = ...
      { ...
          'AirT' ,			'air_t'  ; ...
          'AirT_60' ,			'air_t'  ; ...
          'AirT-60' ,			'air_t'  ; ...
          'AirT_60m' ,			'air_t'  ; ...
          'AirT-60m' ,			'air_t'  ; ...
          'AirT1' ,			'air_t'  ; ...
          'AirT1_60' ,			'air_t'  ; ...
          'AirT1-60' ,			'air_t'  ; ...
          'AirT1_60m' ,			'air_t'  ; ...
          'AirT1-60m' ,			'air_t'  ; ...
          'AirT2' ,			'wxt_air_t'  ; ...
          'AirT2_60' ,			'wxt_air_t'  ; ...
          'AirT2-60' ,			'wxt_air_t'  ; ...
          'AirT2_60m' ,			'wxt_air_t'  ; ...
          'AirT2-60m' ,			'wxt_air_t'  ; ...
          'AT' ,			'air_t'  ; ...
          'Baro' ,			'barom'  ; ...
          'Baro_60' ,			'barom'  ; ...
          'Baro-60' ,			'barom'  ; ...
          'Baro_60m' ,			'barom'  ; ...
          'Baro-60m' ,			'barom'  ; ...
          'Baro1' ,			'barom'  ; ...
          'Baro1_60' ,			'barom'  ; ...
          'Baro1-60' ,			'barom'  ; ...
          'Baro1_60m' ,			'barom'  ; ...
          'Baro1-60m' ,			'barom'  ; ...
          'Baro2' ,			'wxt_barom'  ; ...
          'Baro2_60' ,			'wxt_barom'  ; ...
          'Baro2-60' ,			'wxt_barom'  ; ...
          'Baro2_60m' ,			'wxt_barom'  ; ...
          'Baro2-60m' ,			'wxt_barom'  ; ...
          'BICS305',			'bic_surf_305nm'; ...
          'BICS330',			'bic_surf_330nm'; ...
          'BICS380',			'bic_surf_380nm'; ...
          'BICSCnt',			'bic_surf_count'; ...
          'BICU305',			'bic_shallow_305nm'; ...
          'BICU330',			'bic_shallow_330nm'; ...
          'BICU380',			'bic_shallow_380nm'; ...
          'BICUCnt',			'bic_shallow_count'; ...
          'BICD305',			'bic_deep_305nm'; ...
          'BICD330',			'bic_deep_330nm'; ...
          'BICD380',			'bic_deep_380nm'; ...
          'BICDCnt',			'bic_deep_count'; ...
          'Fluor' ,			'fluorometer_vsignal'  ; ...
          'PAR0' ,			'licor_surf_par' ; ...
          'PAR1' ,			'licor_shallow_par' ; ...
          'PAR-S' ,			'bic_surf_par' ; ...
          'PAR-U' ,			'bic_shallow_par' ; ...
          'PAR-1' ,			'bic_shallow_par' ; ...
          'PAR-1m' ,			'bic_shallow_par' ; ...
          'PAR-D' ,			'bic_deep_par' ; ...
          'SLBaro1' ,			'barom_surf'  ; ...
          'SSC' , 			'ct_shallow_cond'  ; ...
          'SSS' , 			'ct_shallow_salinity'  ; ...
          'SST' , 			'ct_shallow_seatemp'  ; ...
          'Transmis' , 			'xmissometer_vsignal'  ; ...
          'WaterHgt' ,			'tide'  ; ...
          'WaterT' ,			'sea_t'  ; ...
          'WDir' ,			'wind_dir'  ; ...
          'WDir-60' ,			'wind_dir'  ; ...
          'WDir_60' ,			'wind_dir'  ; ...
          'WDc-60' ,			'wind_dir'  ; ...
          'WDc_60' ,			'wind_dir'  ; ...
          'WDc-60m' ,			'wind_dir'  ; ...
          'WDc_60m' ,			'wind_dir'  ; ...
          'WDir-60m' ,			'wind_dir'  ; ...
          'WDir_60m' ,			'wind_dir'  ; ...
          'WDir1' ,			'wind1_dir'  ; ...
          'WDir1-60' ,			'wind1_dir'  ; ...
          'WDir1_60' ,			'wind1_dir'  ; ...
          'WDc1-60' ,			'wind1_dir'  ; ...
          'WDc1_60' ,			'wind1_dir'  ; ...
          'WDc1-60m' ,			'wind1_dir'  ; ...
          'WDc1_60m' ,			'wind1_dir'  ; ...
          'WDir1-60m' ,			'wind1_dir'  ; ...
          'WDir1_60m' ,			'wind1_dir'  ; ...
          'WDir2' ,			'wind2_dir'  ; ...
          'WDir2-60' ,			'wind2_dir'  ; ...
          'WDir2_60' ,			'wind2_dir'  ; ...
          'WDc2-60' ,			'wind2_dir'  ; ...
          'WDc2_60' ,			'wind2_dir'  ; ...
          'WDc2-60m' ,			'wind2_dir'  ; ...
          'WDc2_60m' ,			'wind2_dir'  ; ...
          'WDir2-60m' ,			'wind2_dir'  ; ...
          'WDir2_60m' ,			'wind2_dir'  ; ...
          'WG' ,			'wind_gust'  ; ...
          'WG-01m' ,			'wind_gust'  ; ...
          'WG_01m' ,			'wind_gust'  ; ...
          'WG-01mSp' ,			'wind_gust'  ; ...
          'WG_01mSp' ,			'wind_gust'  ; ...
          'WG1' ,			'wind1_gust'  ; ...
          'WG1-01m' ,			'wind1_gust'  ; ...
          'WG1_01m' ,			'wind1_gust'  ; ...
          'WG1-01mSp' ,			'wind1_gust'  ; ...
          'WG1_01mSp' ,			'wind1_gust'  ; ...
          'WG2' ,			'wind2_gust'  ; ...
          'WG2-01m' ,			'wind2_gust'  ; ...
          'WG2_01m' ,			'wind2_gust'  ; ...
          'WG2-01mSp' ,			'wind2_gust'  ; ...
          'WG2_01mSp' ,			'wind2_gust'  ; ...
          'WS' ,			'wind_speed'  ; ...
          'WS-60' ,			'wind_speed'  ; ...
          'WS_60' ,			'wind_speed'  ; ...
          'WS-60m' ,			'wind_speed'  ; ...
          'WS_60m' ,			'wind_speed'  ; ...
          'WS1' ,			'wind1_speed'  ; ...
          'WS1-60' ,			'wind1_speed'  ; ...
          'WS1_60' ,			'wind1_speed'  ; ...
          'WS1-60m' ,			'wind1_speed'  ; ...
          'WS1_60m' ,			'wind1_speed'  ; ...
          'WS2' ,			'wind2_speed'  ; ...
          'WS2-60' ,			'wind2_speed'  ; ...
          'WS2_60' ,			'wind2_speed'  ; ...
          'WS2-60m' ,			'wind2_speed'  ; ...
          'WS2_60m' ,			'wind2_speed'  ; ...
          ...
          'Cond' ,			'ctd_shallow_cond'  ; ...
          'SeaT' ,			'ctd_shallow_seatemp'  ; ...
          'WT' ,			'ctd_shallow_seatemp'  ; ...
          'SeaPres' ,			'ctd_shallow_seapres'  ; ...
          'IDepth' ,			'ctd_shallow_i_depth'  ; ...
          'I-Depth' ,			'ctd_shallow_i_depth'  ; ...
          'I_Depth' ,			'ctd_shallow_i_depth'  ; ...
          'Sal' ,			'ctd_shallow_salinity'  ; ...
          ...
          'Cond3m' ,			'ctd_deep_cond'  ; ...
          'SeaT3m' ,			'ctd_deep_seatemp'  ; ...
          'SeaPres3m' ,			'ctd_deep_seapres'  ; ...
          'IDepth3m' ,			'ctd_deep_i_depth'  ; ...
          'I-Depth3m' ,			'ctd_deep_i_depth'  ; ...
          'I_Depth3m' ,			'ctd_deep_i_depth'  ; ...
          'Sal3m' ,			'ctd_deep_salinity'  ; ...
          ...
          'SBSeaT' ,			'seabird_seatemp'  ; ...
          'SBPresDB' ,			'seabird_seapres'  ; ...
          'SBDepth' ,			'seabird_i_depth'  ; ...
          ...
          'WXTGust' ,			'wxt_wgust'  ; ...
          'WXTGu_01m' ,			'wxt_wgust'  ; ...
          'WXTSpd' ,			'wxt_wspeed'  ; ...
          'WXTSp-60' ,			'wxt_wspeed'  ; ...
          'WXTSp_60' ,			'wxt_wspeed'  ; ...
          'WXTSp-60m' ,			'wxt_wspeed'  ; ...
          'WXTSp_60m' ,			'wxt_wspeed'  ; ...
          'WXTDir' ,			'wxt_wdir'  ; ...
          'WXTDr-60' ,			'wxt_wdir'  ; ...
          'WXTDr_60' ,			'wxt_wdir'  ; ...
          'WXTDr-60m' ,			'wxt_wdir'  ; ...
          'WXTDr_60m' ,			'wxt_wdir'  ; ...
          'WXTAirT' ,			'wxt_air_t'  ; ...
          'WXTAT-60' ,			'wxt_air_t'  ; ...
          'WXTAT_60' ,			'wxt_air_t'  ; ...
          'WXTAT-60m' ,			'wxt_air_t'  ; ...
          'WXTAT_60m' ,			'wxt_air_t'  ; ...
          'WXTHumid' ,			'wxt_humid'  ; ...
          'Humid-60' ,			'wxt_humid'  ; ...
          'Humid_60' ,			'wxt_humid'  ; ...
          'Humid-60m' ,			'wxt_humid'  ; ...
          'Humid_60m' ,			'wxt_humid'  ; ...
          'WXTBaro' ,			'wxt_barom'  ; ...
          'WXTBa-60' ,			'wxt_barom'  ; ...
          'WXTBa_60' ,			'wxt_barom'  ; ...
          'WXTBa-60m' ,			'wxt_barom'  ; ...
          'WXTBa_60m' ,			'wxt_barom'  ; ...
          'RainAmt' ,			'wxt_rain_amt'  ; ...
          'RainDur' ,			'wxt_rain_dur'  ; ...
          'RnMaxInt' ,			'wxt_rain_max_int'  ; ...
          'RnAvgInt' ,			'wxt_rain_avg_int'  ; ...
          'HailAmt' ,			'wxt_hail_amt'  ; ...
          'HailDur' ,			'wxt_hail_dur'  ; ...
          'HlMaxInt' ,			'wxt_hail_max_int'  ; ...
          'HlAvgInt' ,			'wxt_hail_avg_int'  ; ...
          ...
          'SAMI_Temp',			'sami_seatemp' ; ...
          'SAMI_pCO2',			'sami_pco2' ; ...
          'SAMI_i434',			'sami_i434' ; ...
          'SAMI_i620',			'sami_i620' ; ...
          'SAMI_i740',			'sami_i740' ; ...
          'SAMI_k434',			'sami_k434' ; ...
          'SAMI_k620',			'sami_k620' ; ...
          ...
          'Air_Temp',			'air_t' ; ...
          'Barometer',			'barom' ; ...
          'CT_Voltage',			'ct_shallow_volts' ; ...
          'CT_Voltage_2m',		'ct_deep_volts' ; ...
          'Conductivity',		'ct_shallow_cond' ; ...
          'Conductivity_1m',		'ct_shallow_cond' ; ...
          'Conductivity_2m',		'ct_deep_cond' ; ...
          'DewPt',			'dew_t' ; ...
          'DewPt-60',			'dew_t' ; ...
          'DewPt_60',			'dew_t' ; ...
          'DewPt-60m',			'dew_t' ; ...
          'DewPt_60m',			'dew_t' ; ...
          'Dew_Point',			'dew_t' ; ...
          'Fluorometry',		'fluorometer_vsignal' ; ...
          'IDepth',			'wave_i_depth' ; ...
          'TidalHgt',			'wave_tide_height' ; ...
          'Panel_Temp',			'logger_ptemp' ; ...
          'PAR_1m',			'licor_shallow_par' ; ...
          'PAR_3m',			'licor_deep_par' ; ...
          'PAR_Surface',		'licor_surf_par' ; ...
          'Salinity',			'ct_shallow_salinity' ; ...
          'Salinity_1m',		'ct_shallow_salinity' ; ...
          'Salinity_2m',		'ct_deep_salinity' ; ...
          'Sea_Temp',			'ct_shallow_seatemp' ; ...
          'Sea_Temp_1m',		'ct_shallow_seatemp' ; ...
          'Sea_Temp_2m',		'ct_deep_seatemp' ; ...
          'Sea_Temp_NDBC',		'sea_t' ; ...
          'Signal_Strength',		'logger_psig' ; ...
          'Tide',			'tide' ; ...
          'Transmissometry',		'xmissometer_vsignal' ; ...
          'UVB_Radiation_1m',		'solar_shallow_uvb' ; ...
          'UVB-1',			'solar_shallow_uvb' ; ...
          'UVB_1',			'solar_shallow_uvb' ; ...
          'UVB-1m',			'solar_shallow_uvb' ; ...
          'UVB_1m',			'solar_shallow_uvb' ; ...
          'UVB_Radiation_Surface',	'solar_surf_uvb' ; ...
          'UVB-S',			'solar_surf_uvb' ; ...
          'UVB_S',			'solar_surf_uvb' ; ...
          'WaveHgt',			'wave_height' ; ...
          'Wind_Dir_1',			'wind1_dir' ; ...
          'Wind_Gust_1',		'wind1_gust' ; ...
          'Wind_Speed_1',		'wind1_speed' ; ...
          'Wind_Dir_2',			'wind2_dir' ; ...
          'Wind_Gust_2',		'wind2_gust' ; ...
          'Wind_Speed_2',		'wind2_speed' ; ...
      };

  % If no translation found, use the column heading instead
  var = colhdr;
  idx = find(strcmpi(vars(:, 1), var));
  if ( ~isempty(idx) )
    var = vars{idx, 2};
  end;

return;
