function stn = read_usf_kd_par(stn_or_stnm)
%function stn = read_usf_kd_par(stn_or_stnm)
%
% Load USF A/MODIS Kd_PAR monthly mean file for station (or station named) STN_OR_STNM
%
% Last Saved Time-stamp: <Wed 2012-02-08 12:19:24  Lew.Gramer>

  stn = get_station_from_station_name(stn_or_stnm);
  clear stn_or_stnm;

  datapath = get_ecoforecasts_path('data');
  kdpath = fullfile(datapath,'usf','Kd_PAR');

  %Lat_24.46_Lon_-81.88_month.txt
  fname = fullfile(kdpath,[lower(stn.station_name),'_usf_kd_par.dat']);
  if ( ~exist(fname,'file') )
    error('No data file found: "%s"',fname);
  end;

  rawdat = load(fname);

  % Date for beginning of each month
  dt = rawdat(:,1);
  yr = floor(dt./1e3);
  jd = dt - (yr.*1e3);
  begdts = datenum(yr,1,1) + jd - 1;

  % "Mean Date" of all available data for that month
  avgjd = rawdat(:,4);
  avgdts = datenum(2002,1,1) + avgjd - 1;

  stn.amodis_kd_par.date = avgdts;
  stn.amodis_kd_par.data = rawdat(:,5);
  stn.amodis_kd_par_sd.date = avgdts;
  stn.amodis_kd_par_sd.data = rawdat(:,6);
  stn.amodis_kd_par_min.date = avgdts;
  stn.amodis_kd_par_min.data = rawdat(:,7);
  stn.amodis_kd_par_max.date = avgdts;
  stn.amodis_kd_par_max.data = rawdat(:,8);
  stn.amodis_kd_par_clim.date = begdts;
  stn.amodis_kd_par_clim.data = rawdat(:,9);
  stn.amodis_kd_par_anom.date = avgdts;
  stn.amodis_kd_par_anom.data = rawdat(:,13);

  stn.amodis_kd_par_simple_anom.date = stn.amodis_kd_par.date;
  stn.amodis_kd_par_simple_anom.data = stn.amodis_kd_par.data - nanmean(stn.amodis_kd_par.data);

return;
