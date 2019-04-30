function stn = load_cleaning_dates(stn_or_stnm)
%function stn = load_cleaning_dates(stn_or_stnm)
%
% Load cleaning dates for station struct or station name STN_OR_STNM. Adds a
% DATENUM vector field STN.cleaning_date upon return.
%
% Last Saved Time-stamp: <Wed 2011-12-07 15:59:09  Lew.Gramer>

  ecopath = get_ecoforecasts_path('.');

  stn = get_station_from_station_name(stn_or_stnm);
  stnm = lower(stn.station_name);

  if ( strcmp(stnm,'mlrf1') )
    stnm = 'mlrf2';
  end;

  clngfname = fullfile(ecopath,[stnm,'-cleaning-dates.csv']);
  if ( ~exist(clngfname,'file') )
    warning('Cannot find "%s"!',clngfname);
  else
    x = importdata(clngfname);
    stn.cleaning_date = datenum(x);
    x=[]; clear x
  end;

return;
