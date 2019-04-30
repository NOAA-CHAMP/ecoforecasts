function stn = get_triaxys_waves(fname,stn_or_stnm,PFX)
%function stn = get_triaxys_waves(fname,stn_or_stnm,PFX)
%
% Read Tri-AXYS Wave buoy data summary file FNAME, and store wave measurement
% time series structs in STN.([PFX 'avgwavehgt']), etc. STN_OR_STNM if given
% can specify an existing station struct or known station 5-char code. PFX
% DEFAULTS to 'triaxys_'; if PFX is empty, no fieldname prefix is used.
%
% Last Saved Time-stamp: <Sun 2011-09-25 12:24:52  lew.gramer>

  if ( ~exist('stn_or_stnm','var') || isempty(stn_or_stnm) )
    stn = [];
  else
    stn = get_station_from_station_name(stn_or_stnm);
  end;

  if ( ~exist('PFX','var') )
    PFX = 'triaxys';
  end;
  if ( length(PFX) > 0 && PFX(end) ~= '_' )
    PFX(end+1) = '_';
  end;

  x = importdata(fname);
  dts = datenum(x.textdata(2:end,1));

  stn.([PFX 'avgwavehgt']).date = dts;
  stn.([PFX 'avgwavehgt']).data = x.data(:,4);

  stn.([PFX 'avgwaveper']).date = dts;
  stn.([PFX 'avgwaveper']).data = x.data(:,5);

  stn.([PFX 'maxwavehgt']).date = dts;
  stn.([PFX 'maxwavehgt']).data = x.data(:,6);

  stn.([PFX 'sigwavehgt']).date = dts;
  stn.([PFX 'sigwavehgt']).data = x.data(:,7);

  stn.([PFX 'sigwaveper']).date = dts;
  stn.([PFX 'sigwaveper']).data = x.data(:,8);

  stn.([PFX 'peakwaveper']).date = dts;
  stn.([PFX 'peakwaveper']).data = x.data(:,9);

  stn.([PFX 'peakwavedir']).date = dts;
  stn.([PFX 'peakwavedir']).data = x.data(:,12);

  x=[]; clear x;

return;
