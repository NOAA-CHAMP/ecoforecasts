function stn = station_wind_at_height(stn,wfld,dfld,afld,wz,u10fld,newz)
%function stn = station_wind_at_height(stn,wfld,dfld,afld,wz,u10fld,newz)
%
% Use wind speed time series STN.(WFLD) (Kts), at height WZ [m], to calculate
% wind velocity [m/s] time series STN.(U10FLD) (DEFAULT: [WFLD,'_10m']), at
% height NEWZ [m] (DEFAULT: 10). If optional fieldname for air temperature
% [oC] STN.(AFLD) is given, use it in calculation. If WZ is empty or missing,
% call STATION_INSTRUMENT_HEIGHTS(STN.station_name) to try to look it up. If
% wind direction DFLD is a valid field of STN, add STN.([U10FLD,'_u/v'] also. 
%
% CALLS: INTERSECT_ALL_DATES (Ecoforecasts Toolbox); SPSHFTTC (Air_Sea).
%
% Last Saved Time-stamp: <Thu 2013-12-26 16:42:03 Eastern Standard Time gramer>


  if ( ~exist('wfld','var') || ~isfield(stn,wfld) || ~is_valid_ts(stn.(wfld)) )
    error('Second arg WFLD must name a valid time series field of STN!');
  end;

  if ( ~exist('newz','var') || isempty(newz) )
    newz = 10;
  end;
  if ( ~exist('u10fld','var') || isempty(u10fld) )
    u10fld = [wfld,'_',num2str(newz),'m'];
  end;
  if ( ~isvarname(u10fld) )
    error('U10FLD must be a valid fieldname!');
  end;
  if ( isfield(stn,u10fld) )
    warning('Ecoforecasts:WindAtHeight:Overwrite',...
            'Field STN.(%s) already exists and will be overwritten',u10fld);
    stn.(u10fld) = [];
  end;


  % Height of wind sensor [m]
  if ( ~exist('wz','var') || isempty(wz) || ~isnumeric(wz) )
    [wz,az,pz,stz,dtz,slz,dlz] = station_instrument_heights(stn.station_name);
  end;
  if ( wz == newz )
    warning('Ecoforecasts:WindAtHeight:Confusion',...
            'STATION_INSTRUMENT_HEIGHTS says wind speed already at %g m??',newz);
  end;

  % Wind speed at NEWZ m [ASSUMES WFLD is in KTS, converts M/S]
  if ( exist('afld','var') && isfield(stn,afld) )
    [wix,aix] = intersect_dates(stn.(wfld).date,stn.(afld).date);
    w.date = stn.(wfld).date(wix); w.data = kts2mps(stn.(wfld).data(wix));
    a.date = stn.(afld).date(aix); a.data = stn.(afld).data(aix);

    stn.(u10fld).date = w.date;
    stn.(u10fld).data = spshfttc(w.data,wz,newz,a.data);
  else
    stn.(u10fld).date = stn.(wfld).date;
    stn.(u10fld).data = spshfttc(kts2mps(stn.(wfld).data(wix)),wz,newz);
  end;


  if ( exist('dfld','var') && isfield(stn,dfld) )

    [tix,dix] = intersect_dates(stn.(u10fld).date,stn.(dfld).date);
    t.date = stn.(u10fld).date(tix); t.data = stn.(u10fld).data(tix);
    d.date = stn.(dfld).date(dix); d.data = stn.(dfld).data(dix);

    u10xfld = [u10fld '_u'];
    u10yfld = [u10fld '_v'];
    stn.(u10xfld).date = t.date;
    stn.(u10yfld).date = t.date;
    [stn.(u10xfld).data,stn.(u10yfld).data] = spddir_to_uv(t.data,d.data);

  end;

return;
