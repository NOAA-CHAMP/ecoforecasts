function stn = station_xtide(stn_or_stnm,dts)
%function stn = station_xtide(stn_or_stnm,dts)
%
% Use T_TIDE toolbox function T_XTIDE (v.) to predict tide height in [FEET]
% at station struct STN or station named STNM. If vector DTS is given, only
% predict for those dates. Otherwise, try to find an appropriate date range
% among the other fields in STN; if none found, predict for 1987 to present.
% Struct returned with fields STN.XTIDE and STN.XTIDE_M (in [m]).
%
% WARNING: If STN has no tide gauge, T_XTIDE will only return data for the
% reporting tide station *nearest to* STN: this may be *very approximate*!
%
% Last Saved Time-stamp: <Tue 2010-05-25 11:37:25 Eastern Daylight Time lew.gramer>

  datapath = get_ecoforecasts_path('data');

  if ( isstruct(stn_or_stnm) )
    stn = stn_or_stnm;
    if ( ~isfield(stn,'station_name') )
      error('STN struct must have "station_name" field!');
    end;
    stnm = stn.station_name;
  elseif ( ischar(stn_or_stnm) )
    stnm = stn_or_stnm;
    stn.station_name = stnm;
  else
    error('First arg must be a station struct or station name string!');
  end;

  if ( ~exist('dts','var') || isempty(dts) )
    dts = [];
    tide_fld = 'ctd_deep_i_depth';
    if (~isfield(stn,tide_fld)); tide_fld = 'ctd_shallow_i_depth';	end;
    if (~isfield(stn,tide_fld)); tide_fld = 'tide';			end;
    if (~isfield(stn,tide_fld)); tide_fld = 'ndbc_tide';		end;
    if (~isfield(stn,tide_fld)); tide_fld = 'wave_tide_height';		end;
    if (~isfield(stn,tide_fld)); tide_fld = 'sea_t';			end;
    if (~isfield(stn,tide_fld)); tide_fld = 'ct_shallow_seatemp';	end;
    if (~isfield(stn,tide_fld))
      warning('No DTS given, no comparison field found! Predicting 1987-present...');
      dts = datenum(1987,1,1):(1/24):now;
    else
      dts = stn.(tide_fld).date(1):(1/24):stn.(tide_fld).date(end);
    end;
  end;

  if ( ~isfield(stn,'lon') )
    [stn.lon,stn.lat,stn.depth] = get_station_coords(stn.station_name);
  end;
  if ( ~isfield(stn,'lon') )
    error('Cannot find coordinates for station "%s"!', stn.station_name);
  end;

  stn.xtide.date = dts';
  stn.xtide.data = t_xtide(stn.lon,stn.lat,dts, 'units','feet');

  stn.xtide_m = stn.xtide;
  stn.xtide_m.data = stn.xtide_m.data .* 0.3048;

return;
