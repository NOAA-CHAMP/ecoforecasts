function stn = get_station_from_station_name(stn_or_stnm)
%function stn = get_station_from_station_name(stn_or_stnm)
%
% Return station struct STN with .station_name and .lon,.lat,.depth fields if
% possible, from arg STN_OR_STNM, which is a struct or station name string.
%
% Last Saved Time-stamp: <Wed 2011-05-04 06:59:42  lew.gramer>

  if ( ischar(stn_or_stnm) )
    stn.station_name = stn_or_stnm;
  elseif ( isstruct(stn_or_stnm) )
    stn = stn_or_stnm;
  else
    error('Arg should be either a station STRUCT or name string!');
  end;
  if ( isfield(stn,'station_name') )
    if ( ~isfield(stn,'lon') || ~isfield(stn,'lat') )
      try
        [stn.lon,stn.lat,stn.depth] = get_station_coords(stn.station_name);
      catch
      end;
    end;
  elseif ( ~isfield(stn,'lon') || ~isfield(stn,'lat') )
    warning('Station struct without .station_name or .lon,.lat fields');
  end;

return;
