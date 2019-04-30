function [neareststn,dst,ang] = find_nearest_station(stnm_or_stn_or_loc,excludeself)
%function [neareststn,dst,ang] = find_nearest_station(stnm_or_stn_or_loc,excludeself)
%
% Find the nearest station from GET_ALL_STATION_METADATA, to the given named
% station, station STRUCT with fields .lon and .lat, or [LON,LAT] 2-vector.
% If EXCLUDESELF (DEFAULT: False), the nearest *other* station is returned.
% NOTE: If caller specifies station by name, EXCLUDESELF defaults to True.
%
% Last Saved Time-stamp: <Sun 2016-07-31 16:33:15 Eastern Daylight Time gramer>

  
  if ( ~exist('excludeself','var') )
    excludeself = [];
  end;

  if ( ischar(stnm_or_stn_or_loc) )
    stn = get_station_from_station_name(stnm_or_stn_or_loc);
    if ( isempty(excludeself) )
      excludeself = true;
    end;
  elseif ( isnumeric(stnm_or_stn_or_loc) && numel(stnm_or_stn_or_loc) == 2 )
    stn.lon = stnm_or_stn_or_loc(1);
    stn.lat = stnm_or_stn_or_loc(2);
  else
    stn = stnm_or_stn_or_loc;
  end;
  clear stnm_or_stn_or_loc

  if ( isempty(excludeself) )
    excludeself = false;
  end;

  if ( ~isfield(stn,'lon') || ~isfield(stn,'lat') )
    error('No lon,lat found');
  end;

  STATIONS = get_all_station_metadata;

  [dsts,angs] = distance_wgs84(stn.lat,stn.lon,STATIONS.lats,STATIONS.lons);
  if ( excludeself )
    [ig,ix] = min(dsts);
    dsts(ix) = [];
    angs(ix) = [];
  end;

  [dst,ix] = min(dsts);
  ang = angs(ix);

  neareststn = get_station_from_station_name(STATIONS.codes{ix});

return;
