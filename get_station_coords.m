function [lon,lat,dep] = get_station_coords(stnm)
%function [lon,lat,dep] = get_station_coords(stnm)
%
% Return station coordinates and depth for any known reef monitoring station
%
% Last Saved Time-stamp: <Thu 2010-02-18 13:16:42 Eastern Standard Time gramer>

  STATIONS = get_all_station_metadata();

  ix = find(strcmpi(STATIONS.codes,stnm));
  if ( isempty(ix) )
    error('Do not know coordinates for station "%s"!', stnm);
  else
    lon = STATIONS.lons(ix);
    lat = STATIONS.lats(ix);
    dep = STATIONS.depths(ix);
  end;

return;
