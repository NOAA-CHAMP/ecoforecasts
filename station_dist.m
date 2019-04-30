function [d,degT] = station_dist(stnm1,stnm2)
%function [d,degT] = station_dist(stnm1,stnm2)
%
% Return distance in km between two monitoring stations STNM1 and STNM2: each
% may be 5-char station code strings, 2-vectors [LAT,LON], or STRUCTs with
% fields .lon and .lat. Name strings are passed to GET_STATION_COORDS (v.)
%
% CALLS: SW_DIST (from Sw - the Sea Water Toolkit)
%
% Last Saved Time-stamp: <Sat 2017-08-26 18:09:25 Eastern Daylight Time gramer>

  d = nan;

  if ( isfield(stnm1,'lon') && isfield(stnm1,'lat') )
    lon1 = stnm1.lon; lat1 = stnm1.lat;
  elseif ( isnumeric(stnm1) && numel(stnm1) == 2 )
    lon1 = stnm1(2); lat1 = stnm1(1);
  else
    [lon1,lat1,ig] = get_station_coords(stnm1);
  end;
  if ( isfield(stnm2,'lon') && isfield(stnm2,'lat') )
    lon2 = stnm2.lon; lat2 = stnm2.lat;
  elseif ( isnumeric(stnm2) && numel(stnm2) == 2 )
    lon2 = stnm2(2); lat2 = stnm2(1);
  else
    [lon2,lat2,ig] = get_station_coords(stnm2);
  end;

  [d,cwE] = sw_dist([lat1 lat2], [lon1 lon2], 'km');

  if ( nargout > 1 )
    degT = 90 - cwE;
    degT(degT<0) = 360 + degT(degT<0);
  end;

return;
