function az = azimuth_wgs84(lat1,lon1,lat2,lon2)
%function az = azimuth_wgs84(lat1,lon1,lat2,lon2)
%
% Calls AZIMIUTH (v.) with WGS-84 ellipticity and polar radius to calculate
% heading(s) [degT] between pairs of lat/lon points. Args may be matrices.
% If LAT1 (LON1) scalar but LAT2 (LON2) is not, call REPMAT to match sizes.
% Otherwise, requires SIZE(LAT1) == SIZE(LAT2) and SIZE(LON1) == SIZE(LON2).
%
% Last Saved Time-stamp: <Fri 2017-02-24 16:25:05 Eastern Standard Time gramer>

  if ( isscalar(lat1) );  lat1 = repmat(lat1,size(lat2));  end;
  if ( isscalar(lon1) );  lon1 = repmat(lon1,size(lon2));  end;

  if ( ndims(lat1)~=ndims(lat2) || any(size(lat1)~=size(lat2)) || ...
       ndims(lon1)~=ndims(lon2) || any(size(lon1)~=size(lon2)) )
    error('LAT1,LON1 must be scalar or must match dimensions of LAT2,LON2!');
  end;

  az = azimuth(lat1,lon1, lat2,lon2, wgs84_ellipsoid());

return;
