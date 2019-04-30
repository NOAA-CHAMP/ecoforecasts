function [lons,lats] = transect_wgs84(lon,lat,rng,az)
%function [lons,lats] = transect_wgs84(lon,lat,rng,az)
%
% Return coordinate(s) [LONS,LATS], along a great circle transect starting at
% LON,LAT and defined by azimuth AZ and end-point distance RNG (in [km]). RNG
% may be a vector in which case [LONS,LATS] is a seq. of points on a line. 
%
% Calls RECKON (MAP Toolbox, v.) using the World Geodetic System 84 ellipsoid
% for Earth, as returned by WGS84_ELLIPSOID (v.).
%
% Last Saved Time-stamp: <Tue 2010-10-12 09:40:30 Eastern Daylight Time gramer>

  [lats,lons] = reckon(lat,lon,rng,az,wgs84_ellipsoid());

return;
