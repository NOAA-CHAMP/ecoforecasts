function [LATOUT,LONOUT] = reckon_wgs84(lat,lon,rng,az)
%function [LATOUT,LONOUT] = reckon_wgs84(lat,lon,rng,az)
%
% Calls RECKON (v.) with WGS-84 ellipticity and polar radius to calculate
% a linear track of lat/lon points RNG [km] apart, in heading AZ (degT).
% Args may be scalars or arrays; if any two or more of the args are not
% scalars, then the dimensions of those args must match each other exactly.
% For a mix of scalar and NxM inputs, LATOUT and LONOUT will be NxM also.
%
% Last Saved Time-stamp: <Tue 2011-11-01 14:30:14  lew.gramer>

  [LATOUT,LONOUT] = reckon(lat,lon, rng,az, wgs84_ellipsoid());

return;
