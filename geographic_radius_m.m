function [rx,ry] = geographic_radius_m(lons,lats)
%function [rx,ry] = geographic_radius_m(lons,lats)
%
% Return the radius (range) in m of the rectangle defined by the geographic
% coordinates MIN(LONS(:)),MAX(LONS(:)) and MIN(LATS(:)),MAX(LATS(:)). Uses
% MAX range of rectangle in each direction. Return values RX,RY are suitable
% to be passed to, e.g., READ_HIRES_BATHYMETRY. Calls: DISTANCE_WGS84.
%
% Last Saved Time-stamp: <Mon 2016-10-24 20:51:38 Eastern Daylight Time gramer>

  brx = distance_wgs84(min(lats(:)),min(lons(:)),min(lats(:)),max(lons(:)));
  trx = distance_wgs84(max(lats(:)),min(lons(:)),max(lats(:)),max(lons(:)));
  rx = 1e3 * max(brx,trx) / 2;

  lry = distance_wgs84(min(lats(:)),min(lons(:)),max(lats(:)),min(lons(:)));
  rry = distance_wgs84(min(lats(:)),max(lons(:)),max(lats(:)),max(lons(:)));
  ry = 1e3 * max(lry,rry) / 2;

return;
