function [minlon,maxlon,dlon,minlat,maxlat,dlat,lons,lats,LONS,LATS] = oversample_sst_image_coords(region,subs)
%function [minlon,maxlon,dlon,minlat,maxlat,dlat,lons,lats,LONS,LATS] = oversample_sst_image_coords(region,subs)
%
% Return bounding and grid coordinates and resolution for University of South
% Florida Sea Surface Temperature fields in region REGION, oversampled by
% fractional resolution SUBS (see oversample_sst_image for meaning of SUBS).
%
% Returned values are identical to those from GET_SST_IMAGE_COORDS (v.), but
% at the new oversampled resolution SUBS.
%
% Last Saved Time-stamp: <Sun 2011-03-06 12:48:58  Lew.Gramer>

  [minlon,maxlon,dlon,minlat,maxlat,dlat,lons,lats,LONS,LATS] = get_sst_image_coords(region);

  dlon = dlon*subs;
  dlat = dlat*subs;
  lons = interp1(lons,1:subs:length(lons));
  lats = interp1(lats,1:subs:length(lats));
  [LONS,LATS] = meshgrid(lons,lats);

return;
