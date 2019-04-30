function [X,Y,W,H] = get_ool_pixel_coords(rgn,lats,lons)
%function [X,Y,W,H] = get_ool_pixel_coords(rgn,lats,lons)
%
% Return lower left pixel position X and Y, and pixel width and height W and
% H, for the bounding box of the coordinates LATS, LONS within the region
% named RGN associated with datasets of University of South Florida (USF)
% Optical Oceanography Lab (OOL). See GET_OOL_REGION for valid RGN strings.
%
% Returned values are suitable, e.g., for ImageMagick command geometries:
%
% >> roi = [ -170.7265 -170.6300  -14.3100  -14.2650];
% >> [X,Y,W,H] = get_ool_pixel_coords('SAMOA',roi(3:4),roi(1:2)),
% >> 121 248 42 19
% # display -crop 42x19+121+248 coral/CRCP/Sediment/SAMOA/CI/A2014093*
%
% Last Saved Time-stamp: <Wed 2018-08-29 18:15:33 Eastern Daylight Time gramer>

  [lat,lon,res,bbox] = get_ool_region(rgn);

  [laterr,Y] = min(abs(max(lats(:))-lat(:)));
  %DEBUG:  disp(laterr*111e3);
  if ( laterr > max(diff(unique(lats))) )
    error('MAX(LATS) is outside of image bounding box');
  end;
  Y = numel(lat) - Y;

  [lonerr,X] = min(abs(min(lons(:))-lon(:)));
  %DEBUG:  disp(lonerr*111e3);
  if ( lonerr > max(diff(unique(lons))) )
    error('MIN(LONS) is outside of image bounding box');
  end;

  H = numel(find(min(lats(:))<=lat & lat<=max(lats(:))));
  W = numel(find(min(lons(:))<=lon & lon<=max(lons(:))));

return;
