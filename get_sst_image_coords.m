function [minlon,maxlon,dlon,minlat,maxlat,dlat,lons,lats,LONS,LATS] = get_sst_image_coords(region)
%function [minlon,maxlon,dlon,minlat,maxlat,dlat,lons,lats,LONS,LATS] = get_sst_image_coords(region)
%
% Return bounding box, resolution, and coordinate vectors and (MESHGRID)
% coordinate matrices for University of South Florida SST dataset REGION.
%
% REGION values currently supported: 'florida','keys_sst'
%
% Last Saved Time-stamp: <Mon 2011-02-28 15:04:08 Eastern Standard Time gramer>

  switch (region),
   case 'florida',
    minlon = -91; maxlon = -79; dlon = (1/110);
    minlat =  22; maxlat =  31; dlat = (1/110);
   case 'keys_sst',
    minlon = -83.5; maxlon = -79.5; dlon = (1/110);
    minlat =  24; maxlat =  27.5; dlat = (1/110);
   otherwise,
    error('Do not know coords for region "%s"!',region);
  end;
  
  lons = minlon:dlon:maxlon;
  lats = maxlat:(-dlat):minlat;
  %%%% ??? HACK! Arithmetic above gives a 991x1321 grid for
  %%%% ??? HACK!   'florida', but SST is actually 990x1320
  lons(end) = [];
  lats(end) = [];
  % lons(1) = [];
  % lats(1) = [];

  [LONS, LATS] = meshgrid(lons, lats);

return;
