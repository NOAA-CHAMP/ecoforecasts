function [minlon,maxlon,dlon,minlat,maxlat,dlat,lons,lats,LONS,LATS] = get_avhrr_coords(region)
%function [minlon,maxlon,dlon,minlat,maxlat,dlat,lons,lats,LONS,LATS] = get_avhrr_coords(region)

  if ( ~exist('region','var') || isempty(region) )
    region = 'florida';
  end;

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
