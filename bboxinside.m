function [lonix,latix] = bboxinside(lon,lat,bbox,domesh,tol)
%function [lonix,latix] = bboxinside(lon,lat,bbox,domesh,tol)
% 
% Return indices in the vectors of points LON,LAT, inside or "near" bounding
% box BBOX ([lon1,lon2,lat1,lat2]). If DOMESH is true (DEFAULT), calls INSIDE
% on matrix MESHGRID(LON,LAT). TOL (DEFAULT: 1e-6) is passed to INSIDE: if 0,
% return only *interior* indices. CALLS: INSIDE (Oceans toolbox), MESHGRID.
% 
% Last Saved Time-stamp: <Mon 2017-11-20 13:21:20 Eastern Standard Time gramer>

  if ( ~exist('domesh') || isempty(domesh) )
    domesh = true;
  end;
  if ( ~exist('tol') || isempty(tol) )
    tol = 1e-6;
  end;

  XV = [bbox(1),bbox(2),bbox(2),bbox(1)];
  YV = [bbox(3),bbox(3),bbox(4),bbox(4)];

  if ( domesh )
    [alllon,alllat] = meshgrid(lon,lat);
    ix = find( inside(alllon(:),alllat(:),XV,YV,tol) > 0 );
    [lonix,latix] = ind2sub(size(alllon),ix);
  elseif ( ndims(lon) == ndims(lat) && all(size(lon) == size(lat)) )
    lonix = find( inside(lon(:),lat(:),XV,YV,tol) > 0 );
    latix = lonix;
  else
    error('If DOMESH is not given as True, SIZE(LON) must equal SIZE(LAT)');
  end;

return;
