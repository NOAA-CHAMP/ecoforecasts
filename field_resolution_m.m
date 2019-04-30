function [dy,dx] = field_resolution_m(lons,lats)
%function [dy,dx] = field_resolution_m(lons,lats)
%
% Return the resolution (MINIMUM inter-pixel or gridpoint spacing) in [m] of
% the rectangle defined by geographic coordinates LONS,LATS. A single arg may
% also be passed in which is a STRUCT with fields .lon and .lat, or a cell
% array with vectors {LONS,LATS}. NOTE: Y-resolution (latitude) is returned
% first as a convenience where a scalar is desired. Calls: DISTANCE_WGS84.
%
% Last Saved Time-stamp: <Fri 2017-02-17 14:52:03 Eastern Standard Time gramer>

  if ( nargin < 2 )
    if ( isfield(lons,'lon') && isfield(lons,'lat') )
      lats = lons.lat;
      lons = lons.lon;
    elseif ( iscell(lons) && numel(lons) == 2 )
      lats = lons{2};
      lons = lons{1};
    else
      error('If no second arg, first arg must be STRUCT or cell array');
    end;
  end;
  if ( ~isnumeric(lons) || ~isnumeric(lats) )
    error('LONS and LATS must be numeric arrays');
  end;

  if ( isvector(lons) && isvector(lats) )
    [LONS,LATS] = meshgrid(lons,lats);
  elseif ( ndims(lons) == ndims(lats) && all(size(lons)==size(lats)) )
    LONS = lons;
    LATS = lats;
  else
    error('LONS,LATS must either both be vectors, or must be plaid');
  end;

  [minlon,ix] = min(LONS(:)); [ix,lix] = ind2sub(size(LONS),ix);
  [maxlon,ix] = max(LONS(:)); [ig,rix] = ind2sub(size(LONS),ix);
  [minlat,ix] = min(LATS(:)); [tix,ig] = ind2sub(size(LATS),ix);
  [maxlat,ix] = max(LATS(:)); [bix,ig] = ind2sub(size(LATS),ix);

  ldy = min(distance_wgs84(LATS(1:end-1,lix),minlon,LATS(2:end,rix),minlon));
  rdy = min(distance_wgs84(LATS(1:end-1,rix),maxlon,LATS(2:end,rix),maxlon));
  dy = 1e3 * min(ldy,rdy);

  bdx = min(distance_wgs84(minlat,LONS(bix,1:end-1),minlat,LONS(bix,2:end)));
  tdx = min(distance_wgs84(maxlat,LONS(tix,1:end-1),maxlat,LONS(tix,2:end)));
  dx = 1e3 * min(bdx,tdx);

return;
