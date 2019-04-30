function [lonix, latix] = gridnbhd_km(lons,lats,lon,lat,radiuskm)
%function [lonix, latix] = gridnbhd_km(lons,lats,lon,lat,radiuskm)
%
% Find all grid points in the meshgrid formed by LONS and LATS, that lie
% within RADIUSKM km (DEFAULT: 9.75km) of the point [lat,lon]. Uses WGS-84
% radius and ellipticity, as inputs to MAP Toolbox function DISTANCE (qv).
% If RADIUSKM==0, then just return the gridpoint nearest to [lon,lat].
%
% Last Saved Time-stamp: <Sat 2010-01-30 11:35:37 Eastern Standard Time gramer>

  if (~exist('radiuskm','var') || isempty(radiuskm)); radiuskm = 9.75; end;

  % If given a matrix (ANY matrix), assume input is already gridded!
  if ( size(lons,1) > 1 && size(lons,2) > 1 )
    LONS = lons;
    LATS = lats;
  else
    [LONS,LATS] = meshgrid(lons,lats);
  end;

  LON = repmat(lon, size(LONS));
  LAT = repmat(lat, size(LATS));

  % Assume a spherical earth... Hey, what is this, rocket science??
  %ds = distance(LATS, LONS, LAT, LON, [6371 0]);

  % *OR* use WGS-84 ellipticity and polar radius - why not?
  ds = distance(LATS, LONS, LAT, LON, [6356.752 (1/298.25722356)]);

  if ( radiuskm == 0 )
    % Just return index of one (nearest) gridpoint
    [ig, dix] = min(ds(:));
    [latix, lonix] = ind2sub(size(ds), dix(1));
  else
    % Return indices of all gridpoints within RADIUSKM
    [latix, lonix] = find(ds <= radiuskm);
  end;

return;
