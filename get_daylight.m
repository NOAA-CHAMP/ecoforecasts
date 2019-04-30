function day = get_daylight(dts,lat,lon,minz)
%function day = get_daylight(dts,lat,lon,minz)
% Returns a vector same size as DTS, with 1 for timestamps that occur during
% local daylight, and 0 otherwise. DTS assumed to be in *GMT*. Uses negative
% LON for West longitudes: LON > 180 are converted. Optional MINZ specifies
% minimum solar altitude (degrees) for timestamp to be considered daylit
% (DEFAULT: 0). CALLS: SORADNA1 (from AIR_SEA Toolbox).
%
% Last Saved Time-stamp: <Fri 2011-04-01 07:35:00  Lew.Gramer>

  if (~exist('minz','var')||isempty(minz));
    minz=0;
  end;
  lon(lon>180) = 360 - lon(lon>180);
  [yrs,ig,ig] = datevec(dts);
  yds = dts - datenum(yrs,1,1);
  [z,s] = soradna1(yds,yrs,-lon,lat);
  day = repmat(0,size(dts));
  day(z > minz) = 1;

return;
