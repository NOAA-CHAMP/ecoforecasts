function km = lon2km(lon1,lon2,lat1,lat2)
%function km = lon2km([lon1[,lon2[,lat1[,lat2]]]])
  if (~exist('lon1','var')||isempty(lon1)); lon1=-80; end;
  if (~exist('lon2','var')||isempty(lon2)); lon2=lon1+1; end;
  if (~exist('lat1','var')||isempty(lat1)); lat1=+25; end;
  if (~exist('lat2','var')||isempty(lat2)); lat2=lat1; end;
  if (~exist('az','var')||isempty(az)); az=90; end;
  
  km = distance_wgs84(lat1,lon1,lat2,lon2);
return;
