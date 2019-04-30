function degs=kms2degs(kmxy,lat)
%function degs=kms2degs(kmxy,lat)
%
% Calculate coordinate range box [dLON,dLAT] from distance range box
% kmxy=[kmx,kmy]. CALLS: RECKON_WGS84 (v.)
%
  
  if (~exist('lon','var')||isempty(lon)); lon=-80; end;
  if (~exist('lat','var')||isempty(lat)); lat=+25; end;
  
  [ig,lonout] = reckon_wgs84(lat,lon,kmxy(1),90);
  degs(1)=abs(lonout-lon);
  
  [latout,ig] = reckon_wgs84(lat,lon,kmxy(2),0);
  degs(2)=abs(latout-lat);
  
return;
