function [lons,lats] = ellipse_wgs84(lon,lat,xrad,yrad,npts)
%function [lons,lats] = ellipse_wgs84(lon,lat,xrad,yrad,npts)
%
% Return NPTS points (DEFAULT: 100) with coordinates [LONS,LATS], along the
% ellipse centered at LON,LAT, with semi-x axis XRAD and semi-y axis YRAD.
% N.B.: XRAD,YRAD are assumed to be in [km], but LON,LAT in decimal degrees.
%
% Calls ELLIPSE1 (v.) using the World Geodetic System 84 ellipsoid for Earth,
% as returned by WGS84_ELLIPSOID (v.).  Calls AXES2ECC for geographic ellipse.
%
% Last Saved Time-stamp: <Fri 2010-10-08 13:53:42  Lew.Gramer>

  if ( ~exist('npts','var') || isempty(npts) )
    npts = 100;
  end;

  smaecc(1) = xrad;
  smaecc(2) = axes2ecc(xrad,yrad);

  [lats,lons] = ellipse1(lat,lon,smaecc,90,[],wgs84_ellipsoid(),[],100);

return;
