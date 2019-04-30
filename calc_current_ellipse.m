function [lats,lons] = calc_current_ellipse(lon0,lat0,umj,umn,uinc,uph,ellips)
%function [lats,lons] = calc_current_ellipse(lon0,lat0,umj,umn,uinc,uph,ellips)
%
% Calculate set of LATS and LONS tracing an ellipse centered at location
% LON0,LAT0, with major axis UMJ (NOTE: in cm/s!), minor axis UMN (ditto),
% increment UINC (ditto), and phase UPH (degrees, not currently used). Useful
% for, e.g., outputs from function TMD_Toolbox ELLIPSE (v.). Optional ELLIPS
% is the ELLIPSOID arg passed to ELLIPSE (DEFAULT: WGS84_ELLIPSOID()).
%
% CALLS: AXES2ECC, ELLIPSE1 (Map Toolbox), WGS84_ELLIPSOID (Ecoforecasts).
% 
% Last Saved Time-stamp: <Fri 2017-03-10 15:55:57 Eastern Standard Time gramer>

  if ( ~exist('ellips','var') || isempty(ellips) )
    ellips = wgs84_ellipsoid();
  end;

  % % [LAT,LON] = ellipse1(LAT0,LON0,ELLIPSE) computes ellipse(s) with
  % % center(s) at LAT0, LON0.  ELLIPSE must have the form [semimajor axis,
  % % eccentricity].  LAT0 and LON0 can be scalar or column vectors. The
  % % input and output latitudes and longitudes are in units of degrees.
  % % ELLIPSE must have the same number of rows as the input LAT0 and LON0.
  % % The semimajor axis (ELLIPSE(1)) is in degrees of arc length on a
  % % sphere. All ellipses are oriented so that their major axis runs
  % % north-south.
  % %
  % % [LAT,LON] = ellipse1(LAT0,LON0,ELLIPSE,OFFSET) computes the ellipses
  % % where the major axis is rotated from due north by an azimuth OFFSET.
  % % The offset angle is measure clockwise from due north.
  % [lats,lons] = ellipse1(lat0,lon0,[umj/111,axes2ecc(umj,abs(umn))],uinc);

  % [lat,LON] = ellipse1(LAT0,LON0,ELLIPSE,OFFSET,AZ,ELLIPSOID) computes
  % the ellipse(s) on the reference ellipsoid defined by ELLIPSOID.
  % ELLIPSOID is a reference ellipsoid (oblate spheroid) object, a
  % reference sphere object, or a vector of the form [semimajor_axis,
  % eccentricity].  The semimajor axis of the ellipse must be in the
  % same units as the ellipsoid's semimajor axis, unless ELLIPSOID is [].
  % If ELLIPSOID is [], then the semimajor axis of the ellipse is
  % interpreted as an angle and the ellipse is computed on a sphere,
  % as in the preceding syntax.
  [lats,lons] = ellipse1(lat0,lon0, [umj/111,axes2ecc(umj,abs(umn))], uinc, [], ellips);

return;
