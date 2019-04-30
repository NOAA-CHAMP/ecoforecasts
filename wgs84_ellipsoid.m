function result = wgs84_ellipsoid
% World Geodetic System 84 (WGS84) reference ellipsoid for Earth, in [km].
%
% From http://en.wikipedia.org/wiki/WGS84:
%   As of the latest revision, the WGS 84 datum surface is a pole-flattened
%   (oblate) spheroid, with major (transverse) radius a = 6,378,137 m at the
%   equator, minor (conjugate) radius b = 6,356,752.314 245 m at the poles.
% [Approximated by an ellipsoid with eccentricity 1/298.257223563 ~ 0.335%.]
%
% Last Saved Time-stamp: <Fri 2010-10-08 14:01:31  Lew.Gramer>

  % OOPS!! result = [6356.752 (1/298.25722356)];

  result = [6378.137 (1/298.25722356)];

return;
