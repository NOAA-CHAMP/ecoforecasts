function [u, v] = spddir_to_uv_curr(cspd,cdir)
%function [u, v] = spddir_to_uv_curr(cspd,cdir)
% Convert time series of current speed CSPD (any units), and direction CDIR
% (degrees True), into u and v vector components (in SAME UNITS as CSPD).
%
% NOTE: This version of 'spddir-to-uv' is coded for CURRENTS: "dir" here
% means "target direction" as it would for ocean currents.
% Thus for EXAMPLE:
%  >> [u,v] = spddir_to_uv_curr(10,0)
% ... will return u == 0 and v == 10, *not* v == -10 as it would for winds.
%
% Last Saved Time-stamp: <Mon 2011-05-16 15:58:48 Eastern Daylight Time gramer>

  cdir = cdir - 180;
  negix = find(cdir < 0);
  cdir(negix) = 360 + cdir(negix);

  u = roundn( cspd .* (-sind(cdir)), -8);
  v = roundn( cspd .* (-cosd(cdir)), -8);

return;
