function [u, v] = spddir_to_uv(wspd,wdir)
%function [u, v] = spddir_to_uv(wspd,wdir)
% Convert time series of wind speed WSPD (any units), and wind direction in
% degrees True WDIR, into u and v vector components (in SAME UNITS as WSPD).
%
% NOTE: This version of 'spddir-to-uv' is coded for WINDS: "dir" here means
% "source direction", NOT "target direction" as it would for ocean currents.
% Thus for EXAMPLE:
%  >> [u,v] = spddir_to_uv(10,0)
% ... will return u == 0 and v == -10, *not* v == 10.
%
% Last Saved Time-stamp: <Mon 2011-05-16 15:59:24 Eastern Daylight Time gramer>

  u = roundn( wspd .* (-sind(wdir)), -8);
  v = roundn( wspd .* (-cosd(wdir)), -8);

return;
