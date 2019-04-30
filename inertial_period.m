function [IP_day,IP_hr,IP_sec,f] = inertial_period(lat)
%function [IP_day,IP_hr,IP_sec,f] = inertial_period(lat)
% Convenience function: returns local inertial period at latitude ABS(LAT) in
% days (e.g., for DATENUM calculations), and also in hours and in seconds. If
% desired, the local Coriolis frequency /f/ can also be returned.

  OMEGA   = 7.292e-5;     %s-1   A.E.Gill p.597
  f = 2*OMEGA*sind(lat);  %SW_F: Phil Morgan 93-04-20  (morgan@ml.csiro.au)

  IP_sec = abs( 2*pi./f );
  IP_hr = IP_sec./3600;
  IP_day = IP_sec./3600./24;

return;
