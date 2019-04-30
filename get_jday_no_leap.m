function [jds,yrs] = get_jday_no_leap(dts)
%function [jds,yrs] = get_jday_no_leap(dts)
% Get Julian Day (1-365, integer) for each datenum in DTS. This version
% returns 365 for Julian Day 366, thereby avoiding Leap Year complications.
% Optionally also returns year of each DATENUM in YRS.
% 
% Last Saved Time-stamp: <Fri 2011-04-22 07:21:58  Lew.Gramer>

  [yrs,mos,dys] = datevec(dts);
  jds = floor(dts) - datenum(yrs,1,1) + 1;
  jds(jds==366) = 365;

return;
