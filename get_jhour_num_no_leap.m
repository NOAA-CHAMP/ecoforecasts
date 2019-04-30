function [jhs,yrs] = get_jhour_no_leap(dts)
%function [jhs,yrs] = get_jhour_no_leap(dts)
% Get "Julian hour NUMBER" ([1-366]*24, integer) for each datenum in DTS.
% This version returns JHS-24 for all hours on Julian Day 366, thereby
% avoiding Leap Year complications, i.e., return values are always JHS>=1,
% JHS<=8760. Optionally also returns year of each DATENUM in YRS.
% 
% Last Saved Time-stamp: <Sun 2016-05-01 19:51:09 Eastern Daylight Time gramer>

  [yrs,mos,dys] = datevec(dts);
  jhs = floor((dts-datenum(yrs,1,1))*24)+1;
  jhs(jhs>8760) = jhs(jhs>8760)-24;

return;
