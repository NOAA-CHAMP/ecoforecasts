function [yds,yrs] = get_yearday_no_leap(dts)
%function [yds,yrs] = get_yearday_no_leap(dts)
% 
% Return "year-days" YDS [0.0-364.9999], floating-point, and year number YRS
% (integer) for each date in the vector of DATENUM (qv.) DTS.  Like YEARDAY
% in AIR_SEA Toolbox, but simpler and more robust.  Distinguish from similar
% function GET_JDAY_NO_LEAP, which always returns an INTEGRAL 1-365. This
% version lumps leap-yeardays into day 364 to avoid Leap Year complications.
% 
% Last Saved Time-stamp: <Mon 2013-11-04 12:36:57 Eastern Standard Time gramer>

  [yrs,ig,ig] = datevec(dts);
  yds = dts - datenum(yrs,1,1);
  yds(yds>=365) = yds(yds>=365)-1;

return;
