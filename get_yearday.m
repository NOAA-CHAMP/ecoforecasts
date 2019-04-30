function [yds,yrs] = get_yearday(dts)
%function [yds,yrs] = get_yearday(dts)
%
% Return "year-days" YDS [0.0-365.9999], floating-point, and year number YRS
% (integer) for each date in the vector of DATENUM (qv.) DTS.  Similar to
% YEARDAY in the AIR_SEA Toolbox, but simpler and more robust.  Distinguish
% from similar function GET_JDAY, which always returns an INTEGRAL 1-366.
%
% Last Saved Time-stamp: <Mon 2013-11-04 12:34:08 Eastern Standard Time gramer>

  [yrs,ig,ig] = datevec(dts);
  yds = dts - datenum(yrs,1,1);

return;
