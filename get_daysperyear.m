function ndys = get_daysperyear(dts)
%function ndys = get_daysperyear(dts)
% For each DATENUM in DTS, return the number of days in the year
% corresponding to that date (i.e., returns 366 for leap years).
% 
% Last Saved Time-stamp: <Fri 2011-04-22 08:26:24  Lew.Gramer>

  [yrs,ig,ig] = datevec(dts);
  ndys = datenum(yrs,12,31) - datenum(yrs,1,1) + 1;

return;
