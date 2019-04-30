function ndys = get_dayspermonth(dts)
%function ndys = get_dayspermonth(dts)
% For each DATENUM in DTS, return the number of days in the month
% corresponding to that date (including Feb 29 for leap years).
% 
% Last Saved Time-stamp: <Fri 2011-04-22 08:12:30  Lew.Gramer>

  [yrs,mos,ig] = datevec(dts);
  ndys = eomday(yrs,mos);

return;
