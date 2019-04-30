function ye = get_yearend(dts)
%function ye = get_yearend(dts)
% Return DATENUM of final second of the calendar year, for each DATENUM in DTS.
% 
% Last Saved Time-stamp: <Fri 2018-04-06 14:47:33 Eastern Daylight Time gramer>

  [yrs,mos,dys] = datevec(dts);
  ye = datenum(yrs,12,31,23,59,59);

return;
