function modys = get_monthday(dts)
%function modys = get_monthday(dts)
% Return DATENUM of calendar month and day of month, for each DATENUM in DTS.
% 
% Last Saved Time-stamp: <Sun 2016-03-06 16:25:50 Eastern Standard Time gramer>

  [yrs,mos,dys] = datevec(dts);
  modys = datenum(0,mos,dys);

return;
