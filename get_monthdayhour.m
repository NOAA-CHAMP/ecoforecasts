function modyhrs = get_monthdayhour(dts)
%function modyhrs = get_monthdayhour(dts)
% Return DATENUM of calendar month-day-hour-min-sec, i.e.,
% DATENUM(0,m,d,H,M,S), for each DATENUM in DTS.
% 
% Last Saved Time-stamp: <Sun 2016-03-06 16:28:01 Eastern Standard Time gramer>

  [yrs,mos,dys,hrs,mns,scs] = datevec(dts);
  modyhrs = datenum(0,mos,dys,hrs,mns,scs);

return;
