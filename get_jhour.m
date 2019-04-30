function [jhrs,yrs] = get_jhour(dts)
%function [jhrs,yrs] = get_jhour(dts)
% Get "Julian hour" (time difference vs. DATENUM(YR,1,1,0,0,0)) for each
% datenum in DTS. Optionally also returns year of each DATENUM in YRS.
% 
% Last Saved Time-stamp: <Sun 2016-05-01 20:00:39 Eastern Daylight Time gramer>

  [yrs,mos,dys,hrs,mns,scs] = datevec(dts);
  jhrs = datenum(0,mos,dys-1,hrs,0,0);

return;
