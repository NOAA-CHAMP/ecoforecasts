function [jhrs,yrs] = get_jhour_no_leap(dts)
%function [jhrs,yrs] = get_jhour_no_leap(dts)
% Get "Julian hour" (time difference vs. DATENUM(YR,1,1,0,0,0)) for each
% datenum in DTS. This version returns JHS-1 for all hours on Julian Day
% 366, thereby avoiding Leap Year complications. Optionally, also returns
% year of each DATENUM in YRS. 
% 
% Last Saved Time-stamp: <Sun 2016-05-01 20:01:18 Eastern Daylight Time gramer>

  [yrs,mos,dys,hrs,mns,scs] = datevec(dts);
  jhrs = datenum(0,mos,dys-1,hrs,0,0);
  % Lump day 366 into day 365
  jhrs(jhrs>=365) = jhrs(jhrs>=365) - 1;

return;
