function [jhs,yrs] = get_jhour(dts)
%function [jhs,yrs] = get_jhour(dts)
% Get "Julian hour NUMBER" ([1-366]*24, integer) for each datenum in DTS.
% Optionally also returns year of each DATENUM in YRS.
% 
% Last Saved Time-stamp: <Sun 2016-05-01 19:50:42 Eastern Daylight Time gramer>

  [yrs,mos,dys] = datevec(dts);
  jhs = floor((dts-datenum(yrs,1,1))*24)+1;

return;
