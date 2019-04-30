function [jds,yrs] = get_jday(dts)
%function [jds,yrs] = get_jday(dts)
% Get Julian Day (1-366, integer) for each DATENUM in DTS. Optionally return year too.
% 
% Last Saved Time-stamp: <Fri 2011-04-22 07:18:33  Lew.Gramer>

  [yrs,mos,dys] = datevec(dts);
  jds = floor(dts) - datenum(yrs,1,1) + 1;

return;
