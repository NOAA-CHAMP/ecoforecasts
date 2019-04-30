function [wks,yrs] = get_week(dts)
%function [wks,yrs] = get_week(dts)
% Get Week number (1-52) for each datenum in DTS. Optionally also returns year.
% 
% Last Saved Time-stamp: <Fri 2011-04-22 07:19:16  Lew.Gramer>

  [yrs,mos,dys] = datevec(dts);
  jds = floor(dts) - datenum(yrs,1,1) + 1;
  wks = ceil(jds / 7);
  wks(wks > 52) = 52;

return;
