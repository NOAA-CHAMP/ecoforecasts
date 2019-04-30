function [p,yrs] = get_pentad(dts)
%function [p,yrs] = get_pentad(dts)
% Get "Pentad" number P (1-73) for each datenum in DTS: Pentad is a 5-day
% binning period used for example by Pathfinder v5 satellite SST product.
% Note that in leap years, year-day 366 is forced into Pentad 73. Optionally
% also returns the year for each DATENUM in YRS.
% 
% Last Saved Time-stamp: <Fri 2011-04-22 07:19:02  Lew.Gramer>

  [yrs,mos,dys] = datevec(dts);
  jds = floor(dts) - datenum(yrs,1,1) + 1;
  p = ceil(jds / 5);
  p(p > 73) = 73;

return;
