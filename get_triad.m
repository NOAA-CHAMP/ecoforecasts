function [t,yrs] = get_triad(dts)
%function [t,yrs] = get_triad(dts)
% Get "Triad" number T (1-122) for each datenum in DTS: a Triad is a 3-day
% binning period used for example by some reanalysis and satellite products.
% Note that in non-leap years, Triad 122 only has TWO days in it. Optionally
% also returns the year for each DATENUM in YRS.
% 
% Last Saved Time-stamp: <Wed 2011-08-24 07:52:02  Lew.Gramer>

  [yrs,mos,dys] = datevec(dts);
  jds = floor(dts) - datenum(yrs,1,1) + 1;
  t = ceil(jds / 3);

return;
