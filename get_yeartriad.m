function yrts = get_yeartriad(dts)
%function yrts = get_yeartriad(dts)
% Return the DATENUM of the start of the Triad (3-day period), for each
% DATENUM in DTS, e.g., any date between Jan 1 and Jan 3 of year YR returns
% DATENUM(YR,1,1), and year-day 364 to year-end are always in Triad 122.
% 
% Last Saved Time-stamp: <Wed 2011-08-24 07:54:28  Lew.Gramer>

  [t,yr] = get_triad(dts);

  % GET_TRIAD always returns Triad between 1 and 122
  yrts = datenum(yr,1,1) + ((t-1)*3);

return;
