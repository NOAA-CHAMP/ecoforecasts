function [s,yrs] = get_season(dts)
%function [s,yrs] = get_season(dts)
% Get "season number" for each datenum in DTS. "Season number" is numeric
% equivalent of a triad, one of JFM, AMJ, JAS, OND (integers 1-4, resp.)
% Optionally also returns year for each DATENUM in YRS.
% 
% Last Saved Time-stamp: <Fri 2011-04-22 07:18:17  Lew.Gramer>

  [yrs,mos] = datevec(dts);

  s = repmat(0,size(dts));
  s(mos == 1  | mos ==  2 | mos ==  3) = 1;
  s(mos == 4  | mos ==  5 | mos ==  6) = 2;
  s(mos == 7  | mos ==  8 | mos ==  9) = 3;
  s(mos == 10 | mos == 11 | mos == 12) = 4;

return;
