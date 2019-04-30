function [s,yrs] = get_season(dts)
%function [s,yrs] = get_season(dts)
% Get "season number" (integer 1-4) for each datenum in DTS. NOTE WELL that
% in this context seasons are as defined by CLIPS CHAMP Ecoforecast code of
% Hendee (1998): "season 1" begins on Julian day 356 of one year and ends on
% jday 80 of the FOLLOWING YEAR; "2" is 81-172; "3" is 173-264; "4" 265-355.
% Optionally also returns year of each DATENUM in YRS.
% 
% Last Saved Time-stamp: <Fri 2011-04-22 07:23:46  Lew.Gramer>

  [jds,yrs] = get_jday(dts);

  s = repmat(0,size(dts));

  s(355 < jds | jds <=  80) = 1;
  s( 80 < jds & jds <= 172) = 4;
  s(172 < jds & jds <= 264) = 2;
  s(264 < jds & jds <= 355) = 3;

return;
