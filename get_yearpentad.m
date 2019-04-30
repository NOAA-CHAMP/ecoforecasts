function yrps = get_yearpentad(dts)
%function yrps = get_yearpentad(dts)
% Return the DATENUM of the start of the Pentad (5-day period), for each
% DATENUM in DTS, e.g., any date between Jan 1 and Jan 5 of year YR returns
% DATENUM(YR,1,1), and year-day 361 to year end are always in Pentad 73.
% 
% Last Saved Time-stamp: <Wed 2011-08-24 07:52:25  Lew.Gramer>

  [p,yr] = get_pentad(dts);

  % GET_PENTAD always returns Pentad between 1 and 73
  yrps = datenum(yr,1,1) + ((p-1)*5);

return;
