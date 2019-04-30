function yfs = get_yearfrac(dts)
%function yfs = get_yearfrac(dts)
%
% Return floating-point "factional year" YFS, e.g., 2014.5 for June 1st, 2014
% for each date in the vector of DATENUM (qv.) DTS.
%
% Last Saved Time-stamp: <Tue 2017-12-05 12:12:40 Eastern Standard Time gramer>

  yfs = get_year(dts) + (get_jday(dts)/366);

return;
