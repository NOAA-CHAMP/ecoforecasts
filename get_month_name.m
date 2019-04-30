function monms = get_month_name(dts)
%function monms = get_month_name(dts)
% 
% Get short name of month for current Locale (e.g., 'Jan' through 'Dec' for
% US-English), for each DATENUM in DTS. Returns result of DATESTR(DTS,'mmm').
% 
% Last Saved Time-stamp: <Sat 2017-07-01 13:11:21 Eastern Daylight Time gramer>

  monms = datestr(dts,'mmm');

return;
