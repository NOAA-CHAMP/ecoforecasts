function [anomts,clim,jds] = anomalize_daily_mean_ts(ts)
%function [anomts,clim,jds] = anomalize_daily_mean_ts(ts)
%
% FOR A MORE GENERAL FUNCTION, see: ANOMALIZE_TS.m.
%
% Calculate a year-day climatology CLIM for time series struct TS. Subtract
% the year-day mean from each value of TS to form new time series ANOMTS.
%
% Last Saved Time-stamp: <Wed 2016-08-31 11:26:18 Eastern Daylight Time gramer>

  [clim,jds] = grp_ts(ts.data,ts.data);
  anomts = ts;
  for jdix=1:numel(jds)
    ix = find(get_jday(anomts.date)==jds(jdix));
    anomts.data(ix) = anomts.data(ix) - clim(jdix);
  end;

return;
