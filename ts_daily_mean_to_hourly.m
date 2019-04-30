function hts = ts_daily_mean_to_hourly(ts)
%function hts = ts_daily_mean_to_hourly(ts)
%
% Expand a daily mean time series TS into an hourly time series HTS with the
% same daily mean. Note that the resulting time series is NOT smooth, but it
% will contain a 24-hour gap where ever a day is missing from TS.

  hts.date = [ts.date(1):(1/24):(ts.date(end)+(23/24)+eps)]';
  ts24 = repmat(ts.data,[1,24])';
  hts.data = ts24(:);

return;
