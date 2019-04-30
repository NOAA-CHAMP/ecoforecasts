function idx = ts_boreal_warm(ts)
%function idx = ts_boreal_warm(ts)
%
% Return those indices of "time series" (struct with .DATE and .DATA fields)
% TS, which happen to full during the Boreal warm part of the year - in this
% case, arbitrarily chosen as year-days 120 through 302, inclusive.

  [yr,ig] = datevec(ts.date);
  jd = ts.date - datenum(yr,1,1) + 1;
  idx = find(120 <= jd & jd <= 302);

return;
