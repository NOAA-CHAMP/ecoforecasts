function idx = ts_boreal_cool(ts)
%function idx = ts_boreal_cool(ts)
%
% Return those indices of "time series" (struct with .DATE and .DATA fields)
% TS, which happen to full during the Boreal cool part of the year - in this
% case, arbitrarily chosen as year-days 1 to 119, and 303 to 31 Dec, incl.

  [yr,ig] = datevec(ts.date);
  jd = ts.date - datenum(yr,1,1) + 1;
  idx = find(120 > jd | jd > 302);

return;
