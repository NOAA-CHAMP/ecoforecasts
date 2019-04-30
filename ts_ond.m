function idx = ts_ond(ts)
%function idx = ts_ond(ts)
%
% Return those indices of "time series" (struct with .DATE and .DATA fields)
% TS, which happen to full during the months October, November, or December.

  [yr,mo] = datevec(ts.date);
  idx = find(mo == 10 | mo == 11 | mo == 12);

return;
