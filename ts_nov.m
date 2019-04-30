function idx = ts_nov(ts)
%function idx = ts_nov(ts)
%
% Return those indices of "time series" (struct with .DATE and .DATA fields)
% TS, which happen to full during the given month.

  [yr,mo] = datevec(ts.date);
  idx = find(mo == 11);

return;
