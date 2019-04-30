function idx = ts_dec(ts)
%function idx = ts_dec(ts)
%
% Return those indices of "time series" (struct with .DATE and .DATA fields)
% TS, which happen to full during the given month.

  [yr,mo] = datevec(ts.date);
  idx = find(mo == 12);

return;
