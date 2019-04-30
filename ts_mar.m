function idx = ts_mar(ts)
%function idx = ts_mar(ts)
%
% Return those indices of "time series" (struct with .DATE and .DATA fields)
% TS, which happen to full during the given month.

  [yr,mo] = datevec(ts.date);
  idx = find(mo == 3);

return;
