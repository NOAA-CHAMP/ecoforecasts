function idx = ts_jfm(ts)
%function idx = ts_jfm(ts)
%
% Return those indices of "time series" (struct with .DATE and .DATA fields)
% TS, which happen to full during the months January, February, or March.

  [yr,mo] = datevec(ts.date);
  idx = find(mo == 1 | mo == 2 | mo == 3);

return;
