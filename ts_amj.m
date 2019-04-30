function idx = ts_amj(ts)
%function idx = ts_amj(ts)
%
% Return those indices of "time series" (struct with .DATE and .DATA fields)
% TS, which happen to full during the months April, May, or June.

  [yr,mo] = datevec(ts.date);
  idx = find(mo == 4 | mo == 5 | mo == 6);

return;
