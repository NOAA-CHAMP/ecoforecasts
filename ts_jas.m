function idx = ts_jas(ts)
%function idx = ts_jas(ts)
%
% Return those indices of "time series" (struct with .DATE and .DATA fields)
% TS, which happen to full during the months July, August, or September.

  [yr,mo] = datevec(ts.date);
  idx = find(mo == 7 | mo == 8 | mo == 9);

return;
