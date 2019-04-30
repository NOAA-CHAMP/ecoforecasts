function idx = ts_nonzero(ts)
%function idx = ts_nonzero(ts)
%
% Return all indices of "time series" TS (struct with .DATE and .DATA fields)
% for which .DATA happens to be different from zero by at least a small value

  idx = find(abs(ts.data) > 1e-9);

return;
