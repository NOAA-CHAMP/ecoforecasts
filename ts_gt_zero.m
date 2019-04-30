function idx = ts_gt_zero(ts)
%function idx = ts_gt_zero(ts)
%
% Return those indices of "time series" (struct with .DATE and .DATA fields)
% TS, for which .DATA happens to be Greater Than ZERO.

  idx = find(ts.data > 0);

return;
