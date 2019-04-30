function idx = ts_isfinite(ts)
%function idx = ts_isfinite(ts)
%
% Return those INDICES (not a new TS, like FINITE_TS) of time series TS for which TS.data ISFINITE (v.)
%
% Last Saved Time-stamp: <Fri 2012-06-01 13:01:12  Lew.Gramer>

  idx = find(isfinite(ts.data));

return;
