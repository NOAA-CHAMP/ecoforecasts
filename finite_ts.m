function newts = finite_ts(ts)
%function newts = finite_ts(ts)
%
% Return a new time series struct NEWTS containing only those indices of time
% series TS for which TS.data ISFINITE (v.)  This is just short-hand for the
% call SUBSET_TS(TS,@(x)(find(isfinite(x.data)))) (v.)
%
% Last Saved Time-stamp: <Mon 2012-03-26 11:50:28  Lew.Gramer>

  newts.date = ts.date(isfinite(ts.data));
  newts.data = ts.data(isfinite(ts.data));

return;
