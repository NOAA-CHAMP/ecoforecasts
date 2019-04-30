function newts = good_grp_ts(ts,grpfun,pct,goodfun)
%function newts = good_grp_ts(ts,grpfun,pct,goodfun)
%
% Subset time series TS (STRUCT with .date and .data) to periods with at
% least PCT good values (0<PCT<=1), where "good" is defined as all values
% whose indices are returned by GOODFUN(TS) (DEFAULT: @TS_ISFINITE). Creates
% time series of only good values NEWTS = SUBSET_TS(TS,GOODFUN), then calls
% FIND_PCT_GOOD_DATES(NEWTS.date,GRPFUN,PCT) (v.) to find subset. 
%
% NOTE: GRPFUN must be a function-handle that accepts and returns DATENUM.
%
% Last Saved Time-stamp: <Thu 2017-05-25 13:44:38 Eastern Daylight Time lew.gramer>
%error('UGH! This turns out to be very complex to code in a generic way... For now, see FIND_PCT_GOOD_PERIODS.m instead!');

error('UGH! This turns out to be very complex to code in a generic way... For now, see FIND_PCT_GOOD_PERIODS.m instead!');

  if ( ~exist('ts','var') || ~is_valid_ts(ts) )
    error('First arg must be a valid time series STRUCT');
  end;
  if ( ~exist('goodfun','var') || isempty(goodfun) )
    goodfun = @ts_isfinite;
  end;

  newts = subset_ts(ts,goodfun);
  idx = find_pct_good_dates(newts.date,grpfun,pct);
  newts.date = newts.date(idx);
  newts.data = newts.data(idx);

return;
