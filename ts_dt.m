function result = ts_dt(ts,dt,varargin)
%function result = ts_dt(ts,dt,varargin)
%
% Take the finite difference of successive values in time series TS; if
% optional DT is *non-empty*, only every DT index of TS. Call FILTER_GAPS
% (with VARARGIN) to remove gaps in RESULT corresponding to gaps in TS.
%
% Last Saved Time-stamp: <Tue 2013-02-05 14:47:10 Eastern Standard Time gramer>

  if ( ~is_ts(ts) )
    error('Arg TS must be a time series struct!');
  end;
  if ( ~exist('dt','var') || isempty(dt) )
    dt = 1;
  end;

  x.ts = ts;
  x.result.date = ts.date(2:dt:end);
  x.result.data = diff(ts.data(1:dt:end));

  x = filter_gaps(x,'ts','result',varargin{:});
  result = x.result;
  x=[]; clear x;

return;
