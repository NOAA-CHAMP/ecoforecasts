function result = ts_d_dt(ts,dt,varargin)
%function result = ts_d_dt(ts,dt,varargin)
%
% Take finite difference of successive values in time series TS divided by
% time gaps; if optional DT is *non-empty*, use it as unit time difference
% (DEFAULT: 1.0 == 1 day). Calls FILTER_GAPS (with VARARGIN) to remove gaps
% in RESULT corresponding to gaps in TS.
%
% Last Saved Time-stamp: <Mon 2017-10-30 22:29:41 Eastern Daylight Time gramer>

  if ( ~exist('ts','var') || ~is_ts(ts) )
    error('Arg TS must be a time series struct!');
  end;
  if ( ~exist('dt','var') || isempty(dt) )
    dt = 1;
  end;

  x.ts = ts;
  x.result.date = ts.date(2:end);
  x.result.data = dt .* diff(ts.data) ./ diff(ts.date);

  x = filter_gaps(x,'ts','result',varargin{:});
  result = x.result;
  x=[]; clear x;

return;
