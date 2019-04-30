function newts = nearest_ts(ts1,idx,method)
%function newts = nearest_ts(ts1,idx_or_fcn_or_ts2,method)
%
% Return time series NEWTS interpolated from the values of time series TS1
% that are selected by IDX. IDX may be an index vector, a FUNCTION_HANDLE or
% CHAR function name that accepts and returns vectors of DATENUMs, or another
% time series STRUCT.  If IDX is not specified, NEWTS is all the values
% nearest to each unique whole hour in TS, UNIQUE(ROUND(TS.date*24)/24).
%
% Calls INTERP1 (v.) on TS1.date,data with METHOD (DEFAULT: 'nearest').
%
% Last Saved Time-stamp: <Tue 2016-08-30 16:30:11 Eastern Daylight Time lew.gramer>

  if ( ~exist('ts1','var') || ~is_ts(ts1) )
    error('First agument must be a time series STRUCT');
  end;

  if ( ~exist('idx','var') || isempty(idx) )
    newts.date = unique(round(ts1.date*24)/24);
  elseif ( isa(idx,'function_handle') || ischar(idx) )
    newts.date = feval(idx,ts1.date);
  elseif ( is_ts(idx) )
    newts.date = idx.date;
  else
    newts.date = ts1.date(idx);
  end;

  if ( ~exist('method','var') || isempty(method) )
    method = 'nearest';
  end;

  newts.data = interp1(ts1.date,ts1.data,newts.date,method);

  % Convenience for 2D time-series fields, e.g., ADCP profiles
  if ( isfield(ts1,'prof') )
    newts.prof = interp1(ts1.date,ts1.prof,newts.date,method);
  end;

return;
