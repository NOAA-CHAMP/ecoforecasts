function result = ts_fun(ts,fn,ix,varargin)
%function result = ts_fun(ts,fn,ix,varargin)
%
% Apply FN (a function-handle, name, or unary punctuation operator) to all
% values of the time series TS. If optional IX is *non-empty*, only use those
% points from TS. IX can also be a time series-indexing function_handle,
% e.g., TS_JFM. Any additional args after IX are simply passed on as
% arguments to FN, e.g., FEVAL(FN,TS,VARARGIN{:}).
%
% Last Saved Time-stamp: <Thu 2012-02-09 13:05:41  Lew.Gramer>

  if ( ~is_ts(ts) )
    error('Arg TS must be a time series struct!');
  end;

  if ( ~exist('ix','var') || isempty(ix) )
    ix = 1:length(ts.date);
  end;

  if ( isa(ix,'function_handle') )
    ix = feval(ix,ts);
  end;

  ts.date = ts.date(ix);
  ts.data = ts.data(ix);

  result = ts;
  result.data(:) = feval(fn,ts.data,varargin{:});

return;
