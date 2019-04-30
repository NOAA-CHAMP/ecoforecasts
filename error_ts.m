function result = error_ts(ts1,ts2,method,ix1,ix2)
%function result = error_ts(ts1,ts2,method,ix1,ix2)
%
% Calculate residuals between the two time series structs TS1 and TS2. METHOD
% may be any of 'rss','rmse','nrmse','cvrmse','lad','mem' (DEFAULT: 'rmse').
%
% Last Saved Time-stamp: <Tue 2011-02-08 07:16:58  lew.gramer>

  if ( ~isfield(ts1,'date') || ~isfield(ts2,'date') || ~isfield(ts1,'data') || ~isfield(ts2,'data') )
    error('First and second args TS1 and TS2 must be time series structs!');
  end;
  if ( ~exist('ix1','var') || isempty(ix1) )
    ix1 = 1:length(ts1.date);
  end;
  if ( ~exist('ix2','var') || isempty(ix2) )
    ix2 = 1:length(ts2.date);
  end;
  if ( ~exist('method','var') || isempty(method) )
    method = 'rmse';
  end;

  if ( isa(ix1,'function_handle') )
    ix1 = feval(ix1,ts1);
  end;
  if ( isa(ix2,'function_handle') )
    ix2 = feval(ix2,ts2);
  end;

  ts1.date = ts1.date(ix1);
  ts1.data = ts1.data(ix1);

  ts2.date = ts2.date(ix2);
  ts2.data = ts2.data(ix2);

  [ix1,ix2] = intersect_dates(ts1.date,ts2.date);
  n = length(ix1);
  if (n == 0)
    error('No coincident dates in time series!');
  end;

  resid = abs( ts1.data(ix1) - ts2.data(ix2) );
  ssr = sum(resid.^2);
  rmse = sqrt(ssr/n);

  switch ( lower(method) ),
   case {'resid','residual'},
    % Undocumented, used mainly for debugging
    result = resid;
   case {'ssr','rss'},
    % Sum of Squared Residuals
    result = ssr;
   case {'l2','l-2','rms','rmsd','rmse'},
    % Root Mean Squared Error
    result = rmse;
   case {'nrmsd','nrmse'},
    % Normalized RMSE
    result = rmse / (max(resid) - min(resid));
   case {'cvrmsd','cvrmse'},
    % Covariance(RMSE)
    result = rmse / mean(resid);
   case {'l1','l-1','lad','lav','lae'},
    % Least Absolute Difference error
    result = sum( resid );
   case {'mem'},
    % Maximum Entropy Method error
    result = sum( resid.*log(resid) );
   otherwise,
    error('Unknown error evaluation method "%s"',method);
  end;
  resid = []; ssr = []; rmse = [];

return;
