function station = filter_gaps(station,basevar,var,maxgap,ringlen,tol,maskval)
%function station = filter_gaps(station,basevar,var,maxgap,ringlen,tol,maskval)
%
% Find all gaps in the time series STATION.(BASEVAR) that are longer than MAXGAP
% (default=3, i.e., 72 hrs), and remove any point from time series STATION.(VAR)
% that lies within TOL (default=(29/(24*60), i.e., 29 mins) of these gaps. If
% optional RINGLEN is non-zero, remove that many additional days worth of data
% from each gap: i.e., assume that each gap in BASEVAR causes a "ringing" or
% Gibbs oscillation in time series VAR, and remove the result of that ringing. 
% If optional MASKVAL is non-empty, fill gaps in STATION.(VAR).data with that
% value rather than removing them from the time series.
%
% Last Saved Time-stamp: <Thu 2011-05-19 15:17:09  Lew.Gramer>

  if ( ~exist('maxgap','var') || isempty(maxgap) )
    % Largest gap size to ignore: default is 3 days
    maxgap = 3.0;
  end;
  if ( maxgap <= 0 )
    error('MAXGAP argument should be greater than zero!');
  end;

  if ( ~exist('ringlen','var') || isempty(ringlen) )
    % Time AFTER each gap to remove from VAR: default is 0 days
    ringlen = 0.0;
  end;
  if ( ringlen < 0 )
    error('RINGLEN argument should be >= zero!');
  end;

  if ( ~exist('tol','var') || isempty(tol) )
    % Default tolerance for timestamp matching is 29 minutes
    tol = 29.0/(24.0*60.0);
  end;
  if ( tol <= 0 )
    error('Tolerance argument TOL should be greater than zero!');
  end;

  if ( ~exist('maskval','var') || isempty(maskval) )
    maskval = [];
  end;
  if ( numel(maskval) > 1 )
    error('MASKVAL argument must be empty or a scalar!');
  end;

  % First, remove any data from before and after BASEVAR was available
  gapix = find(station.(var).date < (min(station.(basevar).date) - tol));
  if ( isempty(maskval) )
    station.(var).date(gapix) = [];
    station.(var).data(gapix) = [];
  else
    station.(var).data(gapix) = maskval;
  end;

  gapix = find(station.(var).date > (max(station.(basevar).date) + tol));
  if ( isempty(maskval) )
    station.(var).date(gapix) = [];
    station.(var).data(gapix) = [];
  else
    station.(var).data(gapix) = maskval;
  end;

  % Remove any data in VAR that matches gaps in BASEVAR longer than MAXGAP days
  dtdiff = diff(station.(basevar).date);
  gapix = find( (dtdiff - maxgap) > eps );
  for ixix = length(gapix):-1:1
    ix = gapix(ixix);
    [ig, preix] = min(abs(station.(var).date - station.(basevar).date(ix) + tol));
    [ig, postix] = min(abs(station.(var).date - (station.(basevar).date(ix+1) + ringlen - tol)));
    % Ensure we only remove elements INSIDE the gap (plus RINGLEN minus TOL)...
    if ( (preix+1) <= (postix) )
      if ( isempty(maskval) )
        station.(var).date((preix+1):(postix)) = [];
        station.(var).data((preix+1):(postix)) = [];
      else
        station.(var).data((preix+1):(postix)) = maskval;
      end;
    end;
  end;

return;
