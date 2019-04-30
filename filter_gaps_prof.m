function stn = filter_gaps_prof(stn,basevar,var,maxgap,ringlen,tol,maskval)
%function stn = filter_gaps_prof(stn,basevar,var,maxgap,ringlen,tol,maskval)
%
% Find all gaps in the profile/time series STN.(BASEVAR) longer than MAXGAP
% (DEFAULT: 3 =72 h), and remove any point from profile/time series STN.(VAR)
% lying within TOL (DEFAULT: 29/(24*60) =29 min) of each gap. If optional
% RINGLEN<>0, remove that many additional days worth of data from each gap:
% i.e., assume that each gap in BASEVAR causes Gibbs "ringing" in time series
% VAR, and remove the result of that ringing. If optional MASKVAL non-empty,
% fill gaps in STN.(VAR) with that value instead of removing them from the
% profile/time series. STN.(VAR) must be a struct with fields 'date', 'prof',
% and optionally with field 'data'. STN.(VAR).prof may be either shape.
%
% Last Saved Time-stamp: <Mon 2015-04-13 13:59:18 Eastern Daylight Time gramer>

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

  prof = stn.(var).prof;
  % Ensure we respect the shape of the profile array
  if ( size(prof,2) > size(prof,1) )
    swapped = true;
    prof = prof';
  else
    swapped = false;
  end;

  % First, remove any data from before and after BASEVAR was available
  gapix = find(stn.(var).date < (min(stn.(basevar).date) - tol));
  if ( isempty(maskval) )
    stn.(var).date(gapix) = [];
    prof(gapix,:) = [];
    if ( isfield(stn.(var),'data') )
      stn.(var).data(gapix) = [];
    end;
  else
    prof(gapix,:) = maskval;
    if ( isfield(stn.(var),'data') )
      stn.(var).data(gapix) = maskval;
    end;
  end;

  gapix = find(stn.(var).date > (max(stn.(basevar).date) + tol));
  if ( isempty(maskval) )
    stn.(var).date(gapix) = [];
    prof(gapix,:) = [];
    if ( isfield(stn.(var),'data') )
      stn.(var).data(gapix) = [];
    end;
  else
    prof(gapix,:) = maskval;
    if ( isfield(stn.(var),'data') )
      stn.(var).data(gapix) = maskval;
    end;
  end;

  % Remove any data in VAR that matches gaps in BASEVAR longer than MAXGAP days
  dtdiff = diff(stn.(basevar).date);
  gapix = find( (dtdiff - maxgap) > eps );
  for ixix = length(gapix):-1:1
    ix = gapix(ixix);
    [ig, preix] = min(abs(stn.(var).date - stn.(basevar).date(ix) + tol));
    [ig, postix] = min(abs(stn.(var).date - (stn.(basevar).date(ix+1) + ringlen - tol)));
    % Ensure we only remove elements INSIDE the gap (plus RINGLEN minus TOL)...
    if ( (preix+1) <= (postix) )
      stn.(var).date((preix+1):(postix)) = [];
      if ( isempty(maskval) )
        prof((preix+1):(postix),:) = [];
        if ( isfield(stn.(var),'data') )
          stn.(var).data((preix+1):(postix)) = [];
        end;
      else
        prof((preix+1):(postix),:) = maskval;
        if ( isfield(stn.(var),'data') )
          stn.(var).data((preix+1):(postix)) = maskval;
        end;
      end;
    end;
  end;

  if ( swapped )
    prof = prof';
  end;
  stn.(var).prof = prof;

return;
