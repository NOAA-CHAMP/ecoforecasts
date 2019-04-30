function [ix1,ix2] = intersect_dates_slow(dts1, dts2, tol)
%function [ix1,ix2] = intersect_dates_slow(dts1, dts2, tol)
%
% Return indices of all the elements in two series of timestamps (datenums)
% 'dts1' and 'dts2', that match: a "matching" timestamp is defined as one
% where the difference between the elements of the two series is less than
% some tolerance level 'tol'. DEFAULT 'tol' is 29 minutes, i.e., 29/(24*60).
%
% Last Saved Time-stamp: <Mon 2009-08-17 09:03:04 Eastern Daylight Time Lew.Gramer>

  ix1 = [];
  ix2 = [];

  if ( any(diff(dts1) <= 0) || any(diff(dts2) <= 0) )
    error('Both date series must be monotonically increasing!');
  end;
  if ( ~exist('tol','var') || isempty(tol) )
    % Default tolerance for timestamp matching is 29 minutes
    tol = 29.0/(24.0*60.0);
  end;
  if ( tol <= 0 )
    error('Tolerance argument should be greater than zero!');
  end;

  switched = false;
  if ( numel(dts2) > numel(dts1) )
    switched = true;
    tmp = dts1;
    dts1 = dts2;
    dts2 = tmp;
    tmp = []; clear tmp;
  end;    

  preix1 = find( ((dts2(1)-tol) <= dts1) & (dts1 <= (dts2(end)+tol)) );
  preix2 = find( ((dts1(1)-tol) <= dts2) & (dts2 <= (dts1(end)+tol)) );
  if ( isempty(preix1) || isempty(preix2) )
    return;
  end;

  % SLOOOWOWOWOW! Need a clever MATLAB-ish way to do this...
  for ixi = 1:length(preix1)
    ix = preix1(ixi);
    goodix = find(abs(dts1(ix)-dts2) <= tol);
    if ( ~isempty(goodix) )
      ix1 = [ix1 ix];
      ix2 = [ix2 goodix(1)];
    end;
  end;

  ix1 = unique(ix1);
  ix2 = unique(ix2);

  if ( switched )
    tmp = ix1;
    ix1 = ix2;
    ix2 = tmp;
    tmp = []; clear tmp;
  end;    

return;
