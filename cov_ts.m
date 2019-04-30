function [c,r,p,fix1,fix2] = cov_ts(t1,t2,fix1,fix2)
%function [c,r,p,fix1,fix2] = cov_ts(t1,t2,fix1,fix2)
%
% Return covariance C, correlation coefficient R, and P estimate for the null
% hypothesis (with 95% confidence) of no correlation, among coincident values
% of time series  T1 and T2: T1,T2 may either be structs having fields .date,
% .data, or numeric vectors (in which case they must be identically sized).
% If either of the optional args FIX1 or FIX2 are non-empty, statistics are
% limited to those subindices within T1 and T2, resp. (E.g., if both optional
% args are given, covariance/correlation of coincident values from T1.data(FIX1)
% and T2.data(FIX2), or T1(FIX1) and T2(FIX2), is done.) Either FIX1 or FIX2
% may also be a function handle, in which case FIXn = FEVAL(FIXn,Tn).
%
% Last Saved Time-stamp: <Sun 2012-06-03 20:50:19  Lew.Gramer>

  if ( isfield(t1,'data') )
    if ( ~isfield(t2,'data') )
      error('T1 and T2 must either both be time series structs, or both vectors!');
    end;
    nt1 = numel(t1.data);
    nt2 = numel(t2.data);
  else
    if ( isstruct(t2) )
      error('T1 and T2 must either both be time series structs, or both vectors!');
    end;
    nt1 = numel(t1);
    nt2 = numel(t2);
    if ( nt1 ~= nt2 )
      error('Vectors T1 and T2 must both have the same number of elements!');
    end;
  end;

  if ( ~exist('fix1','var') || isempty(fix1) )
    fix1 = 1:nt1;
  elseif ( isa(fix1,'function_handle') )
    fix1 = fix1(t1);
  end;
  if ( ~exist('fix2','var') || isempty(fix2) )
    fix2 = 1:nt2;
  elseif ( isa(fix2,'function_handle') )
    fix2 = fix2(t2);
  end;

  if ( isfield(t1,'date') )
    [ix1,ix2] = intersect_dates(t1.date(fix1),t2.date(fix2));
    fix1 = fix1(ix1);
    fix2 = fix2(ix2);
    t1 = t1.data(fix1);
    t2 = t2.data(fix2);
  else
    [ig,I1,I2] = intersect(fix1,fix2);
    fix1 = fix1(I1);
    fix2 = fix2(I2);
    t1 = t1(fix1);
    t2 = t2(fix2);
  end;

  if ( isscalar(t1) || isscalar(t2) )
    c = 0;
    r = 0;
    p = NaN;

  else
    C = cov(t1,t2);
    c = C(1,2);

    if ( nargout > 1 )
      [R,P] = corrcoef(t1,t2);
      r = R(1,2);
      p = P(1,2);
    end;
  end;

return;
