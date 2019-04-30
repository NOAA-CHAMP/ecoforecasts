function outvec = findiff(invec,n)
%function outvec = findiff(invec,n)
%
% Calculate N-point *centered* finite difference (DEFAULT: 3) of numeric
% vector INVEC, and return it in OUTVEC. INVEC is assumed to have no gaps.
% If N is imaginary, a *forward* finite difference is used with ABS(N)-point
% stencil; if N<0, use *backward* finite difference with ABS(N) points.
% NOTE: LENGTH(OUTVEC) relative to INVEC depends on the given value of N.
%
% Coefficients below were derived from tables at this Web site:
%   https://en.wikipedia.org/wiki/Finite_difference_coefficient
%
% Last Saved Time-stamp: <Mon 2016-07-18 17:07:21 Eastern Daylight Time gramer>

  % Default length of (centered) finite difference template
  if (~exist('n','var') || isempty(n) )
    n = 3;
  end;
  if ( ~isnumeric(n) || ~isscalar(n) )
    error('If given, N must be a numeric scalar');
  end;

  switch ( n )
   % Centered
   case 3,
    outvec = (-invec(1:(end-2)) + invec(3:end)) / 2;
   case 5,
    outvec = (invec(1:(end-4)) - 8*invec(2:(end-3)) + ...
              8*invec(4:(end-1)) - invec(5:end)) / 12;
   case 7,
    outvec = (-invec(1:(end-6)) + 9*invec(2:(end-5)) - ...
              45*invec(3:(end-4)) + 45*invec(5:(end-2)) - ...
              9*invec(6:(end-1)) + invec(7:end)) / 60;
   case 9,
    outvec = (1/280)*invec(1:end-8) + (-4/105)*invec(2:end-7) + ...
             (1/5)*invec(3:end-6) + (-4/5)*invec(4:end-5) + ...
             (4/5)*invec(6:end-3) + (-1/5)*invec(7:end-2) + ...
             (4/105)*invec(8:end-1) + (-1/280)*invec(9:end);

   % Forward
   case 1i,
    outvec = diff(invec);
   case 2i,
    outvec = (-1)*invec(1:end-1) + (1)*invec(2:end);
   case 3i,
    outvec = (-3/2)*invec(1:end-2) + (2)*invec(2:end-1) + (-1/2)*invec(3:end);
   case 4i,
    outvec = (-11/6)*invec(1:end-3) + (3)*invec(2:end-2) + ...
             (-3/2)*invec(3:end-1) + (1/3)*invec(4:end);
   case 5i,
    outvec = (-25/12)*invec(1:end-4) + (4)*invec(2:end-3) + ...
             (-3)*invec(3:end-2) + (4/3)*invec(4:end-1) + ...
             (-1/4)*invec(5:end);
   case 6i,
    outvec = (-137/60)*invec(1:end-5) + (5)*invec(2:end-4) + ...
             (-5)*invec(3:end-3) + (10/3)*invec(4:end-2) + ...
             (-5/4)*invec(5:end-1) + (1/5)*invec(6:end);
   case 7i,
    outvec = (-49/20)*invec(1:end-6) + (6)*invec(2:end-5) + ...
             (-15/2)*invec(3:end-4) + (20/3)*invec(4:end-3) + ...
             (-15/4)*invec(5:end-2) + (6/5)*invec(6:end-1) + ...
             (-1/6)*invec(7:end);

   % Backward
   case -1,
    outvec = diff(invec(end:-1:1));
   case -2,
    outvec = (1)*invec(2:end) + (-1)*invec(1:end-1);
   case -3,
    outvec = (3/2)*invec(3:end) + (-2)*invec(2:end-1) + (1/2)*invec(1:end-2);
   case -4,
    outvec = (11/6)*invec(4:end) + (-3)*invec(3:end-1) + ...
             (3/2)*invec(2:end-2) + (-1/3)*invec(1:end-3);
   case -5,
    outvec = (25/12)*invec(5:end) + (-4)*invec(4:end-1) + ...
             (3)*invec(3:end-2) + (-4/3)*invec(2:end-3) + ...
             (1/4)*invec(1:end-4);
   case -6,
    outvec = (137/60)*invec(6:end) + (-5)*invec(5:end-1) + ...
             (5)*invec(4:end-2) + (-10/3)*invec(3:end-3) + ...
             (5/4)*invec(2:end-4) + (-1/5)*invec(1:end-5);
   case -7,
    outvec = (49/20)*invec(7:end) + (-6)*invec(6:end-1) + ...
             (15/2)*invec(5:end-2) + (-20/3)*invec(4:end-3) + ...
             (15/4)*invec(3:end-4) + (-6/5)*invec(2:end-5) + ...
             (1/6)*invec(1:end-6);

   otherwise,
    error('Do not know how to finite difference with N=%s!',num2str(n));
  end;

return;
