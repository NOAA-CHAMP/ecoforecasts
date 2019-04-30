function outvec = findiff2(invec,n)
%function outvec = findiff2(invec,n)
%
% Calculate *SECOND-ORDER* finite difference using N points (DEFAULT: 3) on
% numeric vector INVEC, and return it in OUTVEC. INVEC is assumed to have no
% gaps. If N is imaginary, a *forward* finite difference is used with ABS(N)
% points; if N<0, use *backward* finite difference with ABS(N) points. NOTE:
% LENGTH(OUTVEC) relative to INVEC depends on the given value of N.
%
% Coefficients below were derived from tables at this Web site:
%   https://en.wikipedia.org/wiki/Finite_difference_coefficient
%
% Last Saved Time-stamp: <Mon 2016-07-18 18:40:32 Eastern Daylight Time gramer>

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
    outvec = invec(1:end-2) ...
             + (-2)*invec(2:end-1) ...
             + invec(3:end);
   case 5,
    outvec = (-1/12)*invec(1:end-4) + (4/3)*invec(2:end-3) ...
             + (-5/2)*invec(3:end-2) ...
             + (4/3)*invec(4:end-1) + (-1/12)*invec(5:end);
   case 7,
    outvec = (1/90)*invec(1:end-6) + (-3/20)*invec(2:end-5) + (3/2)*invec(3:end-4) ...
             + (-49/18)*invec(4:end-3) ...
             + (3/2)*invec(5:end-2) + (-3/20)*invec(6:end-1) + (1/90)*invec(7:end);
   case 9,
    outvec = (-1/560)*invec(1:end-8) + (8/315)*invec(2:end-7) + (-1/5)*invec(3:end-6) + (8/5)*invec(4:end-5) ...
             + (-205/72)*invec(5:end-4) ...
             + (8/5)*invec(6:end-3) + (-1/5)*invec(7:end-2) + (8/315)*invec(8:end-1) + (-1/560)*invec(9:end);

   % Forward
   case 1i,
    outvec = diff(invec,2);
   case 2i,
    outvec = (1)*invec(1:end-2) + (-2)*invec(2:end-1) + (1)*invec(3:end);
   case 3i,
    outvec = (2)*invec(1:end-3) + (-5)*invec(2:end-2) + (4)*invec(3:end-1) + (-1)*invec(4:end);
   case 4i,
    outvec = (35/12)*invec(1:end-4) + (-26/3)*invec(2:end-3) + (19/2)*invec(3:end-2) + (-14/3)*invec(4:end-1) + (11/12)*invec(5:end);
   case 5i,
    outvec = (15/4)*invec(1:end-5) + (-77/6)*invec(2:end-4) + (107/6)*invec(3:end-3) + (-13)*invec(4:end-2) + (61/12)*invec(5:end-1) + (-5/6)*invec(6:end);
   case 6i,
    outvec = (203/45)*invec(1:end-6) + (-87/5)*invec(2:end-5) + (117/4)*invec(3:end-4) + (-254/9)*invec(4:end-3) + (33/2)*invec(5:end-2) + (-27/5)*invec(6:end-1) + (137/180)*invec(7:end);
   case 7i,
    outvec = (469/90)*invec(1:end-7) + (-223/10)*invec(2:end-6) + (879/20)*invec(3:end-5) + (-949/18)*invec(4:end-4) + (41)*invec(5:end-3) + (-201/10)*invec(6:end-2) + (1019/180)*invec(7:end-1) + (-7/10)*invec(8:end);

   % Backward
   case -1,
    outvec = diff(invec(end:-1:1),2);
   case -2,
    outvec = (1)*invec(1:end-2) + (-2)*invec(2:end-1) + (1)*invec(3:end);
   case -3,
    outvec = (-1)*invec(1:end-3) + (4)*invec(2:end-2) + (-5)*invec(3:end-1) + (2)*invec(4:end);
   case -4,
    outvec = (11/12)*invec(1:end-4) + (-14/3)*invec(2:end-3) + (19/2)*invec(3:end-2) + (-26/3)*invec(4:end-1) + (35/12)*invec(5:end);
   case -5,
    outvec = (-5/6)*invec(1:end-5) + (61/12)*invec(2:end-4) + (-13)*invec(3:end-3) + (107/6)*invec(4:end-2) + (-77/6)*invec(5:end-1) + (15/4)*invec(6:end);
   case -6,
    outvec = (137/180)*invec(1:end-6) + (-27/5)*invec(2:end-5) + (33/2)*invec(3:end-4) + (-254/9)*invec(4:end-3) + (117/4)*invec(5:end-2) + (-87/5)*invec(6:end-1) + (203/45)*invec(7:end);
   case -7,
    outvec = (-7/10)*invec(1:end-7) + (1019/180)*invec(2:end-6) + (-201/10)*invec(3:end-5) + (41)*invec(4:end-4) + (-949/18)*invec(5:end-3) + (879/20)*invec(6:end-2) + (-223/10)*invec(7:end-1) + (469/90)*invec(8:end);

   otherwise,
    error('Do not know how to second-order finite difference with N=%s!',num2str(n));
  end;

return;
