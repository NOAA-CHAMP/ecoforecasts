function y = fixn(x,n) 
% 
% See also ROUNDN, ROUND, FLOOR, CEIL, FIX (built-ins); FLOORN, CEILN (Ecoforecasts).

  if ( exist('narginchk','builtin') )
    narginchk(2,2);
    nargoutchk(0,1);
  else
    error(nargchk(2,2,nargin)) 
    error(nargoutchk(0,1,nargout)) 
  end;

  factors  = 10 .^ (fix(-n));
  y = fix(x .* factors) ./ factors;

return;
