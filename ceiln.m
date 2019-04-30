function y = ceiln(x,n) 
% 
% See also ROUNDN, ROUND, FLOOR, CEIL, FIX (built-ins); FLOORN, FIXN (Ecoforecasts).

  if ( exist('narginchk','builtin') )
    narginchk(2,2);
    nargoutchk(0,1);
  else
    error(nargchk(2,2,nargin));
    error(nargoutchk(0,1,nargout));
  end;

  factors  = 10 .^ (fix(-n));
  y = ceil(x .* factors) ./ factors;

return;
