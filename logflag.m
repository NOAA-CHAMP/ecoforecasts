function cm = logflag(M)
%function cm = logflag(M)
% Return a COLORMAP similar to FLAG(v.) but with logarithmically increasing
% color intensity on both the upper (red) and lower (blue) half.

  if ( nargin > 0 )
    n = floor((M-1)/2);
  else
    n = 5;
  end;

  cm = zeros((n*2)+1,3);

  cm(1:n,3)     = logspace(-3,0,n);
  cm(n+1,:)     = [1 1 1];
  cm(n+2:end,1) = logspace(-3,0,n);

return;
