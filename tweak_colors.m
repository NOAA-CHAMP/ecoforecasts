function clrs = tweak_colors(clrs)
%function clrs = tweak_colors(clrs)
% Convert a list (CHAR or CELLSTR) of PLOT line-style color characters into a
% slightly tweaked (darker) CELL of colors and RGB triplets: 'g'reen becomes
% forest green, 'y'ellow goldenrod, gr'a'y or gr'e'y [.5,.5.,5].

  if ( ~exist('clrs','var') || isempty(clrs) )
    clrs = {'b','m','g','c','r','e','y'};
  elseif ( ischar(clrs) )
    clrs = cellstr(clrs(:));
  end;
  for ix=1:numel(clrs)
    switch ( lower(clrs{ix}) )
     case 'g',
      clrs{ix} = [0.10,0.90,0.10];
     case 'y',
      clrs{ix} = [0.75,0.75,0.00];
     case {'a','e'},
      clrs{ix} = [0.50,0.50,0.50];
     otherwise,
    end;
  end;

return;
