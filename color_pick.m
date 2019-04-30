function clr = color_pick(ix,clrs)
%function clr = color_pick(ix,clrs)
%
% Return the IXth color in the list (CHAR or cell array) of colors CLRS
% (DEFAULT: {'b','c','g','m','r','k'}). If IX > NUMEL(CLRS), cycle through
% the colors repeatedly. See PLOT for the list of colors. NOTE: If CLRS is a
% cell array, any of its elements may also be an RGB triplet. If IX is the
% CHAR string 'NCOLORS', then the number of elements in CLRS is returned.
%
% Last Saved Time-stamp: <Mon 2017-05-01 15:22:00 Eastern Daylight Time gramer>

  if ( ~exist('clrs','var') || isempty(clrs) )
    clrs = {'b','r',[0,0.5,0],'m','c','k'};
  end;
  if ( ~iscell(clrs) )
    if ( ischar(clrs) )
      clrs = cellstr(clrs(:));
    else
      error('Arg CLRS must be a char or cell array of colors (strings or RGB triplets)');
    end;
  end;
  if ( ischar(ix) && strncmpi(ix,'N',1) )
    clr = numel(clrs);
  else
    clr = clrs{mod(ix-1,numel(clrs))+1};
  end;

return;
