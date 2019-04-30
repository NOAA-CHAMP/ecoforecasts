function tighten_axes(h)
%function tighten_axes(h)
%
% Set the 'LooseInset' of an AXES, or all AXES in a given FIGURE, to be the
% same as that AXES' 'TightInset' Property. This has the effect of removing
% almost all bordering white space around the AXES when a figure is plotted.
% If caller passes a handle other than a FIGURE or AXES, search up its chain
% of parents for a FIGURE; if found, tighten all that figure's child AXES.
%
% Last Saved Time-stamp: <Fri 2018-02-09 14:30:23 Eastern Standard Time gramer>

  if ( ~exist('h','var') || isempty(h) )
    h = get(0,'CurrentFigure');
  end;
  if ( ~ishandle(h) )
    error('No current figure and 1st arg was not a valid handle');
  elseif ( strcmpi(get(h,'Type'), 'axes') )
    ahs = h;
  else
    % Find the FIGURE that is the ultimate parent of handle H
    while ( ishandle(h) && ~strcmpi(get(h,'Type'), 'figure') )
      h = get(h,'Parent');
    end;
    if ( ~ishandle(h) )
      error('No FIGURE found associated with handle');
    end;
    ahs = get(h,'Children');
  end;

  for aix = 1:numel(ahs)
    ah = ahs(aix);
    if ( strcmp(get(ah,'Type'), 'axes') )
      set(ah,'LooseInset',get(gca,'TightInset'));
    end;
  end;

return;
