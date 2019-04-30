function ax = subplots_set(varargin)
% axesHandles = SUBPLOTS_SET([figure_handle,]axesPropertyName1,axesPropertyValue1,...)
%
% Set Properties in all (plotting) AXES in the Current Figure (or Figure
% specified in optional first arg), using the specified name-value pairs.
%
% Last Saved Time-stamp: <Fri 2011-06-10 09:40:32  lew.gramer>


  args = varargin;

  fh = args{1};
  if ( ~ischar(fh) )
    if ( isscalar(fh) && ishandle(fh) && strcmpi(get(fh,'Type'),'figure') )
      args = args(2:end);
    else
      error('First arg must be a Figure handle or AXES property name!');
    end;
  else
    fh = gcf;
  end;

  nms = {};
  vals = {};
  for ix=1:2:length(args)
    if ( ~ischar(args{ix}) )
      error('Expected an Axes Property Name in argument %d!',ix);
    end;
    nms{end+1} = args{ix};
    if ( length(args) < ix+1 )
      error('Property name %s was given no value?',nms{end});
    end;
    vals{end+1} = args{ix+1};
  end;

  ax = [];
  kids = get(fh,'Children');
  for ix = 1:length(kids)
    kid = kids(ix);
    if ( ishandle(kid) )
      % We don't want titles, colorbars, legends, etc.!
      if ( strcmpi(get(kid,'Type'),'axes') && isempty(get(kid,'Tag')) )
        ax(end+1) = kids(ix);
      end;
    end;
  end;

  if ( length(ax) < 1 )
    error('No valid AXES in Figure?!');
  elseif ( length(ax) < 2 )
    warning('Only one valid AXES found in Figure');
  end;

  for ix=1:length(ax)
    set(ax(ix),nms,vals);
  end;

  if ( nargout < 1 )
    ax = []; clear ax;
  end;

return;
