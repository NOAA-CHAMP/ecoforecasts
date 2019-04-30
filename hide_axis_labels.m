function hide_axis_labels(ax,xyz)
%function hide_axis_labels([ax[,xyz]])
%

  if ( ~exist('ax','var') )
    ax = gca;
  end;
  if ( ~exist('xyz','var') )
    xyz = 'x';
  end;
  set(get(ax,[xyz,'axis']),'Vis','off');

return;
