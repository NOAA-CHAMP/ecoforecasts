function center_axes(varargin)
%function center_axes([ax,][axis_name])
%
% Ensure that the origin of axis AXIS_NAME ('x', 'y', 'z', or 'all', DEFAULT:
% 'y') in AXES AX (DEFAULT: GCA) lies at the center of the plot.
%
% Last Saved Time-stamp: <Wed 2017-06-07 10:42:06 Eastern Daylight Time gramer>

  ax = gca;
  axis_name = 'y';

  args = varargin;
  if ( numel(args) > 0 )
    if ( ishandle(args{1}) )
      ax = args{1};
      args(1) = [];
    end;
    if ( ischar(args{1}) && numel(args{1}) > 0 )
      axis_name = args{1};
      args(1) = [];
    end;
  end;

  switch ( lower(axis_name) ),
   case 'x',
    lm = xlim(ax);    xlim(ax,[-max(abs(lm)),max(abs(lm))]);
   case 'y',
    lm = ylim(ax);    ylim(ax,[-max(abs(lm)),max(abs(lm))]);
   case 'z',
    lm = zlim(ax);    zlim(ax,[-max(abs(lm)),max(abs(lm))]);
   case 'all',
    lm = xlim(ax);    xlim(ax,[-max(abs(lm)),max(abs(lm))]);
    lm = ylim(ax);    ylim(ax,[-max(abs(lm)),max(abs(lm))]);
    lm = zlim(ax);    zlim(ax,[-max(abs(lm)),max(abs(lm))]);
   otherwise,
    error('Axis name must be ''x'', ''y'', ''z'' or ''all''!');
  end;

return;
