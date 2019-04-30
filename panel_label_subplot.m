function ah=panel_label_subplot(varargin)
%function ah=panel_label_subplot([ax,]lbl[,corner[,ARGS...]])
%
% Add a panel label (ANNOTATION 'textbox') to the corner of the subplot
% specified by AXES handle or SUBPLOT position (v.) AX (DEFAULT: GCA).
% Optional CORNER string may be one of 'lowerleft' (DEFAULT) or 'll',
% 'upperleft' or 'ul', 'upperright' or 'ur', or 'lowerright' or 'lr'.
% If CORNER is given, additional pass-through name-value pairs for
% ANNOTATION may follow it, e.g., 'FontSize',x, 'LineWidth',y, etc.
%
% Last Saved Time-stamp: <Mon 2018-03-26 15:05:53 Eastern Daylight Time gramer>

  args = varargin;
  if ( numel(args) > 0 && ishandle(args{1}) )
    ax = args{1};
    args(1) = [];
  else
    ax = gca;
  end;

  fh = getfig(ax);

  if ( numel(args) > 0 && ischar(args{1}) )
    lbl = args{1};
    args(1) = [];
  else
    error('No label specified for subplot');
  end;
  if ( numel(args) > 0 && ischar(args{1}) )
    crnr = args{1};
    args(1) = [];
  else
    crnr = 'lowerleft';
  end;

  axpos = get(ax, 'position');
  %DEBUG:  disp(axpos);
  switch (lower(crnr)),
   case {'lowerleft','ll'},
    txtpos = axpos .* [1.0 1.0 0.5 0.5];
    ah=annotation(fh,'textbox',txtpos,'String',lbl,'FitBoxToText','on','VerticalAlignment','bo','HorizontalAlignment','l',args{:});

   case {'upperleft','ul'},
    txtpos = axpos .* [1.0 1.0 0.5 1.0];
    ah=annotation(fh,'textbox',txtpos,'String',lbl,'FitBoxToText','on','VerticalAlignment','to','HorizontalAlignment','l',args{:});

   case {'upperright','ur'},
    txtpos = axpos .* [1.0 1.0 1.0 1.0];
    ah=annotation(fh,'textbox',txtpos,'String',lbl,'FitBoxToText','on','VerticalAlignment','to','HorizontalAlignment','r',args{:});

   case {'lowerright','lr'},
    txtpos = axpos .* [1.0 1.0 1.0 1.0];
    ah=annotation(fh,'textbox',txtpos,'String',lbl,'FitBoxToText','on','VerticalAlignment','bo','HorizontalAlignment','r',args{:});

   otherwise,
    error('Have not implemented corner "%s"',crnr);
  end;

return;
