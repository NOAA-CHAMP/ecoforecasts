function th = annotaxis(varargin)
%function th = annotaxis([AX,][X_Y_OR_Z,]X_OR_Y_COORD[,Y_OR_Z_COORD[,X_Y_OR_Z_COORD]],ANNOTSTR[,name1,value1,...])
%
% Add a TEXT annotation string ANNOTSTR outside of the specific axis ('x',
% 'y', or 'z', DEFAULT 'y') of the AXES object AX (DEFAULT: GCA). Required
% first argument is the other coordinate at which to place annotation, e.g.,
% the X coordinate at which to place a 'y'-axis' annotation string. For 3-D
% plots, an optional third coordinate (Z, Y, or X) may also be given. Horiz.
% and Vert. alignment of TEXT is chosen based on calling arguments: any args
% given after ANNOTSTR are passed on to TEXT in sequence following the chosen
% alignment args. Additional args may thus be used to override these default
% text alignments, and/or to control TEXT color, font attrs., rotation, etc.
%
% Optional X_Y_OR_Z may also be "yr" to place Y annotation to right of AXES
% (rather than the default left), or "xt" to place X annotation above AXES.
%
% Last Saved Time-stamp: <Tue 2018-03-27 11:38:20 Eastern Daylight Time gramer>


  %% Process calling arguments
  args = varargin;

  if ( numel(args) < 2 )
    error('At least two args are required: a coordinate and a CHAR annotation string');
  end;

  % Optional: AXES handle
  if ( isscalar(args{1}) && ishandle(args{1}) && strcmpi(get(args{1},'Type'),'AXES') )
    ax = args{1};
    args(1) = [];
  else
    ax = gca;
  end;

  % Optional: axis name ('x', 'y', 'z') to annotate
  if ( ischar(args{1}) )
    x_y_or_z = args{1};
    args(1) = [];
  else
    % DEFAULT: Annotate the Y-axis
    x_y_or_z = 'y';
  end;

  % At this point in the processing of our calling arguments, we need at
  % least two (more) arguments, and the first of them must be numeric
  if ( numel(args) < 2 || ~isnumeric(args{1}) )
    error('First non-optional arg must be a numeric coordinate');
  end;

  % Required: primary coordinate for annotation (e.g., Y for 'y'-axis annotation)
  coord1 = args{1};
  args(1) = [];

  % Optional: secondary coordinate for annotation (e.g., Z for 'y'-axis annotation)
  if ( isnumeric(args{1}) )
    coord2 = args{1};
    args(1) = [];
  else
    coord2 = [];
  end;

  % Optional: third coordinate for annotation (useful in 3-D plots)
  if ( isnumeric(args{1}) )
    coord3 = args{1};
    args(1) = [];
  else
    coord3 = [];
  end;

  if ( numel(args) < 1 || ~ischar(args{1}) )
    error('Second non-optional arg must be a CHAR annotation string');
  end;

  % Required: CHAR annotation string
  txt = args{1};
  args(1) = [];


  %% Calculate coordinates and attributes
  switch ( lower(x_y_or_z) )

   case {'x','xb','xbottom'},
    pos = get(get(ax,'XLabel'),'Pos'); pos = min(pos(1:2));
    lim = min(get(ax,'YLim'));
    xcoord = coord1;
    ycoord = pos(1) - ( (lim(1) - pos(1)) / 2 );
    if ( ~isempty(coord2) )
      zcoord = coord2;
    else
      zcoord = median(get(ax,'ZLim'));
    end;
    align = {'Horizontal','right', 'Vertical','mid', 'Rotation',90};

   case {'xt','xtop'},
    pos = get(get(ax,'XLabel'),'Pos'),
    lim = get(ax,'YLim'),
    yoff = (lim(1) - pos(2)) / 2,
    xcoord = coord1;
    ycoord = lim(2) + yoff,
    if ( ~isempty(coord2) )
      zcoord = coord2;
    else
      zcoord = median(get(ax,'ZLim'));
    end;
    align = {'Horizontal','right', 'Vertical','mid', 'Rotation',270};


   case {'y','yl','yleft',},
    pos = get(get(ax,'YLabel'),'Pos'); pos = min(pos(1:2));
    lim = min(get(ax,'XLim'));
    xcoord = pos(1) - ( (lim(1) - pos(1)) / 2 );
    ycoord = coord1;
    if ( ~isempty(coord2) )
      zcoord = coord2;
    else
      zcoord = median(get(ax,'ZLim'));
    end;
    align = {'Horizontal','right', 'Vertical','mid'};

   case {'yr','yright',},
    pos = get(get(ax,'YLabel'),'Pos');
    lim = get(ax,'XLim');
    xoff = (lim(1) - pos(1)) / 2;
    xcoord = lim(2) + xoff;

    % Careful! There may be a COLORBAR on the right...
    cbh = get(ax,'Colorbar');
    if ( ishandle(cbh) )
      % POSf = DS2NFU(AX, POS) converts 4-element position vector, POS from data space to normalized figure units
      % POSd = NFU2DS(AX, POS) converts 4-element position vector, POS from normalized figure units to data space
      warning('May not do ''yright'' annotation well when there is a COLORBAR');
    end;

    ycoord = coord1;
    if ( ~isempty(coord2) )
      zcoord = coord2;
    else
      zcoord = median(get(ax,'ZLim'));
    end;
    align = {'Horizontal','left', 'Vertical','mid'};


   case 'z',
    if ( ~isempty(coord2) )
      xcoord = coord2;
    else
      xcoord = min(get(ax,'XLim'));
    end;
    if ( ~isempty(coord3) )
      ycoord = coord3;
    else
      ycoord = min(get(ax,'YLim'));
    end;
    zcoord = coord1;
    align = {'Horizontal','center', 'Vertical','mid'};

   otherwise,
    error('X_Y_OR_Z should specify which axis to annotate (''x'',''xt'',''y'',''yr'', or ''z'')');

  end; %switch ( lower(x_y_or_z) )



  %% Add annotation string to the plot - return its handle

  th = text(xcoord,ycoord,zcoord,txt,align{:},args{:});

return;
