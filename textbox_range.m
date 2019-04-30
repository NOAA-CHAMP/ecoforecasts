function anh = textbox_range(varargin)
%function anh = textbox_range([ax,][YnotX,]xcoords[,lbls[,locs]][,args...])
%
% Place very light gray horizontal "textbox" rectangles spaced horizontally
% along axes AX (DEFAULT: GCA) between pairs of x-coordinates XCOORDS(:,1:2).
% Place  optional label strings LBLS at LOCations, each element of which is a
% string concatenation of either 'upper', 'center', or 'lower' followed by
% either 'left', 'center', or 'right'; DEFAULT: 'ul' for 'UpperLeft') in each
% textbox. Optional ARGS are passed through to ANNOTATION (v.), if given.
%
% If optional YnotX is logical True (not numeric), XCOORDS is interpreted as
% an array of Y coordinate ranges for horizontal textboxes instead.
%
% If AX is a FIGURE, e.g., GCF (v.), add a textbox to each AXES child of AX.
%
% CALLS: DS2NFU (convert data space coordinates into Normalized Figure Units)
%
% Last Saved Time-stamp: <Sun 2018-04-15 19:54:26 Eastern Daylight Time gramer>

  args = varargin;

  %lightgray = [.8,.8,.8];
  lightgray = [.6,.6,.6];

  if ( numel(args) > 0 && all(ishandle(args{1})) )
    ax = args{1};
    args(1) = [];
    if ( isscalar(ax) && strcmpi(get(ax,'type'),'figure') )
      ax = findobj(ax,'type','axes');
    elseif ( ~strcmpi(get(ax(1),'type'),'axes') )
      error('Optional first arg neither a FIGURE nor an AXES (array)');
    end;
  else
    ax = gca;
  end;


  if ( numel(args) > 0 && islogical(args{1}) )
    YnotX = args{1};
    args(1) = [];
  else
    YnotX = false;
  end;

  xcoords = args{1};
  args(1) = [];
  if ( numel(xcoords) == 2 )
    % Make sure a single range is a row vector
    xcoords = xcoords(:)';
  end;

  if ( numel(args) > 0 && (ischar(args{1}) || iscellstr(args{1}) || isstring(args{1})) )
    lbls = args{1};
    args(1) = [];
  else
    lbls = '';
  end;
  if ( ischar(lbls) || iscellstr(lbls) )
    lbls = string(lbls);
  end;
  if ( strcmpi(lbls,'ROMAN') )
    lbls = num2roman(1:size(xcoords,1));
    if ( strcmp(lbls,'roman') )
      lbls = lower(lbls);
    end;
  end;
  if ( numel(lbls) == 1 )
    lbls = repmat(lbls,[1,size(xcoords,1)]);
  end;
  if ( numel(lbls) == 1 )
    lbls = {lbls};
  end;

  if ( numel(args) > 0 && (ischar(args{1}) || isstring(args{1}) ) )
    locs = args{1};
    args(1) = [];
  else
    locs = 'ul';
  end;

  switch (lower(locs))
   case {'ul','upperleft'},
    alignargs = {'Vert','top','Horiz','left'};
   case {'uc','uppercenter','um','uppermiddle'},
    alignargs = {'Vert','top','Horiz','center'};
   case {'ur','upperright'},
    alignargs = {'Vert','top','Horiz','right'};

   case {'cl','centerleft','ml','middleleft'},
    alignargs = {'Vert','middle','Horiz','left'};
   case {'cc','center','centercenter','cm','centermiddle','mc','middlecenter','mm','middlemiddle'},
    alignargs = {'Vert','middle','Horiz','center'};
   case {'cr','centerright','mr','middleright'},
    alignargs = {'Vert','middle','Horiz','right'};

   case {'ll','lowerleft'},
    alignargs = {'Vert','bottom','Horiz','left'};
   case {'lc','lowercenter','lm','lowermiddle'},
    alignargs = {'Vert','bottom','Horiz','center'};
   case {'lr','lowerright'},
    alignargs = {'Vert','bottom','Horiz','right'};

   otherwise,
    error('Unknown location "%s"',locs);
  end;

  % Loop through all AXES caller passed in
  for axix = 1:numel(ax)
    axes(ax(axix));
    if ( YnotX )
      yl = xlim;
    else
      yl = ylim;
    end;
    for ix = 1:size(xcoords,1)
      if ( YnotX )
        coords = [yl(1),xcoords(ix,1),(yl(2)-yl(1)),xcoords(ix,2)-xcoords(ix,1)];
      else
        coords = [xcoords(ix,1),yl(1),xcoords(ix,2)-xcoords(ix,1),(yl(2)-yl(1))];
      end;
      rectpos = ds2nfu(coords);
      anh(ix,axix) = ...
          annotation('textbox',rectpos,'Color','k','BackgroundColor',...
                     lightgray,'FaceAlpha',0.3,'EdgeColor','none', ...
                     'String',lbls{ix},'FontSize',24,'FontName','Times', ...
                     alignargs{:},args{:});
    end;
  end;

return;
