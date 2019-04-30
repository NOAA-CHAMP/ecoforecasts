function [isobath,transects,cs,ch,cobj,clhs] = map_freef(boundingbox,isodepth,offsetLabels,axesOpts)
%function [isobath,transects,cs,ch,cobj,clhs] = map_freef(boundingbox,isodepth,offsetLabels,axesOpts)
%
% Draw a map of the Florida Reef Tract coastline, together with one or more
% isobaths at depth(s) ISODEPTH (DEFAULT [-350m]), suitable for plotting
% current vectors or property contour fields. BOUNDINGBOX defines the corners
% of the map drawn. If OFFSETLABELS is given, non-empty, not false (DEFAULT
% False), isobath contour labels are offset from contour lines (v. CONTOURF).
% AXESOPTS (DEFAULT {'same','nohit'}) controls AXES where contours are drawn:
% 'top'/'bot[tom]'=new AXES on top/bottom, 'same'=draw contours in current
% AXES; 'hit'/'nohit' = AXES with contours is/is not selectable (the property
% 'HitTest' in AXES and ContourGroup objects is set to 'on' or 'off', resp.)
%
% Returns coordinates of inner (Florida near-shore) isobath, coordinates of
% periodic TRANSECT lines perpendicular to that isobath, and the matrix and
% CONTOURGROUP handle returned by the CONTOUR of isobaths. If ISODEPTH is all
% NaN values, no isobaths are plotted or returned.
%
% DEFAULT: BOUNDINGBOX = [-80.50 -79.10 +24.80 +25.90];
%
% Last Saved Time-stamp: <Fri 2012-03-02 12:40:52  lew.gramer>

  global sofla_coast;
  global freef_topo;

  isobath=[];
  transects=[];
  cs=[];
  ch=[];
  cobj=[];
  clhs=[];

  if ( ~exist('boundingbox', 'var') || isempty(boundingbox) )
    boundingbox = [-80.50 -79.10 +24.80 +25.90];
  end;
  if ( ~exist('isodepth', 'var') || isempty(isodepth) )
    isodepth = -350;
  end;
  if ( ~isnumeric(isodepth) )
    if ( strcmpi(isodepth,'none') )
      isodepth = [nan nan];
    else
      error('Optional 2nd arg ISODEPTH should be numeric or the string "none"!');
    end;
  elseif ( isscalar(isodepth) && ~isnan(isodepth) )
    isodepth = [ isodepth isodepth ];
  end;
  if ( ~exist('axesOpts','var') || isempty(axesOpts) )
    axesOpts = {'same','nohit'};
  end;

  % AXES Options - see above. Allow caller to specify simple string
  if ( ischar(axesOpts) )
    % But don't be picky - let caller specify comma-separated string as well
    axesOpts = strread(axesOpts,'%s','delimiter',',; ')
  end;
  axesOpts = lower(axesOpts);


  isobath = [];
  transects = [];
  cs = [];
  ch = [];

  nohit = true;
  if ( ismember('hit',axesOpts) )
    nohit = false;
  elseif ( ismember('nohit',axesOpts) )
    nohit = true;
  end;

  cax = gca;
  if ( ismember('same',axesOpts) )
    hold on;
    ax = cax;

  else  % 'top' or 'bottom'
    pos = get(cax,'Position');
    ax = axes('position',pos);
    set(ax,'Color','none');
    if ( ~nohit )
      set(ax,'HitTest','on');
    else
      set(ax,'HitTest','off');
    end;

    if ( ismember('bottom',axesOpts) || ismember('bot',axesOpts) )
      %%%% ??? What does this do exactly?!
    end;
  end;

  axes(ax);
  %DEBUG:  cax,ax,gca,

  % Load and draw (low-) medium- (full-)resolution coastline

  if ( ~exist('sofla_coast', 'var') || isempty(sofla_coast) )
    disp('Reloading coastline');
    % sofla_coast = load('sofla_coast_low.dat');
    sofla_coast = load('sofla_coast_medium.dat');
    % sofla_coast = load('sofla_coast.dat');
  end;

  %line(sofla_coast(:,1), sofla_coast(:,2), ...
  %     'Color', [.4 .3 .2]);
  cobj = fill(sofla_coast(:,1), sofla_coast(:,2), [.4 .3 .2]);
  if ( ~nohit )
    set(cobj,'HitTest','on');
  else
    set(cobj,'HitTest','off');
  end;


  % Load and draw requested isobath from ETOPO1 1-arcmin global relief model:
  % http://www.ngdc.noaa.gov/mgg/global/relief/ 
  % Amante, C. and B. W. Eakins, ETOPO1 1 Arc-Minute Global Relief Model: Procedures, Data Sources and Analysis. NOAA Technical Memorandum NESDIS NGDC-24, 19 pp, March 2009.

  if ( ~exist('freef_topo', 'var') || isempty(freef_topo) )
    disp('Reloading freef_topo...');

    freef_topo = load('freef.topo.mat');

    % This mesh would be HUGE: We don't really seem to need it??
    % [freef_topo.lons,freef_topo.lats] = meshgrid(freef_topo.tlon,freef_topo.tlat);
  end;

  latix = find(boundingbox(3) <= freef_topo.tlat & freef_topo.tlat <= boundingbox(4));
  lonix = find(boundingbox(1) <= freef_topo.tlon & freef_topo.tlon <= boundingbox(2));
  lon = freef_topo.tlon(lonix);
  lat = freef_topo.tlat(latix);
  topo = freef_topo.topo(latix,lonix);

  if ( all(isnan(isodepth)) )
    disp('Isobars not plotted per user request...');
    map_freef_return(cax,ax,boundingbox);
    return;
  end;

  if ( numel(topo) < 4 )
    warning('map_freef:NotEnoughTopo', ...
            'Too few topo points in bounding box: no isobaths rendered');
    map_freef_return(cax,ax,boundingbox);
    return;
  end;

  isobathColor = [.7 .4 .3];

  [cs, ch] = contour(ax, lon, lat, topo, ...
                     isodepth, 'LineColor', isobathColor);
  if ( ~nohit )
    set(ch,'HitTest','on');
  else
    set(ch,'HitTest','off');
  end;
  if ( isempty(cs) )
    warning('map_freef:NoContourLines', ...
            'No isobath contours found in bounding box!');
    map_freef_return(cax,ax,boundingbox);
    return;
  end;

  if ( exist('offsetLabels','var') && ~isempty(offsetLabels) && (offsetLabels ~= false) )
    % Put "+" within contours, and place labels nearby instead...
    clhs = clabel(cs, 'Color',isobathColor);
  else
    clhs = clabel(cs, ch, 'Color',isobathColor);
  end;
  for clh=clhs(:)'
    if ( ~nohit )
      set(clh,'HitTest','on');
    else
      set(clh,'HitTest','off');
    end;
  end;


  if ( nargout > 0 )
    % Find lat/lon positions (nearest Florida) for desired isobath
    isobath = cs(1:2, 2:(cs(2,1) + 1));

    if ( nargout > 1 )
      % Subset gradient vectors at each interior point of isobath
      %[topx, topy] = gradient(sofla_topo);
      %lnix = ismember(sofla_topo_lons, isobath(1,:));
      %ltix = ismember(sofla_topo_lats, isobath(2,:));
      %topx = topx(lnix & ltix);
      %topy = topy(lnix & ltix);
      %transects(1,:,:) = isobath(:,2:end-1) - [topx ; topy];
      %transects(2,:,:) = isobath(:,2:end-1) + [topx ; topy];

      % Calculate secant and normal slopes at each interior point
      sct = isobath(:,3:end) - isobath(:,1:end-2);
      nrm = [sct(2,:) ; -sct(1,:)];
      % Calculate straight transects orthogonal to our isobath
      transects = repmat(nan, [2 2 length(isobath(1,1:end-2))]);
      transects(1,:,:) = isobath(:,2:end-1) - (2.*nrm);
      transects(2,:,:) = isobath(:,2:end-1) + (2.*nrm);

      %scts = [];
      %line(scts(:,1,:), scts(:,2,:));
    end;

  end;

  map_freef_return(cax,ax,boundingbox);

return;



%%%%%%%%%%%%%%%%%%%
% PRIVATE FUNCTIONS
%%%%%%%%%%%%%%%%%%%

function map_freef_return(cax,ax,boundingbox)
  set(ax, 'XLim', [boundingbox(1) boundingbox(2)]);
  set(ax, 'YLim', [boundingbox(3) boundingbox(4)]);
  % linkaxes([ax,cax],'xy');
  hlink = linkprop([cax,ax],{'DataAspectRatio','Position','View','XLim','YLim',});
  set(ax,'UserData',hlink);
  axes(cax);
  %DEBUG:  cax,ax,gca,
return;
