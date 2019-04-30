function [isobath,transects,cs,ch,cobj,clhs] = map_freef_ngdc(boundingbox,isodepth,offsetLabels)
%function [isobath,transects,cs,ch,cobj,clhs] = map_freef_ngdc(boundingbox,isodepth,offsetLabels)
%
% Draw a map of the Florida Reef Tract coastline, together with one or more
% isobaths at depth(s) ISODEPTH (DEFAULT: [-350m]), suitable for plotting
% WERA current vectors or property contour fields.  BOUNDINGBOX defines the
% corners of the map drawn.  If OFFSETLABELS is given, non-empty and does not
% evaluate to false (DEFAULT: False), isobath contour labels are offset from
% contour lines. Returns coordinates of inner (Florida near-shore) isobath,
% and coordinates of periodic TRANSECT lines perpendicular to that isobath.
%
% DEFAULT: BOUNDINGBOX = [-80.50 -79.10 +24.80 +25.90];
%
% Last Saved Time-stamp: <Fri 2012-03-02 12:40:46  lew.gramer>

  global sofla_coast;
  global freef_ngdc_topo;

  isobath=[];
  transects=[];
  cs=[];
  ch=[];
  cobj=[];
  clhs=[];

  if ( ~exist('boundingbox', 'var') || isempty(boundingbox) )
    % boundingbox = [-80.50 -79.10 +24.80 +25.90];
    boundingbox = [-80.99 -80.01 +24.76 +25.64];
  end;
  if ( ~exist('isodepth', 'var') || isempty(isodepth) )
    % isodepth = -350;
    isodepth = -[2 5 10 20 40];
  end;
  if ( isscalar(isodepth) && ~isnan(isodepth) )
    isodepth = [ isodepth isodepth ];
  end;

  isobath = [];
  transects = [];

  hold on;

  % Load and draw (low-) medium- (full-)resolution coastline

  if ( ~exist('sofla_coast', 'var') || isempty(sofla_coast) )
    disp('Reloading coastline');
    % sofla_coast = load('sofla_coast_low.dat');
    % sofla_coast = load('sofla_coast_medium.dat');
    sofla_coast = load('sofla_coast.dat');
  end;

  %line(sofla_coast(:,1), sofla_coast(:,2), ...
  %     'Color', [.4 .3 .2]);
  cobj = fill(sofla_coast(:,1), sofla_coast(:,2), [.4 .3 .2]);


  % Load and draw requested isobath from NGDC Coastal Relief Map topography

  if ( ~exist('freef_ngdc_topo', 'var') || isempty(freef_ngdc_topo) )
    disp('Reloading freef_ngdc_topo...');

    stn.station_name = 'freef_ngdc';
    stn.lon = -80.5; stn.lat = 25.2; stn.depth = 0;
    stn = get_ngdc_bathy_station(stn,50,false);
    freef_ngdc_topo.topo = stn.ngdc_92m_bathy.field';
    freef_ngdc_topo.tlat = unique(stn.ngdc_92m_bathy.lat(:));
    freef_ngdc_topo.tlon = unique(stn.ngdc_92m_bathy.lon(:));
    stn = []; clear stn;
  end;

  latix = find(boundingbox(3) <= freef_ngdc_topo.tlat & freef_ngdc_topo.tlat <= boundingbox(4));
  lonix = find(boundingbox(1) <= freef_ngdc_topo.tlon & freef_ngdc_topo.tlon <= boundingbox(2));
  lon = freef_ngdc_topo.tlon(lonix);
  lat = freef_ngdc_topo.tlat(latix);
  topo = freef_ngdc_topo.topo(latix,lonix);

  if ( numel(topo) < 4 )
    warning('map_freef:NotEnoughTopo', ...
            'Too few topo points in bounding box: no isobaths rendered');
    return;
  end;

  if ( all(isnan(isodepth)) )
    disp('Isobars not plotted per user request...');
    return;
  end;


  [cs, ch] = contour(lon, lat, topo, ...
                     isodepth, 'LineColor', [.6 .5 .4]);
  if ( isempty(cs) )
    warning('map_freef:NoContourLines', ...
            'No isobath contours found in bounding box!');
    return;
  end;

  % Put "+" within contours, and place labels nearby instead...
  if ( exist('offsetLabels','var') && ~isempty(offsetLabels) && (offsetLabels ~= false) )
    clhs = clabel(cs, 'Color',[.6 .5 .4]);
  else
    clhs = clabel(cs, ch, 'Color',[.6 .5 .4]);
  end;

  set(gca, 'xlim', [boundingbox(1) boundingbox(2)]);
  set(gca, 'ylim', [boundingbox(3) boundingbox(4)]);

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

return;
