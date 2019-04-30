function bath = find_shallow_slopes(stn_or_center_point,rad,N);
%function bath = find_shallow_slopes(stn_or_center_point,rad,N);
%
% Adapted from: $GRAMER/Postdoc/FIND_SHALLOW_SLOPES_GOM_CARIB.M
%
% Search seafloor topography (NGDC 10-93 m Digital Elevation Model) around a
% center point, for all gridpoints with "appropriate depth" and "sufficient
% seafloor slope" for horizontal convection (parameters described below).
%
% If DOPRINT, PRINT figure in LOCAL directory; if DOSAVE, save .MAT file.
%
% Last Saved Time-stamp: <Tue 2016-07-19 13:13:35 Eastern Daylight Time lew.gramer>


  % If true, create figure and MAT file (resp.) in LOCAL DIRECTORY
  doPrint = false;
  doSave = false;

  doCoast = true;

  if ( isstruct(stn_or_center_point) && isfield(stn_or_center_point,'lon') )
    stn = stn_or_center_point;
    center_point = [stn.lon,stn.lat];
  elseif ( ischar(stn_or_center_point) )
    stn = get_station_from_station_name(stn_or_center_point);
    center_point = [stn.lon,stn.lat];
  elseif ( isnumeric(stn_or_center_point) && numel(stn_or_center_point) == 2 )
    center_point = stn_or_center_point;
    stn.lon = center_point(1);
    stn.lat = center_point(2);
  else
    error('First arg must be a station STRUCT or name, or a LON,LAT coordinate pair');
  end;

  if ( ~exist('rad','var') || isempty(rad) )
    rad = [13e3,13e3]; % Search for lots of sites! (USVI, E or W PR, SMKF1)
  elseif ( ~isnumeric(rad) || numel(rad) > 2 )
    error('Optional second arg should be a numeric 1- or 2-vector of search radii');
  end;
  if ( numel(rad) == 1 )
    rad = [rad,rad];
  end;

  if ( ~exist('N','var') || isempty(N) )
    N = 5;
  end;

  % Isobath contours and contour labels to plot
  cs = -[0:2:60];
  csl = false;


  %% Suitability criteria for study moorings

  % Min & max depth for *initiating* horizontal convection
  %depth_rng = [-2,-15]; % Search for deeper steeper
  %depth_rng = [-1,-10]; % Search for shallower-shallower
  depth_rng = [-1,-15]; % Conservative w.r.t. HC, liberal with insolation

  % Max depth for *terminating* horizontal convection
  %mld = -50; % Search for deeper steeper
  mld = -45; % Conservative w.r.t. HC, liberal with insolation

  % Minimum seafloor slope for *initiating* horizontal convection
  % Slope = Rise/Run or "beta" (Monismith et al. 2006)
  %min_slope = 0.0200; % Search for deeper steeper
  %min_slope = 0.0100; % Search for shallower-shallower
  %min_slope = 0.0100; % Search for lots of sites
  min_slope = 0.0200; % Conservative w.r.t. HC, liberal with insolation

  % Minimum seafloor slope for *sustaining* horizontal convection
  % Slope = Rise/Run or "beta" (Monismith et al. 2006)
  %min_sustained_slope = 0.0100; % Search for deeper steeper
  %min_sustained_slope = 0.0050; % Search for shallower-shallower
  %min_sustained_slope = 0.0050; % Search for lots of sites
  min_sustained_slope = 0.0075; % Conservative w.r.t. HC, liberal with insolation

  % Maximum gradient deviation (deg) for *sustaining* horizontal convection
  %max_dev_deg = 90; % Search for deeper steeper
  %max_dev_deg = 60; % Search for shallower-shallower
  %max_dev_deg = 45; % Search for lots of sites
  max_dev_deg = 60; % Conservative w.r.t. HC, liberal with insolation


  tic,

  %% Get bathymetry data
  if ( ~isfield(stn,'ngdc_hires_bathy') )
    stn = read_hires_bathymetry(stn,rad,[],true);
  end;

  [bet,ang,stn.ngdc_hires_bathy] = find_ngdc_slope(stn.ngdc_hires_bathy,stn.lon,stn.lat,N);

  % Find all sites with target depth and slope
  [siteyix,sitexix] = find(max(depth_rng) >= stn.ngdc_hires_bathy.field ...
                           & stn.ngdc_hires_bathy.field >= min(depth_rng) ...
                           & stn.ngdc_hires_bathy.beta >= min_slope);
  siteix = sub2ind(size(betas),siteyix,sitexix);
  sites = [sitexix,siteyix,lon(sitexix),lat(siteyix),betas(siteix),dzs(siteix),dzns(siteix),dzdxs(siteix),dzdys(siteix)];

  %% Find "sustained" slopes (see comments below)
  disp('Find sustained-slope sites');

  % OPTION: Look for "N-grid-point" sustained slopes
  %tgtext = 180;
  tgtext = 270;
  %tgtext = 350;
  npts = ceil(tgtext/dx)-1;
  ext = floor((npts+1)*dx);

  mysix = 1:size(sites,1);
  cumbet = sites(mysix,5);
  for nix = 1:npts
    jx = sitexix(mysix);
    ix = siteyix(mysix);
    bet = betas(siteix(mysix));
    dzn = sites(mysix,7);
    dix = nix*round(cosd(dzn));
    djx = nix*round(sind(dzn));

    nbet = betas(sub2ind(size(zs),ix+dix,jx+djx));
    ndzn = dzns(sub2ind(size(zs),ix+dix,jx+djx));
    nmz = mzs(sub2ind(size(zs),ix+dix,jx+djx));
    matchix = find(max(depth_rng)>=nmz & nmz>=mld & ...
                   nbet>=min_sustained_slope & cosd(ndzn-dzn)>cosd(max_dev_deg));
    mysix = mysix(matchix);
    cumbet(mysix) = cumbet(mysix) + nbet(matchix);
  end;
  sustained_sites=sites(mysix,:);
  cumbet(mysix) = cumbet(mysix) / (npts+1);
  sustained_sites(:,5)=cumbet(mysix);

  cumbet=[];
  clear nix jx ix bet dzn dix djx nbet ndzn nmz matchix cumbet


  %% Plot sustained-slope target sites
  disp('Ready to plot');
  %pause;

  % plot_bath_hc_sites(lon,lat,zs,cs,csl,sustained_sites,doCoast);
  % %plot(bath.lon,bath.lat,'wp', 'MarkerEdgeColor','black', 'MarkerFaceColor','white');
  plot_hires_bathymetry(bath,cs,rad,doCoast,@contourf,csl);
  lh=plot_hc_target_sites(sustained_sites(:,3),sustained_sites(:,4),sustained_sites(:,5));
  titlename(sprintf('Sites with (very, somewhat, slight) reduced exposure to thermal stress (%g m)',ext));
  if (doPrint)
    disp(['Printing']);
    print('-dpng',sprintf('shallow_slopes_%g_%g_%g_%gm.png',center_point,min_slope,ext));
  end;


  if (doSave)
    matfname = sprintf('shallow_slopes_%g_%g_%g_%gm.mat',center_point,min_slope,ext);
    disp(['Saving ',matfname]);
    save(matfname);
  end;

  toc,
  size(sites),
  size(sustained_sites),


error('UNEDITED CODE FROM HERE...');



  %% Get bathymetry data
  if ( ~isfield(stn,'ngdc_hires_bathy') )
    stn = read_hires_bathymetry(stn,rad,[],true);
  end;

  zs = stn.ngdc_hires_bathy.field;
  lat = stn.ngdc_hires_bathy.lat(:);
  lon = stn.ngdc_hires_bathy.lon(:);

  %dx = 10;
  dE = median(distance_wgs84(lat(1),lon(1:end-1),lat(1),lon(2:end))*1e3);
  dN = median(distance_wgs84(lat(1:end-1),lon(1),lat(2:end),lon(1))*1e3);
  dx = max(dE,dN);


  %% Calculate seafloor slopes
  disp('Estimate betas');

  % Find all points on the SHELF
  [yix,xix,zix] = find(-200<zs&zs<0);
  % Avoid points on the edge (see smoothing below)
  badix = find(yix==1|yix==size(zs,1)|xix==1|xix==size(zs,2));
  yix(badix)=[];
  xix(badix)=[];

  % Smooth with a 3x3 mean
  zar = repmat(nan,[9,numel(yix)]);
  dix = 1;
  for dyix=-1:1
    for dxix=-1:1
      zar(dix,:) = zs(sub2ind(size(zs),yix+dyix,xix+dxix));
      dix = dix + 1;
    end;
  end;
  mzsvec = nanmean(zar,1);
  zar=[]; clear zar dix dyix dxix

  mzs = repmat(nan,size(zs));
  mzs(sub2ind(size(zs),yix,xix)) = mzsvec;
  mzsvec=[]; clear mzsvec;

  % Use an N-point finite-difference stencil to estimate depth gradients
  nrad = floor(N/2);
  dzdxs = repmat(nan,size(zs));
  dzdys = repmat(nan,size(zs));
  for ixix=nrad+1:size(zs,1)-nrad
    dzdxs(ixix,nrad+1:size(zs,2)-nrad) = findiff(zs(ixix,:),N);
  end;
  for jxix=nrad+1:size(zs,2)-nrad
    dzdys(nrad+1:size(zs,1)-nrad,jxix) = findiff(zs(:,jxix),N);
  end;
  dzs = uv_to_spd(dzdxs,dzdys);
  dzns = uv_to_dir_curr(dzdxs,dzdys);
  betas = dzs./dx;

  % Find all sites with target depth and slope
  [siteyix,sitexix] = find(max(depth_rng)>=mzs & mzs>=min(depth_rng) & betas>=min_slope);
  siteix = sub2ind(size(betas),siteyix,sitexix);
  sites = [sitexix,siteyix,lon(sitexix),lat(siteyix),betas(siteix),dzs(siteix),dzns(siteix),dzdxs(siteix),dzdys(siteix)];



  %% Find "sustained" slopes (see comments below)
  disp('Find sustained-slope sites');

  % OPTION: Look for "N-grid-point" sustained slopes
  %tgtext = 180;
  tgtext = 270;
  %tgtext = 350;
  npts = ceil(tgtext/dx)-1;
  ext = floor((npts+1)*dx);

  mysix = 1:size(sites,1);
  cumbet = sites(mysix,5);
  for nix = 1:npts
    jx = sites(mysix,1);
    ix = sites(mysix,2);
    bet = sites(mysix,5);
    dzn = sites(mysix,7);
    dix = nix*round(cosd(dzn));
    djx = nix*round(sind(dzn));

    nbet = betas(sub2ind(size(zs),ix+dix,jx+djx));
    ndzn = dzns(sub2ind(size(zs),ix+dix,jx+djx));
    nmz = mzs(sub2ind(size(zs),ix+dix,jx+djx));
    matchix = find(max(depth_rng)>=nmz & nmz>=mld & ...
                   nbet>=min_sustained_slope & cosd(ndzn-dzn)>cosd(max_dev_deg));
    mysix = mysix(matchix);
    cumbet(mysix) = cumbet(mysix) + nbet(matchix);
  end;
  sustained_sites=sites(mysix,:);
  cumbet(mysix) = cumbet(mysix) / (npts+1);
  sustained_sites(:,5)=cumbet(mysix);

  cumbet=[];
  clear nix jx ix bet dzn dix djx nbet ndzn nmz matchix cumbet


  %% Plot sustained-slope target sites
  disp('Ready to plot');
  %pause;

  % plot_bath_hc_sites(lon,lat,zs,cs,csl,sustained_sites,doCoast);
  % %plot(bath.lon,bath.lat,'wp', 'MarkerEdgeColor','black', 'MarkerFaceColor','white');
  plot_hires_bathymetry(bath,cs,rad,doCoast,@contourf,csl);
  lh=plot_hc_target_sites(sustained_sites(:,3),sustained_sites(:,4),sustained_sites(:,5));
  titlename(sprintf('Sites with (very, somewhat, slight) reduced exposure to thermal stress (%g m)',ext));
  if (doPrint)
    disp(['Printing']);
    print('-dpng',sprintf('shallow_slopes_%g_%g_%g_%gm.png',center_point,min_slope,ext));
  end;


  if (doSave)
    matfname = sprintf('shallow_slopes_%g_%g_%g_%gm.mat',center_point,min_slope,ext);
    disp(['Saving ',matfname]);
    save(matfname);
  end;

  toc,
  size(sites),
  size(sustained_sites),

return;
