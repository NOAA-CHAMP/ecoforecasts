1;
%% SCRIPT to test Python 'process_nwps' netCDF result file(s)

if ( ~exist('datapath','var') )
  datapath = get_ecoforecasts_path('data');
end;
if ( ~exist('dataset','var') )
  %dataset = 'key_nwps_CG2';
  dataset = 'mfl_nwps_CG1';
end;
if ( ~exist('filedates','var') )
  filedates = [datenum(2016,07,31,6,0,0):0.5:datenum(2016,08,02,18,0,0)];
end;

g = []; clear g
g = struct('lat',[],'lon',[],'t',[],'Ws',[],'Cs',[],'Hs',[],'Sw',[]);
n = []; clear n
n = struct('lat',[],'lon',[],'t',[],'Ws',[],'Cs',[],'Hs',[],'Sw',[],'Lg',[],'Lgm',[]);

for fileix = 1:numel(filedates)
  filedate = filedates(fileix);
  [y,m,d,H,M,S] = datevec(filedate);
  basefile = sprintf('%s_%04d%02d%02d_%02d00',dataset,y,m,d,H);
  clear y m d H M S;

  datafile = fullfile(datapath,[basefile,'.grib2']);
  if ( ~exist(datafile,'file') )
    disp(['Skipping ',datafile]);
    %%% EARLY LOOP TERMINATION %%%
    continue;
  end;
  nc = mDataset(datafile); 
  if ( isempty(g.lat) )
    g.lat = cast(nc{'lat'}(:),'double'); 
    g.lon = cast(nc{'lon'}(:),'double'); 
    % Convert to West/East longitudes
    g.lon(g.lon >= 180) = g.lon(g.lon >= 180) - 360;
    g.nlat = numel(g.lat);
    g.nlon = numel(g.lon);
  end;
  %% We only want the first four forecast periods
  gix = 1:4;
  hrs = cast(nc{'time'}(gix),'double');
  nhrs = numel(hrs);
  g.t(end+1:end+nhrs,1) = filedate + (hrs/24); 
  % Current_direction degrees           Current_direction @ surface
  % Current_speed  m s-1                Current_speed @ surface
  % Deviation_of_sea_level_from_mean m  Deviation_of_sea_level_from_mean @ surface
  % Primary_wave_direction degrees      Primary_wave_direction @ surface
  % Primary_wave_mean_period s          Primary_wave_mean_period @ surface
  % Significant_height_of_combined_wind_waves_and_swell m Significant_height_of_combined_wind_waves_and_swell @ surface
  % Significant_height_of_swell_waves m Significant_height_of_swell_waves @ surface
  % UnknownParameter_14 Unknown         UnknownParameter_14 @ surface
  % UnknownParameter_193 Unknown        UnknownParameter_193 @ surface
  % Wind_direction_from_which_blowing degrees Wind_direction_from_which_blowing @ surface
  % Wind_speed     m s-1                Wind_speed @ surface
  g.Cs(end+1:end+nhrs,1:g.nlat,1:g.nlon) = cast(nc{'Current_speed'}(gix,:,:),'double'); 
  g.Ws(end+1:end+nhrs,1:g.nlat,1:g.nlon) = cast(nc{'Wind_speed'}(gix,:,:),'double'); 
  g.Hs(end+1:end+nhrs,1:g.nlat,1:g.nlon) = cast(nc{'Significant_height_of_combined_wind_waves_and_swell'}(gix,:,:),'double'); 
  g.Sw(end+1:end+nhrs,1:g.nlat,1:g.nlon) = cast(nc{'Significant_height_of_swell_waves'}(gix,:,:),'double'); 
  close(nc); clear nc

  datafile = fullfile(datapath,[basefile,'.nc']);
  nc = mDataset(datafile);
  if ( isempty(n.lat) )
    n.lat = nc{'lat'}(:); 
    n.lon = nc{'lon'}(:); 
    n.nlat = numel(n.lat);
    n.nlon = numel(n.lon);
  end;
  hrs = nc{'time'}(:);
  nhrs = numel(hrs);
  n.t(end+1:end+nhrs,1) = datenum(1,1,1) + (hrs/24); 
  n.Cs(end+1:end+nhrs,1:n.nlat,1:n.nlon) = nc{'currspeed'}(:,:,:);
  n.Ws(end+1:end+nhrs,1:n.nlat,1:n.nlon) = nc{'windspeed'}(:,:,:); 
  n.Hs(end+1:end+nhrs,1:n.nlat,1:n.nlon) = nc{'sigwavehgt'}(:,:,:); 
  n.Sw(end+1:end+nhrs,1:n.nlat,1:n.nlon) = nc{'sigswellhgt'}(:,:,:); 
  n.Lg(end+1:end+nhrs,1:n.nlat,1:n.nlon) = nc{'ardhuin_surface_drift'}(:,:,:); 
  n.Lgm(end+1:end+nhrs,1:n.nlat,1:n.nlon) = nc{'monismith_surface_drift'}(:,:,:); 
  close(nc); clear nc
end;

stn.lon = mean(g.lon);
stn.lat = mean(g.lat);
stn = read_hires_bathymetry(stn,geographic_radius_m(g.lon,g.lat),[],false);

bbox = field_bbox(g);
maxWs = nanmax(nanmax(n.Ws(:)),nanmax(g.Ws(:)));
maxCs = nanmax(nanmax(n.Cs(:)),nanmax(g.Cs(:)));
maxHs = nanmax(nanmax(n.Hs(:)),nanmax(g.Hs(:)));
maxLg = nanmax(nanmax(n.Lg(:)),nanmax(n.Lgm(:)));
stepLg = floorn(maxLg/10,-2);

if 0;
  % Ooops - no longer works: N and G may be on different grids!
  fmg; plot(n.lon,(n.lon-g.lon)*111e3*1e3); titlename('Longitude error??? [m]');
  fmg; plot(n.lat,(n.lat-g.lat)*111e3*1e3); titlename('Latitude error [m]');
end;

if 1;
  %% Time-Index to plot
  %tix = 1;
  %tix = 5;
  % % Date-time with highest (sum of Hs) waves
  % [ig,tix] = nanmax(nansum(n.Hs(:,:),2));
  % Date-time with greatest (sum of Lg) surface transport
  [ig,tix] = nanmax(nansum(n.Lg(:,:),2));

  if 0;
    fmg; contourf(g.lon,g.lat,squeeze(g.Hs(tix,:,:)),0:0.10:maxHs); colorbar; 
    plot_hires_coastline(stn.ngdc_hires_bathy,[],false);
    axis([bbox(:)',0,maxHs,0,maxHs]); daspect([1,cosd(n.lat(1)),1]);
    titlename(['GRiB2 Hs: ',datestr(g.t(tix))]);

    fmg; contourf(g.lon,g.lat,squeeze(g.Sw(tix,:,:)),0:0.10:maxHs); colorbar; 
    plot_hires_coastline(stn.ngdc_hires_bathy,[],false);
    axis([bbox(:)',0,maxHs,0,maxHs]); daspect([1,cosd(n.lat(1)),1]);
    titlename(['GRiB2 Sw: ',datestr(g.t(tix))]);
  end;
  if 1;
    fmg; contourf(n.lon,n.lat,squeeze(n.Ws(tix,:,:)),0:0.10:maxWs); colorbar; 
    plot_hires_coastline(stn.ngdc_hires_bathy,[],false); 
    axis([bbox(:)',0,maxWs,0,maxWs]); daspect([1,cosd(n.lat(1)),1]);
    titlename(['netCDF Ws: ',datestr(n.t(tix))]);

    fmg; contourf(n.lon,n.lat,squeeze(n.Cs(tix,:,:)),0:0.10:maxCs); colorbar; 
    plot_hires_coastline(stn.ngdc_hires_bathy,[],false); 
    axis([bbox(:)',0,maxCs,0,maxCs]); daspect([1,cosd(n.lat(1)),1]);
    titlename(['netCDF Cs: ',datestr(n.t(tix))]);

    fmg; contourf(n.lon,n.lat,squeeze(n.Sw(tix,:,:)),0:0.10:maxHs); colorbar; 
    plot_hires_coastline(stn.ngdc_hires_bathy,[],false); 
    axis([bbox(:)',0,maxHs,0,maxHs]); daspect([1,cosd(n.lat(1)),1]);
    titlename(['netCDF Sw: ',datestr(n.t(tix))]);

    fmg; contourf(n.lon,n.lat,squeeze(n.Lg(tix,:,:)),0:stepLg:maxLg); colorbar; 
    plot_hires_coastline(stn.ngdc_hires_bathy,[],false); 
    axis([bbox(:)',0,maxLg,0,maxLg]); daspect([1,cosd(n.lat(1)),1]);
    titlename(['netCDF Ardhuin Lg: ',datestr(n.t(tix))]);
  end;
end;

if 1;
  fmg; contourf(g.lon,g.lat,squeeze(nanmax(g.Hs)),0:0.10:maxHs); colorbar; 
  plot_hires_coastline(stn.ngdc_hires_bathy,[],false);
  axis([bbox(:)',0,maxHs,0,maxHs]); daspect([1,cosd(n.lat(1)),1]);
  titlename(['GRiB2 MAX Hs: ',datestr(n.t(1)),'-',datestr(n.t(end))]);

  fmg; contourf(n.lon,n.lat,squeeze(nanmax(n.Hs)),0:0.10:maxHs); colorbar; 
  plot_hires_coastline(stn.ngdc_hires_bathy,[],false); 
  axis([bbox(:)',0,maxHs,0,maxHs]); daspect([1,cosd(n.lat(1)),1]);
  titlename(['netCDF MAX Hs: ',datestr(n.t(1)),'-',datestr(n.t(end))]);

  fmg; contourf(n.lon,n.lat,squeeze(nanmax(n.Sw)),0:0.10:maxHs); colorbar; 
  plot_hires_coastline(stn.ngdc_hires_bathy,[],false); 
  axis([bbox(:)',0,maxHs,0,maxHs]); daspect([1,cosd(n.lat(1)),1]);
  titlename(['netCDF MAX Sw: ',datestr(n.t(1)),'-',datestr(n.t(end))]);

  fmg; contourf(n.lon,n.lat,squeeze(nanmax(n.Lg)),0:stepLg:maxLg); colorbar; 
  plot_hires_coastline(stn.ngdc_hires_bathy,[],false); 
  axis([bbox(:)',0,maxLg,0,maxLg]); daspect([1,cosd(n.lat(1)),1]);
  titlename(['netCDF MAX Ardhuin Lg: ',datestr(n.t(1)),'-',datestr(n.t(end))]);

  fmg; contourf(n.lon,n.lat,squeeze(nanmax(n.Lgm)),0:stepLg:maxLg); colorbar; 
  plot_hires_coastline(stn.ngdc_hires_bathy,[],false); 
  axis([bbox(:)',0,maxLg,0,maxLg]); daspect([1,cosd(n.lat(1)),1]);
  titlename(['netCDF MAX Monismith Lg: ',datestr(n.t(1)),'-',datestr(n.t(end))]);
end;
