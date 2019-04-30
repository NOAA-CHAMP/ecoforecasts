1;
%% SCRIPT to test Python 'process_nwps' netCDF result file(s)

if ( ~exist('datapath','var') )
  datapath = get_ecoforecasts_path('data');
end;
if ( ~exist('dataset','var') )
  %dataset = 'key_nwps_CG2';
  dataset = 'mfl_nwps_CG1';
end;
if ( ~exist('filedate','var') )
  filedate = datenum(2016,07,31,6,0,0);
end;
if ( ~exist('basefile','var') )
  [y,m,d,H,M,S] = datevec(filedate);
  basefile = sprintf('%s_%04d%02d%02d_%02d00',dataset,y,m,d,H);
  clear y m d H M S;
end;

if ( ~exist('g','var') || ~isfield(g,'Hs') || ~isfield(g,'basefile') || ~strcmp(g.basefile,basefile) )
  g = []; clear g
  g.basefile = basefile;
  g.datafile = fullfile(datapath,[g.basefile,'.grib2']);
  
  nc = mDataset(g.datafile); 
  g.t = datenum(2016,7,31,6,0,0) + (cast(nc{'time'}(:),'double')/24); 
  g.lat = cast(nc{'lat'}(:),'double'); 
  g.lon = cast(nc{'lon'}(:),'double'); 
  % Convert to West/East longitudes
  g.lon(g.lon >= 180) = g.lon(g.lon >= 180) - 360;
  g.Hs = cast(nc{'Significant_height_of_combined_wind_waves_and_swell'}(:,:,:),'double'); 
  g.Sw = cast(nc{'Significant_height_of_swell_waves'}(:,:,:),'double'); 
  close(nc); clear nc
end;

if ( ~exist('n','var') || ~isfield(n,'Hs') || ~isfield(n,'basefile') || ~strcmp(n.basefile,basefile) )
  n = []; clear n
  n.basefile = basefile;
  n.datafile = fullfile(datapath,[n.basefile,'.nc']);
  
  nc = mDataset(n.datafile);
  n.t = datenum(1,1,1) + (nc{'time'}(:)/24); 
  n.lat = nc{'lat'}(:); 
  n.lon = nc{'lon'}(:); 
  n.Hs = nc{'sigwavehgt'}(:,:,:); 
  n.Sw = nc{'sigswellhgt'}(:,:,:); 
  n.Lg = nc{'ardhuin_surface_drift'}(:,:,:); 
  n.Lgm = nc{'monismith_surface_drift'}(:,:,:); 
  close(nc); clear nc
end;

stn.lon = mean(g.lon);
stn.lat = mean(g.lat);
stn = read_hires_bathymetry(stn,geographic_radius_m(g.lon,g.lat),[],false);

bbox = field_bbox(g);
maxHs = nanmax(nanmax(n.Hs(:)),nanmax(g.Hs(:)));
maxLg = nanmax(nanmax(n.Lg(:)),nanmax(n.Lgm(:)));

%fmg; plot(n.lon,(n.lon-g.lon)*111e3*1e3); titlename('Longitude error??? [m]');
%fmg; plot(n.lat,(n.lat-g.lat)*111e3*1e3); titlename('Latitude error [m]');

fmg; contourf(g.lon,g.lat,squeeze(g.Hs(1,:,:)),0:0.10:3); colorbar; 
plot_hires_coastline(stn.ngdc_hires_bathy,[],false);
axis([bbox(:)',0,maxHs,0,maxHs]); daspect([1,cosd(n.lat(1)),1]);
titlename(['GRiB2 Hs: ',datestr(g.t(1))]);

fmg; contourf(n.lon,n.lat,squeeze(n.Hs(1,:,:)),0:0.10:3); colorbar; 
plot_hires_coastline(stn.ngdc_hires_bathy,[],false); 
axis([bbox(:)',0,maxHs,0,maxHs]); daspect([1,cosd(n.lat(1)),1]);
titlename(['netCDF Hs: ',datestr(n.t(1))]);

fmg; contourf(n.lon,n.lat,squeeze(n.Sw(1,:,:)),0:0.10:3); colorbar; 
plot_hires_coastline(stn.ngdc_hires_bathy,[],false); 
axis([bbox(:)',0,maxHs,0,maxHs]); daspect([1,cosd(n.lat(1)),1]);
titlename(['netCDF Sw: ',datestr(n.t(1))]);

fmg; contourf(n.lon,n.lat,squeeze(nanmax(n.Lg)),0:0.0050:0.0500); colorbar; 
plot_hires_coastline(stn.ngdc_hires_bathy,[],false); 
axis([bbox(:)',0,maxLg,0,maxLg]); daspect([1,cosd(n.lat(1)),1]);
%titlename(['netCDF Lg: ',datestr(n.t(1))]);
titlename(['netCDF MAX Ardhuin Lg: ',datestr(n.t(1)),'-',datestr(n.t(end),'HH:MM:SS')]);

fmg; contourf(n.lon,n.lat,squeeze(nanmax(n.Lgm)),0:0.0050:0.0500); colorbar; 
plot_hires_coastline(stn.ngdc_hires_bathy,[],false); 
axis([bbox(:)',0,maxLg,0,maxLg]); daspect([1,cosd(n.lat(1)),1]);
%titlename(['netCDF Lgm: ',datestr(n.t(1))]);
titlename(['netCDF MAX Monismith Lg: ',datestr(n.t(1)),'-',datestr(n.t(end),'HH:MM:SS')]);
