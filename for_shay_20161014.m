1;
%% SCRIPT to test Python 'process_nwps' netCDF result file(s)

if ( ~exist('datapath','var') )
  datapath = get_ecoforecasts_path('data');
end;
if ( ~exist('dataset','var') )
  %dataset = 'mfl_nwps_CG1';
  if ( ~exist('region','var') )
    region = 'mfl';
  end;
  switch ( region ),
   case 'mfl',
    stage = 'CG1';
   case 'key',
    stage = 'CG2';
   otherwise,
    error('Unknown region %s',region);
  end;
  dataset = [region,'_nwps_',stage];
end;

% if ( ~exist('filedate','var') )
%   filedate = datenum(2016,07,31,6,0,0);
% end;
if ( ~exist('filedates','var') )
  %filedates = datenum(2016,10,01,6,0,0):0.5:datenum(2016,10,13,18,0,0);
  filedates = datenum(2016,10,4,6,0,0);
end;


% GRiB2 time indices to extract
tix = 1:4;

for dtix = 1:numel(filedates)
  filedate = filedates(dtix);
  [y,m,d,H,M,S] = datevec(filedate);
  basefile = sprintf('%s_%04d%02d%02d_%02d00',dataset,y,m,d,H);
  clear y m d H M S;

  g = []; clear g
  g.basefile = basefile;
  g.datafile = fullfile(datapath,[g.basefile,'.grib2']);
  if ( ~exist(g.datafile,'file') )
    disp(['SKIPPING ',g.datafile]);
    continue;
  end;
  nc = mDataset(g.datafile); 
  g.t = datenum(2016,7,31,6,0,0) + (cast(nc{'time'}(tix),'double')/24); 
  g.lat = cast(nc{'lat'}(:),'double'); 
  g.lon = cast(nc{'lon'}(:),'double'); 
  % Convert to West/East longitudes
  g.lon(g.lon >= 180) = g.lon(g.lon >= 180) - 360;
  g.Hs = cast(nc{'Significant_height_of_combined_wind_waves_and_swell'}(tix,:,:),'double'); 
  g.Cs = cast(nc{'Current_speed'}(tix,:,:),'double'); 
  g.Cd = cast(nc{'Current_direction'}(tix,:,:),'double'); 
  [g.Cu,g.Cv] = spddir_to_uv_curr(g.Cs,g.Cd);
  close(nc); clear nc

  n = []; clear n
  n.basefile = basefile;
  n.datafile = fullfile(datapath,[n.basefile,'.nc']);
  if ( exist(n.datafile,'file') )
    doN = true;
  else
    disp(['No netCDF ',n.datafile]);
    doN = false;
  end;
  if ( doN )
    nc = mDataset(n.datafile);
    n.t = datenum(1,1,1) + (nc{'time'}(:)/24); 
    n.lat = nc{'lat'}(:); 
    n.lon = nc{'lon'}(:); 
    n.Hs = nc{'sigwavehgt'}(:,:,:); 
    n.Cs = cast(nc{'currspeed'}(:,:,:),'double'); 
    n.Cd = cast(nc{'currdir'}(:,:,:),'double'); 
    [n.Cu,n.Cv] = spddir_to_uv_curr(n.Cs,n.Cd);
    n.Lg = nc{'gramer_surface_drift'}(:,:,:); 
    n.Ld = nc{'gramer_surface_drift_dir'}(:,:,:); 
    [n.Lu,n.Lv] = spddir_to_uv_curr(n.Lg,n.Ld);
    close(nc); clear nc
  end;

  if ( ~exist('stn','var') || ~isfield(stn,'lon') )
    stn = []; clear stn;
    stn.lon = mean(g.lon);
    stn.lat = mean(g.lat);
    if 0;
      stn = read_hires_bathymetry(stn,geographic_radius_m(g.lon,g.lat),[],false);
    else
      % This relies on the hack that we have a coastline database for South Florida!
      stn.ngdc_hires_bathy.lon = g.lon;
      stn.ngdc_hires_bathy.lat = g.lat;
    end;
  end;

  bbox = field_bbox(g);
  if ( doN )
    maxHs = nanmax(nanmax(n.Hs(:)),nanmax(g.Hs(:)));
    maxCs = nanmax(nanmax(n.Cs(:)),nanmax(g.Cs(:)));
    maxCu = nanmax(nanmax(n.Cu(:)),nanmax(g.Cu(:)));
    maxCv = nanmax(nanmax(n.Cv(:)),nanmax(g.Cv(:)));
    maxLg = nanmax(n.Lg(:));
    maxLu = nanmax(n.Lu(:));
    maxLv = nanmax(n.Lv(:));
  else
    maxHs = nanmax(g.Hs(:));
    maxCs = nanmax(g.Cs(:));
    maxCu = nanmax(g.Cu(:));
    maxCv = nanmax(g.Cv(:));
  end;

  if ( isfield(stn.ngdc_hires_bathy,'field') )
    plotcoast = @(x)(plot_hires_coastline(stn.ngdc_hires_bathy,[],false));
  else
    % This relies on the hack that we have a coastline database for South Florida!
    plotcoast = @(x)(plot_hires_coastline(stn.ngdc_hires_bathy));
  end;

  fmg; contourf(g.lon,g.lat,squeeze(g.Hs(1,:,:)),0:0.05:3); colorbar; 
  plotcoast;
  axis([bbox(:)',0,maxCs,0,maxCs]); daspect([1,cosd(g.lat(1)),1]);
  quiver(g.lon,g.lat,squeeze(g.Cu(1,:,:)),squeeze(g.Cv(1,:,:)));
  titlename(['GRiB Hs vs. Cu/v: ',datestr(g.t(1))]);
  print('-dpng',fullfile(efs_figspath(),sprintf('grib_cs_%s.png',datestr(g.t(1),'yyyymmddHHMM'))));

  if ( doN )
    fmg; contourf(n.lon,n.lat,squeeze(n.Hs(1,:,:)),0:0.05:3); colorbar; 
    plotcoast;
    axis([bbox(:)',0,maxCs,0,maxCs]); daspect([1,cosd(n.lat(1)),1]);
    quiver(n.lon,n.lat,squeeze(n.Lu(1,:,:)),squeeze(n.Lv(1,:,:)));
    titlename(['netCDF Hs vs. L_gu/v: ',datestr(n.t(1))]);
    print('-dpng',fullfile(efs_figspath(),sprintf('nc_cs_%s.png',datestr(n.t(1),'yyyymmddHHMM'))));
  end;

if 0;
  if ( doN )
    fmg; plot(n.lon,(n.lon-g.lon)*111e3*1e3); titlename('Longitude error??? [m]');
    fmg; plot(n.lat,(n.lat-g.lat)*111e3*1e3); titlename('Latitude error [m]');
  end;

  fmg; contourf(g.lon,g.lat,squeeze(g.Hs(1,:,:)),0:0.05:3); colorbar; 
  plotcoast;
  axis([bbox(:)',0,maxHs,0,maxHs]); daspect([1,cosd(g.lat(1)),1]);
  titlename(['GRiB2 Hs: ',datestr(g.t(1))]);
  print('-dpng',fullfile(efs_figspath(),sprintf('grib_hs_%s.png',datestr(g.t(1),'yyyymmddHHMM'))));

  if ( doN )
    fmg; contourf(n.lon,n.lat,squeeze(n.Hs(1,:,:)),0:0.05:3); colorbar; 
    plotcoast;
    axis([bbox(:)',0,maxHs,0,maxHs]); daspect([1,cosd(n.lat(1)),1]);
    titlename(['netCDF Hs: ',datestr(n.t(1))]);
    print('-dpng',fullfile(efs_figspath(),sprintf('nc_hs_%s.png',datestr(g.t(1),'yyyymmddHHMM'))));

    fmg; contourf(n.lon,n.lat,squeeze(n.Cs(1,:,:)),0:0.0025:1.0000); colorbar; 
    plotcoast;
    axis([bbox(:)',0,maxCs,0,maxCs]); daspect([1,cosd(n.lat(1)),1]);
    titlename(['netCDF Cs: ',datestr(n.t(1))]);
    print('-dpng',fullfile(efs_figspath(),sprintf('nc_cs_%s.png',datestr(g.t(1),'yyyymmddHHMM'))));

    fmg; contourf(n.lon,n.lat,squeeze(nanmax(n.Lg)),0:0.0025:0.0500); colorbar; 
    plotcoast;
    axis([bbox(:)',0,maxLg,0,maxLg]); daspect([1,cosd(n.lat(1)),1]);
    titlename(['netCDF MAX Lg: ',datestr(n.t(1))]);
  end;

  fmg; contourf(g.lon,g.lat,squeeze(g.Cs(1,:,:)),0:0.0025:1.0000); colorbar; 
  plotcoast;
  axis([bbox(:)',0,maxCs,0,maxCs]); daspect([1,cosd(g.lat(1)),1]);
  titlename(['netCDF Cs: ',datestr(g.t(1))]);
  print('-dpng',fullfile(efs_figspath(),sprintf('grib_cs_%s.png',datestr(g.t(1),'yyyymmddHHMM'))));
end;

end;
