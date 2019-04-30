1;
%% SCRIPT INTERP_DEM_BATHY.m
% Interleaves a Digitel Elevation Model (DEM, heights above sea level) with a
% bathymetry dataset, both with fields lon,lat,dat,url: both are loaded from
% .MAT files if they exist, or calls MDATASET (v.) to extract and save them.
%
% Last Saved Time-stamp: <Sun 2016-05-22 12:38:08 Eastern Daylight Time gramer>

doSave = false;

% Target resolution in [m]

%target_res = 20;
%% One arcsecond of latitude
%target_res = 30;
%target_res = 45;
%% Three arcseconds of latitude
%target_res = 90;


%region = 'guam';

%region = 'saipan';
%region = 'tinian';
%region = 'rota';

%region = 'tutuila';
%region = 'ofu';
region = 'tau';


switch ( region ),
  %{
  % hmrg_bathytopo_50m_mhi
 case 'hawaii',
  demnm = 'usgs_dem_10m_bigisland';
  batnm = '';
  target_res = 30;
 case 'maui',
  demnm = 'usgs_dem_10m_maui';
  batnm = '';
  target_res = 30;
  %}

 case 'guam',
  demnm = 'usgs_dem_10m_guam';
  batnm = 'pibhmc_bathy_5m_guam';
  target_res = 30;

 case 'saipan',
  demnm = 'usgs_dem_10m_saipan';
  batnm = 'pibhmc_bathy_5m_saipan';
  target_res = 30;
 case 'tinian',
  demnm = 'usgs_dem_10m_tinian';
  % Memory limits forced special handling: LAT(1:END), LON(500:END-500), DAT(1:END,500:END-500)
  batnm = 'pibhmc_bathy_5m_tinian.PARTIAL';
  target_res = 30;
 case 'rota',
  demnm = 'usgs_dem_10m_rota';
  batnm = 'pibhmc_bathy_5m_rota';
  target_res = 30;

 case 'tutuila',
  demnm = 'usgs_dem_10m_tutuila';
  batnm = 'pibhmc_bathy_5m_tutuila';
  target_res = 30;
 case 'ofu',
  demnm = 'usgs_dem_10m_ofuolosega';
  batnm = 'pibhmc_bathy_5m_ofuolosega';
  target_res = 30;
 case 'tau',
  demnm = 'usgs_dem_10m_tau';
  batnm = 'pibhmc_bathy_5m_tau';
  target_res = 30;

end;


cmbpath = fullfile(get_ecoforecasts_path('coast'),[batnm,'.dods.',num2str(target_res),'m.',demnm,'.dods.mat']);

if ( exist(cmbpath,'file') )
  disp(['Loading ',cmbpath]);
  load(cmbpath);
  
else

  dempath = fullfile(get_ecoforecasts_path('coast'),[demnm,'.dods.mat']);
  batpath = fullfile(get_ecoforecasts_path('coast'),[batnm,'.dods.mat']);

  if ( ~exist(dempath,'file') )
    disp(['Extracting ',dempath]);
    url = ['http://oos.soest.hawaii.edu/thredds/dodsC/',demnm];
    nc = mDataset(url);
    dat = cast(nc{'elev'}(:,:),'double');
    lon = cast(nc{'lon'}(:),'double');
    lat = cast(nc{'lat'}(:),'double');
    close(nc); clear nc
    disp(['Saving ',dempath]);
    save(dempath,'url','lon','lat','dat','-v7.3');  
    dat=[]; lat=[]; lon=[]; clear dat lat lon url
  end;
  if ( ~exist(batpath,'file') )
    disp(['Extracting ',batpath]);
    url = ['http://oos.soest.hawaii.edu/thredds/dodsC/',batnm];
    nc = mDataset(url);
    dat = cast(nc{'elev'}(:,:),'double');
    lon = cast(nc{'lon'}(:),'double');
    lat = cast(nc{'lat'}(:),'double');
    close(nc); clear nc
    disp(['Saving ',batpath]);
    save(batpath,'url','lon','lat','dat','-v7.3');  
    dat=[]; lat=[]; lon=[]; clear dat lat lon url
  end;

  dem = load(dempath);
  dem.dat(dem.dat<-eps) = nan;

  bat = load(batpath);
  bat.dat(bat.dat>eps) = nan;

  url = bat.url;
  demurl = dem.url;


  disp(['Interpolating ',batnm,' into ',demnm]);

  % Up- (or down-)sample elevation model to bathymetry resolution and spatial coverage
  rawlon = bat.lon;
  rawlat = bat.lat;
  [RAWLON,RAWLAT]=meshgrid(rawlon,rawlat);
  rawdat = interp2(dem.lon,dem.lat,dem.dat,RAWLON,RAWLAT,'linear',nan);
  RAWLON=[]; RAWLAT=[]; clear RAWLON RAWLAT
  rawdat(rawdat<-eps) = nan;
  rawdat(isnan(rawdat)) = bat.dat(isnan(rawdat));

  bat=[]; clear bat batpath
  dem=[]; clear dem dempath


  disp(['Smoothing to ',num2str(target_res),' m resolution']);

tic, % THIS LOOP IS SLOW!
  % Downsample (using NANEAN) to target horizontal resolution
  dx=1e3*median(distance_wgs84(mean(rawlat(:)),rawlon(1:end-1),mean(rawlat(:)),rawlon(2:end)));
  dy=1e3*median(distance_wgs84(mean(rawlat(:)),rawlon(1:end-1),mean(rawlat(:)),rawlon(2:end)));
  nx = floor(target_res/dx/2);
  ny = floor(target_res/dy/2);
  xixen = nx+1:(nx*2)+1:length(rawlon)-nx;
  yixen = ny+1:(ny*2)+1:length(rawlat)-ny;
  lon = rawlon(xixen);
  lat = rawlat(yixen);
  for xixix=1:numel(xixen);
    for yixix=1:numel(yixen);
      xix = xixen(xixix);
      yix = yixen(yixix);
      datdat = rawdat(yix-ny:yix+ny,xix-nx:xix+nx);
      dat(yixix,xixix) = nanmean(datdat(:));
      datdat=[]; clear datdat
    end;
  end;
  rawdat = [];
  clear dx dy nx ny rawdat rawlat rawlon xix xixen xixix yix yixen yixix
toc,

  % fmg; contourf(lon,lat,dat,-1000:50:100); colorbar;
  % contour(lon,lat,dat,[0,0],'Color','k','LineWidth',2);
  % titlename(['Bathy/DEM at ',num2str(target_res),' m Sans Spline']);

  % Attempt to spline-file all holes at once - probably won't work
  disp(['Spline-fit NaN-filled holes']);
  s = warning('OFF','MATLAB:interp2:NaNstrip');
  nanix = find(isnan(dat));
  [LON,LAT] = meshgrid(lon,lat);
  try,
    dat(nanix) = interp2(lon,lat,dat,LON(nanix),LAT(nanix),'spline',nan)
  catch,
    catchwarn;
    doSave = false;
  end;
  warning(s); clear s
  LON=[]; LAT=[]; clear LON LAT ans nanix

  if ( doSave )
    disp(['Saving ',cmbpath]);
    save(cmbpath,'url','demurl','lon','lat','dat','-v7.3');
  else
    disp(['NOT saving ',cmbpath]);
  end;
end;

% fmg; contourf(lon,lat,dat,-40:2:4); colorbar; 
% contour(lon,lat,dat,[0,0],'Color','k','LineWidth',2);
% daspect([1,cosd(lat(1)),1]); set_surf_cursor;
% titlename(['Bathy/DEM at ',num2str(target_res),' m']);

% fmg; contourf(lon,lat,dat,-1000:50:100); colorbar;
% contour(lon,lat,dat,[0,0],'Color','k','LineWidth',2);
% daspect([1,cosd(lat(1)),1]); set_surf_cursor;
% axis([145.74,145.8,15.12,15.17]);
% titlename(['Bathy/DEM at ',num2str(target_res),' m (Laolao Bay)']);

fmg; contourf(lon,lat,dat,-1000:50:100); colorbar;
contour(lon,lat,dat,[0,0],'Color','k','LineWidth',2);
daspect([1,cosd(lat(1)),1]); set_surf_cursor;
titlename(['Bathy/DEM at ',num2str(target_res),' m wider bathy']);

