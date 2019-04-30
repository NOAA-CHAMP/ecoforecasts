1;

ncfname = fullfile(get_thesis_path('../data/hycom/GOM2'),'archv.2014_001_12.nc');
nc = mDataset(ncfname);
nj_info(nc),
more on
nj_info(nc),
lon = cast(nc{'Longitude'}(:),'double');
lat = cast(nc{'Latitude'}(:),'double');
z = cast(nc{'Depth'}(:),'double');
z
t = cast(nc{'temperature'}(:,:,:,:),'double');
fmg; contourf(lon,lat,squeeze(t(1,1,:,:))); colorbar;
fmg; contourf(lon,lat,squeeze(t(1,:,:))); colorbar;
axis([-83,-79,24.5,26]);
axis([-83,-79,24.5,26]); daspect([1,cosd(25),1]);
axis([-83,-79,24,26]); daspect([1,cosd(25),1]);
fmg; contourf(lon,lat,squeeze(t(2,:,:))); colorbar; titlename('10 m T');
axis([-83,-79,24,26]); daspect([1,cosd(25),1]);
stn = get_station_from_station_name('mlrf1'); stn = load_all_ndbc_data(stn);
fmg; plot_ts(stn.ndbc_air_t,stn.ndbc_sea_t,stn.ndbc_wind1_speed); legend('Air','Sea','Wind');
xlim(datenum(2012,1,[1,8])); datetick3;
fmg; plot_ts(stn.ndbc_air_t,stn.ndbc_sea_t,stn.ndbc_wind1_speed); legend('Air','Sea','Wind');
xlim(datenum(2014,1,[1,8])); datetick3;
ylim([0,35]);



1;

%stn = get_station_from_station_name('mlrf1');
stn = get_station_from_station_name('lonf1');
stn = load_all_ndbc_data(stn);

%nc = mDataset('http://tds.hycom.org/thredds/dodsC/datasets/GOMu0.04/expt_50.1/data/netcdf/1993/hycom_gomu_501_1993010100_t000.nc');
nc = mDataset('http://tds.hycom.org/thredds/dodsC/datasets/GOMu0.04/expt_50.1/data/netcdf/2006/hycom_gomu_501_2006080100_t000.nc');
z = nc{'depth'}(:);
lon = cast(nc{'lon'}(:),'double');
lat = cast(nc{'lat'}(:),'double');
t = cast(nc{'water_temp'}(:),'double');
s = cast(nc{'salinity'}(:),'double');
u = cast(nc{'water_u'}(:),'double');
v = cast(nc{'water_v'}(:),'double');
close(nc); clear nc

nc = mDataset('c:/Users/lew.gramer/Documents/RSMAS/Coastal/thesis/data/hycom/FKEYS/303_archv.2006_213_00_3zt.nc');
flon = cast(nc{'Longitude'}(:),'double');
flat = cast(nc{'Latitude'}(:),'double');
ft = cast(nc{'temperature'}(:),'double');
close(nc); clear nc

nc = mDataset('c:/Users/lew.gramer/Documents/RSMAS/Coastal/thesis/data/hycom/FKEYS/303_archv.2006_213_00_3zu.nc');
flon = cast(nc{'Longitude'}(:),'double');
flat = cast(nc{'Latitude'}(:),'double');
fu = cast(nc{'u'}(:),'double');
close(nc); clear nc

nc = mDataset('c:/Users/lew.gramer/Documents/RSMAS/Coastal/thesis/data/hycom/FKEYS/303_archv.2006_213_00_3zv.nc');
flon = cast(nc{'Longitude'}(:),'double');
flat = cast(nc{'Latitude'}(:),'double');
fv = cast(nc{'v'}(:),'double');
close(nc); clear nc

d = 3;

stn = plot_hires_bathymetry(stn,-[0:5:80],[80e3,40e3],true);
axis(axis);

[ig,lonix] = min(abs(lon-stn.lon));
[ig,latix] = min(abs(lat-stn.lat));
quiver(lon(lonix-20:lonix+20),lat(latix-10:latix+10),squeeze(u(d,latix-10:latix+10,lonix-20:lonix+20)),squeeze(v(d,latix-10:latix+10,lonix-20:lonix+20)),'r');
qh=quiver(flon,flat,squeeze(fu(d,:,:)),squeeze(fv(d,:,:)),4,'k');
axis([-81.242701853251376,-80.516937290426270,24.692041274175342,25.239836350899584]);
titlename(['Currents depth level ',num2str(d)]);
