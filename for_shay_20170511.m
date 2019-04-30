1;

nc = mDataset('D:\python\CURIE\icon\tmp\waves\dl-archive\key_nwps_CG2_20161201_0600.grib2');
t = cast(nc{'time'}(:),'double'); dts = datenum(2016,12,1,6,0,0) + (t/24);
lon = cast(nc{'lon'}(:,:),'double'); lon = lon-360;
lat = cast(nc{'lat'}(:,:),'double');
Hs = cast(nc{'Significant_height_of_combined_wind_waves_and_swell'}(:,:,:),'double');
Tp = cast(nc{'Primary_wave_mean_period'}(:,:,:),'double');
close(nc); clear nc

fmg; contourf(lon,lat,squeeze(nanmean(Hs))); colorbar; titlename([datestr(dts(1)),' to ',datestr(dts(end))]); daspect([1,cosd(25),1]);
bath = read_hires_bathymetry_for_field({lon,lat},false);
plot_hires_coastline(bath.ngdc_hires_bathy);
print('-dpng','figs/for_shay_Hs_2016_12_01_2016_12_05.png');
