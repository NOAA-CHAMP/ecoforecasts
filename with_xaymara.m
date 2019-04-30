1;

if ( ~exist('doPrint','var') )
  doPrint = false;
end;
if ( doPrint )
  if ( ~exist('figspath','var') )
    %figspath = get_ecoforecasts_path('figs');
    figspath = get_postdoc_path;
  end;
end;

if ( ~exist('hacf1','var') || ~isfield(hacf1,'daily_sea_t') )
  %% NOTE: The script called here currently resides in $MYDATA/Postdoc/
  extract_kuffner_hudson;
end;

if 0;
find_date_ranges(ccrf1.daily_sea_t.date,100)
find_date_ranges(cmdf1.daily_sea_t.date,100)
find_date_ranges(cryf1.daily_sea_t.date,100)
find_date_ranges(hacf1.daily_sea_t.date,100)
find_date_ranges(snaf1.daily_sea_t.date,100)
end;

fmg; plot_ts(hacf1.daily_sea_t);


sanf1 = get_station_from_station_name('sanf1'); sanf1 = load_all_ndbc_data(sanf1);
sanf1 = station_optimal_isobath_orientation(sanf1);

smkf1 = get_station_from_station_name('smkf1'); smkf1 = load_all_ndbc_data(smkf1);
smkf1 = station_optimal_isobath_orientation(smkf1);

lonf1 = get_station_from_station_name('lonf1'); lonf1 = load_all_ndbc_data(lonf1);
lonf1 = station_optimal_isobath_orientation(lonf1);


mlrf1 = get_station_from_station_name('mlrf1'); mlrf1 = load_all_ndbc_data(mlrf1);
mlrf1 = station_optimal_isobath_orientation(mlrf1);
mlrf1 = verify_variable(mlrf1,'ndbc_sea_t_1_d_avg');
mlrf1 = verify_variable(mlrf1,'ndbc_sea_t_1_d_decimate');

%fmg; plot_ts(hacf1.daily_sea_t,mlrf1.ndbc_sea_t_1_d_decimate);

fmg; plot_ts(hacf1.daily_sea_t,mlrf1.ndbc_sea_t_1_d_decimate); legend('Hens & Chickens (4 m)','Molasses (2 m)'); titlename('Daily sea temperature [^oC]');
if doPrint; print('-dpng',fullfile(figspath,'mlrf1-vs-hacf1-seatemp-timeseries.png')); end;

scatter_fit_ts_seasons(hacf1.daily_sea_t,mlrf1.ndbc_sea_t_1_d_avg,[],[],'Hens & Chickens (4 m)','Molasses (2 m)',[],[],true); titlename('Daily sea temperature [^oC]');
if doPrint; print('-dpng',fullfile(figspath,'mlrf1-scatter-hacf1-seatemp.png')); end;

% plot_beta_vs_t_season


%mlrf1 = plot_hires_bathymetry(mlrf1,-[0:5:80],[12e3,12e3],true,[],false);
hacf1 = plot_hires_bathymetry(hacf1,-[0:2:40],[25e3,15e3],true,@contour,false);
plot(hacf1.lon,hacf1.lat,'ks','MarkerFaceColor','r'); th = text(hacf1.lon,hacf1.lat,'H & C','Color','r');
plot(sanf1.lon,sanf1.lat,'ks','MarkerFaceColor','r'); th = text(sanf1.lon,sanf1.lat,'SANF1','Color','r');
plot(smkf1.lon,smkf1.lat,'ks','MarkerFaceColor','r'); th = text(smkf1.lon,smkf1.lat,'SMKF1','Color','r');
plot(lonf1.lon,lonf1.lat,'ks','MarkerFaceColor','r'); th = text(lonf1.lon,lonf1.lat,'LONF1','Color','r');
plot(mlrf1.lon,mlrf1.lat,'ks','MarkerFaceColor','r'); th = text(mlrf1.lon,mlrf1.lat,'MLRF1','Color','r');
if doPrint; print('-dpng',fullfile(figspath,'hacf1-bathy.png')); end;

hacf1 = plot_hires_bathymetry(hacf1,-[0:0.5:7],[2e3,2e3],true,@contour,true);
if doPrint; print('-dpng',fullfile(figspath,'hacf1-bathy-blowup.png')); end;

ccrf1 = plot_hires_bathymetry(ccrf1,-[0:2:40],[25e3,15e3],true,@contour,false);
plot(ccrf1.lon,ccrf1.lat,'ks','MarkerFaceColor','r'); th = text(ccrf1.lon,ccrf1.lat,'Crocker','Color','r');
plot(sanf1.lon,sanf1.lat,'ks','MarkerFaceColor','r'); th = text(sanf1.lon,sanf1.lat,'SANF1','Color','r');
plot(smkf1.lon,smkf1.lat,'ks','MarkerFaceColor','r'); th = text(smkf1.lon,smkf1.lat,'SMKF1','Color','r');
plot(lonf1.lon,lonf1.lat,'ks','MarkerFaceColor','r'); th = text(lonf1.lon,lonf1.lat,'LONF1','Color','r');
plot(mlrf1.lon,mlrf1.lat,'ks','MarkerFaceColor','r'); th = text(mlrf1.lon,mlrf1.lat,'MLRF1','Color','r');
if doPrint; print('-dpng',fullfile(figspath,'ccrf1-bathy.png')); end;

ccrf1 = plot_hires_bathymetry(ccrf1,-[0:0.5:7],[2e3,2e3],true,@contour,true);
if doPrint; print('-dpng',fullfile(figspath,'ccrf1-bathy-blowup.png')); end;
