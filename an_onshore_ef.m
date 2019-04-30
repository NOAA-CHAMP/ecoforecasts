1;

mlrf1=read_station_ef_csv('mlrf1','onshore-flux');

fmg; boxplot_ts(mlrf1.onshore_flux.sri); titlename('MLRF1 Onshore Flux S/RI');
print('-djpeg',fullfile(get_ecoforecasts_path('figs'),[lower(mlrf1.station_name),'-monthly-boxplot-sri-onshore-flux.jpg']));

fmg;
cm=cool(numel(unique(roundn(mlrf1.onshore_flux.sri.data,1))));
gscatter(get_year(mlrf1.onshore_flux.sri.date),get_jday(mlrf1.onshore_flux.sri.date),roundn(mlrf1.onshore_flux.sri.data,1),cm); datetick3('y',3);
print('-djpeg',fullfile(get_ecoforecasts_path('figs'),[lower(mlrf1.station_name),'-year-jday-gscatter-sri-onshore-flux.jpg']));
