1;

% 06-14 01:36 station_boxplots.m
% 06-14 01:36 scatter_multifit.m
% 06-14 01:36 scatter_fit.m
% 06-14 01:36 rstool_station.m
% 06-14 01:36 get_doy.m
% 06-14 01:35 station_subset_anova.m
% 06-14 01:35 station_plot_fields_x.m
% 06-14 01:35 station_anova_multcompare.m
% 06-14 01:35 get_daylight.m
% 06-14 01:35 get_all_station_metadata.m
% 06-14 01:35 fix_ndbc.m
% 06-14 01:35 calc_kd_daily.m
% 06-14 01:35 ancold.m
% 05-31 23:16 multiplot_station.m
% 05-31 23:16 multiplot_datetick.m
% 05-27 13:44 station_plot_fields.m

% stnms = {'fwyf1','mlrf1','lonf1'};
stnms = {'mlrf1'};

for ix=1:length(stnms)
  stnm = stnms{ix};
  stn = load_all_ndbc_data([],stnm);
  station_subset_anova(stn,'ndbc_sea_t',@(x)(find(get_month(x.date)==1)),'January');
  station_subset_anova(stn,'ndbc_sea_t',@(x)(find(get_month(x.date)==2)),'February');
  station_subset_anova(stn,'ndbc_sea_t',@(x)(find(get_month(x.date)==3)),'March');
  station_subset_anova(stn,'ndbc_sea_t',@(x)(find(get_month(x.date)==4)),'April');
  station_subset_anova(stn,'ndbc_sea_t',@(x)(find(get_month(x.date)==5)),'May');
  station_subset_anova(stn,'ndbc_sea_t',@(x)(find(get_month(x.date)==6)),'June');
%   for wk = 1:15
%     wknm = sprintf('Week #%d',wk);
%     station_subset_anova(stn,'ndbc_sea_t',@(x)(find(get_week(x.date)==wk)),wknm);
%   end;
  stn = []; clear stn;
end;
