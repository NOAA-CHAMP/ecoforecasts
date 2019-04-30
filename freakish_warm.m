1;

figspath = get_ecoforecasts_path('figs');

if 0;
  lkwf1=[]; pvgf1=[]; fwyf1=[]; vakf1=[]; mlrf1=[];
  clear lkwf1 pvgf1 fwyf1 vakf1 mlrf1
end;

lkwf1 = get_station_from_station_name('lkwf1'); lkwf1 = load_all_ndbc_data(lkwf1);
pvgf1 = get_station_from_station_name('pvgf1'); pvgf1 = load_all_ndbc_data(pvgf1);
fwyf1 = get_station_from_station_name('fwyf1'); fwyf1 = load_all_ndbc_data(fwyf1);
vakf1 = get_station_from_station_name('vakf1'); fwyf1 = load_all_ndbc_data(fwyf1);
mlrf1 = get_station_from_station_name('mlrf1'); mlrf1 = load_all_ndbc_data(mlrf1);

lkwf1 = verify_variable(lkwf1,'ndbc_sea_t_1_d_avg'); lkwf1 = verify_variable(lkwf1,'ndbc_air_t_1_d_avg'); 
fwyf1 = verify_variable(fwyf1,'ndbc_sea_t_1_d_avg'); fwyf1 = verify_variable(fwyf1,'ndbc_air_t_1_d_avg');
vakf1 = verify_variable(vakf1,'ndbc_sea_t_1_d_avg'); vakf1 = verify_variable(vakf1,'ndbc_air_t_1_d_avg');
mlrf1 = verify_variable(mlrf1,'ndbc_sea_t_1_d_avg'); mlrf1 = verify_variable(mlrf1,'ndbc_air_t_1_d_avg');

anomts=[]; climts=[]; hiasdts=[]; hipctts=[]; loasdts=[]; lopctts=[];
a_anomts=[]; a_climts=[]; a_hiasdts=[]; a_hipctts=[]; a_loasdts=[]; a_lopctts=[];
clear a_anomts a_climts a_hiasdts a_hipctts a_loasdts a_lopctts anomts ans asd clim climts cumfun hiasd hiasdts hipct hipctts loasd loasdts lopct lopctts per tid

%[anomts,climts,loasdts,hiasdts,lopctts,hipctts,cumfun,per] = anomalize_ts_to_ts(fwyf1.ndbc_sea_t,@get_jhour_no_leap);
%[anomts,climts,loasdts,hiasdts,lopctts,hipctts,cumfun,per] = anomalize_ts_to_ts(fwyf1.ndbc_sea_t,1,@get_jhour_no_leap);
[anomts,climts,loasdts,hiasdts,lopctts,hipctts,cumfun,per] = anomalize_ts_to_ts(fwyf1.ndbc_sea_t_1_d_avg,3);
%[anomts,climts,loasdts,hiasdts,lopctts,hipctts,cumfun,per] = anomalize_ts_to_ts(vakf1.ndbc_sea_t_1_d_avg,3);

[a_anomts,a_climts,a_loasdts,a_hiasdts,a_lopctts,a_hipctts,a_cumfun,a_per] = anomalize_ts_to_ts(fwyf1.ndbc_air_t_1_d_avg,3);

%% SEA TEMP COMPARISON
fmg;
plot_ts(lkwf1.ndbc_sea_t_1_d_avg,'.',...
        fwyf1.ndbc_sea_t_1_d_avg,'.',...
        vakf1.ndbc_sea_t_1_d_avg,'.',...
        mlrf1.ndbc_sea_t_1_d_avg,'.',...
        climts,'k','LineWidth',1.5,...
        lopctts,'--','Color',[.5,.5,.5],'LineWidth',1.5,...
        hipctts,'--','Color',[.5,.5,.5],'LineWidth',1.5);
legend('Lake Worth pier (LKWF1)','Off Key Biscayne (FWYF1)',...
       'Virginia Key (VAKF1)','Off Upper Keys (MLRF1)','FWYF1 Climatology \pm 47%');
%xlim([datenum(2016,12,1),now]); datetick3;
xlim([datenum(2016,6,1),now]); datetick3;
print('-dpng',fullfile(figspath,'South_Florida_water_temps_Jun 2016-Jan_2017.png'));
%print('-dpng',fullfile(figspath,'South_Florida_water_temps_Jun 2016-Jan_2017_VAKF1_clim.png'));


%% AIR-SEA TEMP COMPARISON
fmg;
plot_ts(fwyf1.ndbc_air_t_1_d_avg,'.',...
        fwyf1.ndbc_sea_t_1_d_avg,'.',...
        vakf1.ndbc_sea_t_1_d_avg,'.',...
        a_climts,'k','LineWidth',1.5,...
        a_lopctts,'--','Color',[.5,.5,.5],'LineWidth',1.5,...
        a_hipctts,'--','Color',[.5,.5,.5],'LineWidth',1.5);
legend('Air temp. off Key Biscayne (FWYF1)','Sea temp. off Key Biscayne (FWYF1)',...
       'Virginia Key sea temp. (VAKF1)','FWYF1 Air temp. Climatology \pm 47%');
%xlim([datenum(2016,12,1),now]); datetick3;
xlim([datenum(2016,6,1),now]); datetick3;
print('-dpng',fullfile(figspath,'South_Florida_AIR_AND_WATER_temps_Jun 2016-Jan_2017.png'));


%% AIR TEMP COMPARISON
fmg;
plot_ts(lkwf1.ndbc_air_t_1_d_avg,'.',...
        fwyf1.ndbc_air_t_1_d_avg,'.',...
        vakf1.ndbc_air_t_1_d_avg,'.',...
        a_climts,'Color',[0.8,0.8,1.0],'LineWidth',1.5,...
        a_lopctts,'--','Color',[.5,.5,.75],'LineWidth',1.5,...
        a_hipctts,'--','Color',[.5,.5,.75],'LineWidth',1.5);
legend('Air temp. Lake Worth pier (LKWF1)','Air temp. off Key Biscayne (FWYF1)',...
       'Virginia Key air temp. (VAKF1)','FWYF1 air temp. Climatology \pm 47%');
%xlim([datenum(2016,12,1),now]); datetick3;
xlim([datenum(2016,6,1),now]); datetick3;
print('-dpng',fullfile(figspath,'South_Florida_AIR_temps_Jun 2016-Jan_2017.png'));



%% REPRODUCE Brian McNoldy Twitter graph
lkwf1.ndbc_sea_t_anom = anomalize_ts_to_ts(lkwf1.ndbc_sea_t_1_d_avg,3);
fwyf1.ndbc_sea_t_anom = anomalize_ts_to_ts(fwyf1.ndbc_sea_t_1_d_avg,3);
mlrf1.ndbc_sea_t_anom = anomalize_ts_to_ts(mlrf1.ndbc_sea_t_1_d_avg,3);
