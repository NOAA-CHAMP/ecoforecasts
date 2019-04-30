1;

figspath = get_ecoforecasts_path('figs');

lciy2 = get_station_from_station_name('lciy2');

lciy2 = load_station_data(lciy2);

lciy2 = get_erai_station(lciy2);

%%
% Surface PAR dose

% Manual QA/QC
lciy2.bic_surf_par.data(ismember(floor(lciy2.bic_surf_par.date),datenum(2011,1,[26,27,28])))=[];
lciy2.bic_surf_par.date(ismember(floor(lciy2.bic_surf_par.date),datenum(2011,1,[26,27,28])))=[];

% ERA-Interim parW [in W/m^2] is cr@p...? Calc directly from insolation
lciy2.erai_par.date=lciy2.erai_dsrf.date;
lciy2.erai_par.data=insol_to_par(lciy2.erai_dsrf.data);

lciy2.erai_par.data=lciy2.erai_par.data*1.05*1.05;

lciy2.bic_surf_par_daily_dose = par_dose(lciy2.bic_surf_par);
lciy2.erai_par_daily_dose = par_dose(lciy2.erai_par);

fmg; plot_ts(lciy2.bic_surf_par_daily_dose,lciy2.erai_par_daily_dose);
ylabel('PAR dose [mol-quanta^.m^-^2^.day^-^1]');
legend('In situ surface BIC','ERA-Interim', 'Location','Best');
titlename('Little Cayman CREWS/ICON Station LCIY2');

lciy2.bic_surf_par_30_day_dose = par_dose(lciy2.bic_surf_par,30*24);
lciy2.erai_par_30_day_dose = par_dose(lciy2.erai_par,30*24);

fmg; plot_ts(lciy2.bic_surf_par_30_day_dose,lciy2.erai_par_30_day_dose);
ylabel('PAR dose [mol-quanta^.m^-^2^.month^-^1]');
legend('In situ surface BIC','ERA-Interim', 'Location','Best');
titlename('Little Cayman CREWS/ICON Station LCIY2');



%%
% Attenuated benthic PAR dose

lciy2 = station_tmd_tide(lciy2);

lciy2 = calc_kd_daily(lciy2,'bic_surf_par','bic_deep_par','ctd_deep_i_depth');
fmg; plot_ts(lciy2.kd_1_day_bic_surf_par_bic_deep_par);
ylabel('K_d^P^A^R [m^-^1]');
titlename('Little Cayman CREWS/ICON Station LCIY2');


lciy2.bic_deep_par_daily_dose_calc.date = lciy2.bic_surf_par_daily_dose.date;
lciy2.bic_deep_par_daily_dose_calc.data = lciy2.bic_surf_par_daily_dose.data.*exp(-0.2.*5.5);
fmg; plot_ts(lciy2.bic_deep_par_daily_dose_calc);
ylabel('[mol-quanta^.m^-^2^.day^-^1]');
titlename('CREWS/ICON Station LCIY2 estimated benthic daily PAR dose @5.5m');


% fname = fullfile(get_ecoforecasts_path('data'),'lciy2_par_analysis.csv');
% csv_station_data(lciy2,{'bic_surf_par','bic_deep_par','kd_1_day_bic_surf_par_bic_deep_par','bic_surf_par_daily_dose','bic_deep_par_daily_dose_calc'},fname);


cleandts = unique(floor(lciy2.gSeaT.date(lciy2.gSeaT.data>15))) - 1;


lciy2.bic_deep_par_daily_dose = par_dose(lciy2.bic_deep_par);
lciy2.bic_deep_380nm_daily_dose = par_dose(lciy2.bic_deep_380nm);
lciy2.bic_deep_330nm_daily_dose = par_dose(lciy2.bic_deep_330nm);
lciy2.bic_deep_305nm_daily_dose = par_dose(lciy2.bic_deep_305nm);

lciy2 = verify_variable(lciy2,'bic_deep_par_daily_dose_3_day_average');
lciy2 = verify_variable(lciy2,'bic_deep_380nm_daily_dose_3_day_average');
lciy2 = verify_variable(lciy2,'bic_deep_330nm_daily_dose_3_day_average');
lciy2 = verify_variable(lciy2,'bic_deep_305nm_daily_dose_3_day_average');
lciy2 = verify_variable(lciy2,'wind2_speed_3_day_average');


fmg; plot_ts(lciy2.bic_surf_par,lciy2.bic_deep_par);
ylabel('PAR [mol-quanta^.m^-^2^.s^-^1]');
legend('Surface PAR','4.5m u/w PAR', 'Location','SouthWest');
for ix=1:length(cleandts)
  arrow([cleandts(ix),2500],[cleandts(ix),2300]);
end;
titlename('CREWS/ICON Station LCIY2 raw PAR data');
print('-dpng',fullfile(figspath,'lciy2-surface-vs-deep-par.png'));


fmg;
[ax,h1,h2] = ...
    plotyy(...
        [lciy2.bic_deep_380nm_daily_dose_3_day_average.date,...
         lciy2.bic_deep_330nm_daily_dose_3_day_average.date],...
        [lciy2.bic_deep_380nm_daily_dose_3_day_average.data,...
         lciy2.bic_deep_330nm_daily_dose_3_day_average.data],...
        lciy2.bic_deep_305nm_daily_dose_3_day_average.date,...
        lciy2.bic_deep_305nm_daily_dose_3_day_average.data);
datetick3('x',17,'keeplimits');
ylabel(ax(1),'UV-330 and UV-380 Dose [MJ^.m^-^2^.day^-^1]');
ylabel(ax(2),'UV-305 Dose [MJ^.m^-^2^.day^-^1]');
titlename('CREWS/ICON Station LCIY2 3-day mean 4.5m UV_3_8_0_,_3_3_0_,_3_0_5 daily dose');
for ix=1:length(cleandts)
  arrow([cleandts(ix),2],[cleandts(ix),1.8]);
end;
print('-dpng',fullfile(figspath,'lciy2-deep-uv-doses.png'));


fmg;
[ax,h1,h2] = ...
    plotyy(lciy2.bic_deep_par_daily_dose_3_day_average.date,...
           lciy2.bic_deep_par_daily_dose_3_day_average.data,...
           lciy2.wind2_speed_3_day_average.date,...
           lciy2.wind2_speed_3_day_average.data);
datetick3('x',17,'keeplimits');
ylabel(ax(1),'PAR Dose [mol-quanta^.m^-^2^.day^-^1]');
ylabel(ax(2),'Wind speed [m^.s^-^1]');
titlename('CREWS/ICON Station LCIY2 3-day mean wind speed and 4.5m PAR daily dose');
for ix=1:length(cleandts)
  arrow([cleandts(ix),40],[cleandts(ix),35]);
end;
print('-dpng',fullfile(figspath,'lciy2-deep-par-dose-vs-wind.png'));


clear ix ax h1 h2
clear figspath
