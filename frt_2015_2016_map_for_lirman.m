1;

%basenm = 'frt_2015_2016_map_for_lirman';
basenm = mfilename;

if ( ~exist('doPrint','var') || isempty(doPrint) )
  doPrint = false;
end;
if ( ~exist('figspath','var') )
  %figspath = get_coral_path('Bleach');
  figspath = get_ecoforecasts_path('figs');
end;

lonf1 = get_station_from_station_name('lonf1'); lonf1 = load_all_ndbc_data(lonf1);
mlrf1 = get_station_from_station_name('mlrf1'); mlrf1 = load_all_ndbc_data(mlrf1);
fwyf1 = get_station_from_station_name('fwyf1'); fwyf1 = load_all_ndbc_data(fwyf1);

lonf1 = verify_variable(lonf1,'ndbc_sea_t_24_h_avg');
mlrf1 = verify_variable(mlrf1,'ndbc_sea_t_24_h_avg');
fwyf1 = verify_variable(fwyf1,'ndbc_sea_t_24_h_avg');

lonf1 = verify_variable(lonf1,'ndbc_sea_t_7_d_avg');
mlrf1 = verify_variable(mlrf1,'ndbc_sea_t_7_d_avg');
fwyf1 = verify_variable(fwyf1,'ndbc_sea_t_7_d_avg');


%{
disp('2014 is very like to 1997');
fmg;
subplot_tight(2,1,1);
plot_ts(lonf1.ndbc_sea_t_7_d_avg,mlrf1.ndbc_sea_t_7_d_avg,fwyf1.ndbc_sea_t_7_d_avg);
legend('LONF1','MLRF1','FWYF1');
axis([datenum(1997,4,15),datenum(1997,12,15),26.0,33.0]); datetick('x',12,'keeplimits');
grid on;
titlename('7-day average sea temperature');

subplot_tight(2,1,2);
plot_ts(lonf1.ndbc_sea_t_7_d_avg,mlrf1.ndbc_sea_t_7_d_avg,fwyf1.ndbc_sea_t_7_d_avg);
axis([datenum(2014,4,15),datenum(2014,12,15),26.0,33.0]); datetick('x',12,'keeplimits');
grid on;
if doPrint; print('-dpng',fullfile(figspath,[basenm,'-lonf1-mlrf1-fwyf1-ndbc_sea_t_7_d_avg_1997_vs_2014.png'])); end;

disp('2015 heating may have been more prolonged than 1998');
fmg;
subplot_tight(2,1,1);
plot_ts(lonf1.ndbc_sea_t_7_d_avg,mlrf1.ndbc_sea_t_7_d_avg,fwyf1.ndbc_sea_t_7_d_avg);
legend('LONF1','MLRF1','FWYF1');
axis([datenum(1998,4,15),datenum(1998,12,15),26.0,33.0]); datetick('x',12,'keeplimits');
grid on;
titlename('7-day average sea temperature');

subplot_tight(2,1,2);
plot_ts(lonf1.ndbc_sea_t_7_d_avg,mlrf1.ndbc_sea_t_7_d_avg,fwyf1.ndbc_sea_t_7_d_avg);
axis([datenum(2015,4,15),datenum(2015,12,15),26.0,33.0]); datetick('x',12,'keeplimits');
grid on;
if doPrint; print('-dpng',fullfile(figspath,[basenm,'-lonf1-mlrf1-fwyf1-ndbc_sea_t_7_d_avg_1998_vs_2015.png'])); end;
%}
%{
%}


if ( ~exist('hb','var') || ~isfield(hb,'lonf1') )
  x = load('d:/thesis/data/lonf1-heat_budget-erai-avhrr_weekly-ndbc-tpxo_tide-erai.mat');
  hb.lonf1 = x.stn; x=[]; clear x
  x = load('d:/thesis/data/mlrf1-heat_budget-erai-avhrr_weekly-ndbc-tpxo_tide-erai_DISSERT.mat');
  hb.mlrf1 = x.stn; x=[]; clear x

  hb.lonf1 = verify_variable(hb.lonf1,'simple_ndbc_erai_erai_30a_net_flux_1_d_sum');
  hb.mlrf1 = verify_variable(hb.mlrf1,'simple_ndbc_erai_erai_30a_net_flux_1_d_sum');
  hb.lonf1 = verify_variable(hb.lonf1,'simple_ndbc_erai_erai_30a_net_flux_term_1_d_sum');
  hb.mlrf1 = verify_variable(hb.mlrf1,'simple_ndbc_erai_erai_30a_net_flux_term_1_d_sum');

  hb.lonf1 = verify_variable(hb.lonf1,'simple_ndbc_erai_erai_30a_net_flux_7_d_sum');
  hb.mlrf1 = verify_variable(hb.mlrf1,'simple_ndbc_erai_erai_30a_net_flux_7_d_sum');
  hb.lonf1 = verify_variable(hb.lonf1,'simple_ndbc_erai_erai_30a_net_flux_term_7_d_sum');
  hb.mlrf1 = verify_variable(hb.mlrf1,'simple_ndbc_erai_erai_30a_net_flux_term_7_d_sum');
end;

disp('We use "simple_" for sea-surface, instead of accounting for *absorbed* shortwave radiation with "bottom" ("b_")');

disp('MLRF1 experienced more rapid heating than LONF1!');
%{
fmg; plot_ts(hb.lonf1.simple_ndbc_erai_erai_30a_net_flux,hb.mlrf1.simple_ndbc_erai_erai_30a_net_flux);
axis([datenum(1997,3,15),datenum(1997,12,15),-1000,1000]); datetick3;
legend('LONF1','MLRF1');
titlename('Hourly heating');
if doPrint; print('-dpng',fullfile(figspath,[basenm,'-lonf1-mlrf1-simple_ndbc_erai_erai_30a_net_flux_1_d_sum.png'])); end;

fmg; plot_ts(hb.lonf1.simple_ndbc_erai_erai_30a_net_flux_1_d_sum,hb.mlrf1.simple_ndbc_erai_erai_30a_net_flux_1_d_sum);
axis([datenum(1997,3,15),datenum(1997,12,15),-5000,5000]); datetick3;
legend('LONF1','MLRF1');
titlename('One-day heating');
if doPrint; print('-dpng',fullfile(figspath,[basenm,'-lonf1-mlrf1-simple_ndbc_erai_erai_30a_net_flux_1_d_sum.png'])); end;
%}

%{
fmg; plot_ts(hb.lonf1.simple_ndbc_erai_erai_30a_net_flux_term,hb.mlrf1.simple_ndbc_erai_erai_30a_net_flux_term);
axis([datenum(1997,3,15),datenum(1997,12,15),-0.35,+0.35]); datetick3;
legend('LONF1','MLRF1');
titlename('Hourly heating');
if doPrint; print('-dpng',fullfile(figspath,[basenm,'-lonf1-mlrf1-simple_ndbc_erai_erai_30a_net_flux_1_d_sum.png'])); end;

fmg; plot_ts(hb.lonf1.simple_ndbc_erai_erai_30a_net_flux_term_1_d_sum,hb.mlrf1.simple_ndbc_erai_erai_30a_net_flux_term_1_d_sum);
axis([datenum(1997,3,15),datenum(1997,12,15),-0.35,+0.35]); datetick3;
legend('LONF1','MLRF1');
titlename('Seven-day heating');
if doPrint; print('-dpng',fullfile(figspath,[basenm,'-lonf1-mlrf1-simple_ndbc_erai_erai_30a_net_flux_term_1_d_sum.png'])); end;
%}

%{
%}
% %t0 = datenum(1997,6,1); tN = datenum(1997,11,1);
% t0 = datenum(1997,7,1); tN = datenum(1997,10,1);
t0 = datenum(1994,7,1); tN = datenum(1994,10,1);
error('2004 looks best for "normal"!');
for yr=1993:2009
t0 = datenum(yr,7,1); tN = datenum(yr,10,1);
%%fld = 'simple_ndbc_erai_erai_30a_net_flux_term';
%fld = 'ndbc_erai_erai_30a_net_flux_term';
fld = 'ndbc_erai_erai_30a_avhrr_hc_dTdt';
lon = date_range_ts(hb.lonf1.(fld),[t0,tN]);
mlr = date_range_ts(hb.mlrf1.(fld),[t0,tN]);
[ig,lont0ix] = min(abs(hb.lonf1.ndbc_sea_t.date-t0));
lont0 = hb.lonf1.ndbc_sea_t.data(lont0ix);
[ig,mlrt0ix] = min(abs(hb.mlrf1.ndbc_sea_t.date-t0));
mlrt0 = hb.mlrf1.ndbc_sea_t.data(mlrt0ix);

fmg;
plot(lon.date,lont0+cumsum(lon.data),mlr.date,mlrt0+cumsum(mlr.data)); datetick3;
axis(axis);
plot_ts(hb.lonf1.ndbc_sea_t,'k:',hb.mlrf1.ndbc_sea_t,'k-');
legend('LONF1','MLRF1');
titlename('Modeled vs. Measured Sea Temperature');
pause;
end;
if doPrint; print('-dpng',fullfile(figspath,[basenm,'-lonf1-mlrf1-',fld,'_',datestr(t0,'yyyy-mm-dd'),'_',datestr(tN,'yyyy-mm-dd'),'.png'])); end;


if ( ~exist('subrgn','var') )
  subrgn='FRT';
  use_habitat_map=true;
  allow_hard_bottom=true; 
end;
if ( ~exist('dTdt_SS','var') )
  doFigs=false; 
  calc_spatial_dt_hc_thermal_stress;
end;

if ( ~exist('bath','var') || ~isfield(bath,'field') )
  x = read_hires_bathymetry_for_field({lon,lat},false);
  bath = x.ngdc_hires_bathy; x=[]; clear x
end;

%scenarios = 4; % Only interested in the bleaching summer (e.g., 1997 which is like 2014)
scenarios = 3:4; % Interested in the bleaching summer (e.g., 1997, like 2014)
                 % and a "max monthly mean" year for comparison (e.g., 1994)
%doPrint=false; plot_sites=true; plot_coastline=true; plot_spatial_dt_hc
plot_sites=false;
plot_coastline=true;
plot_spatial_dt_hc_thermal_stress
