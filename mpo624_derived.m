function station = mpo624_derived(stanm_or_station)

  if ( isstruct(stanm_or_station) )
    station = stanm_or_station;
    stanm = station.station_name;
  elseif ( ischar(stanm_or_station) )
    stanm = stanm_or_station;
    load(['data/' stanm '-ndbc.mat'], 'station');
  else
    error('First arg must be a station-name string, or station STRUCT!');
  end;

  station = verify_variable(station, 'sea_t_1_day_deviation_3_day_average');
  station = verify_variable(station, 'sea_t_3_day_average_0_day_asof_diff_sea_t_1_day_minimum');
  station = verify_variable(station, 'wind1_u');
  station = verify_variable(station, 'wind1_v');

  %%% Do some monthly and annual boxplots - Tsea, Tair, Wu, Wv
%   do_boxplots(station, 'air_t', 'T_a_i_r', [5 35]);
%   do_boxplots(station, 'sea_t', 'T_s_e_a', [5 35]);
%   do_boxplots(station, 'wind1_u', 'W_U', [-30 30]);
%   do_boxplots(station, 'wind1_v', 'W_V', [-30 30]);

  fidx = 0;
  clear multidat;
  clear fldabbr;

  fidx = fidx + 1;
  multidts{fidx} = station.sea_t.date;
  multidat{fidx} = station.sea_t.data;
  fldabbr{fidx} = 'T_s_e_a';

  fidx = fidx + 1;
  multidts{fidx} = station.sea_t_1_day_deviation_3_day_average.date;
  multidat{fidx} = station.sea_t_1_day_deviation_3_day_average.data;
  fldabbr{fidx} = '\mu_3_d\sigma_1_d(T_s_e_a)';

  fidx = fidx + 1;
  multidts{fidx} = station.sea_t_3_day_average.date;
  multidat{fidx} = station.sea_t_3_day_average.data;
  fldabbr{fidx} = '\mu_3_d(T_s_e_a)';

  fidx = fidx + 1;
  multidts{fidx} = station.sea_t_3_day_average_0_day_asof_diff_sea_t_1_day_minimum.date;
  multidat{fidx} = station.sea_t_3_day_average_0_day_asof_diff_sea_t_1_day_minimum.data;
  fldabbr{fidx} = 'anom^m^i^n_3_d(T_s_e_a)';

  multiplot(multidts, multidat, ...
            'YLabel', fldabbr, ...
            'LineSpec', {'b-','g-','r-','c-'}, ...
            'Title', [upper(stanm) ': variables derived from in situ data']);
  set(gcf, 'units','normalized', 'outerposition',[0 0 1 1]);
  datetick2('keepticks','keeplimits'); set_datetick_cursor; drawnow;
%   print('-dpng', [stanm '-mpo624-derived.png']);

return;


function do_boxplots(station, fld, ylbl, ylm)

  stnm = station.station_name;
  [yrs,mos,dys] = datevec(station.(fld).date);

  figure;
  set(gcf, 'units','normalized', 'outerposition',[0 0 1 1]);
  boxplot(station.(fld).data, mos, 'notch','on', 'whisker',1.5);
  xlabel('Year Month'); title([upper(stnm) ' ' ylbl]);
  ylim(ylm);
  print('-dpng', sprintf('%s-mpo624-month-%s-boxplot.png', stnm, fld));

  figure;
  set(gcf, 'units','normalized', 'outerposition',[0 0 1 1]);
  boxplot(station.(fld).data, yrs, 'notch','on', 'whisker',1.5);
  xlabel('Year'); title([upper(stnm) ' ' ylbl]);
  ylim(ylm);
  print('-dpng', sprintf('%s-mpo624-year-%s-boxplot.png', stnm, fld));

return;
