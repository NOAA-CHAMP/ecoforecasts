1;

more_status = get(0, 'More'); more('off');

figspath = 'figs';

flds = {'sea_t','sea_t_1_day_deviation_3_day_average','air_t', ...
        'air_t_1_day_deviation_3_day_average','wind1_u','wind1_v', ...
        'wind1_speed_3_day_average','wind1_u_7_day_deviation_sum_wind1_v'};

ylbls = {'T_s_e_a', '\mu_3_d(\sigma_1_d(T_s_e_a))', 'T_a_i_r', ...
         '\mu_3_d(\sigma_1_d(T_a_i_r))', 'W_U', 'W_V', '\mu_3_d(W)', ...
         '\sigma_7_d(W_U) + \sigma_7_d(W_V)'};

% First data from MLRF1 was archived for this date...
startdt = datenum(1987,12,04);

stanms = {'lkwf1','fwyf1','mlrf1','smkf1','sanf1','dryf1','42003'};

for ix = 1:length(stanms)
  stanm = stanms{ix};
  stn=load_all_ndbc_data([],stanm);
  multiplot_station(stn,flds,'','',ylbls,[startdt now]);
  print('-dtiff','-r300',fullfile(figspath,[stanm '-multiplots.tiff']));
  stn = [];
  clear stn;
end;

more(more_status);
clear more_status;
