1;

figspath = 'figs';

if (~exist('fwyf1','var'))
  fwyf1=load_all_ndbc_data([],'fwyf1');
end;
if (~exist('mlrf1','var'))
  mlrf1=load_all_ndbc_data([],'mlrf1');
end;
if (~exist('smkf1','var'))
  smkf1=load_all_ndbc_data([],'smkf1');
end;
if (~exist('sanf1','var'))
  sanf1=load_all_ndbc_data([],'sanf1');
end;
if (~exist('dryf1','var'))
  dryf1=load_all_ndbc_data([],'dryf1');
end;

flds = {'sea_t','sea_t_1_day_deviation_3_day_average','air_t', ...
        'air_t_1_day_deviation_3_day_average','wind1_u','wind1_v', ...
        'wind1_speed_3_day_average','wind1_u_7_day_deviation_sum_wind1_v'};

ylbls = {'T_s_e_a', '\mu_3_d(\sigma_1_d(T_s_e_a))', 'T_a_i_r', ...
         '\mu_3_d(\sigma_1_d(T_a_i_r))', 'U_W', 'V_W', '\mu_3_d(SPD_W)', ...
         '\mu_7_d(U_W) + \mu_7_d(V_W)'};

fwyf1=multiplot_station(fwyf1,flds,'','',ylbls);
print('-dtiff','-r300',fullfile(figspath,'fwyf1-multiplots.tiff'));

mlrf1=multiplot_station(mlrf1,flds,'','',ylbls);
print('-dtiff','-r300',fullfile(figspath,'mlrf1-multiplots.tiff'));

smkf1=multiplot_station(smkf1,flds,'','',ylbls);
print('-dtiff','-r300',fullfile(figspath,'smkf1-multiplots.tiff'));

sanf1=multiplot_station(sanf1,flds,'','',ylbls);
print('-dtiff','-r300',fullfile(figspath,'sanf1-multiplots.tiff'));

dryf1=multiplot_station(dryf1,flds,'','',ylbls);
print('-dtiff','-r300',fullfile(figspath,'dryf1-multiplots.tiff'));
