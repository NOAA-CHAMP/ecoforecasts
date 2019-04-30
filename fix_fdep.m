1;

more off

fdepk = read_fdep_stevens_data('k',[],true,true);
fmg; plotyy_ts(fdepk.fdep_seatemp_deep,'-','Color',[0,0.5,0],fdepk.fdep_battv,'r-');
fmg; plot_ts(fdepk.fdep_battv,'.',fdepk.fdep_airtemp,'.',fdepk.fdep_seatemp_deep,'.',fdepk.fdep_wind_speed,'.'); legend('batt','Ta','Ts','U');

more on
