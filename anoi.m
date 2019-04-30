function stn = anoi(stn_or_stnm)

stn = get_station_from_station_name(stn_or_stnm);

stn = load_all_ndbc_data(stn);
stn = get_misst_station(stn);
stn = read_oisst2(stn);

stn = verify_variable(stn,'ndbc_sea_t_30_day_average');
stn = verify_variable(stn,'misst_sst_30_day_average');
stn = verify_variable(stn,'oisst2_sst_30_day_average');

mmmi=prctile(stn.ndbc_sea_t_30_day_average.data(ismember(get_year(stn.ndbc_sea_t_30_day_average.date),[1992:2006]),98.3));
disp(['In situ MMM = ' num2str(mmmi)]);
mmmm=prctile(stn.misst_sst_30_day_average.data(ismember(get_year(stn.misst_sst_30_day_average.date),[1992:2006]),98.3));
mmmo=prctile(stn.oisst2_sst_30_day_average.data(ismember(get_year(stn.oisst2_sst_30_day_average.date),[1992:2006]),98.3));

fmg;
plot_ts(stn.ndbc_sea_t_30_day_average,'-',stn.misst_sst_30_day_average,'-',stn.oisst2_sst_30_day_average,'-');
legend('In situ','MISST','OISST v2');
annotline([],mmmi,'In situ','blue');
annotline([],mmmm,'MISST',[.2 .6 .2]);
annotline([],mmmo,'OISSTv2','red');
titlename([upper(stn.station_name) ' monthly mean sea temperatures']);

return;
