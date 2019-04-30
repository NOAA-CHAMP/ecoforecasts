1;

stn = load_station_data('mlrf1');
[stn.lon,stn.lat,stn.depth] = get_station_coords(stn.station_name);
stn = calc_kd(stn,'bic_surf_par','bic_shallow_par',3);

%x = importdata('mlrf2-cleaning-dates.csv'); stn.cleaning_date = datenum(x); x=[]; clear x
stn = load_cleaning_dates(stn);

fmg;
plot_ts(stn.bic_surf_par,stn.bic_shallow_par); titlename('Molasses Reef MLRF2'); xlim([min(stn.bic_surf_par.date),max(stn.bic_surf_par.date)]); datetick3('x',2,'keeplimits');
arrow([stn.cleaning_date(1),2500],[stn.cleaning_date(1),2300]); 
legend('Surface PAR','U/W PAR','Cleaning', 'Location','Best');
for ix=2:length(stn.cleaning_date); arrow([stn.cleaning_date(ix),2500],[stn.cleaning_date(ix),2300]); end;
% print('-dpng','figs/mlrf2-clean-par.png');

kd = stn.kd_bic_surf_par_bic_shallow_par;
ix = find( get_daylight(kd.date,stn.lat,stn.lon,40) & ...
           ismember(get_hour(kd.date),14:20) & ...
           (kd.data>0.0001) );
kd.date = kd.date(ix); kd.data = kd.data(ix);
keepix = [];
for ix=1:length(stn.cleaning_date)
  dtdif = kd.date - stn.cleaning_date(ix);
  keepix = union(keepix,find( (18/24) < dtdif & dtdif < 7 ));
end;

[cum,tid] = grp_ts(kd.data(keepix),kd.date(keepix),[],[],2);
fmg; plot(tid,cum,'*'); xlim([0,366]); datetick3('x',3,'keeplimits');
titlename('Seasonal Mean K_d^P^A^R');
print('-dpng','figs/mlrf2-clean-kd-par-climatology.png');

[cum,tid] = grp_ts(kd.data(keepix),kd.date(keepix),[],'numel',2);
fmg; plot(tid,cum,'*'); xlim([0,366]); datetick3('x',3,'keeplimits');
titlename('Total # Data Points');
print('-dpng','figs/mlrf2-clean-kd-par-climatology-N.png');
