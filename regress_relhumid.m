1;

stnm1 = 'smkf1';
stnm2 = 'lonf1';
% stnm2 = 'lkwf1';
% stnm2 = '42003';

disp('load');
if ( ~exist('stn1','var') || ~strcmp(stn1.station_name,stnm1) )
  stn1 = []; clear stn1;
  % stn1 = load_all_ndbc_data([], stnm1);
  load(['data/' stnm1 '-ndbc.mat'], 'station');
  stn1 = station; station = []; clear station;
end;
if ( ~exist('stn2','var') || ~strcmp(stn2.station_name,stnm2) )
  stn2 = []; clear stn2;
  stn2 = load_all_ndbc_data([], stnm2);
  station = stn2; save(['data/' stnm2 '-ndbc.mat'], 'station'); station = []; clear station;
  % load(['data/' stnm2 '-ndbc.mat'], 'station');
  % stn2 = station; station = []; clear station;
end;

disp('relhumid');
if ( ~isfield(stn1,'relhumid') )
  stn1 = station_dewp_to_relhumid(stn1,'air_t','dew_t','relhumid');
end;
if ( ~isfield(stn2,'relhumid') )
  stn2 = station_dewp_to_relhumid(stn2,'air_t','dew_t','relhumid');
end;

disp('intersect_dates');
[ix1,ix2] = intersect_dates(stn1.relhumid.date,stn2.relhumid.date);

dv = datevec(stn1.relhumid.date(ix1([1 end])));

disp('scatter_fit');
[B,Stats] = scatter_fit(stn1.relhumid.data(ix1), ...
                        stn2.relhumid.data(ix2), ...
                        [upper(stnm1) ' RH'], ...
                        [upper(stnm2) ' RH']);
title(sprintf('Relative Humidities %s v. %s: %d-%d',stnm1,stnm2,dv(1,1),dv(end,1)));
print('-dpng',sprintf('figs/scatter_fit-%s-%s-relhumid.png',stnm1,stnm2));

disp('mscohere');
[Cxy,Fn] = mscohere(stn1.relhumid.data(ix1), ...
                    stn2.relhumid.data(ix2));
figure;
set(gcf, 'units','normalized', 'outerposition',[0 0 1 1]);
% loglog(Fn,Cxy);
semilogx(Fn,Cxy);
title(sprintf('Relative Humidity coherence %s v. %s: %d-%d',stnm1,stnm2,dv(1,1),dv(end,1)));
print('-dpng',sprintf('figs/coherence-%s-%s-relhumid.png',stnm1,stnm2));
