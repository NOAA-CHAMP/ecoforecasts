1;

load('data/smkf1.mat');

station = qa_ts(station, 'tide', [], 1);
figure;
plot(station.tide_qc.date, station.tide_qc.data);
datetick; set_datetick_cursor; title('Tide_QC');


tide = station.tide.data;
tide = tide(~isnan(tide));
tide = (tide - mean(tide));
[tidepsd,fr] = pwelch(tide);
figure; loglog((fr*24),tidepsd); title('PWelch (raw)');

x=max(station.tide.date(isnan(station.tide.data)));
tide = station.tide.data(station.tide.date > x);
tide = (tide - mean(tide));
[tidepsd,fr] = pwelch(tide);
figure; loglog((fr*24),tidepsd); title('PWelch (>nan)');

tide = station.tide_qc.data;
tide=(tide-mean(tide));
[tidepsd,fr] = pwelch(tide);
figure; loglog((fr*24),tidepsd); title('PWelch (QC)');

[tidepsd,fr] = pcov(tide,20);
figure; loglog((fr*24),tidepsd); title('PCov (QC)');

[tidepsd,fr] = peig(tide,10);
figure; loglog((fr*24),tidepsd); title('PEig (QC)');



[name,fr,tcon,xout] = t_tide(station.tide.data, 'output', 'none');
size(xout)
figure; plot(station.tide.date, (station.tide.data - xout));
datetick; set_datetick_cursor;
title('Tide - T\_Tide');

dts = station.tide.date(station.tide.date >= datenum(2006,1,1));
tide = station.tide.data(station.tide.date >= datenum(2006,1,1));
[name,fr,tcon,xout] = t_tide(tide, 'output', 'none');
figure; plot(dts, (tide - xout));
datetick; set_datetick_cursor;
title('Tide - T\_Tide (2006-2008)');


[name,fr,tcon,xout] = t_tide(station.tide_qc.data, 'output', 'none');
figure; plot(station.tide_qc.date, (station.tide_qc.data - xout));
datetick; set_datetick_cursor;
title('Tide - T\_Tide (QC)');

dts = station.tide_qc.date(station.tide_qc.date >= datenum(2006,1,1));
tide = station.tide_qc.data(station.tide_qc.date >= datenum(2006,1,1));
[name,fr,tcon,xout] = t_tide(tide, 'output', 'none');
figure; plot(dts, (tide - xout)); datetick;
set_datetick_cursor;
title('Tide - T\_Tide (QC, 2006-2008)');



x=max(station.tide.date(isnan(station.tide.data)));
dts = station.tide.date(station.tide.date > x);
tide = station.tide.data(station.tide.date > x);

dts = station.tide.date(~isnan(station.tide.data));
tide = station.tide.data(~isnan(station.tide.data));
length(find(diff(dts) > 0.09)),

mu_tide = mean(tide);
dtide = tide - mu_tide;
dts(-2.5 > dtide | dtide > 2.5) = [];
tide(-2.5 > dtide | dtide > 2.5) = [];

mu_tide = mean(tide);
tide = tide - mu_tide;

fdts = dts(1):(1/24):dts(end);
ftide = interp1(dts, tide, fdts);

% fts = filter_ts(ftide, 3);
% figure; plot(fdts(40:end-40), mu_tide + (fts(40:end-40)));
% datetick('x', 'mmmyy', 'keepticks', 'keeplimits'); set_datetick_cursor;
% title('Tide 3HLP');

fts = filter_ts(ftide, 40);
figure; plot(fdts(40:end-40), mu_tide + (fts(40:end-40)));
datetick('x', 'mmmyy', 'keepticks', 'keeplimits'); set_datetick_cursor;
title('Tide 40HLP');

station = verify_variable(station, 'sea_t_1_day_deviation_3_day_average');
station = qa_ts(station, 'sea_t', [], 1);
figure;
plot(station.sea_t_1_day_deviation_3_day_average.date, ...
     station.sea_t_1_day_deviation_3_day_average.data);
datetick; set_datetick_cursor;
title('Sea T - \mu_3_d(\sigma_1_d(T))');







ymd2jday(2006,6,14)
[LONS, LATS, sst_iw, chl_iw, tsm_iw] = anmodis('2006165.18', [1 1 1]);
set(gca, 'clim', [22 29]);

sw_dist([24.63 , 24.70], [-81.11 , -80.71], 'km')
sw_dist([24.70, 24.98], [-80.71, -80.41], 'km')
sw_dist([24.57, 24.6], [-81.31, -81.22], 'km')

1/10e3
(2*pi)/10e3

k=(2*pi)/10e3
om = k*0.30,
om = k*0.50
om = k*0.10

k=(2*pi)/40e3
om = k*0.30,
om = k*0.50,

m=1/200; N=1e-3; f=7e-5;
om=sqrt(((N*k)^2)/(m^2))
om<f
N=1e-2;
om=sqrt(((N*k)^2)/(m^2))
om<f
N=5e-3;
om=sqrt(((N*k)^2)/(m^2))
N=3e-3;
om=sqrt(((N*k)^2)/(m^2))
N=2e-3;
om=sqrt(((N*k)^2)/(m^2))

sw_dist([24.57, 24.77], [-81.27, -80.76], 'km')

pd ../../RSMAS/Coastal/thesis
% [f,z] = animwera(num2str([2006164:2006166]'), [1 0 0]);
[f,z] = animwera(num2str([2006166:2006167]'), [1 0 0]);
for fi=1:length(f); figure(f(fi)); set(gca, 'clim', [3e-3 9e-3]); end;
reviewanim(f);

for fi=f; close(fi); end;



station = load_station_data('smkf1');

dts = station.tide.date(~isnan(station.tide.data));
tide = station.tide.data(~isnan(station.tide.data));

dts(abs(tide-mean(tide)) > 2.5) = [];
tide(abs(tide-mean(tide)) > 2.5) = [];

fts = filter_ts(tide, 40);
figure; plot(dts, fts); datetick('x', 'mmmyy', 'keepticks', 'keeplimits');

set_datetick_cursor; title('40HLP Tide - SMKF1');
[tidepsd,fr] = pwelch(fts(80:end-80));
figure; loglog(fr,tidepsd); title('PWelch (40HLP)');


x = max(station.tide.date(isnan(station.tide.data)));
dts = station.tide.date(station.tide.date > x);
tide = station.tide.data(station.tide.date > x);

dts(abs(tide-mean(tide)) > 2.5) = [];
tide(abs(tide-mean(tide)) > 2.5) = [];

fts = filter_ts(tide, 40);
[tidepsd,fr] = pwelch(fts(80:end-80));
figure; loglog(fr,tidepsd); title('PWelch (40HLP; 2004-2008)');

[tidepsd,fr] = pmtm(fts(80:end-80));
figure; loglog(fr,tidepsd); title('PSD - MTM (40HLP; 2004-2008)');

fts = filter_ts(tide, 72);
[tidepsd,fr] = pmtm(fts(144:end-144));
figure; loglog(fr,tidepsd); title('PSD - MTM (72HLP; 2004-2008)');
