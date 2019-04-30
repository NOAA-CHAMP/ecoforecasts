%-- 5/25/2016 1:53 PM --%
quit
%-- 5/25/2016 1:58 PM --%
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
help extract_wave_model
extract_wave_model('hi',datenum(2012,07,1):now)
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
extract_wave_model('nwhi')
%-- 5/25/2016 4:36 PM --%
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
dir *.m
dir *slop*
dir *slop*m
pd
pd +2
dir *slop*m
dir *fig*m
pd
dir *fig*.m
pd +9
cd ../../Gramer_1
dirs
pd
cd ../../Gramer_1
cd ../../Gramer_2
cd ../Gramer_2
dirs
pd
po
pd +11
po
pd ../in_prep/Gramer_Mariano_2
dir *.m
popd
pd ../in_review/Gramer_Mariano_1
dir *.m
po
pd +2
pd +16
dirs
help popd
type popd
dirs
popd +16
pd
pd ../Walker_et_al
dir Flanary
dir *.m
pd
pd +2
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
quit
%-- 5/25/2016 5:16 PM --%
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
help rename
help move
help mv
help ren
help rem
stn = get_wave_model('manm1',[
])
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
quit
%-- 5/26/2016 1:04 PM --%
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
quit
%-- 5/26/2016 7:26 PM --%
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
d
dirs
pd +4
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
cmpefkeys
%stnm = 'mlrf1';
stnm = 'sanf1';
stn = []; clear stn
stn = get_station_from_station_name(stnm);
stn = station_optimal_isobath_orientation(stn);
stn = plot_hires_bathymetry(stn,-[0:5:80],[12e3,12e3],true,[],false);
fh = gcf;
axis(axis);
[c,h] = contour(mdl.lon,mdl.lat,squeeze(mdl.t(1,1,:,:)),[floor(min_t):0.25:ceil(max_t)]);
clabel(c,h);
c=[]; clear c h
[res,stn] = plot_bathy_transect(stn,[],stn.isobath_orientation+90,'ngdc_hires_bathy');
[dx,az] = distance_wgs84(stn.lat,stn.lon,res.lat,res.lon);
dx(round(az) ~= stn.isobath_orientation+90) = -dx(round(az) ~= stn.isobath_orientation+90);
res.dx = dx;
lonix = interp1(mdl.lon,1:mdl.nlon,res.lon,'nearest');
latix = interp1(mdl.lat,1:mdl.nlat,res.lat,'nearest');
ix = sub2ind([size(mdl.t,3),size(mdl.t,4)],latix,lonix);
t = squeeze(mdl.t(1,:,ix));
[c,h] = contourf(res.dx,-mdl.z,t);
clabel(c,h);
c=[]; clear c h
plot(res.dx,res.field,'k-','LineWidth',2);
titlename([upper(stnm),' bathymetry vs. cross-shore T']);
figure(fh);
plot(res.lon,res.lat,'k-','LineWidth',2);
clear ans cstnms stix stnm
clear az dx ig ix latix lonix t
%clear res stn min_t max_t min_curr max_curr max_spd
close all
stn=[]; clear stn
clear datpath dt dtix dts fh hr hycpath jd matfname ncname ncpath
clear datpath dt dtix dts fh hr hycpath jd matfname ncname ncpath yr
fmg; contourf(mdl.lon,mdl.lat,squeeze(nanmean(mdl.t(:,2,:,:)))); colorbar;
fmg; contourf(mdl.lon,mdl.lat,squeeze(nanmean(mdl.t(:,2,:,:)))); colorbar; daspect([1,cosd(mdl.lat(1)),1]);
fmg; contourf(mdl.lon,mdl.lat,squeeze(nanmean(mdl.t(:,1,:,:)))); colorbar; daspect([1,cosd(mdl.lat(1)),1]);
clear ans
stn = get_station_from_station_name(stnm);
stnm = 'fwyf1';
stn = get_station_from_station_name(stnm);
stn = load_all_ndbc_data(stn);
scatter_fit_ts(mdl.(stnm).t,stn.ndbc_sea_t);
stnm = 'mlrf1';
stn = get_station_from_station_name(stnm); stn = load_all_ndbc_data(stn);
scatter_fit_ts(mdl.(stnm).t,stn.ndbc_sea_t);
mdl.(stnm)
mdl.(stnm).t
stnm = 'sanf1';
stn = get_station_from_station_name(stnm); stn = load_all_ndbc_data(stn);
scatter_fit_ts(mdl.(stnm).t,stn.ndbc_sea_t);
stnm = 'smkf1';
stn = get_station_from_station_name(stnm); stn = load_all_ndbc_data(stn);
scatter_fit_ts(mdl.(stnm).t,stn.ndbc_sea_t);
stnm = 'lonf1';
stn = get_station_from_station_name(stnm); stn = load_all_ndbc_data(stn);
scatter_fit_ts(mdl.(stnm).t,stn.ndbc_sea_t);
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
cmpefkeys
mdl.z
stn
stn = load_all_ndbc_data(stn);
scatter_fit_ts(mdl.(stnm).t,stn.ndbc_sea_t);
stnm='fwyf1';
scatter_fit_ts(mdl.(stnm).t,stn.ndbc_sea_t);
scatter_fit_ts(mdl.(stnm).t1,stn.ndbc_sea_t);
scatter_fit_ts(mdl.(stnm).t10,stn.ndbc_sea_t);
mdl.(stnm)
mdl.(stnm).t10
mdl.(stnm).t10.data
mdl.(stnm).t1.data
stnm='mlrf1';
stn = get_station_from_station_name(stnm); stn = load_all_ndbc_data(stn);
scatter_fit_ts(mdl.(stnm).t10,stn.ndbc_sea_t);
scatter_fit_ts(mdl.(stnm).t1,stn.ndbc_sea_t);
scatter_fit_ts(mdl.(stnm).t20,stn.ndbc_sea_t);
mdl.(stnm)
scatter_fit_ts(mdl.(stnm).t50,stn.ndbc_sea_t);
scatter_fit_ts(mdl.(stnm).t40,stn.ndbc_sea_t);
scatter_fit_ts(mdl.(stnm).t30,stn.ndbc_sea_t);
close all
scatter_fit_ts(mdl.(stnm).t20,stn.ndbc_sea_t);
scatter_fit_ts(mdl.(stnm).t10,stn.ndbc_sea_t);
close all
scatter_fit_ts(mdl.(stnm).t10,stn.ndbc_sea_t);
scatter_fit_ts(mdl.(stnm).t20,stn.ndbc_sea_t);
close all
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
stn = get_station_from_station_name('lciy7'); stn = load_station_data(stn);
stn
stn = load_ftpxls_data('data/lciy2-2010.xls');
sth
stn
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
stn = get_station_from_station_name('lciy2'); stn = load_station_data(stn);
stn
stn.ctd_deep_seatemp
find_date_ranges(stn.ctd_deep_seatemp.date)
fmg; boxplot_ts(stn.ctd_deep_seatemp,@get_year)
fmg; boxplot_ts(stn.ctd_deep_seatemp,@get_yearmonth)
datetick3;
fmg; boxplot_ts(stn.ctd_deep_seatemp,@(x)(datestr(get_yearmonth(x.date))))
fmg; boxplot_ts(stn.ctd_deep_seatemp,@(x)(datestr(get_yearmonth(x))))
help datestr
fmg; boxplot_ts(stn.ctd_deep_seatemp,@(x)(datestr(get_yearmonth(x),'mmm-yy')));
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
stn = get_station_from_station_name('llbp7'); stn = load_station_data(stn);
fmg; boxplot_ts(stn.ctd_deep_seatemp,@(x)(datestr(get_yearmonth(x),'mmm-yy')));
help boxplot_ts
fmg; boxplot_ts(stn.ctd_deep_seatemp,@(x)(datestr(get_yearmonth(x),'mmm-yy')));
ix=find(get_yearmonth(stn.ctd_deep_seatemp.date)~=datenum(2012,8,1));
stn.ctd_deep_seatemp
fmg; boxplot_ts(stn.ctd_deep_seatemp,@(x)(datestr(get_yearmonth(x),'mmm-yy')),'indices',ix);
help ismember
help get_yearmonth
unitsratio('ft','m')*100
unitsratio('m','ft')*100
good
calendar
datestr(datenum(2016,5,31)+10)
calendar(2016,6)
fmg; boxplot_ts(stn.ctd_deep_seatemp,@(x)(datestr(get_yearmonth(x),'mmm-yy')),'indices',ix);
ix=find(~ismember(get_yearmonth(stn.ctd_deep_seatemp.date),[datenum([2011,2011,2012,2012,2013],[8,9,3,8,4],1)));
ix=find(~ismember(get_yearmonth(stn.ctd_deep_seatemp.date),datenum([2011,2011,2012,2012,2013],[8,9,3,8,4],1)));
fmg; boxplot_ts(stn.ctd_deep_seatemp,@(x)(datestr(get_yearmonth(x),'mmm-yy')),'indices',ix);
all
which all
help isdatetime
help datetime
fmg; boxplot_ts(stn.ctd_deep_seatemp,@(x)(datestr(get_yearmonth(x),'mmm-yy')),'indices',@(x)(find_pct_good_dates(x.date,@get_yearmonth,0.75)));
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
stn = get_station_from_station_name('srvi2'); stn = load_station_data(stn);
find_date_ranges(stn.ctd_deep_seatemp.date)
find_date_ranges(stn.ctd_shallow_seatemp.date)
fmg; boxplot_ts(stn.ctd_shallow_seatemp,@(x)(datestr(get_yearmonth(x),'mmm-yy')),'indices',@(x)(find_pct_good_dates(x.date,@get_yearmonth,0.75)));
numel(find(ismember(get_yearmonth(stn.ctd_shallow_seatemp.date),datenum(2006,1,1))))
ans/24
numel(find(ismember(get_yearmonth(stn.ctd_shallow_seatemp.date),datenum(2006,1,1))))
min(diff(unique(stn.ctd_shallow_seatemp.date)))
min(diff(unique(stn.ctd_shallow_seatemp.date)))*24
median(diff(unique(stn.ctd_shallow_seatemp.date)))
min(diff(unique(stn.ctd_shallow_seatemp.date)))*24*60
help ts_dt
dir *_ts*
help interp_ts
t = interp_ts(stn.ctd_shallow_seatemp);
t = interp_ts(stn.ctd_shallow_seatemp,1/24);
min(diff(unique(t.date)))*24*60
fmg; plot_ts(t,stn.ctd_shallow_seatemp);
fmg; plot_ts(stn.ctd_shallow_seatemp,t);
legend('CTD','Hourly');
datetick3;
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
anefkeys
fmg; contourf(mdl.lon,mdl.lat,squeeze(nanmean(mdl.t(:,1,:,:)))); colorbar; daspect([1,cosd(mdl.lat(1)),1]);
set_surf_cursor
495-199
394-163
231*296*7*7*8*4
231*296*7*7*8*4/1e9
27*231*296*7*7*8*4/1e9
27*4*231*296*4*8/1e9
27*4*300*296*4*8/1e9
mdl.lat
mdl
27*4*589*296*4*8/1e9
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
anefkeys
whos
2^16
2^32
2^31
2^30
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
anefkeys
fmg; contourf(mdl.lon,mdl.lat,squeeze(nanmean(mdl.t(:,1,:,:)))); colorbar; daspect([1,cosd(mdl.lat(1)),1]);
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
cmpefkeys
mdl.true_z
mdl.z
close all
get_station_from_station_name('dryf1')
get_station_from_station_name('plsf1')
if ( 'abc' ); disp('y'); end;
clear ans cstnms stix stnm
clear az dx ig ix latix lonix s t
stn
stn=[]; clear res stn min_t max_t min_s max_s min_curr max_curr max_spd
clear datpath hycpath matfname
clear fh
dl
mdl
mdl.z
1.06e9*(10/9)
mdl.true_z
clear ans
cmpefkeys
reviewanim([],0,0,0)
get(0,'child')
reviewanim(1:12,0,0,0)
reviewanim(figsnamed('cross-shore'),0,0,0)
mdl.lonf1.t
mdl.lonf1.t.data
mdl.lonf1.t.prof(:,1)
mdl.lonf1.t.prof(:,2)
reviewanim(figsnamed('cross-shore'),0,0,0)
ylim([-100,0])
ylim([-120,0])
ylim tight
ylim default
ylim([-10,0])
stn=[]; clear res stn min_t max_t min_s max_s min_curr max_curr max_spd
stn=[]; clear ans res stn min_t max_t min_s max_s min_curr max_curr max_spd
fmg; contourf(mdl.lon,mdl.lat,squeeze(nanmean(mdl.t(:,1,:,:)))); colorbar; daspect([1,cosd(mdl.lat(1)),1]);
load('c:\users\gramer\Documents\MATLAB\ecoforecasts\coast\key_west_fl_mhw.mat
load('c:\users\gramer\Documents\MATLAB\ecoforecasts\coast\key_west_fl_mhw.mat')
min(diff(unique(lat)))
min(diff(unique(lat)))*111e3
min(diff(unique(lon)))*111e3
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
anefkeys
mdl.true_z
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
anefkeys
dbcont
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
cmpefkeys
reviewanim(1:12,0,0,0)
reviewanim(figsnamed('cross-shore'),0,0,0)
close all
mdl.z
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
1.31e9*(4/3)
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
anefkeys
whos
mdl.s
mdl
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
anefkeys
dbcont
whos
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
cmpefkeys
stn=[]; clear ans res stn min_t max_t min_s max_s min_curr max_curr max_spd
clear datpath hycpath matfname
mdl
mdl.mlrf1
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
load('coast/eFKEYS_01-1-7-2012_ORIG.mat');
load('data/eFKEYS_01-1-7-2012_ORIG.mat');
mdl.z
mdl
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
load('data/eFKEYS_01-1-7-2012_FULL_REGION.mat')
mdl
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
12/sqrt(2)
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
cmpefkeys
mdl.z
reviewanim(figsnamed('cross-shore'),0,0,0)
max_t
min_t
mdl.smkf1
mdl.smkf1.t120
nansummary(mdl.smkf1.t120.data)
nansummary(mdl.smkf1.t1-0.data)
nansummary(mdl.smkf1.t100.data)
nansummary(mdl.smkf1.t80.data)
nansummary(mdl.smkf1.t30.data)
nansummary(mdl.smkf1.t40.data)
nansummary(mdl.smkf1.t50.data)
nansummary(mdl.smkf1.t60.data)
reviewanim(figsnamed('cross-shore'),0,0,0)
ylim([-200,0])
reviewanim(figsnamed('cross-shore'),0,0,0)
reviewanim([],0,0,0)
ylim([-200,0])
reviewanim(figsnamed('cross-shore'),0,0,0)
figsnamed('cross-shore')
for fh=figsnamed('cross-shore'); figure(fh); xlim([-150,0]); end;
fh
for fh=figsnamed('cross-shore')'; figure(fh); xlim([-150,0]); end;
xlim default
for fh=figsnamed('cross-shore')'; figure(fh); xlim([-10,+10]); ylim([-150,0]); end;
timenow
ylim default
ylim([-10,0])
.
ylim([-10,0])
whos
mdl.true_z
mdl.z
ylim([-150,0]);
calendar(2016,6)
ylim([-140,0]);
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
anefkeys
mdl
mdl.z
close(nc); clear nc
dbcont
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
anefkeys
dbcont
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
anefkeys
clear datpath hycpath matfname ans
cmpefkeys
reviewanim(figsnamed('cross-shore'),0,0,0)
ylim default
ylim([-150,0]);
ylim([-140,0]);
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
stn = get_station_from_station_name('looe1'); stn = get_fkeys_hycom(stn);
stn
stn.fkeys_hycom_u
stn.fkeys_hycom_u_field
clear ans
stn
stn.fkeys_hycom_u
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
dirs
pd +11
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
stn = get_station_from_station_name(stnm); stn = load_all_ndbc_data(stn);
stnm = 'mlrf1;
stnm = 'mlrf1';
stn = get_station_from_station_name(stnm); stn = load_all_ndbc_data(stn);
fmg; plot_ts(stn.ndbc_sea_t,stn.ndbc_air_t);
fmg; plot_ts(stn.ndbc_air_t,stn.ndbc_sea_t);
xlim(datenum([2004,2008],[1,2],1)); datetick3;
fmg; plot_ts(stn.ndbc_air_t,stn.ndbc_sea_t);
xlim(datenum(2012,1,[1,8])); datetick3;
for yr=2004:2008; xlim(datenum(yr,1,[1,8])); datetick3; pause; end;
for yr=2004:2008; xlim(datenum(yr,1,[1,8])); ylim([8,26]); datetick3; pause; end;
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
which anfkeys
hycpath = get_thesis_path('../data/hycom/FKEYS');
dts = datenum(2008,1,1,[6:6:(7*24)-1],0,0);
for dtix=1:numel(dts)
dt = dts(dtix);
yr = get_year(dt);
jd = get_jday(dt);
hr = get_hour(dt);
disp(datestr(dt));
ncname = sprintf('902_archv.%04d_%03d_%02d_3zu.nc',yr,jd,hr);
ncpath = fullfile(hycpath,ncname);
if ( ~exist(ncpath,'file') )
warning('Skipping %s',ncpath);
else
mdl.d(dtix,1) = dt;
nc = mDataset(ncpath);
break;
end;
end;
nj_info(nc),
close(nc); clear nc
ncname = sprintf('902_archv.%04d_%03d_%02d_3zt.nc',yr,jd,hr);
ncpath = fullfile(hycpath,ncname);
nc = mDataset(ncpath);
nj_info(nc),
close(nc); clear nc
ncname = sprintf('902_archv.%04d_%03d_%02d_3zs.nc',yr,jd,hr);
ncpath = fullfile(hycpath,ncname);
nc = mDataset(ncpath);
nj_info(nc),
close(nc); clear nc
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
anfkeys
mdl
whos
mdl.z
mdl.true_z
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
anfkeys
cmpefkeys
reviewanim(figsnamed('cross-shore'),0,0,0)
min_t
max_t
stn=[]; clear ans res stn min_t max_t min_s max_s min_curr max_curr max_spd
clear datpath hycpath matfname ans
mdl
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
doPrint = true;
doFKEYS = true;
cmpefkeys
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
doFKEYS = true;
cmpefkeys
mdl=[]; clear all
pack
cmpefkeys
max_t
min_t
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
doPrint = true;
cmpefkeys
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
cmpefkeys
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
cmpefkeys
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
doPrint = true;
doFKEYS = true;
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
doFKEYS = true;
cmpefkeys
reviewanim(figsnamed('cross-shore'),0,0,0)
reviewanim([],0,0,0)
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
doPrint = true;
cmpefkeys
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
doPrint = true;
doFKEYS = true;
cmpefkeys
mdl=[]; clear all
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
cmpefkeys
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
help read_hires_bathymetry
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
help plot_bathy_transect
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
doPrint = true;
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
cmpefkeys
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
cmpefkeys
close all
stn=[]; clear ans res stn min_t max_t min_s max_s min_curr max_curr max_spd
clear datpath hycpath matfname ans
clear datpath hycpath matfname ncols ncstnms ans
cstnms = {'fwyf1','mlrf1','lonf1','smkf1','looe1','sanf1'};
ncstnms = numel(cstnms);
fmg;
ncols = floor(ncstnms/2);
for stix=1:2:ncstnms
stnm = cstnms{stix};
spt(2,floor(ncstnms/2),stix);
hist(mdl.(stnm).l.data);
xlim([-max_spd,+max_spd]);
xlabel([upper(stnm),' alongshore (m/s)']);
spt(2,ncols,stix+1);
hist(mdl.(stnm).x.data);
xlim([-max_spd,+max_spd]);
xlabel([upper(stnm),' cross-shore (m/s)']);
end;
min_t = +inf;
max_t = -inf;
min_s = +inf;
max_s = -inf;
min_curr = +inf;
max_curr = -inf;
for stix=1:ncstnms
stnm = cstnms{stix};
min_t = nanmin([min_t,nanmin(mdl.(stnm).t.prof(:))]);
max_t = nanmax([max_t,nanmax(mdl.(stnm).t.prof(:))]);
if ( isfield(mdl.(stnm),'s') )
min_s = nanmin([min_s,nanmin(mdl.(stnm).s.prof(:))]);
max_s = nanmax([max_s,nanmax(mdl.(stnm).s.prof(:))]);
end;
min_curr = nanmin([min_curr,nanmin(mdl.(stnm).u.prof(:)),nanmin(mdl.(stnm).v.prof(:))]);
max_curr = nanmax([max_curr,nanmax(mdl.(stnm).u.prof(:)),nanmax(mdl.(stnm).v.prof(:))]);
end;
max_spd = max(abs(max_curr),abs(min_curr));
clear ans
fmg;
ncols = floor(ncstnms/2);
for stix=1:2:ncstnms
stnm = cstnms{stix};
spt(2,floor(ncstnms/2),stix);
hist(mdl.(stnm).l.data);
xlim([-max_spd,+max_spd]);
xlabel([upper(stnm),' alongshore (m/s)']);
spt(2,ncols,stix+1);
hist(mdl.(stnm).x.data);
xlim([-max_spd,+max_spd]);
xlabel([upper(stnm),' cross-shore (m/s)']);
end;
fmg;
ncols = floor(ncstnms/2);
for stix=1:2:ncstnms
stnm = cstnms{stix};
spt(2,floor(ncstnms/2),stix);
hist(mdl.(stnm).l.data,floor(-max_spd):0.05:ceil(+max_spd));
xlim([-max_spd,+max_spd]);
xlabel([upper(stnm),' alongshore (m/s)']);
spt(2,ncols,stix+1);
hist(mdl.(stnm).x.data,floor(-max_spd):0.05:ceil(+max_spd));
xlim([-max_spd,+max_spd]);
xlabel([upper(stnm),' cross-shore (m/s)']);
end;
fmg;
ncols = floor(ncstnms/2);
for stix=1:2:ncstnms
stnm = cstnms{stix};
%spt(2,floor(ncstnms/2),((stix-1)*ncols)+1);
spt(2,ncols,stix);
hist(mdl.(stnm).l.data,floor(-max_spd):0.05:ceil(+max_spd));
xlim([-max_spd,+max_spd]);
xlabel([upper(stnm),' alongshore (m/s)']);
spt(2,ncols,stix+1);
hist(mdl.(stnm).x.data,floor(-max_spd):0.05:ceil(+max_spd));
xlim([-max_spd,+max_spd]);
xlabel([upper(stnm),' cross-shore (m/s)']);
end;
fmg;
ncols = floor(ncstnms/2);
for stix=1:2:ncstnms
stnm = cstnms{stix};
spt(2,ncols,((stix-1)*ncols)+1);
hist(mdl.(stnm).l.data,floor(-max_spd):0.05:ceil(+max_spd));
xlim([-max_spd,+max_spd]);
xlabel([upper(stnm),' alongshore (m/s)']);
spt(2,ncols,((stix-1)*ncols)+2);
hist(mdl.(stnm).x.data,floor(-max_spd):0.05:ceil(+max_spd));
xlim([-max_spd,+max_spd]);
xlabel([upper(stnm),' cross-shore (m/s)']);
end;
fmg;
ncols = floor(ncstnms/2);
for stix=1:2:ncstnms
stnm = cstnms{stix};
spt(2,ncols,floor(stix/2)+1);
hist(mdl.(stnm).l.data,floor(-max_spd):0.05:ceil(+max_spd));
xlim([-max_spd,+max_spd]);
xlabel([upper(stnm),' alongshore (m/s)']);
spt(2,ncols,ncols+floor(stix/2)+1);
hist(mdl.(stnm).x.data,floor(-max_spd):0.05:ceil(+max_spd));
xlim([-max_spd,+max_spd]);
xlabel([upper(stnm),' cross-shore (m/s)']);
end;
fmg;
ncols = floor(ncstnms/2);
for stix=1:2:ncstnms
stnm = cstnms{stix};
spt(2,ncols,floor(stix/2)+1);
hist(mdl.(stnm).l.data,floor(-max_spd):0.025:ceil(+max_spd));
xlim([-max_spd,+max_spd]);
xlabel([upper(stnm),' alongshore (m/s)']);
spt(2,ncols,ncols+floor(stix/2)+1);
hist(mdl.(stnm).x.data,floor(-max_spd):0.025:ceil(+max_spd));
xlim([-max_spd,+max_spd]);
xlabel([upper(stnm),' cross-shore (m/s)']);
end;
max_curr
max_spd
min_curr
cstnms
fmg;
stnixen = [1,2,4,5];
ncols = numel(stnixen);
colix = 1;
for stix=stnixen(:)';
stnm = cstnms{stix};
spt(2,ncols,colix);
hist(mdl.(stnm).l.data,floor(-max_spd):0.025:ceil(+max_spd));
xlim([-max_spd,+max_spd]);
xlabel([upper(stnm),' alongshore (m/s)']);
spt(2,ncols,colix+ncols);
hist(mdl.(stnm).x.data,floor(-max_spd):0.025:ceil(+max_spd));
xlim([-max_spd,+max_spd]);
xlabel([upper(stnm),' cross-shore (m/s)']);
colix = colix + 1;
end;
fmg;
stnixen = [1,2,4,5];
ncols = numel(stnixen);
colix = 1;
%min_hist = floor(-max_spd); max_hist = ceil(+max_spd);
min_hist = -0.2; max_hist = +0.2;
for stix=stnixen(:)';
stnm = cstnms{stix};
spt(2,ncols,colix);
hist(mdl.(stnm).l.data,min_hist:0.02:max_hist);
xlim([-max_spd,+max_spd]);
xlabel([upper(stnm),' alongshore (m/s)']);
spt(2,ncols,colix+ncols);
hist(mdl.(stnm).x.data,min_hist:0.02:max_hist);
xlim([-max_spd,+max_spd]);
xlabel([upper(stnm),' cross-shore (m/s)']);
colix = colix + 1;
end;
fmg;
stnixen = [1,2,4,5];
ncols = numel(stnixen);
colix = 1;
%min_hist = floor(-max_spd); max_hist = ceil(+max_spd);
min_hist = -0.2; max_hist = +0.2;
for stix=stnixen(:)';
stnm = cstnms{stix};
spt(2,ncols,colix);
hist(mdl.(stnm).l.data,min_hist:0.02:max_hist);
xlabel([upper(stnm),' alongshore (m/s)']);
spt(2,ncols,colix+ncols);
hist(mdl.(stnm).x.data,min_hist:0.02:max_hist);
xlabel([upper(stnm),' cross-shore (m/s)']);
colix = colix + 1;
end;
fmg;
stnixen = [1,2,4,5];
ncols = numel(stnixen);
colix = 1;
%min_hist = floor(-max_spd); max_hist = ceil(+max_spd);
min_hist = -0.2; max_hist = +0.2;
for stix=stnixen(:)';
stnm = cstnms{stix};
spt(2,ncols,colix);
hist(mdl.(stnm).l.data,min_hist:0.02:max_hist);
xlim([min_hist,max_hist]);
xlabel([upper(stnm),' alongshore (m/s)']);
spt(2,ncols,colix+ncols);
hist(mdl.(stnm).x.data,min_hist:0.02:max_hist);
xlim([min_hist,max_hist]);
xlabel([upper(stnm),' cross-shore (m/s)']);
colix = colix + 1;
end;
xlim default
xlim([min_hist,max_hist]);
xlim default
xlim([min_hist,max_hist]);
fmg;
stnixen = [1,2,4,5];
ncols = numel(stnixen);
colix = 1;
%min_hist = floor(-max_spd); max_hist = ceil(+max_spd); d_hist = 0.025;
min_hist = -0.2; max_hist = +0.2; d_hist = 0.02;
for stix=stnixen(:)';
stnm = cstnms{stix};
spt(2,ncols,colix);
hist(mdl.(stnm).l.data,min_hist:d_hist:max_hist);
xlim([min_hist-d_hist,max_hist+d_hist]);
xlabel([upper(stnm),' alongshore (m/s)']);
spt(2,ncols,colix+ncols);
hist(mdl.(stnm).x.data,min_hist:d_hist:max_hist);
xlim([min_hist-d_hist,max_hist+d_hist]);
xlabel([upper(stnm),' cross-shore (m/s)']);
colix = colix + 1;
end;
mdl
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
doPrint = true;
cmpefkeys
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
doPrint = true;
doFKEYS = true;
cmpefkeys
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
cmpefkeys
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
help plot_bathy_transect
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
cmpefkeys
set(gca,'clim',[0,1])
help station_reorient_field
help station_reorient_vectors
mdl
help station_reorient_vectors
help station_reorient_field
x = station_reorient_vectors(mdl,stn.isobath_orientation,mdl.u,mdl.v)
x = station_reorient_vectors(mdl,stn.isobath_orientation,'u','v')
mdl.x
mdl.l
foo = station_reorient_vectors(mdl,stn.isobath_orientation,'u','v','x','l')
help reorient_vectors
[foo.x,foo.l] = reorient_vectors(stn.isobath_orientation,mdl.u,mdl.v)
foo
size(foo.x)
mdl.u
size(mdl.u)
size(foo.x)
foo=[]; clear foo
res
u = squeeze(nanmean(mdl.u(:,:,ix)));
v = squeeze(nanmean(mdl.v(:,:,ix)));
disp([nanmin(u(:)),nanmax(u(:))]);
disp([nanmin(v(:)),nanmax(v(:))]);
[x,l] = reorient_vectors(stn.isobath_orientation,u,v)
disp([nanmin(l(:)),nanmax(l(:))]);
disp([nanmin(x(:)),nanmax(x(:))]);
size(latix)
lonix = interp1(mdl.lon,1:mdl.nlon,res.lon,'nearest');
latix = interp1(mdl.lat,1:mdl.nlat,res.lat,'nearest');
ix = sub2ind([size(mdl.t,3),size(mdl.t,4)],latix,lonix);
%t = squeeze(mdl.t(1,:,ix));
t = squeeze(nanmean(mdl.t(:,:,ix)));
%disp([nanmin(t(:)),nanmax(t(:))]);
u = squeeze(nanmean(mdl.u(:,:,ix)));
v = squeeze(nanmean(mdl.v(:,:,ix)));
disp([nanmin(u(:)),nanmax(u(:))]);
disp([nanmin(v(:)),nanmax(v(:))]);
[x,l] = reorient_vectors(stn.isobath_orientation,u,v)
disp([nanmin(l(:)),nanmax(l(:))]);
disp([nanmin(x(:)),nanmax(x(:))]);
lonix = interp1(mdl.lon,1:mdl.nlon,res.lon,'nearest');
latix = interp1(mdl.lat,1:mdl.nlat,res.lat,'nearest');
ix = sub2ind([size(mdl.t,3),size(mdl.t,4)],latix,lonix);
%t = squeeze(mdl.t(1,:,ix));
t = squeeze(nanmean(mdl.t(:,:,ix)));
%disp([nanmin(t(:)),nanmax(t(:))]);
u = squeeze(nanmean(mdl.u(:,:,ix)));
v = squeeze(nanmean(mdl.v(:,:,ix)));
disp([nanmin(u(:)),nanmax(u(:))]);
disp([nanmin(v(:)),nanmax(v(:))]);
[x,l] = reorient_vectors(stn.isobath_orientation,u,v);
disp([nanmin(l(:)),nanmax(l(:))]);
disp([nanmin(x(:)),nanmax(x(:))]);
whos
stn
mdl.fwyf1
whos
help whos
x
foo = mdl.fwyf1;
foo
whos
foo=[]; clear foo
who
whos
6*676e3
6*676e3/1e6
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
cmpefkeys
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
doPrint = true;
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
cmpefkeys
reviewanim([],0,0,0)
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
doFKEYS = true;
cmpefkeys
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
cmpefkeys
reviewanim([],0,0,0)
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
dirs
pd
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
pd +9
doPrint = false;
analyze_face_nutrients
rev
rev
reviewanim([],0,0,0)
%-- 5/30/2016 9:33 PM --%
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
pwd
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
cN = 'NF09';
face.(cN).fname = 'c:\users\gramer\Documents\Coral\FACE\Master_Data_Sheet_01Feb10.xls';
xl = importdata(face.(cN).fname);
sht = 'Sheet3';
endoff = 0; dtendoff = endoff+21;
face.(cN).dt = datenum(xl.textdata.(sht)(4:end-dtendoff,1))+xl.data.(sht)(1:end-endoff,1);
face.(cN).lat = xl.data.(sht)(1:end-endoff,3)+(xl.data.(sht)(1:end-endoff,4)./60);
face.(cN).lon = -( xl.data.(sht)(1:end-endoff,5)+(xl.data.(sht)(1:end-endoff,6)./60) );
face.(cN).z = -xl.data.(sht)(1:end-endoff,7);
face.(cN).t = xl.data.(sht)(1:end-endoff,8);
face.(cN).nh4 = xl.data.(sht)(1:end-endoff,15);
badix = find(isnan(face.(cN).z) | isnan(face.(cN).t));
face.(cN).dt(badix) = [];
face.(cN).lat(badix) = [];
face.(cN).lon(badix) = [];
face.(cN).z(badix) = [];
face.(cN).t(badix) = [];
face.(cN).nh4(badix) = [];
face.(cN).den = sw_dens(face.(cN).s,face.(cN).t,-face.(cN).z);
face.(cN).dena = face.(cN).den - 1000;
xl=[]; clear xl sht endoff dtendoff
zmin   = -290;             zmax   = 0;
tmin   = 7.9;              tmax   = 29.1;
nh4min = 0;                nh4max = 20;
dmin   = 0;                dmax   = 27;
scatter_fit(face.(cN).t,face.(cN).nh4,[cN,' T'],'NH_4 \mumol');
axis([tmin,tmax,nmin,nmax]); legend('Location','NorthEast');
scatter_fit(face.(cN).t,face.(cN).no2,[cN,' T'],'NO_2 \mumol');
axis([tmin,tmax,no2min,no2max]); legend('Location','NorthEast');
xl=[]; clear xl sht endoff dtendoff
zmin   = -290;             zmax   = 0;
tmin   = 7.9;              tmax   = 29.1;
nh4min = 0;                nh4max = 20;
dmin   = 0;                dmax   = 27;
scatter_fit(face.(cN).t,face.(cN).nh4,[cN,' T'],'NH_4 \mumol');
axis([tmin,tmax,nmin,nmax]); legend('Location','NorthEast');
scatter_fit(face.(cN).t,face.(cN).no2,[cN,' T'],'NO_2 \mumol');
axis([tmin,tmax,no2min,no2max]); legend('Location','NorthEast');
ylim([0,2])
ylim([0,1.4])
print('-dpng','NF09-T24-NH4-scatter.png');
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
dirs
pd
pd +7
dir *ef_*m
type ef_eval_ool_250m.m
stn = ef_eval_ool_250m('manm1','CI',false,true,false,false,true);
stn
stn.CI.ef_events
stn.CI
stn.ww3
stn.CI.ef_events
dataset='CI';
fmg; contourf(stn.(dataset).lon,stn.(dataset).lat,stn.(dataset).ef_events.field_n);
plot_hires_coastline(stn.ngdc_hires_bathy); daspect([1,cosd(stn.lat),1]);
plot(stn.lon,stn.lat,'kp','MarkerFaceColor','w');
titlename('Event days (# per pixel)'); colorbar; axis(ef_roi);
stn.(dataset).ef_events
stn.(dataset)
scatter_fit_ts(stn.(dataset).depth(:),stn.(dataset).ef_events.field_n(:))
scatter_fit(stn.(dataset).depth(:),stn.(dataset).ef_events.field_n(:))
stn=[]; clear stn ans
stn = ef_eval_ool_250m('pomf1','CI',false,true,false,false,true);
scatter_fit(stn.(dataset).depth(:),stn.(dataset).ef_events.field_n(:))
fmg; contourf(stn.(dataset).lon,stn.(dataset).lat,stn.(dataset).ef_events.field_n);
plot_hires_coastline(stn.ngdc_hires_bathy); daspect([1,cosd(stn.lat),1]);
plot(stn.lon,stn.lat,'kp','MarkerFaceColor','w');
titlename('Event days (# per pixel)'); colorbar; axis(ef_roi);
stn=[]; clear stn ans
stn = ef_eval_ool_250m('pvgf1','CI',false,true,false,false,true);
fmg; contourf(stn.(dataset).lon,stn.(dataset).lat,stn.(dataset).ef_events.field_n);
plot_hires_coastline(stn.ngdc_hires_bathy); daspect([1,cosd(stn.lat),1]);
plot(stn.lon,stn.lat,'kp','MarkerFaceColor','w');
titlename('Event days (# per pixel)'); colorbar;
scatter_fit(stn.(dataset).depth(:),stn.(dataset).ef_events.field_n(:))
stn=[]; clear stn ans
axis([-300,0,0,100]); xlabel('Depth [m]'); ylabel('S/RI'); print('-dpng','pvgf1-depth-scatter-sri.png');
axis([-300,20,0,100]); xlabel('Depth [m]'); ylabel('S/RI'); print('-dpng','pvgf1-depth-scatter-sri.png');
axis([-300,20,0,100]); xlabel('Depth [m]'); ylabel('S/RI'); print('-dpng','pomf1-depth-scatter-sri.png');
axis([-300,20,0,100]); xlabel('Depth [m]'); ylabel('S/RI'); print('-dpng','manm1-depth-scatter-sri.png');
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
quit
%-- 6/1/2016 1:53 PM --%
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
stn = ef_eval_ool_250m('manm1','CI',false,true,false,false,true);
ef_eval_ool_250m_do_figs
close all
ef_eval_ool_250m_do_figs
close all
clear ans ix ixen maxci maxgraphs minci selectixen
ef_eval_ool_250m_do_figs
stn.CI.ef_events
clear ans ix ixen maxci maxgraphs minci selectixen
ef_eval_ool_250m_run_efs
ef_eval_ool_250m_disp_figs
stn.CI.ef_events
fmg; contourf(stn.(dataset).lon,stn.(dataset).lat,squeeze(stn.(dataset).ef_events(835,:,:)));
fmg; contourf(stn.(dataset).lon,stn.(dataset).lat,squeeze(stn.(dataset).ef_events.field(835,:,:)));
stn
timenow
close all
clear ans ix ixen maxci maxgraphs minci selectixen
clear evtix latix lonix ndts sat_extreme_pct sat_min_pct sz2d sz3d wave_extreme_pct wave_max_pct wix
ef_eval_ool_250m_run_efs
ef_eval_ool_250m_disp_figs
reviewanim([],0,0,0)
close all
clear ans cmdFile dt evtix ix ixen maxci maxgraphs minci selectixen
cmdFileName='default';
help fopen
islogical(1)
isnumeric(true)
numel(true)
isscalar(true)
isscalar('a')
isscalar('ab')
cmdFileName=true;
ef_eval_ool_250m_disp_figs
reviewanim([],0,0,0)
close all
stn
close all
ef_eval_ool_250m_disp_figs
close all
ef_eval_ool_250m_disp_figs
close all
ef_eval_ool_250m_disp_figs
close all
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
quit
%-- 6/3/2016 5:55 PM --%
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
showSynoptic
showSynoptic=false;
stn = ef_eval_ool_250m('manm1','CI',false,true,false,false,true);
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
quit
%-- 6/3/2016 7:26 PM --%
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
pd ../../../LJGramer_papers\Gramer_2
x = importdata('Xaymara FL Keys sites for Lew.xlsx');
x
x.textdata
x.textdata.FLKeysSitesMcav
edit map_xaymara_sites.m
3392-2711
2938-ans
x.textdata.FLKeysSitesMcav
x.data.FLKeysSitesMcav
mstns = repmat(struct,[nstns,1]);
nmstns = size(x.data.FLKeysSitesMcav,1);
mstns = repmat(struct,[nstns,1]);
mstns = repmat(struct,[nmstns,1]);
mstns(3)
npstns = size(x.data.FLKeysSitesPast,1);
pstns = repmat(struct,[nstns,1]);
nmstns = size(x.data.FLKeysSitesMcav,1);
mstns = repmat(struct,[nmstns,1]);
npstns = size(x.data.FLKeysSitesPast,1);
pstns = repmat(struct,[npstns,1]);
clear ans
mstns(:).station_name = x.textdata.FLKeysSitesMcav{2:end,4};
mstns(:).station_desc = x.textdata.FLKeysSitesMcav{2:end,3};
mstns(:).lon = x.data.FLKeysSitesMcav{:,3};
mstns(:).lat = x.data.FLKeysSitesMcav{:,2};
mstns(:).depth = x.data.FLKeysSitesMcav{:,1};
for stix = 1:numel(mstns)
mstns(stix).station_name = x.textdata.FLKeysSitesMcav{stix+1,4};
mstns(stix).station_desc = x.textdata.FLKeysSitesMcav{stix+1,3};
mstns(stix).lon = x.data.FLKeysSitesMcav{stix,3};
mstns(stix).lat = x.data.FLKeysSitesMcav{stix,2};
mstns(stix).depth = x.data.FLKeysSitesMcav{stix,1};
end;
mstns(stix)
for stix = 1:numel(mstns)
mstns(stix).station_name = x.textdata.FLKeysSitesMcav{stix+1,4};
mstns(stix).station_desc = x.textdata.FLKeysSitesMcav{stix+1,3};
mstns(stix).lon = x.data.FLKeysSitesMcav(stix,3);
mstns(stix).lat = x.data.FLKeysSitesMcav(stix,2);
mstns(stix).depth = x.data.FLKeysSitesMcav(stix,1);
end;
mstns
mstns(9)
for stix = 1:numel(pstns)
pstns(stix).station_name = x.textdata.FLKeysSitesPast{stix+1,4};
pstns(stix).station_desc = x.textdata.FLKeysSitesPast{stix+1,3};
pstns(stix).lon = x.data.FLKeysSitesPast(stix,3);
pstns(stix).lat = x.data.FLKeysSitesPast(stix,2);
pstns(stix).depth = x.data.FLKeysSitesPast(stix,1);
end;
pstns(stix)
nmstns = size(x.data.FLKeysSitesMcav,1);
mstns = repmat(struct,[nmstns,1]);
npstns = size(x.data.FLKeysSitesPast,1);
pstns = repmat(struct,[npstns,1]);
for stix = 1:nmstns
mstns(stix).station_name = x.textdata.FLKeysSitesMcav{stix+1,4};
mstns(stix).station_desc = x.textdata.FLKeysSitesMcav{stix+1,3};
mstns(stix).lon = x.data.FLKeysSitesMcav(stix,3);
mstns(stix).lat = x.data.FLKeysSitesMcav(stix,2);
mstns(stix).depth = x.data.FLKeysSitesMcav(stix,1);
lon(stix,1) = mstns(stix).lon;
lat(stix,1) = mstns(stix).lat;
end;
for stix = 1:npstns
pstns(stix).station_name = x.textdata.FLKeysSitesPast{stix+1,4};
pstns(stix).station_desc = x.textdata.FLKeysSitesPast{stix+1,3};
pstns(stix).lon = x.data.FLKeysSitesPast(stix,3);
pstns(stix).lat = x.data.FLKeysSitesPast(stix,2);
pstns(stix).depth = x.data.FLKeysSitesPast(stix,1);
lon(nmstns+stix,1) = pstns(stix).lon;
lat(nmstns+stix,1) = pstns(stix).lat;
end;
x=[]; clear x
clear ans
distance_wgs84(min(lat),min(lon),max(lat),max(lon))
distance_wgs84(lat(1:end-1),lon(1:end-1),lat(2:end),lon(2:end))
fmg; hist(distance_wgs84(lat(1:end-1),lon(1:end-1),lat(2:end),lon(2:end)))
fmg; hist(distance_wgs84(lat(1:end-1),lon(1:end-1),lat(2:end),lon(2:end)),100)
min(diff(lon))
max(diff(lon))
min(diff(lat))
max(diff(lat))
dx = distance_wgs84(lat(1:end-1),lon(1:end-1),lat(2:end),lon(2:end))
slat
slon
lon
[slat,slatix] = sort(lat); slon = lon(slatix); clear slatix
dx = distance_wgs84(slat(1:end-1),slon(1:end-1),slat(2:end),slon(2:end))
slat([1,2]),slon([1,2])
help sort
[slat,slatix] = sort(lat); slon = lon(slatix);
slatix
dx = distance_wgs84(slat(1:end-1),slon(1:end-1),slat(2:end),slon(2:end))
slat([1,2]),slon([1,2])
[slat,slatix] = sort(lat); slon = lon(slatix);
lat
slat([1,2]),slon([1,2])
slatix
lat([17,42])
dx = distance_wgs84(slat(1:end-1),slon(1:end-1),slat(2:end),slon(2:end))
distance_wgs84(min(lat),min(lon),max(lat),max(lon))
distance_wgs84(min(lat),min(lon),max(lat),min(lon))
distance_wgs84(min(lat),min(lon),min(lat),max(lon))
sqrt( (276^2)+(61^2) )
help plot_hires_bathymetry
stn.lon = mean(lon); stn.lat = mean(lat),
stn = plot_hires_bathymetry(stn,-[0:5:50],[145e3,32e3])
plot(lon,lat,'kp','MarkerFaceColor','w');
+-
stn=[]; clear stn ans
stn.lon = mean(lon); stn.lat = mean(lat);
stn = plot_hires_bathymetry(stn,-[0:5:50],[155e3,35e3])
plot(lon,lat,'kp','MarkerFaceColor','w');
[slon,slonix] = sort(lon); slat = lat(slonix);
dx = distance_wgs84(slat(1:end-1),slon(1:end-1),slat(2:end),slon(2:end))
[ig,ix] = min(lon)
distance_wgs84(lat(ix),lon(ix),stn.lat,stn.lon)
[ig,ix] = max(lon)
distance_wgs84(lat(ix),lon(ix),stn.lat,stn.lon)
stn=[]; clear stn ans
pack
stn.lon = mean([min(lon),max(lon)]); stn.lat = mean([min(lat),max(lat)]);
distance_wgs84(lat(ix),lon(ix),stn.lat,stn.lon)
[ig,ix] = min(lon)
distance_wgs84(lat(ix),lon(ix),stn.lat,stn.lon)
stn.lon = mean([min(lon),max(lon)]);
stn.lat = mean([min(lat),max(lat)]);
[ig,minix] = min(lon);
[ig,maxix] = max(lon);
dx = ceil(max([distance_wgs84(lat(minix),lon(minix),stn.lat,stn.lon),...
distance_wgs84(lat(maxix),lon(maxix),stn.lat,stn.lon)]));
[ig,minix] = min(lat);
[ig,maxix] = max(lat);
dy = ceil(max([distance_wgs84(lat(minix),lon(minix),stn.lat,stn.lon),...
distance_wgs84(lat(maxix),lon(maxix),stn.lat,stn.lon)]));
dx
dy
stn.lon = mean([min(lon),max(lon)]);
stn.lat = mean([min(lat),max(lat)]);
[ig,minix] = min(lon);
[ig,maxix] = max(lon);
dx = ceil(max([distance_wgs84(lat(minix),lon(minix),stn.lat,lon(minix)),...
distance_wgs84(lat(maxix),lon(maxix),stn.lat,lon(maxix))]));
[ig,minix] = min(lat);
[ig,maxix] = max(lat);
dy = ceil(max([distance_wgs84(lat(minix),lon(minix),lat(minix),stn.lon),...
distance_wgs84(lat(maxix),lon(maxix),lat(maxix),stn.lon)]));
dx,dy,
stn.lon = mean([min(lon),max(lon)]);
stn.lat = mean([min(lat),max(lat)]);
[ig,minix] = min(lon);
[ig,maxix] = max(lon);
dy = ceil(max([distance_wgs84(lat(minix),lon(minix),stn.lat,lon(minix)),...
distance_wgs84(lat(maxix),lon(maxix),stn.lat,lon(maxix))]));
[ig,minix] = min(lat);
[ig,maxix] = max(lat);
dx = ceil(max([distance_wgs84(lat(minix),lon(minix),lat(minix),stn.lon),...
distance_wgs84(lat(maxix),lon(maxix),lat(maxix),stn.lon)]));
dx,dy,
stn = plot_hires_bathymetry(stn,-[0:5:50],[dx,dy])
help plot_hires_bathymetry
stn = plot_hires_bathymetry(stn,-[0:5:50],[dx,dy],true,[],[],[],false)
stn=[]; clear stn ans
stn.lon = mean([min(lon),max(lon)]);
stn.lat = mean([min(lat),max(lat)]);
[ig,minix] = min(lon);
[ig,maxix] = max(lon);
dy = ceil(max([distance_wgs84(lat(minix),lon(minix),stn.lat,lon(minix)),...
distance_wgs84(lat(maxix),lon(maxix),stn.lat,lon(maxix))]))*1e3;
[ig,minix] = min(lat);
[ig,maxix] = max(lat);
dx = ceil(max([distance_wgs84(lat(minix),lon(minix),lat(minix),stn.lon),...
distance_wgs84(lat(maxix),lon(maxix),lat(maxix),stn.lon)]))*1e3;
dx,dy,
stn = plot_hires_bathymetry(stn,-[0:5:50],[dx,dy],true,[],[],[],false)
plot(lon,lat,'kp','MarkerFaceColor','w');
mstns(1) = plot_hires_bathymetry(mstns(1),-[0:1:30],[1.5e3,1.5e3],true,[],[],[],true)
x
x = plot_hires_bathymetry(mstns(1),-[0:1:30],[1.5e3,1.5e3],true,[],[],[],true)
whos
x = plot_hires_bathymetry(mstns(1),-[0:1:30],[2e3,2e3],true,[],[],[],true)
x = plot_hires_bathymetry(mstns(1),-[0:0.5:30],[2e3,2e3],true,[],[],[],true)
x = plot_hires_bathymetry(mstns(1),-[0:0.5:30],[2e3,2e3],true,[],[],[],true); set(gca,'CLim',-[0,30]);
x = plot_hires_bathymetry(mstns(1),-[0:0.5:30],[2e3,2e3],true,[],[],[],true); set(gca,'CLim',[-30,0]);
whos
22668*48
ans/1e3
dx = distance_wgs84(slat(1:end-1),slon(1:end-1),slat(2:end),slon(2:end))
dx = distance_wgs84(slat(1:end-1),slon(1:end-1),slat(2:end),slon(2:end));
median(dx(dx>0))
x = plot_hires_bathymetry(mstns(1),-[0:0.5:30],[4e3,4e3],true,[],[],[],true); set(gca,'CLim',[-30,0]);
x = plot_hires_bathymetry(mstns(1),-[0:0.5:30],[4e3,4e3],true,[],[],[],true); set(gca,'CLim',[-30,0]); axis(axis);
plot(lon,lat,'kp','MarkerFaceColor','w');
plot(lon,lat,'kp','MarkerFaceColor','r');
axis(axis);
plot(lon,lat,'kp','MarkerFaceColor','r');
plot(lon,lat,'rp','MarkerFaceColor','r');
plot(lon,lat,'wp','MarkerFaceColor','w');
x = plot_hires_bathymetry(mstns(1),-[0:0.5:30],[4e3,4e3],true,[],[],[],true); set(gca,'CLim',[-30,0]); axis(axis);
plot(lon,lat,'wp','MarkerFaceColor','w');
help read_hires_bathymetry
12.27/2.84
9.259e-5*111e3
11557*10/111e3
(11557*10/111e3) + 25.25
(10909*10/111e3) + (-80.41)
mstns
stn
mstns.ngdc_hires_bathy = [];
mstns(:).ngdc_hires_bathy = [];
mstns(1).ngdc_hires_bathy = [];
stn
stn = plot_hires_bathymetry(stn,-[0:5:50,100:100:800],[dx,dy],true,[],[],[],false);
axis(axis);
dx
dx = ceil(max([distance_wgs84(lat(minix),lon(minix),lat(minix),stn.lon),...
distance_wgs84(lat(maxix),lon(maxix),lat(maxix),stn.lon)]))*1e3;
stn = plot_hires_bathymetry(stn,-[0:5:50,100:100:800],[dx,dy],true,[],[],[],false);
axis(axis);
plot(lon,lat,'ws','MarkerFaceColor','w');
stn=[]; clear stn ans
stn = plot_hires_bathymetry(stn,-[0:5:50,100:100:800],[dx*1.05,dy*1.10],true,[],[],[],false);
axis(axis);
plot(lon,lat,'ws','MarkerFaceColor','w');
stn.lon = mean([min(lon),max(lon)]);
stn.lat = mean([min(lat),max(lat)]);
[ig,minix] = min(lon);
[ig,maxix] = max(lon);
dy = ceil(max([distance_wgs84(lat(minix),lon(minix),stn.lat,lon(minix)),...
distance_wgs84(lat(maxix),lon(maxix),stn.lat,lon(maxix))]))*1e3;
[ig,minix] = min(lat);
[ig,maxix] = max(lat);
dx = ceil(max([distance_wgs84(lat(minix),lon(minix),lat(minix),stn.lon),...
distance_wgs84(lat(maxix),lon(maxix),lat(maxix),stn.lon)]))*1e3;
% Wide-view map
stn = plot_hires_bathymetry(stn,-[0:5:50,100:100:800],[dx*1.05,dy*1.10],true,[],[],[],false);
axis(axis);
plot(lon,lat,'ws','MarkerFaceColor','w');
stn = plot_hires_bathymetry(stn,-[0:5:50,100:100:1000],[dx*1.05,dy*1.10],true,[],[],[],false);
axis(axis);
plot(lon,lat,'ws','MarkerFaceColor','w');
plot(lon,lat,'k.');
stn = plot_hires_bathymetry(stn,-[0:5:50,100:100:800],[dx*1.05,dy*1.10],true,[],[],[],false);
axis(axis);
plot(lon,lat,'ws','MarkerFaceColor','w');
plot(lon,lat,'k.','MarkerSize',0.5);
mstns(1).ngdc_depth = [];
mstns(1).ngdc_hires_bathy = [];
for stix = 1:nmstns
mstns(stix).ngdc_hires_bathy = plot_hires_bathymetry(mstns(1),-[0:1:30],[4e3,4e3],true,[],[],[],true);
set(gca,'CLim',[-30,0]);
axis(axis);
plot(lon,lat,'kp','MarkerFaceColor','w');
end;
stix
mstns = rmfield(mstns,'ngdc_hires_bathy')
mstns = rmfield(mstns,'ngdc_depth')
stix
mstns
x
x=[]; clear x
x = plot_hires_bathymetry(mstns(1),-[0:1:30],[4e3,4e3],true,[],[],[],true);
set(gca,'CLim',[-30,0]);
axis(axis);
plot(lon,lat,'kp','MarkerFaceColor','w');
mstns(stix).ngdc_hires_bathy = x.ngdc_hires_bathy;
mstns
mstns(1).ngdc_hires_bathy
mstns(2).ngdc_hires_bathy
mstns(stix).ngdc_depth = interp2(x.ngdc_hires_bathy.lon,x.ngdc_hires_bathy.lat,x.ngdc_hires_bathy.field,mstns(stix).lon,mstns(stix).lat);
x=[]; clear x
mstns(stix)
for stix = 1:nmstns
x = plot_hires_bathymetry(mstns(1),-[0:1:30],[4e3,4e3],true,[],[],[],true);
set(gca,'CLim',[-30,0]);
axis(axis);
plot(lon,lat,'kp','MarkerFaceColor','w');
mstns(stix).ngdc_hires_bathy = x.ngdc_hires_bathy;
mstns(stix).ngdc_depth = interp2(x.ngdc_hires_bathy.lon,x.ngdc_hires_bathy.lat,x.ngdc_hires_bathy.field,mstns(stix).lon,mstns(stix).lat);
x=[]; clear x
end;
close all
mstns = rmfield(mstns,'ngdc_hires_bathy')
mstns = rmfield(mstns,'ngdc_depth')
for stix = 1:nmstns
x = plot_hires_bathymetry(mstns(stix),-[0:1:30],[4e3,4e3],true,[],[],[],true);
set(gca,'CLim',[-30,0]);
axis(axis);
plot(lon,lat,'kp','MarkerFaceColor','w');
mstns(stix).ngdc_hires_bathy = x.ngdc_hires_bathy;
ndeps(stix) = interp2(x.ngdc_hires_bathy.lon,x.ngdc_hires_bathy.lat,x.ngdc_hires_bathy.field,mstns(stix).lon,mstns(stix).lat);
mstns(stix).ngdc_depth = ndeps(stix);
deps(stix) = mstns(stix).depth;
x=[]; clear x
end;
mstns
for stix = 1:nmstns
x.lon = mstns(stix).lon; x.lat = mstns(stix).lat;
x = plot_hires_bathymetry(x,-[0:1:30],[4e3,4e3],true,[],[],[],true);
set(gca,'CLim',[-30,0]);
axis(axis);
plot(lon,lat,'kp','MarkerFaceColor','w');
mstns(stix).ngdc_hires_bathy = x.ngdc_hires_bathy;
ndeps(stix) = interp2(x.ngdc_hires_bathy.lon,x.ngdc_hires_bathy.lat,x.ngdc_hires_bathy.field,mstns(stix).lon,mstns(stix).lat);
mstns(stix).ngdc_depth = ndeps(stix);
deps(stix) = mstns(stix).depth;
x=[]; clear x
end;
reviewanim([],0,0,0)
scatter_fit(deps,ndeps)
size(deps),size(ndeps)
deps
ndeps
scatter_fit(deps,ndeps')
scatter_fit(deps,ndeps)
fmg; hist(deps-ndeps)
close all
mstns
whos
for stix = 1:nmstns
x.lon = mstns(stix).lon; x.lat = mstns(stix).lat;
%x = plot_hires_bathymetry(x,-[0:1:30],[4e3,4e3],true,[],[],[],true);
x = plot_hires_bathymetry(x,-[0:1:30],[4e3,4e3],true,[],[],[],false);
set(gca,'CLim',[-30,0]);
axis(axis);
plot(lon,lat,'kp','MarkerFaceColor','w');
mstns(stix).ngdc_hires_bathy = x.ngdc_hires_bathy;
ndeps(stix) = interp2(x.ngdc_hires_bathy.lon,x.ngdc_hires_bathy.lat,x.ngdc_hires_bathy.field,mstns(stix).lon,mstns(stix).lat);
mstns(stix).ngdc_depth = ndeps(stix);
deps(stix) = mstns(stix).depth;
x=[]; clear x
end;
close all
lookfor slope
type find_station_ngdc_offshore_slope.m
type station_ngdc_offshore_slope.m
help gradientm
[alp,bet,dY,dX] = gradientm(stn.ngdc_hires_bathy.lat,stn.ngdc_hires_bathy.lon,stn.ngdc_hires_bathy.field,wgs84_ellipsoid);
wgs84_ellipsoid
[alp,bet,dY,dX] = gradientm(stn.ngdc_hires_bathy.lat,stn.ngdc_hires_bathy.lon,stn.ngdc_hires_bathy.field);
dbstop if error\
dbstop if error
[alp,bet,dY,dX] = gradientm(stn.ngdc_hires_bathy.lat,stn.ngdc_hires_bathy.lon,stn.ngdc_hires_bathy.field);
spheroid
dbup
dbquit
stns(1)
mstns(1)
[LAT,LON] = meshgrid(mstns(1).lat,mstns(1).lon);
[alp,bet,dY,dX] = gradientm(LAT,LON,stn.ngdc_hires_bathy.field);
dbquit
[alp,bet,dY,dX] = gradientm(LAT,LON,stn.ngdc_hires_bathy.field);
dbup
nargin,size1, size2, size(varargin{3})
size1
dbup
nargin
varargin
dbup
dbquit
dbstop gradientm
[alp,bet,dY,dX] = gradientm(LAT,LON,stn.ngdc_hires_bathy.field);
size1
varargin
dbquit
[LAT,LON] = meshgrid(mstns(1).ngdc_hires_bathy.lat,mstns(1).ngdc_hires_bathy.lon);
[alp,bet,dY,dX] = gradientm(LAT,LON,mstns(1).ngdc_hires_bathy.field);
size
size1
size2
size(varargin{3})
nargin
dbquit
[LON,LAT] = meshgrid(mstns(1).ngdc_hires_bathy.lon,mstns(1).ngdc_hires_bathy.lat);
[alp,bet,dY,dX] = gradientm(LAT,LON,mstns(1).ngdc_hires_bathy.field);
fmg; contour(bet);
fmg; contourf(lon,lat,bet); colorbar;
dbquit
dbclear all
fmg; contourf(LON,LAT,bet); colorbar;
axis(axis);
plot(lon,lat,'kp','MarkerFaceColor','w');
fmg; contourf(LON,LAT,alp); colorbar;
help gradientm
nansummary(alp)
help wgs84Ellipsoid
help wgs84_ellipsoid
wgs84Ellipsoid
wgs84_ellipsoid
format long
wgs84_ellipsoid
x=wgs84_ellipsoid; x(2)
format short compact
format short
format compact
help wgs84_ellipsoid
wgs84Ellipsoid
[alp,bet,dY,dX] = gradientm(LAT,LON,mstns(1).ngdc_hires_bathy.field,wgs84Ellipsoid);
fmg; contourf(LON,LAT,alp); colorbar;
fmg; contourf(LON,LAT,sqrt((dY^2)+(dX^2))); colorbar;
fmg; contourf(LON,LAT,sqrt((dY.^2)+(dX.^2))); colorbar;
get(gca)
plot_hires_coastline(mstns(1))
mstns(1)
mstns(1).ngdc_hires_bathy
plot_hires_coastline(mstns(1).ngdc_hires_bathy)
fmg; plot_hires_coastline(mstns(1).ngdc_hires_bathy);
fmg; plot_hires_bathymetry(mstns(1),-[0:1:50],[],true);
plot_hires_bathymetry(mstns(1),-[0:1:50],[],true);
help plot_hires_bathymetry
plot_hires_bathymetry(mstns(1),-[0:1:50],[],true,@contour);
contourf(LON,LAT,sqrt((dY.^2)+(dX.^2))); colorbar;
set(gca)
help alphamap
set(gca)
set(gcf)
help alphamap
which alphamap
help alphamap
help alpha
plot_hires_bathymetry(mstns(1),-[0:1:50],[],true,@contour);
[c,h]=contourf(LON,LAT,sqrt((dY.^2)+(dX.^2))); colorbar;
alpha(c,0.5)
size(c)
alpha(gca,0.5)
alpha(gcf,0.5)
size(h)
alpha(h,0.5)
help alpha
help alim
al = alim
al = alim([0,0.5])
alim([0,0.5])
plot_hires_bathymetry(mstns(1),-[0:1:50],[],true,@contour);
sind(6)
sind(2)
close all
[alp,bet,dY,dX] = gradientm(LAT,LON,stn.ngdc_hires_bathy.field,wgs84Ellipsoid);
[LON,LAT] = meshgrid(stn.ngdc_hires_bathy.lon,stn.ngdc_hires_bathy.lat);
fmg; contourf(LON,LAT,bet); colorbar;
size(bet),size(LON)
[alp,bet,dY,dX] = gradientm(LAT,LON,stn.ngdc_hires_bathy.field,wgs84Ellipsoid);
size(bet),size(LON)
fmg; contourf(LON,LAT,bet); colorbar;
set(gca,'clim',[0,5])
set(gca,'clim',[0,1])
set(gca,'clim',[0,0.10])
plot_hires_coastline(stn.ngdc_hires_bathy)
daspect([1,cosd(25),1])
set(gca,'clim',[0,0.05])
set(gca,'clim',[0.05,0.10])
fmg; contourf(LON,LAT,alp); colorbar;
daspect([1,cosd(25),1])
nansummary(dX)
fmg; contourf(LON,LAT,sqrt((dY.^2)+(dX.^2))); colorbar;
set(gca,'clim',[0.00,0.10])
set(gca,'clim',[0.00,0.05])
set(gca,'clim',[0.00,0.02])
beta
slop
help sloop
slop = sqrt((dY.^2)+(dX.^2));
slop(slop<0.005)=nan;
fmg; contourf(LON,LAT,slop); colorbar;
fmg; contourf(LON,LAT,slop,[0.005:0.005:0.050]); colorbar;
plot_hires_coastline(stn.ngdc_hires_bathy)
daspect([1,cosd(25),1])
set(gca,'clim',[0.005,0.050])
daspect([1,cosd(25),1])
set(gca,'clim',[0.005,0.050])
nansummary(stn.ngdc_hires_bathy.field)
numel( find(-2>=stn.ngdc_hires_bathy.field & stn.ngdc_hires_bathy.field>=-30)
numel( find(-2>=stn.ngdc_hires_bathy.field & stn.ngdc_hires_bathy.field>=-30) )
slop( find(-2<stn.ngdc_hires_bathy.field | stn.ngdc_hires_bathy.field<-30) )=nan;
numel( find(~isnan(slop)) )
fmg; contour(LON,LAT,slop); colorbar;
plot_hires_coastline(stn.ngdc_hires_bathy)
help gradientm
[alp,bet,dY,dX] = gradientm(LAT,LON,stn.ngdc_hires_bathy.field,wgs84Ellipsoid);
help asp
asp
alp=[]; bet=[]; dY=[]; dX=[]; clear al alp ans bet c dx dX dY h ig ix
alp=[]; bet=[]; dY=[]; dX=[]; clear al alp ans bet c dx dX dY h ig ix maxix minix slat slatix slon
alp=[]; bet=[]; slop=[]; dY=[]; dX=[]; clear al alp ans bet c dx dX dY h ig ix maxix minix slat slatix slon slonix slop stix x
scatter_fit(deps,ndeps)
dbstop if error
scatter_fit(deps,ndeps)
size(y)),size(X)
size(y),size(X)
dbquit
size(deps),size(ndeps)
scatter_fit(deps,ndeps')
size(y),n
p1
dbup
size(y),Stats
size(y),(Stats.yhat)
size(y),size(Stats.yhat)
dbquit
scatter_fit(deps',ndeps')
scatter_fit(deps,ndeps)
dbquit
dbclear all
size(deps)
size(deps(:))
scatter_fit(deps,ndeps)
size(deps)
scatter_fit(deps,ndeps)
scatter_fit(deps',ndeps')
scatter_fit(deps',ndeps)
help read_hires_bathymetry
help plot_hires_bathymetry
whos
for stix = 1:npstns
x.lon = pstns(stix).lon; x.lat = pstns(stix).lat;
%x = plot_hires_bathymetry(x,-[0:1:30],[4e3,4e3],true,[],[],[],true);
%x = plot_hires_bathymetry(x,-[0:1:30],[4e3,4e3],true,[],[],[],false);
x = read_hires_bathymetry(x,[4e3,4e3],[],false);
set(gca,'CLim',[-30,0]);
axis(axis);
plot(lon,lat,'kp','MarkerFaceColor','w');
ndeps(nmstns+stix) = interp2(x.ngdc_hires_bathy.lon,x.ngdc_hires_bathy.lat,x.ngdc_hires_bathy.field,pstns(stix).lon,pstns(stix).lat);
pstns(stix).ngdc_hires_bathy = x.ngdc_hires_bathy;
pstns(stix).ngdc_depth = ndeps(stix);
deps(nmstns+stix) = pstns(stix).depth;
x=[]; clear x
end;
scatter_fit(deps,ndeps)
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
help pdep
help mdep
help pdeps
help mdeps
mlon
mlat
mlats
mlons
dir
po
which extract_kuffner_hudson
map_xaymara_sites
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
rawin
help save
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
map_xaymara_sites
rawin
rawin.textdata
rawin.textdata.FLKeysSitesMcav
rawin.textdata.FLKeysSitesPast
dbquit
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
help load
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
map_xaymara_sites
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
map_xaymara_sites
deps
epsmdeps
mdeps
pdeps
reviewanim([],0,0,0)
help gradientm
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
map_xaymara_sites
stn
mlrf1 = get_station_from_station_name('mlrf1');
interp2(stn.ngdc_hires_bathy.lon,stn.ngdc_hires_bathy.lat,stn.beta,mlrf1.lon,mlrf1.lat)
mlrf1 = get_station_from_station_name('mlrf1'); mlrf1 = station_optimal_isobath_orientation(mlrf1);
mlrf1 = get_station_from_station_name('mlrf1'); mlrf1 = station_optimal_isobath_orientation(mlrf1); mlrf1 = station_ngdc_offshore_slope(mlrf1),
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
map_xaymara_sites
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
axis(axis);
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
map_xaymara_sites
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
map_xaymara_sites
close all
x
x.ngdc_hires_bathy
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
map_xaymara_sites
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
map_xaymara_sites
help text
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
map_xaymara_sites
reviewanim([],0,0,0)
close all
help plot_hires_bathymetry
help text
plot_hires_bathymetry(stn.ms(end),-[0:1:50]);
set(gca,'CLim',[-30,0]);
set(gca,'CLim',[-50,0]);
plot_hires_bathymetry(stn.ms(end-2),-[0:1:50]);
plot_hires_bathymetry(stn.ms(end-4),-[0:1:50]);
set(gca,'CLim',[-50,0]);
set(gca,'CLim',[-30,0]);
set(gca,'CLim',[-50,0]);
set(gca,'CLim',[-30,0]);
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
fmg; th=text(1,1,'foo','Color','red');
get(th)
fmg; th=text(.5,.5,'bar','Color','k','HorizontalAlignment','right');
text(.5,.5,'bif','Color','b','HorizontalAlignment','left');
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
fmg; th=text(.5,.5,'bar','Color','k','HorizontalAlignment','right');
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
fmg; th=text(.5,.5,'bar','Color','k','HorizontalAlignment','right');
fmg; th=text(.5,.5,'bar','Color','k','HorizontalAlignment','left');
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
map_xaymara_sites
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
map_xaymara_sites
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
map_xaymara_sites
stns.ms(2)
stn.ms(2)
stns
stn
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
map_xaymara_sites
fmg; th=text(.5,.5,'bar','Color','k','HorizontalAlignment','left'); set(th)
fmg; th=text(.5,.5,'bar','Color','k','HorizontalAlignment','left'); get(th)
360/25
close all
strip(' space')
strips(' space')
help strstack
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
stn = extract_xaymara_sites;
map_xaymara_sites
close all
stn
stn.ps
help despace
lookfor strip
map_xaymara_sites
close all
map_xaymara_sites
fmg; th=text(.5,.5,'bar','Color','k','HorizontalAlignment','left'); get(th)
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
map_xaymara_sites
xxxprint('-dpng','foo.png')
dir('foo.png')
print('-dpng','foo.png')
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
help plot_bath_hc_sites
which plot_bath_hc_sites
delete('foo.png');
dirs
pd
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
stn = extract_xaymara_sites;
for ix=1:numel(stn.plons); if ( isempty(find(stn.mlons==stn.plons(ix) & stn.mlats==stn.plats(ix))) ); disp(stn.ps(ix).station_name); end; end;
for ix=1:numel(stn.plons); if ( isempty(find(abs(stn.mlons-stn.plons(ix))<1e-5 & abs(stn.mlats==stn.plats(ix))<1e-5)) ); disp(stn.ps(ix).station_name); end; end;
for ix=1:numel(stn.plons); if ( isempty(find(abs(stn.mlons-stn.plons(ix))>1e-5 & abs(stn.mlats==stn.plats(ix))>1e-5)) ); disp(stn.ps(ix).station_name); end; end;
for ix=1:numel(stn.plons); if ( ~isempty(find(abs(stn.mlons-stn.plons(ix))>1e-5 & abs(stn.mlats==stn.plats(ix))>1e-5)) ); disp(stn.ps(ix).station_name); end; end;
stn.ms
stn.ms.station_name
for ix=1:numel(stn.plons); if ( ~isempty(find(abs(stn.mlons-stn.plons(ix))>1e-5 & abs(stn.mlats==stn.plats(ix))>1e-5)) ); disp(stn.ps(ix).station_name); end; end;
[stn.ms.station_name]
[stn.ms.station_name']
{stn.ms.station_name'}
{stn.ms.station_name}
for ix=1:numel(stn.plons); if ( ~strfind(stn.ps(ix).station_name,{stn.ms.station_name}) ); disp(stn.ps(ix).station_name); end; end;
help strfind
for ix=1:numel(stn.plons); if ( isempty(strfind({stn.ms.station_name},stn.ps(ix).station_name)) ); disp(stn.ps(ix).station_name); end; end;
for ix=1:numel(stn.ms); if ( isempty(strfind({stn.ps.station_name},stn.ms(ix).station_name)) ); disp(stn.ms(ix).station_name); end; end;
{stn.ms.station_name}
{stn.ps.station_name}
for ix=1:numel(stn.ps); if ( isempty(strfind({stn.ms.station_name},stn.ps(ix).station_name)) ); disp(stn.ps(ix).station_name); end; end;
help strfind
ix
{stn.ps.station_name}
{stn.ms.station_name}
stn.ps(15)
for ix=15:16; if ( isempty(strfind({stn.ms.station_name},stn.ps(ix).station_name)) ); disp(stn.ps(ix).station_name); end; end;
stn.ps(ix)
strfind({stn.ms.station_name},stn.ps(ix).station_name)
isempty(strfind({stn.ms.station_name},stn.ps(ix).station_name))
cellfun(@isempty(strfind({stn.ms.station_name},stn.ps(ix).station_name)))
cellfun(@isempty,strfind({stn.ms.station_name},stn.ps(ix).station_name)))
cellfun(@isempty,strfind({stn.ms.station_name},stn.ps(ix).station_name))
~all(cellfun(@isempty,strfind({stn.ms.station_name},stn.ps(ix).station_name)))
all(cellfun(@isempty,strfind({stn.ms.station_name},stn.ps(ix).station_name)))
all(cellfun(@isempty,strfind({stn.ms.station_name},stn.ps(ix+1).station_name)))
for ix=15:16; if ( all(cellfun(@isempty,strfind({stn.ms.station_name},stn.ps(ix).station_name))) ); disp(stn.ps(ix).station_name); end; end;
timenow
sort({stn.ps.station_name})
help sort
for ix=1:numel(stn.plons); if ( ~isempty(find(abs(stn.mlons-stn.plons(ix))>1e-5 & abs(stn.mlats==stn.plats(ix))>1e-5)) ); disp(stn.ps(ix).station_name); end; end;
stn.ps(23)
{stn.ps.station_name}
stn.ps(20)
for ix=1:numel(stn.plons); msix=find(strcmp(stn.ps(ix).station_name,{stn.ms.station_name})); if ( ~isempty(find(abs(stn.mlons-stn.plons(ix))>1e-5 & abs(stn.mlats==stn.plats(ix))>1e-5)) ); disp(stn.ps(ix).station_name); end; end;
msix=find(ismember({stn.ps.station_name},{stn.ms.station_name}))
{stn.ps(msix).station_name}
{stn.ms.station_name}
stns
stn
timenow
{stn.ms.station_name}
{stn.ps.station_name}
stn.ps.ngdc_beta
[stn.ps.ngdc_beta]
nansummary([stn.ps.ngdc_beta])
nansummary([stn.ms.ngdc_beta])
for ix=1:numel(stn.mlons); if ( ~isempty(find(abs(stn.plons-stn.mlons(ix))>1e-5 & abs(stn.plats-stn.mlats(ix))>1e-5)) ); disp(stn.ps(ix).station_name); end; end;
for ix=1:numel(stn.mlons); if ( ~isempty(find(abs(stn.plons(msix)-stn.mlons(ix))>1e-5 & abs(stn.plats(msix)-stn.mlats(ix))>1e-5)) ); disp(stn.ps(ix).station_name); end; end;
{stn.ps.station_name}
size(stn.mlons),size(stn.plons(msix))
stn.mlons-stn.plons(msix)
stn.mlats-stn.plats(msix)
for ix=1:numel(stn.mlons); if ( ~isempty(find(abs(stn.plons(msix)-stn.mlons(ix))>1e-5 & abs(stn.plats(msix)-stn.mlats(ix))>1e-5)) ); disp(stn.ps(ix).station_name); end; end;
msix=find(ismember({stn.ps.station_name},{stn.ms.station_name}));
for ix=1:numel(stn.mlons);
pix = find( abs(stn.plons(msix)-stn.mlons(ix))>1e-5 & abs(stn.plats(msix)-stn.mlats(ix))>1e-5))
if ( isempty(pix) );
disp(stn.ps(ix).station_name);
end;
end;
pix
msix=find(ismember({stn.ps.station_name},{stn.ms.station_name}));
for ix=1:numel(stn.mlons);
pix = find( abs(stn.plons(msix)-stn.mlons(ix))>1e-5 & ...
abs(stn.plats(msix)-stn.mlats(ix))>1e-5 )
if ( isempty(pix) )
disp(stn.ps(ix).station_name);
end;
end;
pix = find( abs(stn.plons(msix)-stn.mlons(ix))>1e-5 & ...
abs(stn.plats(msix)-stn.mlats(ix))>1e-5 );
msix=find(ismember({stn.ps.station_name},{stn.ms.station_name}));
for ix=1:numel(stn.mlons);
pix = find( abs(stn.plons(msix)-stn.mlons(ix))>1e-5 & ...
abs(stn.plats(msix)-stn.mlats(ix))>1e-5 );
if ( isempty(pix) )
disp(stn.ps(ix).station_name);
end;
end;
pix
msix=find(ismember({stn.ps.station_name},{stn.ms.station_name}));
for ix=1:numel(stn.mlons);
pix = find( abs(stn.plons(msix)-stn.mlons(ix))>1e-6 & ...
abs(stn.plats(msix)-stn.mlats(ix))>1e-6 );
if ( isempty(pix) )
disp(stn.ps(ix).station_name);
end;
end;
msix=find(ismember({stn.ps.station_name},{stn.ms.station_name}));
for ix=1:numel(stn.mlons);
pix = find( abs(stn.plons(msix)-stn.mlons(ix))<1e-5 & ...
abs(stn.plats(msix)-stn.mlats(ix))<1e-5 );
if ( isempty(pix) )
disp(stn.ps(ix).station_name);
end;
end;
onem
help onem
msix=find(ismember({stn.ps.station_name},{stn.ms.station_name}));
min_d = 1/111e3;
for ix=1:numel(stn.mlons);
pix = find( abs(stn.plons(msix)-stn.mlons(ix))<min_d & ...
abs(stn.plats(msix)-stn.mlats(ix))<min_d );
if ( isempty(pix) )
disp(stn.ps(ix).station_name);
end;
end;
msix=find(ismember({stn.ps.station_name},{stn.ms.station_name}));
min_d = 5/111e3;
for ix=1:numel(stn.mlons);
pix = find( abs(stn.plons(msix)-stn.mlons(ix))<min_d & ...
abs(stn.plats(msix)-stn.mlats(ix))<min_d );
if ( isempty(pix) )
disp(stn.ps(ix).station_name);
end;
end;
perr
help gperr
msix=find(ismember({stn.ps.station_name},{stn.ms.station_name}));
min_d = 5/111e3;
for ix=1:numel(stn.mlons);
pix = find( abs(stn.plons(msix)-stn.mlons(ix))<min_d & ...
abs(stn.plats(msix)-stn.mlats(ix))<min_d );
if ( isempty(pix) )
disp(stn.ps(ix).station_name);
end;
end;
msix=find(ismember({stn.ps.station_name},{stn.ms.station_name}));
min_d = 50/111e3;
for ix=1:numel(stn.mlons);
pix = find( abs(stn.plons(msix)-stn.mlons(ix))<min_d & ...
abs(stn.plats(msix)-stn.mlats(ix))<min_d );
if ( isempty(pix) )
disp(stn.ps(ix).station_name);
end;
end;
msix=find(ismember({stn.ps.station_name},{stn.ms.station_name}));
min_d = 500/111e3;
for ix=1:numel(stn.mlons);
pix = find( abs(stn.plons(msix)-stn.mlons(ix))<min_d & ...
abs(stn.plats(msix)-stn.mlats(ix))<min_d );
if ( isempty(pix) )
disp(stn.ps(ix).station_name);
end;
end;
msix=find(ismember({stn.ps.station_name},{stn.ms.station_name}));
min_d = 500/111e3;
for ix=1:numel(stn.mlons);
d = distance_wgs84(stn.mlons(ix).lat,stn.mlons(ix).lon,stn.plons(msix).lat,stn.plons(msix).lon);
if ( all(d > min_d) )
disp({stn.ms(ix).station_name);
end;
end;
msix=find(ismember({stn.ps.station_name},{stn.ms.station_name}));
min_d = 500/111e3;
for ix=1:numel(stn.mlons);
d = distance_wgs84(stn.mlons(ix).lat,stn.mlons(ix).lon,stn.plons(msix).lat,stn.plons(msix).lon);
if ( all(d > min_d) )
disp({stn.ms(ix).station_name,min(d)});
end;
end;
size(d)
d
help distance_wgs84
msix=find(ismember({stn.ps.station_name},{stn.ms.station_name}));
min_d = 500/111e3;
for ix=1:numel(stn.mlons);
d = distance_wgs84(stn.mlats(ix),stn.mlons(ix),stn.plats(msix),stn.plons(msix));
if ( all(d > min_d) )
disp({stn.ms(ix).station_name,min(d)});
end;
end;
{stn.ps.station_name}
d = distance_wgs84(stn.plats(15),stn.plons(15),stn.mlats,stn.mlons);
min(d)
stn.ps(15)
[m,ix]=min(d)
stn.ms(12)
stn
scatter_fit(stn.mbetas,stn.pbetas(msix))
x = get_station_from_station_name('mlrf1'); x=station_optimal_isobath_orientation(x); x=station_ngdc_offshore_slope(x),
interp2(stn.ngdc_hires_bathy.lon,stn.ngdc_hires_bathy.lat,stn.beta,x.lon,x.lat)
x = get_station_from_station_name('mlrf1'); x=station_optimal_isobath_orientation(x); x=station_ngdc_offshore_slope(x), interp2(stn.ngdc_hires_bathy.lon,stn.ngdc_hires_bathy.lat,stn.beta,x.lon,x.lat)
x = get_station_from_station_name('fwyf1'); x=station_optimal_isobath_orientation(x); x=station_ngdc_offshore_slope(x), interp2(stn.ngdc_hires_bathy.lon,stn.ngdc_hires_bathy.lat,stn.beta,x.lon,x.lat)
stn.ngdc_hires_bathy.lat(1)
stn.ngdc_hires_bathy.lat(end)
stn.ngdc_hires_bathy.lon(end)
stn.ngdc_hires_bathy.lon(1)
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
msix=find(ismember({stn.ps.station_name},{stn.ms.station_name}));
min_d = 500/111e3;
for ix=1:numel(stn.mlons);
d = distance_wgs84(stn.mlats(ix),stn.mlons(ix),stn.plats(msix),stn.plons(msix));
% if ( all(d > min_d) ); disp({stn.ms(ix).station_name,min(d)}); end;
[m,mx] = min(d);
disp({stn.ms(ix).station_name,stn.ps(msix(mx)).station_name,m});
end;
stn = extract_xaymara_sites;
msix=find(ismember({stn.ps.station_name},{stn.ms.station_name}));
min_d = 500/111e3;
for ix=1:numel(stn.mlons);
d = distance_wgs84(stn.mlats(ix),stn.mlons(ix),stn.plats(msix),stn.plons(msix));
% if ( all(d > min_d) ); disp({stn.ms(ix).station_name,min(d)}); end;
[m,mx] = min(d);
disp({stn.ms(ix).station_name,stn.ps(msix(mx)).station_name,m});
end;
clear min_d
msix=find(ismember({stn.ps.station_name},{stn.ms.station_name}));
for ix=1:numel(stn.mlons);
d = distance_wgs84(stn.mlats(ix),stn.mlons(ix),stn.plats(msix),stn.plons(msix));
% if ( all(d > 500/111e3) ); disp({stn.ms(ix).station_name,min(d)}); end;
[m,mx] = min(d);
str=''; if ( ~strcmpi(stn.ms(ix).station_name,stn.ps(msix(mx)).station_name) ); str='  ***'; end;
disp({stn.ms(ix).station_name,stn.ps(msix(mx)).station_name,m,str});
end;
nansummary(stn.pbetas)
nansummary(stn.mbetas)
stn
stn.mdeps(stix) = stn.ms(stix).depth;
% clear d ix m mx str
clear d ix m mx str
map_xaymara_sites
close all
25*7
map_xaymara_sites
fmg; th=text(.5,.5,'bar','Color','k','HorizontalAlignment','left'); get(th)
set(gco,'VerticalAlignment','top')
set(gco,'VerticalAlignment','bottom')
set(gco,'VerticalAlignment','middle')
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
stn = extract_xaymara_sites;
x = get_station_from_station_name('fwyf1'); x=station_optimal_isobath_orientation(x); x=station_ngdc_offshore_slope(x), interp2(stn.ngdc_hires_bathy.lon,stn.ngdc_hires_bathy.lat,stn.beta,x.lon,x.lat)
min(distance_wgs84(x.lat,x.lon,stn.ngdc_hires_bathy.lat,stn.ngdc_hires_bathy.lon))
[LON,LAT] = meshgrid(stn.ngdc_hires_bathy.lon,stn.ngdc_hires_bathy.lat);
min(distance_wgs84(x.lat,x.lon,LAT,LON))
min(distance_wgs84(x.lat,x.lon,LAT(:),LON(:))
min(distance_wgs84(x.lat,x.lon,LAT(:),LON(:)))
x
max(LAT(:)),max(LON(:))
min(distance_wgs84(x.lat,x.lon,LAT(:),max(LON(:)))
min(distance_wgs84(x.lat,x.lon,LAT(:),max(LON(:))))
min(distance_wgs84(x.lat,x.lon,max(LAT(:)),LON(:)))
[ig,ix]=min(abs(x.lon-LON(:)))
timenow
min(diff(unique(LON(:))))
x = get_station_from_station_name('smkf1'); x=station_optimal_isobath_orientation(x); x=station_ngdc_offshore_slope(x), interp2(stn.ngdc_hires_bathy.lon,stn.ngdc_hires_bathy.lat,stn.beta,x.lon,x.lat)
x = get_station_from_station_name('smkf1'); x=station_optimal_isobath_orientation(x); x=station_ngdc_offshore_slope(x), interp2(stn.ngdc_hires_bathy.lon,stn.ngdc_hires_bathy.lat,stn.beta,x.lon,x.lat,'nearest')
x = get_station_from_station_name('smkf1'); x=station_optimal_isobath_orientation(x); x=station_ngdc_offshore_slope(x), interp2(stn.ngdc_hires_bathy.lon,stn.ngdc_hires_bathy.lat,stn.beta,x.lon,x.lat,'spline')
x = get_station_from_station_name('sanf1'); x=station_optimal_isobath_orientation(x); x=station_ngdc_offshore_slope(x), interp2(stn.ngdc_hires_bathy.lon,stn.ngdc_hires_bathy.lat,stn.beta,x.lon,x.lat)
help interp2
x = get_station_from_station_name('sanf1'); x=station_optimal_isobath_orientation(x); x=station_ngdc_offshore_slope(x), interp2(stn.ngdc_hires_bathy.lon,stn.ngdc_hires_bathy.lat,stn.beta,x.lon,x.lat)
[ig.xix] = min(abs(stn.ngdc_hires_bathy.lon-x.lon); [ig.yix] = min(abs(stn.ngdc_hires_bathy.lat-x.lat);
[ig.xix] = min(abs(stn.ngdc_hires_bathy.lon-x.lon)); [ig.yix] = min(abs(stn.ngdc_hires_bathy.lat-x.lat))
x
[ig,xix] = min(abs(stn.ngdc_hires_bathy.lon-x.lon)); [ig,yix] = min(abs(stn.ngdc_hires_bathy.lat-x.lat))
[ig,xix] = min(abs(stn.ngdc_hires_bathy.lon-x.lon)), [ig,yix] = min(abs(stn.ngdc_hires_bathy.lat-x.lat))
stn
stn.ngdc_hires_bathy
stn.beta(yix-2:yix+2,xix-2:xix+2)
stn.beta_x(yix-2:yix+2,xix-2:xix+2)
stn.beta_y(yix-2:yix+2,xix-2:xix+2)
stn.beta(yix-2:yix+2,xix-2:xix+2)
x
stn.ngdc_hires_bathy.field(yix-2:yix+2,xix-2:xix+2)
stn.ngdc_hires_bathy.field(xix-2:xix+2,yix-2:yix+2)
fmg; contourf(stn.ngdc_hires_bathy.lon(xix-40:xix+40),stn.ngdc_hires_bathy.lat(yix-30:yix+30),stn.ngdc_hires_bathy.field(yix-30:yix+30,xix-40:xix+40)
fmg; contourf(stn.ngdc_hires_bathy.lon(xix-40:xix+40),stn.ngdc_hires_bathy.lat(yix-30:yix+30),stn.ngdc_hires_bathy.field(yix-30:yix+30,xix-40:xix+40))
xix,yix
fmg; contourf(stn.ngdc_hires_bathy.lon(xix-40:xix+40),stn.ngdc_hires_bathy.lat(yix-25:yix+25),stn.ngdc_hires_bathy.field(yix-25:yix+25,xix-40:xix+40))
fmg; contourf(stn.ngdc_hires_bathy.lon(xix-40:xix+40),stn.ngdc_hires_bathy.lat(yix-25:yix+25),stn.ngdc_hires_bathy.field(yix-25:yix+25,xix-40:xix+40)); colorbar;
plot(x.lon,x.lat,'wp')
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
stn = load_portal_data('srvi2','ALL')
dbstop load_portal_data
stn = load_portal_data('srvi2','ALL')
dbcont
pd ../../coral/ICON
stn = load_portal_data('srvi2','ALL')
fmg; plot_ts(stn.air_t_2_degc,stn.seatemp_shallow_degc,stn.seatemp_deep);
fmg; plot_ts(stn.air_t_2_degc,stn.seatemp_shallow_degc,stn.seatemp_deep_degc);
fmg; plot_ts(stn.air_t_2_degc,stn.seatemp_shallow_degc,stn.seatemp_deep_degc); legend('Air','Shallow','Deep');
fmg; plot_ts(stn.air_t_2_degc,stn.seatemp_deep_degc); legend('Air','Deep');
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
pd
cmpefkeys
close(2:30)
mdl=[]; stn=[];
clear all
pack
doFKEYS=true;
cmpefkeys
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
quit
%-- 6/7/2016 9:47 PM --%
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
pd
stn = load_portal_data('lppr1','ALL');
dbstop if error
stn = load_portal_data('lppr1','ALL');
stn
dir('lppr1*')
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
stn = load_portal_data('lppr1','ALL');
dbstop load_portal_data
stn = load_portal_data('lppr1','ALL');
dbquit
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
stn = load_portal_data('lppr1','ALL');
stn
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
dir('lppr1*')
delete('lppr1_portal_ALL.mat')
stn = load_portal_data('lppr1','ALL');
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
stn = load_portal_data('lppr1','ALL');
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
stn = load_portal_data('LPPR1','ALL');
dbstop load_portal_data
stn = load_portal_data('LPPR1','ALL');
x
dbquit
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
stn = load_portal_data('xLPPR1','ALL');
stn = load_portal_data('LPPR1','ALL');
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
stn = load_portal_data('LPPR1','ALL');
fmg; plot_ts(stn.air_t_2_degc,stn.seatemp_deep_degc); legend('Air','Deep');
fmg; hist(diff(stn.air_t_2_degc.data),100);
fmg; hist(diff(stn.air_t_2_degc.data)./diff(stn.air_t_2_degc.date),100);
dpdt = diff(stn.air_t_2_degc.data)./diff(stn.air_t_2_degc.date);
x = stn.air_t_2_degc; x.date(abs(dpdt)>5)=[]; x.data(abs(dpdt)>5)=[];
fmg; plot_ts(x,stn.seatemp_deep_degc); legend('Air','Deep');
fmg; plot_ts(x,stn.air_t_2_degc);
fmg; plot_ts(stn.air_t_2_degc,x);
x = stn.air_t_2_degc; x.date(abs(dpdt)>5)=[]; x.data(abs(dpdt)>100)=[];
fmg; plot_ts(stn.air_t_2_degc,x);
x = stn.air_t_2_degc; x.date(abs(dpdt)>100)=[]; x.data(abs(dpdt)>100)=[];
fmg; plot_ts(stn.air_t_2_degc,x);
fmg; hist(diff(stn.seatemp_deep_degc.data)./diff(stn.seatemp_deep_degc.date),100);
fmg; hist(diff(stn.seatemp_deep_degc.data)./diff(stn.seatemp_deep_degc.date),1000);
fmg; hist(diff(stn.seatemp_deep_degc.data)./diff(stn.seatemp_deep_degc.date),1000); xlim([-1000,1000]);
fmg; hist(diff(stn.seatemp_deep_degc.data)./diff(stn.seatemp_deep_degc.date),1000); xlim([-200,200]);
fmg; hist(diff(stn.seatemp_deep_degc.data)./diff(stn.seatemp_deep_degc.date),10000); xlim([-200,200]);
fmg; hist(diff(stn.seatemp_deep_degc.data)./diff(stn.seatemp_deep_degc.date),10000); xlim([-20,20]);
dpdt = diff(stn.seatemp_deep_degc.data)./diff(stn.seatemp_deep_degc.date);
x = stn.seatemp_deep_degc; x.date(abs(dpdt)>20)=[]; x.data(abs(dpdt)>20)=[];
fmg; plot_ts(stn.seatemp_deep_degc,x);
get(gca)
lookfor select
getSelectedAnnotations
help getselectobjects
getselectobjects
type getselectobjects
h = getselectobjects(gcf)
help datacursormode
h = datacursormode(gcf)
get(h)
ans.Figure
get(h,'Figure')
get(h)
help datacursormode
help getCursorInfo
help datacursormode
getCursorInfo(h)
hi = getCursorInfo(h)
get(hi,'Target')
hi.Target
edit qcbegdel.m
getCursorInfo(h)
help ginput
[X,Y] = ginput(1)
help gtext
gtext('foo')
gtext('bar')
pause
qcbegdel
hi
h
help datacursormode
qcbegdel
hi
qcbegdel
hi
hi.Position
hi.Position(1)
qcbegdel
begdt
x=[]; dpdt=[]; clear ans begdt dpdt enddt h hi x X Y
x{2} = nan
x(1) = []
clear x
help gcf
get(0,'CurrentFigure')
help isvalid
help ishandle
help ndims
[pi,0]
[pi,0]'
dts = [];
dts = selectpointrange(dts)
ishandle([])
get(0,'CurrentFigure')
fh = get(0,'CurrentFigure')
ishandle(fh)
if ( ishandle(fh) ); disp('N'); end;
dts
dbstop selectpointrange
dts = selectpointrange(dts)
nargs
args
varargin
dbquit
dbclear all
clear ans fh
dts = selectpointrange(dts)
dbstop selectpointrange
dts = selectpointrange(dts)
fh
args
nargs
ismatrix(args{1})
dbquit
dbclear all
dts = selectpointrange(dts)
dbstop selectpointrange
dts = selectpointrange(dts)
nargs
args
xes
args
nargs
fh
ishandle(fh)
~ishandle(fh)
dbquit
dts = selectpointrange(dts)
dbquit
dbclear all
dts = selectpointrange(dts)
fmg; plot_ts(stn.seatemp_deep_degc);
dts = selectpointrange(dts)
dts = selectpointrange(dts); datestr(dts(:,end))
dts
datestr(dts)
dts=[];
help brush
which brush
help brush
info bruch
info brush
doc brush
rms = [732691.913889	22.240000
732691.955556	26.660000
732691.997222	26.670000
732692.038889	26.670000
732692.080556	26.690000
732692.122222	26.690000
732692.163889	26.670000
732692.205556	26.700000
732692.247222	26.700000
732692.288889	26.700000
732692.330556	26.660000
732692.372222	26.550000
732692.413889	26.490000
732692.455556	26.480000
732692.497222	26.470000
732692.538889	26.480000
732692.580556	26.510000
732692.622222	26.510000
732692.663889	26.570000
732692.705556	26.640000
732692.747222	26.670000
732692.788889	14.270000
732692.830556	26.700000
732692.872222	26.720000
732692.913889	26.770000
732692.955556	26.750000
732692.997222	26.720000
732693.038889	26.740000
732693.080556	26.830000
732693.122222	26.900000
732693.163889	26.890000
732693.205556	26.840000
732693.247222	26.810000
732693.288889	26.770000
732693.330556	26.680000
732693.372222	26.580000
732693.413889	26.540000
732693.455556	26.540000
732693.497222	26.490000
732693.538889	26.490000
732693.580556	26.570000
732693.622222	26.630000
732693.663889	26.700000
732693.705556	26.770000
732693.747222	26.780000
732693.788889	26.810000
732693.830556	26.790000
732693.872222	26.790000
732693.913889	26.770000
732693.955556	26.760000
732693.997222	26.820000
732694.038889	27.140000
732694.080556	27.130000
732694.122222	27.060000
732694.163889	27.020000
732694.205556	26.990000
732694.247222	26.890000
732694.288889	26.720000
732694.330556	26.610000
732694.372222	26.570000
732694.413889	26.570000
732694.455556	26.580000
732694.497222	26.570000
732694.538889	26.580000
732694.580556	26.610000
732694.622222	26.630000
732694.663889	26.700000
732694.705556	26.770000
732694.747222	26.780000
732694.788889	26.790000
732694.830556	26.780000
732694.872222	26.790000
732694.913889	26.800000
732694.955556	26.810000
732694.997222	26.810000
732695.038889	26.820000
732695.080556	26.910000
732695.122222	26.900000
732695.163889	26.840000
732695.205556	26.730000
732695.247222	26.650000
732695.288889	26.550000
732695.330556	26.550000
732695.372222	26.580000
732695.413889	26.580000
732695.455556	26.570000
732695.497222	26.600000
732695.538889	26.620000
732695.580556	26.690000
732695.622222	26.760000
732695.663889	26.810000
732695.705556	26.810000
732695.747222	26.800000
732695.788889	26.820000
732695.830556	26.820000
732695.872222	26.840000
732695.913889	26.860000
732695.955556	26.830000
732695.997222	26.830000
732696.038889	26.840000
732696.080556	26.880000
732696.122222	26.890000
732696.163889	26.870000
732696.205556	26.850000
732696.247222	26.850000
732696.288889	26.850000
732696.330556	26.820000
732696.372222	26.800000
732696.413889	26.810000
732696.455556	26.840000
732696.497222	26.870000
732696.538889	26.880000
732696.580556	26.890000
732696.622222	26.930000
732696.663889	26.980000
732696.705556	27.010000
732696.747222	27.030000
732696.788889	27.010000
732696.830556	27.020000
732696.872222	27.050000
732696.913889	27.040000
732696.955556	27.110000
732696.997222	27.170000];
clear rms
help jumpax
jumpax
q
clear rms
jumpax([],[],@(ax)(xlim('default'),datetick3(ax)));
jumpax([],[],@(ax)(xlim('default');datetick3(ax)));
jumpax([],[],@(ax)(xlim('default') && datetick3(ax)));
q
ylim default
jumpax([],[],@(ax)(ylim('default') && datetick3(ax)));
b
help ylim
ylim('default')
help inlineeval
help brushing
help brush
help str2cell
help cellstr
help brush
h = brush
get(h)
hi = get(h)
hi.FigureHandle
get(h.FigureHandle)
help brush
h = brush(gcf,'on')
help jumpax
x = [pi,0,3]
pi(end-1)
x(end-1)
x(end/2)
x = [pi,1,3,exp(1)]
x(end/2)
x(floor(end/2))
help isdatenum(pi)
help isdatetime
isdatetime(now)
help datetime
isvector(pi)
numel([])
isnumeric([])
isvector([])
help plot
clear x hi h dts ans
fmg; plot_ts(stn.seatemp_deep_degc);
xlim([datenum(2005,7,7),inf])
fmg; plot_ts(stn.seatemp_deep_degc);
xlim([-inf,+inf])
xlim default
help ver
help version
v = version
help verLessThan
isvector([pi])
isvector([pi;pi])
isvector([])
inf==inf
inf>0
isinf(inf)
inf >= 1
clear v ans
stn = qcstation(stn)
help jumpax
help ylim
help datetick3
type datetick3
stn = qcstation(stn)
help ylim
ylim('default')
ylim(gca,'default')
help strrep
jumpax([],[],@qcstation_jumpaction)
jumpax([],[],@qcstation>qcstation_jumpaction)
jumpax([],[],@qcstation_jumpaction)
b
q
dbstop if error
stn = qcstation(stn)
ax
exist('ax','var')
dbquit
who
whos
clear ans
help input
strncmpi('','q')
strncmpi(' ','q')
strncmpi('','q',1)
skip
clear ans
stn = qcstation(stn)
b
q
s
skip
stn.(fld)
size(dts)
all(dat==stn.(fld).dat)
all(dat==stn.(fld).data)
all(dat==stn.(fld).data | isnan(stn.(fld).data))
numel(find(dat~=stn.(fld).data))
nansummary(dat)
numel(find(dat~=stn.(fld).date))
dbquit
newstn = qcstation(stn)
q
stn.(fld)
size(dts),size(dat)
all(dat==stn.(fld).data | isnan(stn.(fld).data))
all(dat==stn.(fld).data | (isnan(dat) & isnan(stn.(fld).data)))
dbquit
newstn = qcstation(stn)
q
skipo
skip
size(dts),size(dat)
stn.(fld)
dbquit
inp = input('Enter to save and move to next field, or ''q''uit or ''s''kip\n','s');
newstn = qcstation(stn)
q
size(dts),size(dat)
all(dat==stn.(fld).data | (isnan(dat) & isnan(stn.(fld).data)))
dat(1:10)
dts(1:10)
dbquit
clear inp
newstn = qcstation(stn)
q
ax
get(ax)
get(lh)
ax=lh;
all(dat==stn.(fld).data | (isnan(dat) & isnan(stn.(fld).data)))
stn.(fld)
dts
stn.(fld)
size(dts),size(dat)
all(dat==stn.(fld).data | (isnan(dat) & isnan(stn.(fld).data)))
dbquit
lh=plot(1:10,sin(1:10));
get(lh,'XData')
get(lh,'XData')'
clear lh ans
lh=plot(1:10,sin(1:10));
get(lh,'XData')
ans(:)
isnan([1,nan,3])
(isnan([1,nan,3]) & isnan([nan;nan;3]))
clear lh ans
newstn = qcstation(stn)
q
stn.air_t_2_degc
newstn.air_t_2_degc
newstn=[]; clear newstn
clear ans
x={'a','b'}
x{2:end}
x{3:end}
disp(['X is ',x])
disp(['X is ',x{:}])
disp(['X is ',char(x)])
char(x)
disp(x)
newstn = qcstation(stn,{'linked','speed_2_kt','dir_2_deg'});
q
s
clear ans x
close all
pack
strncmpi([],'s',1)
strncmpi('s',[],1)
clear ans
x{2} = nan
x{2}
x{[]}
numel(x{[]})
isempty(x{[]})
newstn = qcstation(stn,{'linked','speed_2_kt','dir_2_deg'});
x
x = {{'foo','bar'},{'bif','bam'}};
find(strcmpi('bif',x))
strcmpi('bif',x)
help strmatch
strmatch('bif',x)
dbquit
dbclear all
cellfun(@(x)(strmatch(x,'bif')),x)
strmatch('bif',x{1})
strmatch('bif',x{2})
cellfun(@(q)(strmatch(q,'bif')),x)
help cellfun
cellfun(@(q)(strmatch(q,'bif')),x,'UniformOutput',false)
x
cellfun(@(q)(strmatch(q,'bif')),x,'UniformOutput',false)
x
x{2}
help strmatch
strvcat(x)
strvcat(x{:})
strcat(x{:})
x = {{'foo','bar'},{'bif','bam'},{'boo'}};
strcat(x{:})
cellfun(@(q)(strmatch(q,'bif')),x,'UniformOutput',false)
strmatch('bif',x{2})
strmatch('bif',x{3})
cellfun(@(q)(strmatch(q,'bif'),),x,'UniformOutput',false)
cellfun(@(q)(strmatch(q,'bif')),x,'UniformOutput',false)
cellfun(@(q)(~isempty(strmatch(q,'bif'))),x,'UniformOutput',false)
x{:}
x{1}{1}
x(:){1}
x{:}(1)
char(x)
char(x{:})
clear ans x
newstn = qcstation(stn,{'linked','speed_2_kt','dir_2_deg'});
dbstop if error
newstn = qcstation(stn,{'linked','speed_2_kt','dir_2_deg'});
args{1}
timenow
linkedflds
linkedflds(end+1) = [pi];
linkedflds(end+1) = args{1}(2:end);
linkedflds(end), args{1}(2:end),
dbquit
whos
newstn = qcstation(stn,{'linked','speed_2_kt','dir_2_deg'});
linkedflds
linkedflds(end+1) = {'pi'}
linkedflds={};
linkedflds(end+1) = { args{1}(2:end) };
linkedflds{1}
linkedflds{1}{1}
dbcont
args{1}
args{1}(2:end)
dbquit
newstn = qcstation(stn,{'linked','speed_2_kt','dir_2_deg'});
q
s
q
b
h
b
h
r
q
fmg; plot_ts(newstn.i_depth_deep_m,newstn.salinity_deep_psu,newstn.seatemp_deep_degc);
stn = newstn; newstn=[]; clear newstn
stn = qcstation(stn,{'seatemp_deep_degc','seatemp_shallow_degc'});
q
fmg; plot_ts(stn.i_depth_deep_m,stn.salinity_deep_psu,stn.seatemp_deep_degc);
fmg; spt(3,1,1); plot_ts(stn.i_depth_deep_m); spt(3,1,2); plot_ts(stn.salinity_deep_psu); spt(3,1,3); plot_ts(stn.seatemp_deep_degc);
fmg; a(1)=spt(3,1,1); plot_ts(stn.i_depth_deep_m); a(2)=spt(3,1,2); plot_ts(stn.salinity_deep_psu); a(3)=spt(3,1,3); plot_ts(stn.seatemp_deep_degc);
linkaxes(a)
help linkaxes
linkaxes(a,'x')
ylim default
stn
dir *.mat
whos('LPPR1_portal_ALL.mat')
whos('-mat','LPPR1_portal_ALL.mat')
help whos
whos('-file','LPPR1_portal_ALL.mat')
save('LPPR1_portal_ALL_qc.mat','stn')
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
load('LPPR1_portal_ALL_qc.mat','stn')
fmg; boxplot_ts(stn.seatemp_deep_degc,@get_year,'indices',@get_jas);
fmg; boxplot_ts(stn.seatemp_deep_degc,@get_year,'indices',@ts_jas);
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
stn = load_portal_data('SRVI2','ALL');
stn = qcstation(stn,{'linked','speed_2_kt','dir_2_deg'});
q
s
q
s
b
q
save('SRVI2_portal_ALL_qc.mat','stn')
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
fmg; boxplot_ts(stn.seatemp_deep_degc,@get_year,'indices',@ts_jas);
load('LPPR1_portal_ALL_qc.mat','stn')
fmg; boxplot_ts(stn.seatemp_deep_degc,@get_year,'indices',@ts_jas);
fmg; boxplot_ts(stn.seatemp_deep_degc,@get_year,'indices',@ts_jas); titlename([upper(stn.station_name),' summer (JAS) deep sea temperatures'); print('-dpng',[lower(stn.station_name),'_seatemp_deep_years.png']);
fmg; boxplot_ts(stn.seatemp_deep_degc,@get_year,'indices',@ts_jas); titlename([upper(stn.station_name),' summer (JAS) deep sea temperatures']); print('-dpng',[lower(stn.station_name),'_seatemp_deep_years.png']);
fmg; boxplot_ts(stn.seatemp_deep_degc,@get_year,'indices',@ts_jas,'mean',true); titlename([upper(stn.station_name),' summer (JAS) deep sea temperatures']); print('-dpng',[lower(stn.station_name),'_seatemp_deep_years.png']);
fmg; boxplot_ts(stn.seatemp_shallow_degc,@get_year,'indices',@ts_jas,'mean',true); titlename([upper(stn.station_name),' summer (JAS) shallow sea temperatures']); print('-dpng',[lower(stn.station_name),'_seatemp_shallow_years.png']);
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
load('SRVI2_portal_ALL_qc.mat','stn')
fmg; boxplot_ts(stn.seatemp_deep_degc,@get_year,'indices',@ts_jas,'mean',true); titlename([upper(stn.station_name),' summer (JAS) deep sea temperatures']); print('-dpng',[lower(stn.station_name),'_seatemp_deep_years.png']);
fmg; boxplot_ts(stn.seatemp_shallow_degc,@get_year,'indices',@ts_jas,'mean',true); titlename([upper(stn.station_name),' summer (JAS) shallow sea temperatures']); print('-dpng',[lower(stn.station_name),'_seatemp_shallow_years.png']);
fmg; boxplot_ts(stn.seatemp_deep_degc,@get_year,'indices',@ts_jfm,'mean',true); titlename([upper(stn.station_name),' winter (JFM) deep sea temperatures']); print('-dpng',[lower(stn.station_name),'_winter_seatemp_deep_years.png']);
fmg; boxplot_ts(stn.seatemp_shallow_degc,@get_year,'indices',@ts_jfm,'mean',true); titlename([upper(stn.station_name),' winter (JFM) shallow sea temperatures']); print('-dpng',[lower(stn.station_name),'_winter_seatemp_shallow_years.png']);
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
load('LPPR1_portal_ALL_qc.mat','stn')
fmg; boxplot_ts(stn.seatemp_deep_degc,@get_year,'indices',@ts_jfm,'mean',true); titlename([upper(stn.station_name),' winter (JFM) deep sea temperatures']); print('-dpng',[lower(stn.station_name),'_winter_seatemp_deep_years.png']);
fmg; boxplot_ts(stn.seatemp_shallow_degc,@get_year,'indices',@ts_jfm,'mean',true); titlename([upper(stn.station_name),' winter (JFM) shallow sea temperatures']); print('-dpng',[lower(stn.station_name),'_winter_seatemp_shallow_years.png']);
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
quit
%-- 6/10/2016 9:17 AM --%
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
stn = load_portal_data('LPPR1','ALL');
fmg; plot_ts(stn.seatemp_shallow_degc,stn.seatemp_deep_degc); legend('Shallow','Deep');
ylabel('^oC'); datetick3('x',3);
ylabel('^oC'); datetick3('x',2);
stnm = stn.station_name;
if (doPrint); print('-dpng',[lower(stnm),'_seatemps_raw.png']); end;
doPrint = true;
if (doPrint); print('-dpng',[lower(stnm),'_seatemps_raw.png']); end;
fmg; plot_ts(stn.seatemp_shallow_degc,stn.air_t_2_degc); legend('Shallow','Air'); ylabel('^oC'); datetick3('x',2);
fmg; plot_ts(stn.seatemp_shallow_degc,stn.air_t_2_degc);
legend('Shallow','Air'); ylabel('^oC'); datetick3('x',2);
titlename('Raw Sea and Air Temperatures');
if (doPrint); print('-dpng',[lower(stnm),'_sea_air_temps_raw.png']); end;
fmg; plot_ts(stn.seatemp_shallow_degc,stn.seatemp_deep_degc);
legend('Shallow','Deep'); ylabel('^oC'); datetick3('x',2);
titlename('Raw Sea Temperature');
if (doPrint); print('-dpng',[lower(stnm),'_seatemps_raw.png']); end;
fmg;
boxplot_ts(qc.stn.air_t_2_degc,@get_year,'mean',true,'indices',@ts_jas,...
'title',[upper(stnm),' Summer (JAS) Air Temperature']);
if (doPrint); print('-dpng',[lower(stnm),'_summer_airtemp_years.png']); end;
qc = load([upper(stnm),'_portal_ALL_qc.mat']);
fmg;
boxplot_ts(qc.stn.air_t_2_degc,@get_year,'mean',true,'indices',@ts_jas,...
'title',[upper(stnm),' Summer (JAS) Air Temperature']);
if (doPrint); print('-dpng',[lower(stnm),'_summer_airtemp_years.png']); end;
qc.stn = qcstation(qc.stn,'air_temp_2_degc');
stn.air_t_2_degcfmg;
boxplot_ts(qc.stn.air_t_2_degc,@get_year,'mean',true,'indices',@ts_jas,...
'title',[upper(stnm),' Summer (JAS) Air Temperature']);
if (doPrint); print('-dpng',[lower(stnm),'_summer_airtemp_years.png']); end;
fmg;
boxplot_ts(qc.stn.air_t_2_degc,@get_year,'mean',true,'indices',@ts_jas,...
'title',[upper(stnm),' Summer (JAS) Air Temperature']);
if (doPrint); print('-dpng',[lower(stnm),'_summer_airtemp_years.png']); end;
qc.stn = qcstation(qc.stn,'air_t_2_degc');
h
q
save('LPPR1_portal_ALL_qc.mat','-struct','qc','-v7.3')
fmg;
boxplot_ts(qc.stn.air_t_2_degc,@get_year,'mean',true,'indices',@ts_jas,...
'title',[upper(stnm),' Summer (JAS) Air Temperature']);
if (doPrint); print('-dpng',[lower(stnm),'_summer_airtemp_years.png']); end;
min([nanmin(stn.air_t_2_degc.data),nanmin(stn.seatemp_shallow_degc.data),nanmin(stn.seatemp_deep_degc.data)]);
min([nanmin(stn.air_t_2_degc.data),nanmin(stn.seatemp_shallow_degc.data),nanmin(stn.seatemp_deep_degc.data)])
min([nanmin(qc.stn.air_t_2_degc.data),nanmin(qc.stn.seatemp_shallow_degc.data),nanmin(qc.stn.seatemp_deep_degc.data)])
min([prctile(qc.stn.air_t_2_degc.data,3),prctile(qc.stn.seatemp_shallow_degc.data,3),prctile(qc.stn.seatemp_deep_degc.data,3)])
min([prctile(qc.stn.air_t_2_degc.data,1),prctile(qc.stn.seatemp_shallow_degc.data,1),prctile(qc.stn.seatemp_deep_degc.data,1)])
floor(min([prctile(qc.stn.air_t_2_degc.data,1),prctile(qc.stn.seatemp_shallow_degc.data,1),prctile(qc.stn.seatemp_deep_degc.data,1)]))
boxmin = floor(min([prctile(qc.stn.air_t_2_degc.data,1),prctile(qc.stn.seatemp_shallow_degc.data,1),prctile(qc.stn.seatemp_deep_degc.data,1)]));
boxmax = ceil(max([prctile(qc.stn.air_t_2_degc.data,99),prctile(qc.stn.seatemp_shallow_degc.data,99),prctile(qc.stn.seatemp_deep_degc.data,99)]));
boxmax
for_jim
close all
help good_grp_ts
type good_grp_ts.m
help find_pct_good_dates
type find_pct_good_dates
stn
fmg; plot_ts(stn.i_depth_deep_m,stn.i_depth_shallow_m);
fmg; plot_ts(qc.stn.i_depth_deep_m,qc.stn.i_depth_shallow_m);
fmg; plot_ts(stn.i_depth_deep_m,stn.i_depth_shallow_m);
qc.stn = qcstation(qc.stn,'i_depth_shallow_m');
q
save('LPPR1_portal_ALL_qc.mat','-struct','qc','-v7.3')
whos('-file','LPPR1_portal_ALL_qc.mat')
qc.stn = qcstation(qc.stn,'i_depth_deep_m');
s
q
s
fmg; plot_ts(stn.i_depth_deep_m,stn.i_depth_shallow_m); legend('Deep','Shallow');
titlename([upper(stn.station_name),' Raw Instrument Depths']);
if (doPrint); print('-dpng',[lower(stnm),'_depths_raw.png']); end;
fmg; plot_ts(qc.stn.i_depth_deep_m,qc.stn.i_depth_shallow_m); legend('Deep','Shallow');
titlename([upper(stn.station_name),' PARTIAL QC - Instrument Depths']);
if (doPrint); print('-dpng',[lower(stnm),'_depths_vs_qc.png']); end;
scatter_fit_ts(qc.stn.seatemp_shallow_degc,qc.stn.seatemp_deep_degc)
scatter_fit_ts_seasons(qc.stn.seatemp_shallow_degc,qc.stn.seatemp_deep_degc)
legend('Location
legend('Location','SouthEast')
scatter_fit_ts_seasons(qc.stn.seatemp_shallow_degc,qc.stn.seatemp_deep_degc,[],[],'Shallow Seatemp','Deep Seatemp',[],[],true)
scatter_fit_ts_seasons(qc.stn.seatemp_shallow_degc,qc.stn.seatemp_deep_degc,[],[],'Shallow Seatemp','Deep Seatemp')
subplots_set('xlim',[25,31.5],'ylim',[25,31.5]);
legend('Location','SouthEast')
if (doPrint); print('-dpng',[lower(stnm),'_seatemp_shallow_scatter_deep.png']); end;
legend('Location','SouthEast')
if (doPrint); print('-dpng',[lower(stnm),'_seatemp_shallow_scatter_deep.png']); end;
legend('Location','NorthEast')
if (doPrint); print('-dpng',[lower(stnm),'_seatemp_shallow_scatter_deep.png']); end;
qc.stn
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
for_jim
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
for_jim
l = load_portal_data('LPPR1','ALL');
help boxplot_tses
fmg; boxplot_tses({qc.stn.seatemp_deep_degc,l.seatemp_deep_degc},@get_year,'mean',true,'indices',@ts_jfm);
fmg; boxplot_tses({qc.stn.seatemp_deep_degc,l.seatemp_deep_degc},@get_year,'indices',@ts_jfm);
ylim([24,31])
y
y = load_portal_data('LPPR1','ALL');
y = load_portal_data('LCIY2','ALL');
fmg; boxplot_tses({qc.stn.seatemp_deep_degc,y.seatemp_deep_degc},@get_year,'indices',@ts_jfm);
ylim([24,31])
fmg; boxplot_tses({qc.stn.seatemp_deep_degc,l.seatemp_deep_degc,y.seatemp_deep_degc},@get_year,'indices',@ts_jfm); ylim([24,31]);
xlim
ylim([24,31])
stnm
if (doPrint); print('-dpng',[lower(stnm),'_winter_seatemp_deep_years.png']); end;
fmg; boxplot_tses({qc.stn.seatemp_deep_degc,l.seatemp_deep_degc,y.seatemp_deep_degc},@get_year,'indices',@ts_jfm);
fmg; boxplot_tses({qc.stn.seatemp_deep_degc,y.seatemp_deep_degc},@get_year,'indices',@ts_jfm);
fmg; boxplot_tses({qc.stn.seatemp_deep_degc,l.seatemp_deep_degc},@get_year,'indices',@ts_jfm);
ylim([24,31])
titlename('SRVI2 (black) vs. LPPR1 (blue) - Deep Seatemp');
doPrint = true;
if (doPrint); print('-dpng',[lower(stnm),'_vs_',lower(l.station_name),'_winter_seatemp_deep_years.png']); end;
if (doPrint); print('-dpng',[lower(stnm),'_winter_seatemp_deep_years.png']); end;
fmg;
boxplot_ts(qc.stn.air_t_2_degc,@get_year,'mean',true,'indices',@ts_jfm,...
'title',[upper(stnm),' Winter (JAS) Air Temperature']);
ylim([boxmin,boxmax]);
if (doPrint); print('-dpng',[lower(stnm),'_winter_airtemp_years.png']); end;
fmg;
boxplot_ts(qc.stn.speed_2_kt,@get_year,'mean',true,'indices',@ts_jfm,...
'title',[upper(stnm),' Winter (JAS) Wind Speeds']);
qc.stn.speed_2_kt
min(diff(qc.stn.speed_2_kt.date))
dbstop if error
fmg;
boxplot_ts(qc.stn.speed_2_kt,@get_year,'mean',true,'indices',@ts_jfm,...
'title',[upper(stnm),' Winter (JAS) Wind Speeds']);
ts
is_valid_ts(ts)
is_ts(ts)
type is_valid_ts.m
size(ts.date),size(ts.data)
nansummary(ts.date)
all(diff(ts.date)>0)
all(diff(ts.date)>eps)
min(diff(ts.date))
numel(find(isnan(ts.date)))
dbquit
badix = find(isnan(qc.stn.speed_2_kt.date));
qc.stn.speed_2_kt.date(badix) = [];
qc.stn.speed_2_kt.data(badix) = [];
is_valid_ts(qc.stn.speed_2_kt)
fmg;
boxplot_ts(qc.stn.speed_2_kt,@get_year,'mean',true,'indices',@ts_jfm,...
'title',[upper(stnm),' Winter (JAS) Wind Speeds']);
if (doPrint); print('-dpng',[lower(stnm),'_winter_speed_years.png']); end;
help multcompare
pd
dir *mult*
pd
help station_anova_multcompare
boxmin,boxmax
help ts_isfinite
is_valid_ts([])
help finite_ts
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
stnm = 'LPPR1';
for_jim
stn=[]; qc=[]; clear stn qc ans
clear all
pack
help good_grp_ts
help finite_ts
dir *portal*
stn = load_portal_data('LPPR1','ALL');
flds=fieldnames(stn); for fldix=1:numel(flds); fld=flds{fldix}; if ( is_ts(stn.(fld)) && ~is_valid_ts(stn.(fld)) ); disp fld; end; end;
flds=fieldnames(stn); for fldix=1:numel(flds); fld=flds{fldix}; disp fld; if ( is_ts(stn.(fld)) && ~is_valid_ts(stn.(fld)) ); disp BAD; end; end;
flds=fieldnames(stn); for fldix=1:numel(flds); fld=flds{fldix}; disp(fld); if ( is_ts(stn.(fld)) && ~is_valid_ts(stn.(fld)) ); disp BAD; end; end;
flds=fieldnames(stn); for fldix=1:numel(flds); fld=flds{fldix}; disp(fld); if ( ~is_valid_ts(stn.(fld)) ); disp BAD; end; end;
all(diff(stn.air_t_2_degc.date)>0)
all(diff(stn.speed_2_kt.date)>0)
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
stn = load_portal_data('LCIY2','ALL');
flds=fieldnames(stn); for fldix=1:numel(flds); fld=flds{fldix}; disp(fld); if ( ~is_valid_ts(stn.(fld)) ); disp BAD; end; end;
flds=fieldnames(stn); for fldix=1:numel(flds); fld=flds{fldix}; if ( ~is_valid_ts(stn.(fld)) ); disp(fld); end; end;
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
load('LPPR1_portal_ALL_qc.mat','stn')
flds=fieldnames(stn); for fldix=1:numel(flds); fld=flds{fldix}; if ( ~is_valid_ts(stn.(fld)) ); disp(fld); end; end;
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
load('LCIY2_portal_ALL_qc.mat','stn')
load('LPPR1_portal_ALL_qc_SAVE.mat','stn')
flds=fieldnames(stn); for fldix=1:numel(flds); fld=flds{fldix}; if ( ~is_valid_ts(stn.(fld)) ); disp(fld); end; end;
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
load('SRVI2_portal_ALL_qc_SAVE.mat','stn')
flds=fieldnames(stn); for fldix=1:numel(flds); fld=flds{fldix}; if ( ~is_valid_ts(stn.(fld)) ); disp(fld); end; end;
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
load('SRVI2_portal_ALL_qc.mat','stn')
flds=fieldnames(stn); for fldix=1:numel(flds); fld=flds{fldix}; if ( ~is_valid_ts(stn.(fld)) ); disp(fld); end; end;
x = finite_ts(stn.dir_2_deg);
x
stn.dir_2_deg
is_valid_ts(x)
is_valid_ts(stn.dir_2_deg)
numel(find(isnan(x.date)))
numel(find(isnan(stn.dir_2_deg.date)))
stn.dir_2_deg = finite_ts(stn.dir_2_deg);
stn.speed_2_kt = finite_ts(stn.speed_2_kt);
numel(find(isnan(stn.speed_2_kt.date)))
is_valid_ts(stn.speed_2_kt)
save('SRVI2_portal_ALL_qc.mat','stn','-v7.3')
flds=fieldnames(stn); for fldix=1:numel(flds); fld=flds{fldix}; if ( ~is_valid_ts(stn.(fld)) ); disp(fld); end; end;
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
stn = load_portal_data('LCIY2','ALL');
stn = qcstation(stn,{'linked','speed_2_kt','dir_2_deg'});
q
get(gca,'XData')
get(gco,'XData')
get(gca,'xlim')
clear ans
help brush
type set_pcolor_cursor.m
help brush
help brushing
help select2d
help brushing>select2d
help brushpoints
doc brush
xl = xlim
xlim(datenum(2012,[1,3],1))
xlim(datenum(2012,[1,7],1))
xlim(datenum(2012,[1,7],1)); datetick3;
xlim(datenum(2014,[1,7],1)); datetick3;
xlim(datenum(2013,[1,7],1)); datetick3;
xlim(datenum(2012,[1,7],1)); datetick3;
jumpax
q
get(gca,'child')
help ord
get(gca)
get(gca,'xaxis')
get(gca)
ax=gca
ploth = get(ax,'Children');
ploth
get(ploth,'foo')
help get
plothatt=get(ploth)
plothatt
isnumeric(plothatt.XData)
ismatrix(plothatt.XData)
isvector(plothatt.XData)
xlim(datenum(2012,[1,7],1)); datetick3;
xlim(datenum(2011,[1,7],1)); datetick3;
xlim(datenum(2010,[1,7],1)); datetick3;
fmg; plot_ts(stn.air_t_2_degc);
xlim(datenum(2010,[1,7],1)); datetick3;
xlim(datenum(2012,[1,7],1)); datetick3;
xlim(datenum(2013,[1,7],1)); datetick3;
xlim(datenum(2012,[1,7],1)); datetick3;
ax=gca
ploth = get(ax,'Children');
for hix = 1:numel(ploth)
h = ploth(hix);
att = get(h);
% Stop as soon as we find some data
if ( isfield(att,ordattnm) )
orddat = att.(ordattnm);
ordmin = nanmin(orddat(:));
ordmax = nanmax(orddat(:));
end; %if ( isfield(att,ordattnm) )
end; %for hix = 1:numel(ploth)
ordattnm='XData';
ploth = get(ax,'Children');
for hix = 1:numel(ploth)
h = ploth(hix);
att = get(h);
% Stop as soon as we find some data
if ( isfield(att,ordattnm) )
orddat = att.(ordattnm);
ordmin = nanmin(orddat(:));
ordmax = nanmax(orddat(:));
end; %if ( isfield(att,ordattnm) )
end; %for hix = 1:numel(ploth)
datestr([ordmin,ordmax])
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
for_jim
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
doPrint = true;
for_jim
qc.stn
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
doPrint = true;
for_jim
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
datestr(datenum(2012,1,0)+303)
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
68e3*365
68e3*365/1e9
68e6*365/1e9
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
datestr(datenum(2014,1,360))
datestr(datenum(2014,1,365))
mod(7,3)
mod(8,3)
mod(3,8)
mod(3,7)
mod(100,9)
mod(9,100)
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
hycpath = get_thesis_path('../data/hycom/GOM2');
yr=2014; jd=1;
hycpath = get_thesis_path('../data/hycom/GOM2');
ncname = sprintf('archv.%04d_%03d_12.nc',yr,jd);
ncpath = fullfile(hycpath,ncname);
nc = mDataset(ncpath);
nj_info(nc),
lon = cast(nc{'Longitude'}(:),'double');
lat = cast(nc{'Latitude'}(:),'double');
t = cast(nc{'temperature'}(:,:,:,:),'double');
19e6*365
z = cast(nc{'Depth'}(:),'double');
z
mld = cast(nc{'mld'}(:,:,:,:),'double');
mlp = cast(nc{'mlp'}(:,:,:,:),'double');
fmg; contourf(lon,lat,mld); colorbar;
fmg; contourf(lon,lat,mlp); colorbar;
set(gca,'clim',[0,180])
nansummary(mld)
nansummary(mlp)
nansummary(mlp-mld)
nansummary(mlp./mld)
close(nc); clear nc
min(diff(unique(lon))
min(diff(unique(lon)))
medain(diff(unique(lon)))
median(diff(unique(lon)))
median(diff(unique(lat)))
min(diff(unique(lat)))
max(diff(unique(lat)))
timenow
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
pd
angom2
mdl
who
whos
5847880*365
mdl
size(mdl.t)
mdl.z
datestr(mdl.d)
timenow
close(nc); clear nc
help ftp
f = ftp('ftp.aoml.noaa.gov','anonymous','lew.gramer@noaa.gov');
f
help cd
help ftp
doc ftp
close(f);
help cd
help ftp>cd
help mget
help catchwarn
365/100
ans*6
close(f); clear f
close(nc); clear nc
mdl
whos
6e6*365
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
help try
try, foobie; catch, catchwarn; try, disp('bar'); catch, end; end;
try, foobie; catch, catchwarn; try, bifbar; catch, end; end;
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
nc = mDataset(ncpath);
ncpath = fullfile(hycpath,ncname);
hycpath = get_thesis_path('../data/hycom/GOM2');
ncname = sprintf('archv.%04d_%03d_12.nc',yr,jd);
yr=2014; jd=1;
yr=2014; jd=99;
ncname = sprintf('archv.%04d_%03d_12.nc',yr,jd);
ncpath = fullfile(hycpath,ncname);
nc = mDataset(ncpath);
nj_info(nc),
close(nc); clear nc
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
d
dirs
pd +4
dir ../data/*heat*
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
load('mlrf1-heat_budget-erai-avhrr_weekly-ndbc-tpxo_tide-erai.mat');
load('../data/mlrf1-heat_budget-erai-avhrr_weekly-ndbc-tpxo_tide-erai.mat');
stn
plot_daily_clim_heat_budget(stn)
stn.optim
stn.optim.climsq
stn.optim.climsq.data
stn.optim.climq
nansummary(stn.optim.climq.data)
nansummary(stn.optim.climqt.data)
nansummary(stn.optim.climbqt)
nansummary(stn.optim.climbqt.data)
stn.optim.climsq
stn.optim.climsq(1)
stn.optim.climsq(1,1,1,1)
stn.optim.climsq(:,1,1,1)
stn.optim
nansummary(stn.optim.q.data)
nansummary(stn.optim.sq.data)
stn.optim.dayq
stn.optim.daysq
nansummary(stn.optim.daysq.data)
nansummary(stn.optim.dayq.data)
nansummary(stn.optim.dayq)
nansummary(stn.optim.dayq.data)
help genvarname
stn.optim
stn
grepstruct(stn,'_dly')
nansummary(stn.ndbc_erai_erai_30a_avhrr_hc_dTdt_flux_dly.data)
10*1000
-19/3/1.03e3/sw_cp(35,30,0)
help sw_cp
-19/3/1.03e3/sw_cp(35,25,0)
19*365*-19/3/1.03e3/sw_cp(35,25,0)
19*365*(-19/3/1.03e3/sw_cp(35,25,0))
19*365*(-19/3/sw_dens0(35,25)/sw_cp(35,25,0))
clear ans
nansummary(stn.ndbc_erai_erai_30a_avhrr_hc_dTdt_flux_dly.data)
nansummary(stn.ndbc_erai_erai_30a_avhrr_hc_dTdt_dly.data)
stn.ndbc_erai_erai_30a_avhrr_hc_dTdt_dly
fmg; boxplot_ts(stn.ndbc_erai_erai_30a_avhrr_hc_dTdt_dly);
find_date_ranges(stn.ndbc_erai_erai_30a_avhrr_hc_dTdt_dly.date)
grepstruct(stn,'_dly')
find_date_ranges(stn.ndbc_erai_erai_30a_avhrr_hc_dTdt_dly.date)
find_date_ranges(stn.ndbc_sea_t_dly.date)
find_date_ranges(stn.ndbc_erai_erai_30a_avhrr_hc_dTdt_err.date)
find_date_ranges(stn.ndbc_erai_erai_30a_avhrr_hc_dTdt_dly.date)
265/40
ans*3e9
find_date_ranges(stn.ndbc_erai_erai_30a_avhrr_hc_dTdthc_err.date)
find_date_ranges(stn.ndbc_erai_erai_30a_avhrr_hc_dTdt_err_1_d_sum.date)
reviewanim([],0,0,0)
print('-dpng','../../../../coral/ICRS/mlrf1_heat_budget.png')
clear ans
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
dirs
pd +5
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
quit
%-- 6/17/2016 9:44 PM --%
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
quit
%-- 6/18/2016 8:37 AM --%
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
plot_dt_vs_beta_with_hc
set(gca,'xscale','linear')
xl=xlim
xlim([0,0.301])
xlim([0,0.0301])
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
plot_dt_vs_beta_with_hc
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
plot_dt_vs_beta_with_hc
help num2str
num2str(pi,2)
num2str(pi/10,2)
num2str(pi/100,2)
num2str(pi*100,2)
num2str(pi*10,2)
num2str(35.1,2)
num2str(39.1,2)
num2str(35.9,2)
num2str(-1.26,2)
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
doPrint = true;
plot_dt_vs_beta_with_hc
exp(-1.3)
1e-1.3
10^(-1.3)
10^(-1.1)
10^(-1.5)
10^(-1.4)
10^(-1.45)
10^(-1.35)
10^(-1.33)
10^(-1.36)
log(0.4)
log(0.04)
log10(0.04)
log10(0.0401)
log10(0.0405)
log10(0.0010)
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
doPrint = true;
plot_dt_vs_beta_with_hc
text(min(bets),0,'Flat bottom');
get(gco)
text(bets(floor(end/4)),0,'Flat bottom','Rotation',90);
text(bets(floor(end/4)),0,'Flat bottom','Rotation',90,'HorizontalAlign','center');
text(bets(floor(end/4)),0,'Flat bottom','Rotation',90,'HorizontalAlign','center','FontSize',12);
text(bets(floor(end/4)),0,'Reef flats','Rotation',90,'HorizontalAlign','center','FontSize',12);
text(bets(ceil(3*end/4)),0,'Platform slope','Rotation',90,'HorizontalAlign','center','FontSize',12);
text(bets(round(end/2)),maxdT*0.75,'Diurnal warming','Rotation',0,'HorizontalAlign','center','FontSize',12);
text(bets(round(end/2)),mindT*0.75,'Cold-front passage','Rotation',0,'HorizontalAlign','center','FontSize',12);
mindT
maxdT
text(bets(floor(end/4)),0,'Reef flats','Rotation',90,'HorizontalAlign','center','FontSize',12);
text(bets(ceil(3*end/4)),0,'Platform slope','Rotation',90,'HorizontalAlign','center','FontSize',12);
text(bets(round(end/2)),maxdT*0.75,'Diurnal warming','Rotation',0,'HorizontalAlign','center','FontSize',12);
text(bets(round(end/2)),mindT*0.75,'Cold-front passage','Rotation',0,'HorizontalAlign','center','FontSize',12);
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
load(get_thesis_path('../data/mlrf1-heat_budget-erai-avhrr_weekly-ndbc-tpxo_tide-erai.mat'));
nansummary(stn.ndbc_erai_erai_30a_avhrr_hc_dTdt_dly.data)
nansummary(stn.ndbc_erai_erai_30a_avhrr_hc_dTdt_flux_dly.data)
fmg; hist(stn.ndbc_erai_erai_30a_avhrr_hc_dTdt_flux_dly.data,100)
fmg; hist(stn.ndbc_erai_erai_30a_avhrr_hc_dTdt_flux_dly.data,100); xlim([-5000,5000]);
fmg; hist(stn.ndbc_erai_erai_30a_avhrr_hc_dTdt_flux_1_d_avg.data,100); xlim([-5000,5000]);
fmg; hist(stn.ndbc_erai_erai_30a_avhrr_hc_dTdt_flux_1_d_avg.data,100); xlim([-500,500]);
fmg; hist(stn.ndbc_erai_erai_30a_avhrr_hc_dTdt_flux_dly.data,100); xlim([-5000,5000]);
fmg; hist(stn.ndbc_erai_erai_30a_avhrr_hc_dTdt_flux_dly.data/24,100); xlim([-500,500]);
fmg; hist(stn.ndbc_erai_erai_30a_avhrr_hc_dTdt_flux_dly.data*3/24,100); xlim([-500,500]);
sumix = find(get_season(stn.ndbc_erai_erai_30a_avhrr_hc_dTdt_flux_dly.date)==3);
winix = find(get_season(stn.ndbc_erai_erai_30a_avhrr_hc_dTdt_flux_dly.date)==1);
fmg; hist(stn.ndbc_erai_erai_30a_avhrr_hc_dTdt_flux_dly.data(sumix)/24,100); xlim([-500,500]); titlename('Summer');
fmg; hist(stn.ndbc_erai_erai_30a_avhrr_hc_dTdt_flux_dly.data(winix)/24,100); xlim([-500,500]); titlename('Winter');
fmg; hist(stn.ndbc_erai_erai_30a_avhrr_hc_dTdt_flux_dly.data(sumix)/24,100); xlim([-500,500]); titlename('Summer');
fmg; hist(stn.ndbc_erai_erai_30a_avhrr_hc_dTdt_flux_dly.data(sumix)*3/24,100); xlim([-500,500]); titlename('Summer');
sumix = find(get_season(stn.ndbc_erai_erai_30a_avhrr_hc_dTdt_flux.date)==3);
winix = find(get_season(stn.ndbc_erai_erai_30a_avhrr_hc_dTdt_flux.date)==1);
sumix = find(get_season(stn.ndbc_erai_erai_30a_avhrr_qtadv_flux.date)==3);
winix = find(get_season(stn.ndbc_erai_erai_30a_avhrr_qtadv_flux.date)==1);
fmg; hist(stn.ndbc_erai_erai_30a_avhrr_qtadv_flux.data(sumix),100); xlim([-500,500]); titlename('Summer');
fmg; hist(stn.ndbc_erai_erai_30a_avhrr_qtadv_flux.data(winix),100); xlim([-500,500]); titlename('Winter');
fmg; hist(stn.ndbc_erai_erai_30a_avhrr_qtadv_flux.data(sumix),100); titlename('Summer');
fmg; hist(stn.ndbc_erai_erai_30a_avhrr_qtadv_flux.data(winix),100); titlename('Winter');
stn = verify_variable(stn,'ndbc_erai_erai_30a_avhrr_qtadv_flux_3_h_sum');
sumix = find(get_season(stn.ndbc_erai_erai_30a_avhrr_qtadv_flux_3_h_sum.date)==3);
winix = find(get_season(stn.ndbc_erai_erai_30a_avhrr_qtadv_flux_3_h_sum.date)==1);
fmg; hist(stn.ndbc_erai_erai_30a_avhrr_qtadv_flux_3_h_sum.data(sumix),100); titlename('Summer');
fmg; hist(stn.ndbc_erai_erai_30a_avhrr_qtadv_flux_3_h_sum.data(winix),100); titlename('Winter');
close all
stn = verify_variable(stn,'ndbc_erai_erai_30a_avhrr_qtadv_flux_3_h_avg');
sumix = find(get_season(stn.ndbc_erai_erai_30a_avhrr_qtadv_flux_3_h_avg.date)==3);
winix = find(get_season(stn.ndbc_erai_erai_30a_avhrr_qtadv_flux_3_h_avg.date)==1);
fmg; hist(stn.ndbc_erai_erai_30a_avhrr_qtadv_flux_3_h_avg.data(sumix),100); titlename('Summer');
fmg; hist(stn.ndbc_erai_erai_30a_avhrr_qtadv_flux_3_h_avg.data(winix),100); titlename('Winter');
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
doPrint = true;
plot_dt_vs_beta_with_hc
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
doPrint = true;
plot_dt_vs_beta_with_hc
set(gca,'xscale','log')
prctile(bets,7)
std(bets,7)
std(bets)
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
nbets = 40;
bets = logspace(log10(0.0010),log10(0.0405),nbets);
prctile(bets,25)
prctile(bets,75)
prctile(bets,93)
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
load(get_thesis_path('../data/mlrf1-heat_budget-erai-avhrr_weekly-ndbc-tpxo_tide-erai.mat'));
fmg; boxplot_ts(stn.ndbc_erai_erai_30a_avhrr_qtadv_flux);
stn = verify_variable(stn,'ndbc_erai_erai_30a_avhrr_qtadv_flux_3_h_avg');
fmg; boxplot_ts(stn.ndbc_erai_erai_30a_avhrr_qtadv_flux_3_h_avg);
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
load(get_thesis_path('../data/lonf1-heat_budget-erai-avhrr_weekly-ndbc-tpxo_tide-erai.mat'));
stn = verify_variable(stn,'ndbc_erai_erai_30a_avhrr_qtadv_flux_3_h_avg');
fmg; boxplot_ts(stn.ndbc_erai_erai_30a_avhrr_qtadv_flux_3_h_avg);
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
load(get_thesis_path('../data/mlrf1-heat_budget-erai-avhrr_weekly-ndbc-tpxo_tide-erai.mat'));
stn = verify_variable(stn,'ndbc_erai_erai_30a_avhrr_qtadv_flux_24_h_min');
stn = verify_variable(stn,'ndbc_erai_erai_30a_avhrr_qtadv_flux_24_h_max');
fmg; plot_ts(stn.ndbc_erai_erai_30a_avhrr_qtadv_flux_24_h_min,stn.ndbc_erai_erai_30a_avhrr_qtadv_flux_24_h_max);
stn = verify_variable(stn,'ndbc_erai_erai_30a_avhrr_qtadv_flux_24_h_iqr');
fmg; plot_ts(stn.ndbc_erai_erai_30a_avhrr_qtadv_flux_24_h_iqr);
fmg; boxplot_ts(stn.ndbc_erai_erai_30a_avhrr_qtadv_flux_24_h_iqr);
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
doPrint = true;
plot_dt_vs_beta_with_hc
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
doPrint = true;
plot_dt_vs_beta_with_hc
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
angom2
doGOM2 = true;
doPrint = true;
cmpefkeys
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
quit
%-- 6/18/2016 2:26 PM --%
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
angom2
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
angom2
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
angom2
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
doGOM2 = true;
cmpefkeys
reviewanim([],0,0,0)
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
angom2
mdl
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
anfkeys
find_date_ranges(mdl.d)
datestr(mdl.d)
dts = datenum(2008,1,1,6,0,0):(6/24):datenum(2008,1,7,18,0,0);
datestr(dts)
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
anefkeys
datestr(mdl.d)
mdl.mlrf1
mdl
dts = datenum(2012,1,1,6,0,0):(6/24):datenum(2012,1,7,18,0,0);
dix = find(ismember(mdl.d,dts));
cmpefkeys
reviewanim([],0,0,0)
mdl.mlrf1
mdl.mlrf1.t
stnm
stnm = 'mlrf1';
min_t = nanmin([min_t,nanmin(mdl.(stnm).t.prof(:))]);
max_t = nanmax([max_t,nanmax(mdl.(stnm).t.prof(:))]);
min_t,max_t
dix = find(ismember(mdl.d,dts));
min_t = nanmin([min_t,nanmin(mdl.(stnm).t.prof(dix,:))]);
max_t = nanmax([max_t,nanmax(mdl.(stnm).t.prof(dix,:))]);
max_t = -inf;
max_t = nanmax([max_t,nanmax(mdl.(stnm).t.prof(dix,:))]);
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
cmpefkeys
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
doPrint = true;
cmpefkeys
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
doPrint = true;
doFKEYS = true;
cmpefkeys
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
doPrint = true;
doGOM2 = true;
cmpefkeys
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
quit
%-- 6/19/2016 8:58 PM --%
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
which gradientm
type gradientm
stn = plot_hires_bathymetry(stn,-[0:5:50,100:100:800],[dx*1.05,dy*1.10],true,[],[],[],false);
stn = get_station_from_station_name(stnm); stn = load_all_ndbc_data(stn);
stnm = 'mlrf1';
stn = get_station_from_station_name(stnm); stn = load_all_ndbc_data(stn);
stn = plot_hires_bathymetry(stn,-[0:5:50,100:100:800],[dx*1.05,dy*1.10],true,[],[],[],false);
dx=1e3; dy=1e3;
stn = plot_hires_bathymetry(stn,-[0:5:50,100:100:800],[dx*1.05,dy*1.10],true,[],[],[],false);
stn.ngdc_hires_bathy
stn.ngdc_hires_bathy.field(floor(end/2),floor(end/2))
stn.ngdc_hires_bathy.field(floor(end/2)+1,floor(end/2)+1)
stn.ngdc_hires_bathy.field(floor(end/2)+[0,1,0,1],floor(end/2)+[0,0,1,1])
stn=[]; qc=[]; clear stn qc ans
stnm = 'sanf1';
stn = get_station_from_station_name(stnm); stn = load_all_ndbc_data(stn);
stn = plot_hires_bathymetry(stn,-[0:5:50,100:100:800],[dx*1.05,dy*1.10],true,[],[],[],false);
stn.ngdc_hires_bathy.field(floor(end/2)+[0,1,0,1],floor(end/2)+[0,0,1,1])
stn=[]; qc=[]; clear stn qc ans
stn = plot_hires_bathymetry('sanf1',-[0:5:50,100:100:800],[dx*1.05,dy*1.10],true,[],[],[],true);
stn.ngdc_hires_bathy.field(floor(end/2)+[0,1,0,1],floor(end/2)+[0,0,1,1])
stn=[]; qc=[]; clear stn qc ans
stn = plot_hires_bathymetry('sanf1',-[0:5:50,100:100:800],[dx*1.05,dy*1.10],true,[],[],[],false);
stn.ngdc_hires_bathy.field(floor(end/2)+[0,1,0,1],floor(end/2)+[0,0,1,1])
stn=[]; qc=[]; clear stn qc ans
stn = plot_hires_bathymetry('smkf1',-[0:5:50,100:100:800],[dx*1.05,dy*1.10],true,[],[],[],false);
stn.ngdc_hires_bathy.field(floor(end/2)+[0,1,0,1],floor(end/2)+[0,0,1,1])
stn=[]; qc=[]; clear stn qc ans
stn = plot_hires_bathymetry('smkf1',-[0:5:50,100:100:800],[dx*1.05,dy*1.10],true,[],[],[],true);
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
dir *xaym*
dir *xay*
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
pd
pd +@
pd +2
dir *xay*
extract_xaymara_sites
stn = ans;
ans=[]; clear ans
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
stn = extract_xaymara_sites
stn.lons
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
stn = extract_xaymara_sites
stn.ngdc_hires_bathy
clear ans
[loner,lonix] = min(abs(stn.lons(1)-stn.ngdc_hires_bathy.lon));
[later,latix] = min(abs(stn.lats(1)-stn.ngdc_hires_bathy.lat));
stn.ngdc_hires_bathy.field(latix+[0,1,0,1],lonix+[0,0,1,1])
x.lon=stn.lons(1); x.lat=stn.lats(1);
x = plot_hires_bathymetry(x,-[0:5:50,100:100:800],[dx*1.05,dy*1.10],true,[],[],[],true);
dx=1e3; dy=1e3;
x = plot_hires_bathymetry(x,-[0:5:50,100:100:800],[dx*1.05,dy*1.10],true,[],[],[],true);
[ig,ix]=min(stn.lons)
x.lon=stn.lons(ix); x.lat=stn.lats(ix);
x = plot_hires_bathymetry(x,-[0:5:50,100:100:800],[dx*1.05,dy*1.10],true,[],[],[],true);
x = plot_hires_bathymetry(x,-[0:5:50,100:100:800],[dx*1.05,dy*1.10],true,[],[],[],false);
x=[]; clear x
x.lon=stn.lons(ix); x.lat=stn.lats(ix);
x = plot_hires_bathymetry(x,-[0:5:50,100:100:800],[dx*1.05,dy*1.10],true,[],[],[],true);
x.ngdc_hires_bathy.field(floor(end/2)+[0,1,0,1],floor(end/2)+[0,0,1,1])
x=[]; x.lon=stn.lons(ix); x.lat=stn.lats(ix);
x = plot_hires_bathymetry(x,-[0:5:50,100:100:800],[dx*1.05,dy*1.10],true,[],[],[],false);
help license
