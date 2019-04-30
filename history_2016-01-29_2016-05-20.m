%-- 1/29/2016 12:59 AM --%
dir
stn = get_station_from_station_name('mlrf1')
pwd
cd ecoforecasts\
stn = get_station_from_station_name('mlrf1')
edit startup.m
which load_all_ndbc_data.m
which get_fkeys_hycom
dirs
help fileparts
ver
help ver
help verLessThan
verLessThan('matlab','R2008a')
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
verLessThan('matlab','7.5')
version
help version
version -release
stn = get_station_from_station_name('mlrf1')
stn = load_all_ndbc_data(stn);
help fileparts
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
stn = get_station_from_station_name('mlrf1')
stn = load_all_ndbc_data(stn);
[pathroot, ig, ig] = fileparts(mfilename('fullpath'));
clear pathroot ig
fmg; plot_ts(stn.ndbc_sea_t);
scatter_fit_ts(stn.ndbc_wind1_speed,stn.ndbc_sea_t)
scatter_fit_ts_seasons(stn.ndbc_wind1_speed,stn.ndbc_sea_t)
subplots_set('ylim',[16,32])
subplots_set('ylim',[17,32])
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
path
which startup
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
quit
%-- 1/29/2016 4:57 PM --%
help iscolumn
which iscolumn
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
quit
%-- 1/31/2016 1:23 PM --%
dbstop if error
matlabrc
dbstop startup
matlabrc
dbquit
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
edit my_startup
matlabrc
dbstop initdirs
matlabrc
dbquit
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
matlabrc
[m,d]=lastwarn
help warning
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
quit
%-- 2/1/2016 1:34 PM --%
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
dirs
pd +5
type extract_sfomc_etc.m
extract_sfomc_etc
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
extract_sfomc_etc
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
extract_sfomc_etc
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
dbstop if error
extract_sfomc_etc
errmsg
dbup
size(res.c_buoy.(fld).date(xlix)),size(res.c_buoy.adcp_heights)
size(res.c_buoy.(fld).date(xlix)),size(res.c_buoy.adcp_heights),size(res.c_buoy.(fld).prof(xlix,:)')
xlix
xl
datestr(xl)
res.c_buoy
res.c_buoy.(fld)
datestr(xl)
find_date_ranges(res.c_buoy.(fld).date,1)
dbquit
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
extract_sfomc_etc
q
res
res.nw_w_btm
edit export_to_coris.m
export_to_coris
dir CoRIS
close all
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
quit
%-- 2/3/2016 2:46 PM --%
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
extract_sfomc_etc
res
res.nw_w_btm = verify_variable(res.nw_w_btm,'adcp_seatemp_diff_24_h_sum');
res.nw_w_btm = verify_variable(res.nw_w_btm,'adcp_x_btm_24_h_sum');
res.nw_w_btm = verify_variable(res.nw_w_btm,'adcp_baroclinic_x_btm_24_h_sum');
res.nw_w_btm = verify_variable(res.nw_w_btm,'adcp_seatemp_diff_48_h_sum');
res.nw_w_btm = verify_variable(res.nw_w_btm,'adcp_x_btm_48_h_sum');
res.nw_w_btm = verify_variable(res.nw_w_btm,'adcp_baroclinic_x_btm_48_h_sum');
res.nw_w_btm = verify_variable(res.nw_w_btm,'adcp_seatemp_diff_72_h_sum');
res.nw_w_btm = verify_variable(res.nw_w_btm,'adcp_x_btm_72_h_sum');
res.nw_w_btm = verify_variable(res.nw_w_btm,'adcp_baroclinic_x_btm_72_h_sum');
res.pier_cc = verify_variable(res.pier_cc,'wind_stress_48_h_sum');
res.pier_cc = verify_variable(res.pier_cc,'air_t_diff_48_h_sum');
res
help fieldnames
res.nw_w_btm
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
pd
pd +8
pd +2
quit
%-- 2/4/2016 7:48 PM --%
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
quit
%-- 2/5/2016 12:01 PM --%
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
help fopen
export_to_coris
res.(stnm)
res.nw_w_btm
res.(stnm).adcp_depth
nansummary(res.(stnm).adcp_depth.data)
numel(find(res.(stnm).adcp_depth.data<10.9))
numel(find(res.(stnm).adcp_depth.data<10))
numel(find(res.(stnm).adcp_depth.data<9))
numel(find(res.(stnm).adcp_depth.data<4))
numel(find(res.(stnm).adcp_depth.data<2))
res.nw_w_btm.adcp_heights
res.(stnm).adcp_depths = res.(stnm).depth res.(stnm).adcp_heights;
res.(stnm).adcp_depths = res.(stnm).depth - res.(stnm).adcp_heights;
nansummary(res.(stnm).adcp_depths)
res.(stnm).adcp_depths
res
res.c_buoy
res.c_buoy.adcp_heights
res.
res
res.sw_buoy.adcp_heights
res.se_buoy.adcp_heights
res.se_buoy
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
export_to_coris
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
export_to_coris
% Calculate multi-bin current averages to match temperature profiles
for cstnm = fieldnames(res)';
stnm = cstnm{:};
if ( isfield(res.(stnm),'seatemp') )
deps = res.(stnm).seatemp.depths;
res.(stnm).xshore.date = res.(stnm).adcp_x.date;
res.(stnm).lshore.date = res.(stnm).adcp_l.date;
for dix = 1:numel(deps)
[derr,hix] = min(abs(deps(dix) - res.(stnm).adcp_depths));
if ( derr > 5 )
error('No ADCP bins for %s depth %g',stnm,deps(dix));
end;
hix = hix-1:hix+1;
res.(stnm).xshore.prof(:,dix) = nanmean(res.(stnm).adcp_x.prof(:,hix),1);
res.(stnm).lshore.prof(:,dix) = nanmean(res.(stnm).adcp_l.prof(:,hix),1);
end;
res.(stnm).xshore.data = nanmean(res.(stnm).xshore.prof,2);
res.(stnm).lshore.data = nanmean(res.(stnm).lshore.prof,2);
res.(stnm).xshore.depths = deps;
res.(stnm).lshore.depths = deps;
end;
end;
hix
res.(stnm).adcp_x
% Calculate multi-bin current averages to match temperature profiles
for cstnm = fieldnames(res)';
stnm = cstnm{:};
if ( isfield(res.(stnm),'seatemp') )
deps = res.(stnm).seatemp.depths;
res.(stnm).xshore.date = res.(stnm).adcp_x.date;
res.(stnm).lshore.date = res.(stnm).adcp_l.date;
for dix = 1:numel(deps)
[derr,hix] = min(abs(deps(dix) - res.(stnm).adcp_depths));
if ( derr > 5 )
error('No ADCP bins for %s depth %g',stnm,deps(dix));
end;
hix = hix-1:hix+1;
res.(stnm).xshore.prof(:,dix) = nanmean(res.(stnm).adcp_x.prof(:,hix),1);
res.(stnm).lshore.prof(:,dix) = nanmean(res.(stnm).adcp_l.prof(:,hix),1);
end;
res.(stnm).xshore.data = nanmean(res.(stnm).xshore.prof,2);
res.(stnm).lshore.data = nanmean(res.(stnm).lshore.prof,2);
res.(stnm).xshore.depths = deps;
res.(stnm).lshore.depths = deps;
end;
end;
% Calculate multi-bin current averages to match temperature profiles
for cstnm = fieldnames(res)';
stnm = cstnm{:};
if ( isfield(res.(stnm),'seatemp') )
deps = res.(stnm).seatemp.depths;
res.(stnm).xshore.date = res.(stnm).adcp_x.date;
res.(stnm).lshore.date = res.(stnm).adcp_l.date;
for dix = 1:numel(deps)
[derr,center_ix(dix)] = min(abs(deps(dix) - res.(stnm).adcp_depths));
if ( derr > 5 )
error('No ADCP bins for %s depth %g',stnm,deps(dix));
end;
end;
for dix = 1:numel(deps)
if ( dix == 1 )
loix = 1;
hiix = mean([center_ix(dix),center_ix(dix+1)]);
elseif ( dix == numel(deps) )
loix = mean([center_ix(dix-1),center_ix(dix)]);
hiix = length(res.(stnm).adcp_depths);
else
loix = mean([center_ix(dix-1),center_ix(dix)]);
hiix = mean([center_ix(dix),center_ix(dix+1)]);
end;
disp({stnm,deps(dix),loix,hiix,});
res.(stnm).xshore.prof(:,dix) = nanmean(res.(stnm).adcp_x.prof(:,loix:hiix),1);
res.(stnm).lshore.prof(:,dix) = nanmean(res.(stnm).adcp_l.prof(:,loix:hiix),1);
end;
res.(stnm).xshore.data = nanmean(res.(stnm).xshore.prof,2);
res.(stnm).lshore.data = nanmean(res.(stnm).lshore.prof,2);
res.(stnm).xshore.depths = deps;
res.(stnm).lshore.depths = deps;
end;
end;
derr
res.(stnm).adcp_depths
% For more readable Hoevmullers
res.nw_w_btm.adcp_depths = -(res.nw_w_btm.depth - res.nw_w_btm.adcp_heights);
res.c_buoy.adcp_depths = -(res.c_buoy.depth - res.c_buoy.adcp_heights);
res.e_buoy.adcp_depths = -(res.e_buoy.depth - res.e_buoy.adcp_heights);
res.ne_buoy.adcp_depths = -(res.ne_buoy.depth - res.ne_buoy.adcp_heights);
res.sw_buoy.adcp_depths = -(res.sw_buoy.depth - res.sw_buoy.adcp_heights);
% Calculate multi-bin current averages to match temperature profiles
for cstnm = fieldnames(res)';
stnm = cstnm{:};
if ( isfield(res.(stnm),'seatemp') )
deps = res.(stnm).seatemp.depths;
res.(stnm).xshore.date = res.(stnm).adcp_x.date;
res.(stnm).lshore.date = res.(stnm).adcp_l.date;
for dix = 1:numel(deps)
[derr,center_ix(dix)] = min(abs(deps(dix) - res.(stnm).adcp_depths));
if ( derr > 5 )
error('No ADCP bins for %s depth %g',stnm,deps(dix));
end;
end;
for dix = 1:numel(deps)
if ( dix == 1 )
loix = 1;
hiix = mean([center_ix(dix),center_ix(dix+1)]);
elseif ( dix == numel(deps) )
loix = mean([center_ix(dix-1),center_ix(dix)]);
hiix = length(res.(stnm).adcp_depths);
else
loix = mean([center_ix(dix-1),center_ix(dix)]);
hiix = mean([center_ix(dix),center_ix(dix+1)]);
end;
disp({stnm,deps(dix),loix,hiix,});
res.(stnm).xshore.prof(:,dix) = nanmean(res.(stnm).adcp_x.prof(:,loix:hiix),1);
res.(stnm).lshore.prof(:,dix) = nanmean(res.(stnm).adcp_l.prof(:,loix:hiix),1);
end;
res.(stnm).xshore.data = nanmean(res.(stnm).xshore.prof,2);
res.(stnm).lshore.data = nanmean(res.(stnm).lshore.prof,2);
res.(stnm).xshore.depths = deps;
res.(stnm).lshore.depths = deps;
end;
end;
round(0.5)
% Calculate multi-bin current averages to match temperature profiles
for cstnm = fieldnames(res)';
stnm = cstnm{:};
if ( isfield(res.(stnm),'seatemp') )
deps = res.(stnm).seatemp.depths;
res.(stnm).xshore.date = res.(stnm).adcp_x.date;
res.(stnm).lshore.date = res.(stnm).adcp_l.date;
for dix = 1:numel(deps)
[derr,center_ix(dix)] = min(abs(deps(dix) - res.(stnm).adcp_depths));
if ( derr > 5 )
error('No ADCP bins for %s depth %g',stnm,deps(dix));
end;
end;
for dix = 1:numel(deps)
if ( dix == 1 )
loix = 1;
hiix = floor(mean([center_ix(dix),center_ix(dix+1)]));
elseif ( dix == numel(deps) )
loix = ceil(mean([center_ix(dix-1),center_ix(dix)]));
hiix = length(res.(stnm).adcp_depths);
else
loix = ceil(mean([center_ix(dix-1),center_ix(dix)]));
hiix = floor(mean([center_ix(dix),center_ix(dix+1)]));
end;
disp({stnm,deps(dix),loix,hiix,});
res.(stnm).xshore.prof(:,dix) = nanmean(res.(stnm).adcp_x.prof(:,loix:hiix),1);
res.(stnm).lshore.prof(:,dix) = nanmean(res.(stnm).adcp_l.prof(:,loix:hiix),1);
end;
res.(stnm).xshore.data = nanmean(res.(stnm).xshore.prof,2);
res.(stnm).lshore.data = nanmean(res.(stnm).lshore.prof,2);
res.(stnm).xshore.depths = deps;
res.(stnm).lshore.depths = deps;
end;
end;
center_ix
deps
res.nw_w_btm.adcp_heights
res.nw_w_btm.adcp_depths
% Calculate multi-bin current averages to match temperature profiles
for cstnm = fieldnames(res)';
stnm = cstnm{:};
if ( isfield(res.(stnm),'seatemp') )
deps = res.(stnm).seatemp.depths;
res.(stnm).xshore.date = res.(stnm).adcp_x.date;
res.(stnm).lshore.date = res.(stnm).adcp_l.date;
for dix = 1:numel(deps)
[derr,center_ix(dix)] = min(abs(deps(dix) - res.(stnm).adcp_depths));
if ( derr > 5 )
error('No ADCP bins for %s depth %g',stnm,deps(dix));
end;
end;
for dix = 1:numel(deps)
if ( dix == 1 )
loix = ceil(mean([center_ix(dix-1),center_ix(dix)]));
hiix = length(res.(stnm).adcp_depths);
elseif ( dix == numel(deps) )
loix = 1;
hiix = floor(mean([center_ix(dix-1),center_ix(dix)]));
else
loix = ceil(mean([center_ix(dix),center_ix(dix+1)]));
hiix = floor(mean([center_ix(dix-1),center_ix(dix)]));
end;
disp({stnm,deps(dix),loix,hiix,});
res.(stnm).xshore.prof(:,dix) = nanmean(res.(stnm).adcp_x.prof(:,loix:hiix),1);
res.(stnm).lshore.prof(:,dix) = nanmean(res.(stnm).adcp_l.prof(:,loix:hiix),1);
end;
res.(stnm).xshore.data = nanmean(res.(stnm).xshore.prof,2);
res.(stnm).lshore.data = nanmean(res.(stnm).lshore.prof,2);
res.(stnm).xshore.depths = deps;
res.(stnm).lshore.depths = deps;
end;
end;
center_ix
loix
hiix
deps
% Calculate multi-bin current averages to match temperature profiles
for cstnm = fieldnames(res)';
stnm = cstnm{:};
if ( isfield(res.(stnm),'seatemp') )
deps = res.(stnm).seatemp.depths;
res.(stnm).xshore.date = res.(stnm).adcp_x.date;
res.(stnm).lshore.date = res.(stnm).adcp_l.date;
for dix = 1:numel(deps)
[derr,center_ix(dix)] = min(abs(deps(dix) - res.(stnm).adcp_depths));
if ( derr > 5 )
error('No ADCP bins for %s depth %g',stnm,deps(dix));
end;
end;
for dix = 1:numel(deps)
if ( dix == 1 )
loix = 1;
hiix = floor(mean([center_ix(dix-1),center_ix(dix)]));
elseif ( dix == numel(deps) )
loix = ceil(mean([center_ix(dix-1),center_ix(dix)]));
hiix = length(res.(stnm).adcp_depths);
else
loix = ceil(mean([center_ix(dix),center_ix(dix+1)]));
hiix = floor(mean([center_ix(dix-1),center_ix(dix)]));
end;
disp({stnm,deps(dix),loix,hiix,});
res.(stnm).xshore.prof(:,dix) = nanmean(res.(stnm).adcp_x.prof(:,loix:hiix),1);
res.(stnm).lshore.prof(:,dix) = nanmean(res.(stnm).adcp_l.prof(:,loix:hiix),1);
end;
res.(stnm).xshore.data = nanmean(res.(stnm).xshore.prof,2);
res.(stnm).lshore.data = nanmean(res.(stnm).lshore.prof,2);
res.(stnm).xshore.depths = deps;
res.(stnm).lshore.depths = deps;
end;
end;
loix
hiix
res.(stnm).adcp_x
dix
% Calculate multi-bin current averages to match temperature profiles
for cstnm = fieldnames(res)';
stnm = cstnm{:};
if ( isfield(res.(stnm),'seatemp') )
deps = res.(stnm).seatemp.depths;
res.(stnm).xshore.date = res.(stnm).adcp_x.date;
res.(stnm).lshore.date = res.(stnm).adcp_l.date;
for dix = 1:numel(deps)
[derr,center_ix(numel(deps)-dix+1)] = min(abs(deps(dix) - res.(stnm).adcp_depths));
if ( derr > 5 )
error('No ADCP bins for %s depth %g',stnm,deps(dix));
end;
end;
for dix = 1:numel(deps)
if ( dix == 1 )
loix = 1;
hiix = floor(mean([center_ix(dix-1),center_ix(dix)]));
elseif ( dix == numel(deps) )
loix = ceil(mean([center_ix(dix-1),center_ix(dix)]));
hiix = length(res.(stnm).adcp_depths);
else
loix = ceil(mean([center_ix(dix),center_ix(dix+1)]));
hiix = floor(mean([center_ix(dix-1),center_ix(dix)]));
end;
disp({stnm,deps(dix),loix,hiix,});
res.(stnm).xshore.prof(:,dix) = nanmean(res.(stnm).adcp_x.prof(:,loix:hiix),1);
res.(stnm).lshore.prof(:,dix) = nanmean(res.(stnm).adcp_l.prof(:,loix:hiix),1);
end;
res.(stnm).xshore.data = nanmean(res.(stnm).xshore.prof,2);
res.(stnm).lshore.data = nanmean(res.(stnm).lshore.prof,2);
res.(stnm).xshore.depths = deps;
res.(stnm).lshore.depths = deps;
end;
end;
dix
loix
hiix
center_ix
% Calculate multi-bin current averages to match temperature profiles
for cstnm = fieldnames(res)';
stnm = cstnm{:};
if ( isfield(res.(stnm),'seatemp') )
deps = res.(stnm).seatemp.depths;
res.(stnm).xshore.date = res.(stnm).adcp_x.date;
res.(stnm).lshore.date = res.(stnm).adcp_l.date;
for dix = 1:numel(deps)
[derr,center_ix(numel(deps)-dix+1)] = min(abs(deps(dix) - res.(stnm).adcp_depths));
if ( derr > 5 )
error('No ADCP bins for %s depth %g',stnm,deps(dix));
end;
end;
for dix = 1:numel(deps)
if ( dix == 1 )
loix = 1;
hiix = floor(mean([center_ix(dix),center_ix(dix+1)]));
elseif ( dix == numel(deps) )
loix = ceil(mean([center_ix(dix-1),center_ix(dix)]));
hiix = length(res.(stnm).adcp_depths);
else
loix = ceil(mean([center_ix(dix-1),center_ix(dix)]));
hiix = floor(mean([center_ix(dix),center_ix(dix+1)]));
end;
disp({stnm,deps(dix),loix,hiix,});
res.(stnm).xshore.prof(:,dix) = nanmean(res.(stnm).adcp_x.prof(:,loix:hiix),1);
res.(stnm).lshore.prof(:,dix) = nanmean(res.(stnm).adcp_l.prof(:,loix:hiix),1);
end;
res.(stnm).xshore.data = nanmean(res.(stnm).xshore.prof,2);
res.(stnm).lshore.data = nanmean(res.(stnm).lshore.prof,2);
res.(stnm).xshore.depths = deps;
res.(stnm).lshore.depths = deps;
end;
end;
dix
loix
hiix
center_ix
res.(stnm).adcp_x
res.(stnm)
res.(stnm).xshore
res.(stnm).adcp_x
size(nanmean(res.(stnm).adcp_x.prof(:,loix:hiix),1))
size(nanmean(res.(stnm).adcp_x.prof(:,loix:hiix),2))
% Calculate multi-bin current averages to match temperature profiles
for cstnm = fieldnames(res)';
stnm = cstnm{:};
if ( isfield(res.(stnm),'seatemp') )
deps = res.(stnm).seatemp.depths;
res.(stnm).xshore.date = res.(stnm).adcp_x.date;
res.(stnm).lshore.date = res.(stnm).adcp_l.date;
for dix = 1:numel(deps)
[derr,center_ix(numel(deps)-dix+1)] = min(abs(deps(dix) - res.(stnm).adcp_depths));
if ( derr > 5 )
error('No ADCP bins for %s depth %g',stnm,deps(dix));
end;
end;
for dix = 1:numel(deps)
if ( dix == 1 )
loix = 1;
hiix = floor(mean([center_ix(dix),center_ix(dix+1)]));
elseif ( dix == numel(deps) )
loix = ceil(mean([center_ix(dix-1),center_ix(dix)]));
hiix = length(res.(stnm).adcp_depths);
else
loix = ceil(mean([center_ix(dix-1),center_ix(dix)]));
hiix = floor(mean([center_ix(dix),center_ix(dix+1)]));
end;
disp({stnm,deps(dix),loix,hiix,});
res.(stnm).xshore.prof(:,dix) = nanmean(res.(stnm).adcp_x.prof(:,loix:hiix),2);
res.(stnm).lshore.prof(:,dix) = nanmean(res.(stnm).adcp_l.prof(:,loix:hiix),2);
end;
res.(stnm).xshore.data = nanmean(res.(stnm).xshore.prof,2);
res.(stnm).lshore.data = nanmean(res.(stnm).lshore.prof,2);
res.(stnm).xshore.depths = deps;
res.(stnm).lshore.depths = deps;
end;
end;
dis
dix
loix,hiix
res.(stnm).adcp_x
res.(stnm).xshore
% Calculate multi-bin current averages to match temperature profiles
for cstnm = fieldnames(res)';
stnm = cstnm{:};
if ( isfield(res.(stnm),'seatemp') )
deps = res.(stnm).seatemp.depths;
res.(stnm).xshore.date = res.(stnm).adcp_x.date;
res.(stnm).lshore.date = res.(stnm).adcp_l.date;
res.(stnm).xshore.prof = repmat(nan,[numel(res.(stnm).adcp_x.date),numel(deps)]);
res.(stnm).lshore.prof = repmat(nan,[numel(res.(stnm).adcp_l.date),numel(deps)]);
for dix = 1:numel(deps)
[derr,center_ix(numel(deps)-dix+1)] = min(abs(deps(dix) - res.(stnm).adcp_depths));
if ( derr > 5 )
error('No ADCP bins for %s depth %g',stnm,deps(dix));
end;
end;
for dix = 1:numel(deps)
if ( dix == 1 )
loix = 1;
hiix = floor(mean([center_ix(dix),center_ix(dix+1)]));
elseif ( dix == numel(deps) )
loix = ceil(mean([center_ix(dix-1),center_ix(dix)]));
hiix = length(res.(stnm).adcp_depths);
else
loix = ceil(mean([center_ix(dix-1),center_ix(dix)]));
hiix = floor(mean([center_ix(dix),center_ix(dix+1)]));
end;
disp({stnm,deps(dix),loix,hiix,});
res.(stnm).xshore.prof(:,dix) = nanmean(res.(stnm).adcp_x.prof(:,loix:hiix),2);
res.(stnm).lshore.prof(:,dix) = nanmean(res.(stnm).adcp_l.prof(:,loix:hiix),2);
end;
res.(stnm).xshore.data = nanmean(res.(stnm).xshore.prof,2);
res.(stnm).lshore.data = nanmean(res.(stnm).lshore.prof,2);
res.(stnm).xshore.depths = deps;
res.(stnm).lshore.depths = deps;
end;
end;
% Calculate multi-bin current averages to match temperature profiles
for cstnm = fieldnames(res)';
stnm = cstnm{:};
if ( isfield(res.(stnm),'seatemp') )
deps = res.(stnm).seatemp.depths;
res.(stnm).xshore.date = res.(stnm).adcp_x.date;
res.(stnm).lshore.date = res.(stnm).adcp_l.date;
res.(stnm).xshore.prof = repmat(nan,[numel(res.(stnm).adcp_x.date),numel(deps)]);
res.(stnm).lshore.prof = repmat(nan,[numel(res.(stnm).adcp_l.date),numel(deps)]);
for dix = 1:numel(deps)
[derr,center_ix(numel(deps)-dix+1)] = min(abs(deps(dix) - res.(stnm).adcp_depths));
if ( derr > 5 )
error('No ADCP bins for %s depth %g',stnm,deps(dix));
end;
end;
for dix = 1:numel(deps)
if ( dix == 1 )
loix = 1;
hiix = floor(mean([center_ix(dix),center_ix(dix+1)]));
elseif ( dix == numel(deps) )
%loix = ceil(mean([center_ix(dix-1),center_ix(dix)]));
loix = hiix+1;
hiix = length(res.(stnm).adcp_depths);
else
%loix = ceil(mean([center_ix(dix-1),center_ix(dix)]));
loix = hiix+1;
hiix = floor(mean([center_ix(dix),center_ix(dix+1)]));
end;
disp({stnm,deps(dix),loix,hiix,});
res.(stnm).xshore.prof(:,dix) = nanmean(res.(stnm).adcp_x.prof(:,loix:hiix),2);
res.(stnm).lshore.prof(:,dix) = nanmean(res.(stnm).adcp_l.prof(:,loix:hiix),2);
end;
res.(stnm).xshore.data = nanmean(res.(stnm).xshore.prof,2);
res.(stnm).lshore.data = nanmean(res.(stnm).lshore.prof,2);
res.(stnm).xshore.depths = deps;
res.(stnm).lshore.depths = deps;
end;
end;
res.nw_w_btm
res.pier_cc
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
export_to_coris
stnm
res.(stnm)
[se26,se35,se98] = ...
intersect_tses([],res.se_buoy.at_26m.sbe_seatemp,res.se_buoy.at_35m.sbe_seatemp,res.se_buoy.at_98m.sbe_seatemp);
res.se_buoy.seatemp.date = se26.date;
res.se_buoy.seatemp.prof = [se26.data,se35.data,se98.data];
res.se_buoy.seatemp.data = nanmean(res.se_buoy.seatemp.prof,2);
res.se_buoy.seatemp.depths = -[26,35,98];
fclose(fid);
fname = fullfile(corispath,'CRCP_789_Product_1435_metadata_NSUOC_sites.csv');
fid = fopen(fname,'w+');
if fid < 3; error('Unable to open %s!',fname); end;
fprintf(fid,'METADATA: NOAA CRCP Project #789: Upwelling off southeast Florida\n');
fprintf(fid,'Project POC: Dr. Lew Gramer, lew.gramer@noaa.gov, 305-562-1735\n');
fprintf(fid,'NOAA Fed POC: Dr. Jim Hendee, jim.hendee@noaa.gov, 305-361-4396\n');
fprintf(fid,'Station,Depth (m),Lon,Lat,Instrument depth (m),Sea temperature start date,End date,Ocean currents start date,End date,Other measurements\n');
cstnms = fieldnames(res);
for cix = 1:numel(cstnms)
stnm = cstnms{cix};
if ( ~isfield(res.(stnm),'seatemp') )
deps = [res.(stnm).depth];
else
deps = res.(stnm).seatemp.depths;
end;
for dix = 1:numel(deps)
dep = deps(dix);
tbeg='-';
tend='-';
cbeg='-';
cend='-';
otherflds = '-';
switch (stnm),
case 'pier_cc',
otherflds = 'Meteorology';
case 'nw_w_btm',
tbeg = datestr(res.(stnm).sbe_seatemp.date(1));
tend = datestr(res.(stnm).sbe_seatemp.date(end));
cbeg = datestr(res.(stnm).adcp_x.date(1));
cend = datestr(res.(stnm).adcp_x.date(end));
otherwise,
tbeg = datestr(res.(stnm).seatemp.date(1));
tend = datestr(res.(stnm).seatemp.date(end));
cbeg = datestr(res.(stnm).xshore.date(1));
cend = datestr(res.(stnm).xshore.date(end));
end;
fprintf(fid,'%s,%g,%g,%g,%g,%s,%s,%s,%s,%s\n',upper(stnm),res.(stnm).depth,...
res.(stnm).lon,res.(stnm).lat,dep,tbeg,tend,cbeg,cend,otherflds);
end;
end;
fclose(fid);
% Calculate multi-bin current averages to match temperature profiles
for cstnm = fieldnames(res)';
stnm = cstnm{:};
if ( isfield(res.(stnm),'seatemp') )
deps = res.(stnm).seatemp.depths;
res.(stnm).xshore.date = res.(stnm).adcp_x.date;
res.(stnm).lshore.date = res.(stnm).adcp_l.date;
res.(stnm).xshore.prof = repmat(nan,[numel(res.(stnm).adcp_x.date),numel(deps)]);
res.(stnm).lshore.prof = repmat(nan,[numel(res.(stnm).adcp_l.date),numel(deps)]);
for dix = 1:numel(deps)
[derr,center_ix(numel(deps)-dix+1)] = min(abs(deps(dix) - res.(stnm).adcp_depths));
if ( derr > 5 )
error('No ADCP bins for %s depth %g',stnm,deps(dix));
end;
end;
for dix = 1:numel(deps)
if ( dix == 1 )
loix = 1;
hiix = floor(mean([center_ix(dix),center_ix(dix+1)]));
elseif ( dix == numel(deps) )
%loix = ceil(mean([center_ix(dix-1),center_ix(dix)]));
loix = hiix+1;
hiix = length(res.(stnm).adcp_depths);
else
%loix = ceil(mean([center_ix(dix-1),center_ix(dix)]));
loix = hiix+1;
hiix = floor(mean([center_ix(dix),center_ix(dix+1)]));
end;
disp({stnm,deps(dix),loix,hiix,});
res.(stnm).xshore.prof(:,dix) = nanmean(res.(stnm).adcp_x.prof(:,loix:hiix),2);
res.(stnm).lshore.prof(:,dix) = nanmean(res.(stnm).adcp_l.prof(:,loix:hiix),2);
end;
res.(stnm).xshore.data = nanmean(res.(stnm).xshore.prof,2);
res.(stnm).lshore.data = nanmean(res.(stnm).lshore.prof,2);
res.(stnm).xshore.depths = deps;
res.(stnm).lshore.depths = deps;
end;
end;
stnm
% Calculate multi-bin current averages to match temperature profiles
for cstnm = fieldnames(res)';
stnm = cstnm{:};
if ( isfield(res.(stnm),'seatemp') && isfield(res.(stnm),'adcp_x') )
deps = res.(stnm).seatemp.depths;
res.(stnm).xshore.date = res.(stnm).adcp_x.date;
res.(stnm).lshore.date = res.(stnm).adcp_l.date;
res.(stnm).xshore.prof = repmat(nan,[numel(res.(stnm).adcp_x.date),numel(deps)]);
res.(stnm).lshore.prof = repmat(nan,[numel(res.(stnm).adcp_l.date),numel(deps)]);
for dix = 1:numel(deps)
[derr,center_ix(numel(deps)-dix+1)] = min(abs(deps(dix) - res.(stnm).adcp_depths));
if ( derr > 5 )
error('No ADCP bins for %s depth %g',stnm,deps(dix));
end;
end;
for dix = 1:numel(deps)
if ( dix == 1 )
loix = 1;
hiix = floor(mean([center_ix(dix),center_ix(dix+1)]));
elseif ( dix == numel(deps) )
%loix = ceil(mean([center_ix(dix-1),center_ix(dix)]));
loix = hiix+1;
hiix = length(res.(stnm).adcp_depths);
else
%loix = ceil(mean([center_ix(dix-1),center_ix(dix)]));
loix = hiix+1;
hiix = floor(mean([center_ix(dix),center_ix(dix+1)]));
end;
disp({stnm,deps(dix),loix,hiix,});
res.(stnm).xshore.prof(:,dix) = nanmean(res.(stnm).adcp_x.prof(:,loix:hiix),2);
res.(stnm).lshore.prof(:,dix) = nanmean(res.(stnm).adcp_l.prof(:,loix:hiix),2);
end;
res.(stnm).xshore.data = nanmean(res.(stnm).xshore.prof,2);
res.(stnm).lshore.data = nanmean(res.(stnm).lshore.prof,2);
res.(stnm).xshore.depths = deps;
res.(stnm).lshore.depths = deps;
end;
end;
res.se_buoy
res.nw_w_btm
res.ne_buoy
find_date_ranges(res.ne_buoy.wind_speed.date)
res.se_buoy
res.ne_buoy
res.e_buoy
res.sw_buoy
fclose(fid);
fname = fullfile(corispath,'CRCP_789_Product_1435_metadata_NSUOC_sites.csv');
fid = fopen(fname,'w+');
if fid < 3; error('Unable to open %s!',fname); end;
fprintf(fid,'METADATA: NOAA CRCP Project #789: Upwelling off southeast Florida\n');
fprintf(fid,'Project POC: Dr. Lew Gramer, lew.gramer@noaa.gov, 305-562-1735\n');
fprintf(fid,'NOAA Fed POC: Dr. Jim Hendee, jim.hendee@noaa.gov, 305-361-4396\n');
fprintf(fid,'Station,Depth (m),Lon,Lat,Instrument depth (m),Sea temperature start date,End date,Ocean currents start date,End date,Other measurements\n');
cstnms = fieldnames(res);
for cix = 1:numel(cstnms)
stnm = cstnms{cix};
if ( ~isfield(res.(stnm),'seatemp') )
deps = [res.(stnm).depth];
else
deps = res.(stnm).seatemp.depths;
end;
for dix = 1:numel(deps)
dep = deps(dix);
tbeg='-';
tend='-';
cbeg='-';
cend='-';
otherflds = '-';
if ( strcmpi(stnm,'pier_cc') || strcmpi(stnm,'ne_buoy') )
otherflds = 'Meteorology';
end;
switch (stnm),
case 'nw_w_btm',
tbeg = datestr(res.(stnm).sbe_seatemp.date(1));
tend = datestr(res.(stnm).sbe_seatemp.date(end));
cbeg = datestr(res.(stnm).adcp_x.date(1));
cend = datestr(res.(stnm).adcp_x.date(end));
otherflds = 'Waves';
case 'se_buoy',
tbeg = datestr(res.(stnm).seatemp.date(1));
tend = datestr(res.(stnm).seatemp.date(end));
otherwise,
tbeg = datestr(res.(stnm).seatemp.date(1));
tend = datestr(res.(stnm).seatemp.date(end));
cbeg = datestr(res.(stnm).xshore.date(1));
cend = datestr(res.(stnm).xshore.date(end));
end;
fprintf(fid,'%s,%g,%g,%g,%g,%s,%s,%s,%s,%s\n',upper(stnm),res.(stnm).depth,...
res.(stnm).lon,res.(stnm).lat,dep,tbeg,tend,cbeg,cend,otherflds);
end;
end;
fclose(fid);
stnm
res.pier_cc
fclose(fid);
fname = fullfile(corispath,'CRCP_789_Product_1435_metadata_NSUOC_sites.csv');
fid = fopen(fname,'w+');
if fid < 3; error('Unable to open %s!',fname); end;
fprintf(fid,'METADATA: NOAA CRCP Project #789: Upwelling off southeast Florida\n');
fprintf(fid,'Project POC: Dr. Lew Gramer, lew.gramer@noaa.gov, 305-562-1735\n');
fprintf(fid,'NOAA Fed POC: Dr. Jim Hendee, jim.hendee@noaa.gov, 305-361-4396\n');
fprintf(fid,'Station,Depth (m),Lon,Lat,Instrument depth (m),Sea temperature start date,End date,Ocean currents start date,End date,Other measurements\n');
cstnms = fieldnames(res);
for cix = 1:numel(cstnms)
stnm = cstnms{cix};
if ( ~isfield(res.(stnm),'seatemp') )
deps = [res.(stnm).depth];
else
deps = res.(stnm).seatemp.depths;
end;
for dix = 1:numel(deps)
dep = deps(dix);
tbeg='-';
tend='-';
cbeg='-';
cend='-';
otherflds = '-';
if ( strcmpi(stnm,'pier_cc') || strcmpi(stnm,'ne_buoy') )
otherflds = 'Meteorology';
end;
switch (stnm),
case 'pier_cc',
case 'nw_w_btm',
tbeg = datestr(res.(stnm).sbe_seatemp.date(1));
tend = datestr(res.(stnm).sbe_seatemp.date(end));
cbeg = datestr(res.(stnm).adcp_x.date(1));
cend = datestr(res.(stnm).adcp_x.date(end));
otherflds = 'Waves';
case 'se_buoy',
tbeg = datestr(res.(stnm).seatemp.date(1));
tend = datestr(res.(stnm).seatemp.date(end));
otherwise,
tbeg = datestr(res.(stnm).seatemp.date(1));
tend = datestr(res.(stnm).seatemp.date(end));
cbeg = datestr(res.(stnm).xshore.date(1));
cend = datestr(res.(stnm).xshore.date(end));
end;
fprintf(fid,'%s,%g,%g,%g,%g,%s,%s,%s,%s,%s\n',upper(stnm),res.(stnm).depth,...
res.(stnm).lon,res.(stnm).lat,dep,tbeg,tend,cbeg,cend,otherflds);
end;
end;
fclose(fid);
export_to_coris
help intersect_ts
help intersect_ts_data
fprintf(2,'%g',nan)
fprintf(2,'%g\n',nan)
fprintf(2,'%s,\n',datestr(res.nw_w_btm.sbe_seatemp.date))
fprintf(2,'%s,\n',datestr(res.nw_w_btm.sbe_seatemp.date'))
datestr(res.nw_w_btm.sbe_seatemp.date')
datestr(res.nw_w_btm.sbe_seatemp.date)
strvcat(datestr(res.nw_w_btm.sbe_seatemp.date))
fprintf(2,'%s,\n',strvcat(datestr(res.nw_w_btm.sbe_seatemp.date)))
fprintf(2,'%s,\n',strcat(datestr(res.nw_w_btm.sbe_seatemp.date)))
x=strvcat(datestr(res.nw_w_btm.sbe_seatemp.date))
x
size(x)
fprintf(2,'%20c,\n',x)
q
fprintf(2,'%20c,\n',x')
q
info
doc
datestr(res.nw_w_btm.sbe_seatemp.date,30)
datestr(res.nw_w_btm.sbe_seatemp.date,31)
datestr(res.nw_w_btm.sbe_seatemp.date,32)
datestr(res.nw_w_btm.sbe_seatemp.date,30)
datestr(res.nw_w_btm.sbe_seatemp.date,31)
size(ans)
datestr(res.nw_w_btm.sbe_seatemp.date)
datestr(res.nw_w_btm.sbe_seatemp.date,31)
min(diff(res.nw_w_btm.sbe_seatemp.date))
max(diff(res.nw_w_btm.sbe_seatemp.date))
help interspace
help gap_fill
help gap_filter_ts
help gap_fill
help gap_expand
[newt,newd] = gap_expand(res.nw_w_btm.sbe_seatemp.date,res.nw_w_btm.sbe_seatemp.data);
size(newt)
res.nw_w_btm.sbe_seatemp
help interp1
interp1(dts,date
interp1(res.nw_w_btm.sbe_seatemp.date,res.nw_w_btm.sbe_seatemp.data,dts,'linear',nan)
help gap_expand
help filter_gaps
help gap_filter_ts
x.origt = res.nw_w_btm.sbe_seatemp;
x.origc = res.nw_w_btm.adcp_x;
begdt = min(x.origt.date(1),x.origc.date(1));
enddt = max(x.origt.date(end),x.origc.date(end));
deldt = min([min(diff(x.origt.date)),min(diff(x.origc.date))]);
x.t.date = begdt:deldt:enddt;
x.t.data = interp1(x.origt.date,x.origt.data,x.t.date,[],nan);
x = filter_gaps(x,'origt','t',(1/24),0,[],nan);
x.c.date = begdt:deldt:enddt;
x.c.data = interp1(x.origc.date,x.origc.data,x.c.date,[],nan);
x = filter_gaps(x,'origc','c',(1/24),0,[],nan);
x
x=[]; clear x
help filter_gaps
x=[]; clear x
x.origt = res.nw_w_btm.sbe_seatemp;
x.origc = res.nw_w_btm.adcp_x;
begdt = min(x.origt.date(1),x.origc.date(1));
enddt = max(x.origt.date(end),x.origc.date(end));
deldt = min([min(diff(x.origt.date)),min(diff(x.origc.date))]);
x.t.date = begdt:deldt:enddt;
x.t.data = interp1(x.origt.date,x.origt.data,x.t.date,[],nan);
x = filter_gaps(x,'origt','t',(1/24),0,[],nan);
x.c.date = begdt:deldt:enddt;
x.c.data = interp1(x.origc.date,x.origc.data,x.c.date,[],nan);
x = filter_gaps(x,'origc','c',(1/24),0,[],nan);
dbstop filter_aps 66
dbstop filter_gaps 66
x = filter_gaps(x,'origt','t',(1/24),0,[],nan);
dbquit
x = filter_gaps(x,'origt','t',(2/24),0,[],nan);
dbquit
x = filter_gaps(x,'origt','t',(1.5/24),0,[],nan);
dbquit
x=[]; clear x
x.origt = res.nw_w_btm.sbe_seatemp;
x.origc = res.nw_w_btm.adcp_x;
begdt = min(x.origt.date(1),x.origc.date(1));
enddt = max(x.origt.date(end),x.origc.date(end));
deldt = min([min(diff(x.origt.date)),min(diff(x.origc.date))]);
x.t.date = begdt:deldt:enddt;
x.t.data = interp1(x.origt.date,x.origt.data,x.t.date,[],nan);
x = filter_gaps(x,'origt','t',(1.5/24),0,[],nan);
x.c.date = begdt:deldt:enddt;
x.c.data = interp1(x.origc.date,x.origc.data,x.c.date,[],nan);
x = filter_gaps(x,'origc','c',(1.5/24),0,[],nan);
x
x.origt
x.t
delt
deldt
deldt*24
x=[]; clear x
x.origt = res.nw_w_btm.sbe_seatemp;
x.origc = res.nw_w_btm.adcp_x;
begdt = min(x.origt.date(1),x.origc.date(1));
enddt = max(x.origt.date(end),x.origc.date(end));
x.t.date = [begdt:(1/24):enddt]';
x.t.data = interp1(x.origt.date,x.origt.data,x.t.date,[],nan);
x = filter_gaps(x,'origt','t',(1.5/24),0,[],nan);
x.c.date = [begdt:(1/24):enddt]';
x.c.data = interp1(x.origc.date,x.origc.data,x.c.date,[],nan);
x = filter_gaps(x,'origc','c',(1.5/24),0,[],nan);
x.origt
x.t
min(diff(x.t.date))
find_date_ranges(res.nw_w_btm.sbe_seatemp.date,1)
fprintf(1,'%g',nan)
res.nw_w_btm
scatter_fit_ts(res.nw_w_btm.adcp_seatemp,res.nw_w_btm.sbe_seatemp)
fmg; plot_ts(ts_op(res.nw_w_btm.adcp_seatemp,res.nw_w_btm.sbe_seatemp,'-'));
find_date_ranges(res.nw_w_btm.sbe_seatemp.date,1)
find_date_ranges(res.e_buoy.at_5m.sbe_seatemp.date,1)
find_date_ranges(res.ne_buoy.at_5m.sbe_seatemp.date,1)
tic, x=[]; clear x
x.origt = res.nw_w_btm.adcp_seatemp;
x.origc = res.nw_w_btm.adcp_x;
begdt = min(x.origt.date(1),x.origc.date(1));
enddt = max(x.origt.date(end),x.origc.date(end));
x.t.date = [begdt:(1/24):enddt]';
x.t.data = interp1(x.origt.date,x.origt.data,x.t.date,[],nan);
x = filter_gaps(x,'origt','t',(1.5/24),0,[],nan);
x.c.date = [begdt:(1/24):enddt]';
x.c.data = interp1(x.origc.date,x.origc.data,x.c.date,[],nan);
x = filter_gaps(x,'origc','c',(1.5/24),0,[],nan);
toc,
x
x.t
x.c
numel
numel(find(isnan(x.c.data)))
numel(find(isnan(x.t.data)))
numel(find(isnan(x.t.data)&isnan(x.c.data)))
min(diff(x.t.date))
min(diff(x.c.date))
fprintf(1,'%d,%d,%d,%d,%d,%d,%g,%g,%g\n',datevec(x.c.date),x.t.data,x.l.date,x.c.data)
x=[]; clear x
x.origt = res.nw_w_btm.adcp_seatemp;
x.origc = res.nw_w_btm.adcp_x;
x.origl = res.nw_w_btm.adcp_l;
begdt = min([x.origt.date(1),  x.origc.date(1),  x.origl.date(1)]);
enddt = max([x.origt.date(end),x.origc.date(end),x.origl.date(end)]);
x.t.date = [begdt:(1/24):enddt]';
x.t.data = interp1(x.origt.date,x.origt.data,x.t.date,[],nan);
x = filter_gaps(x,'origt','t',(1.5/24),0,[],nan);
x.c.date = [begdt:(1/24):enddt]';
x.c.data = interp1(x.origc.date,x.origc.data,x.c.date,[],nan);
x = filter_gaps(x,'origc','c',(1.5/24),0,[],nan);
x.l.date = [begdt:(1/24):enddt]';
x.l.data = interp1(x.origl.date,x.origl.data,x.l.date,[],nan);
x = filter_gaps(x,'origl','l',(1.5/24),0,[],nan);
badix = find(isnan(x.t.data) & isnan(x.c.data) & isnan(x.l.data));
x.t.date(badix)=[]; x.t.data(badix)=[]; x.c.date(badix)=[]; x.c.data(badix)=[]; x.l.date(badix)=[]; x.l.data(badix)=[];
fprintf(1,'%d,%d,%d,%d,%d,%d,%g,%g,%g\n',datevec(x.c.date),x.t.data,x.l.date,x.c.data)
%-- 2/8/2016 12:46 PM --%
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
export_to_coris
fprintf(1,'%d,%d,%d,%d,%d,%d,-11,%g,%g,%g\n',datevec(x.c.date),x.t.data,x.c.data,x.l.date)
%-- 2/8/2016 1:27 PM --%
more on
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
extract_sfomc_etc
q
close all
[m,d]=lastwarn
x
x=[]; clear x
x.origt = res.nw_w_btm.adcp_seatemp;
x.origc = res.nw_w_btm.adcp_x;
x.origl = res.nw_w_btm.adcp_l;
begdt = min([x.origt.date(1),  x.origc.date(1),  x.origl.date(1)]);
enddt = max([x.origt.date(end),x.origc.date(end),x.origl.date(end)]);
x.t.date = [begdt:(1/24):enddt]';
x.t.data = interp1(x.origt.date,x.origt.data,x.t.date,[],nan);
x = filter_gaps(x,'origt','t',(1.5/24),0,[],nan);
x.c.date = [begdt:(1/24):enddt]';
x.c.data = interp1(x.origc.date,x.origc.data,x.c.date,[],nan);
x = filter_gaps(x,'origc','c',(1.5/24),0,[],nan);
x.l.date = [begdt:(1/24):enddt]';
x.l.data = interp1(x.origl.date,x.origl.data,x.l.date,[],nan);
x = filter_gaps(x,'origl','l',(1.5/24),0,[],nan);
badix = find(isnan(x.t.data) & isnan(x.c.data) & isnan(x.l.data));
x.t.date(badix)=[]; x.t.data(badix)=[];
x.c.date(badix)=[]; x.c.data(badix)=[];
x.l.date(badix)=[]; x.l.data(badix)=[];
x
x.c
x.t
x.origt
datestr(begdt),datestr(enddt)
find_date_ranges(x.t.date,10)
3289*24
find_date_ranges(x.t.date,1)
fprintf(2,'%d,%d,%d,%d,%d,%d,-11,%g,%g,%g\n',datevec(x.c.date(1:3)),x.t.data(1:3),x.c.data(1:3),x.l.date(1:3))
fprintf(2,'%d,%d,%d,%d,%d,%d,-11,%g,%g,%g\n',datevec(x.c.date(1:3)'),x.t.data(1:3),x.c.data(1:3),x.l.date(1:3))
fprintf(2,'%d,%d,%d,%d,%d,%d,-11,%g,%g,%g\n',datevec(x.c.date(1:3)'),x.t.data(1:3)',x.c.data(1:3)',x.l.date(1:3)')
fprintf(2,'%d,%d,%d,%d,%d,%d,-11,%g,%g,%g\n',y(1:3),m(1:3),d(1:3),H(1:3),M(1:3),S(1:3),x.t.data(1:3),x.c.data(1:3),x.l.date(1:3))
[y,m,d,H,M,S] = datevec(x.c.date);
fprintf(2,'%d,%d,%d,%d,%d,%d,-11,%g,%g,%g\n',y(1:3),m(1:3),d(1:3),H(1:3),M(1:3),S(1:3),x.t.data(1:3),x.c.data(1:3),x.l.date(1:3))
fprintf(2,'%d,%d,%d,%d,%d,%d,-11,%g,%g,%g\n',y(1:3)',m(1:3)',d(1:3)',H(1:3)',M(1:3)',S(1:3)',x.t.data(1:3)',x.c.data(1:3)',x.l.date(1:3)')
fprintf(2,'%d,%d,%d,%d,%d,%d,-11,%g,%g,%g\n',y(1:3)',m(1:3)',d(1:3)',H(1:3)',M(1:3)',S(1:3)',x.t.data(1:3)',x.c.data(1:3)',x.l.data(1:3)')
[y,m,d,H,M,S] = datevec(x.t.date);
vec = [y,m,d,H,M,S,x.t.data,x.c.data,x.l.data];
size(vec)
fprintf(2,'%d,%d,%d,%d,%d,%d,-11,%g,%g,%g\n',vec(1:3,:))
fprintf(2,'%d,%d,%d,%d,%d,%d,-11,%g,%g,%g\n',vec(1:3,:)')
x=[]; clear x
x.origt = res.nw_w_btm.adcp_seatemp;
x.origc = res.nw_w_btm.adcp_x;
x.origl = res.nw_w_btm.adcp_l;
begdt = min([x.origt.date(1),  x.origc.date(1),  x.origl.date(1)]);
enddt = max([x.origt.date(end),x.origc.date(end),x.origl.date(end)]);
x.t.date = [begdt:(1/24):enddt]';
x.t.data = interp1(x.origt.date,x.origt.data,x.t.date,[],nan);
x = filter_gaps(x,'origt','t',(1.5/24),0,[],nan);
x.c.date = [begdt:(1/24):enddt]';
x.c.data = interp1(x.origc.date,x.origc.data,x.c.date,[],nan);
x = filter_gaps(x,'origc','c',(1.5/24),0,[],nan);
x.l.date = [begdt:(1/24):enddt]';
x.l.data = interp1(x.origl.date,x.origl.data,x.l.date,[],nan);
x = filter_gaps(x,'origl','l',(1.5/24),0,[],nan);
badix = find(isnan(x.t.data) & isnan(x.c.data) & isnan(x.l.data));
x.t.date(badix)=[]; x.t.data(badix)=[];
x.c.date(badix)=[]; x.c.data(badix)=[];
x.l.date(badix)=[]; x.l.data(badix)=[];
[y,m,d,H,M,S] = datevec(x.t.date);
vec = [y,m,d,H,M,S,x.t.data,x.c.data,x.l.data];
fname = fullfile(corispath,'CRCP_789_Product_1435_data_NSUOC_NW_W_btm.csv');
fid = fopen(fname,'w+');
if fid < 3; error('Unable to open %s!',fname); end;
fprintf(fid,'Date-Time,Instrument depth (m),Sea temperature,Cross-shore (m/s),Along-shore (m/s)\n');
fprintf(fid,'%d,%d,%d,%d,%d,%d,-11,%g,%g,%g\n',vec');
fclose(fid);
vec=[]; x=[]; clear y m d H M S vec x
corispath = get_upwelling_path('CoRIS');
fid
x=[]; clear x
x.origt = res.nw_w_btm.adcp_seatemp;
x.origc = res.nw_w_btm.adcp_x;
x.origl = res.nw_w_btm.adcp_l;
begdt = min([x.origt.date(1),  x.origc.date(1),  x.origl.date(1)]);
enddt = max([x.origt.date(end),x.origc.date(end),x.origl.date(end)]);
x.t.date = [begdt:(1/24):enddt]';
x.t.data = interp1(x.origt.date,x.origt.data,x.t.date,[],nan);
x = filter_gaps(x,'origt','t',(1.5/24),0,[],nan);
x.c.date = [begdt:(1/24):enddt]';
x.c.data = interp1(x.origc.date,x.origc.data,x.c.date,[],nan);
x = filter_gaps(x,'origc','c',(1.5/24),0,[],nan);
x.l.date = [begdt:(1/24):enddt]';
x.l.data = interp1(x.origl.date,x.origl.data,x.l.date,[],nan);
x = filter_gaps(x,'origl','l',(1.5/24),0,[],nan);
badix = find(isnan(x.t.data) & isnan(x.c.data) & isnan(x.l.data));
x.t.date(badix)=[]; x.t.data(badix)=[];
x.c.date(badix)=[]; x.c.data(badix)=[];
x.l.date(badix)=[]; x.l.data(badix)=[];
[y,m,d,H,M,S] = datevec(x.t.date);
vec = [y,m,d,H,M,S,x.t.data,x.c.data,x.l.data];
fname = fullfile(corispath,'CRCP_789_Product_1435_data_NSUOC_NW_W_btm.csv');
fid = fopen(fname,'w+');
if fid < 3; error('Unable to open %s!',fname); end;
fprintf(fid,'Date-Time,Instrument depth (m),Sea temperature,Cross-shore (m/s),Along-shore (m/s)\n');
fprintf(fid,'%d,%d,%d,%d,%d,%d,-11,%g,%g,%g\n',vec');
fclose(fid);
vec=[]; x=[]; clear y m d H M S vec x
res.ne_buoy.seatemp
printf('\n')
fprintf(1,'%3.6g',pi)
fprintf(1,'%3.6g\n',pi)
fprintf(1,'%3.3g\n',pi)
fprintf(1,'%3.4g\n',pi)
fprintf(1,'%3.4f\n',pi)
fprintf(1,'%3.4f\n',pi*10)
fprintf(1,'%3.4f\n',nan)
res.pier_cc.air_t
res.pier_cc.wind_u
export_to_coris
res.pier_cc.wind_u
export_to_coris
stnm
fclose(fid);
res.se_buoy
x.origt
nansummary(res.pier_cc.wind_speed.data)
fwyf1 = get_station_from_station_name('fwyf1'); fwyf1 = load_all_ndbc_data(fwyf1);
scatter_fit_ts(res.pier_cc.wind_speed,fwyf1.ndbc_wind1_speed)
scatter_fit_ts_seasons(res.pier_cc.wind_speed,fwyf1.ndbc_wind1_speed)
scatter_fit_ts_seasons(res.pier_cc.wind_speed,fwyf1.ndbc_wind1_speed,[],[],'Pier','FWY',[],[],true)
res.pier_cc.wind_speed_mps
fwyf1=[]; clear fwyf1
help station_bulk_windstress
fprintf(1,'%.4f\n',nan)
fprintf(1,'%.4f\n',pi)
fprintf(1,'%.4f\n',pi*10)
fprintf(1,'%.3f\n',pi*10)
fprintf(1,'%.3f\n',pi*100)
fprintf(1,'%d\n',pi*100)
fprintf(1,'%d\n',pi)
fprintf(1,'%g\n',pi)
fprintf(1,'%f\n',3)
fprintf(1,'%.0f\n',3)
export_to_coris
fclose(fid);
stnm
x
export_to_coris
fprintf(1,'%g\n',pi)
fprintf(1,'%.0g\n',pi)
fprintf(1,'%.0g\n',pi+0.5)
fprintf(1,'%.0f\n',pi+0.5)
pi+0.5
fprintf(1,'%.1f\n',pi+0.5)
fprintf(1,'%.1f\n',pi+0.05)
fprintf(1,'%.1f\n',pi+0.0)
fprintf(1,'%.1f\n',pi+0.01)
fprintf(1,'%.1f\n',pi+0.005)
fprintf(1,'%.1f\n',pi+0.0075)
fprintf(1,'%.1f\n',pi+0.008)
fprintf(1,'%.1f\n',pi+0.009)
fprintf(1,'%f\n',pi+0.009)
export_to_coris
res.c_buoy.dept
res.c_buoy.depth
export_to_coris
res.c_buoy.air_t
res.ne_buoy.air_t
fmg; plot_ts(res.c_buoy.air_t,res.ne_buoy.air_t,res.pier_cc.air_t);
fmg; plot_ts(res.pier_cc.air_t,res.c_buoy.air_t,res.ne_buoy.air_t);
fmg; plot_ts(res.pier_cc.air_t,res.c_buoy.air_t,res.ne_buoy.air_t); legend('Pier','C','NE')
fprintf(1,'%f\n',-pi+0.009)
fprintf(1,'%f\n',-0)
fprintf(1,'%\n',-0)
fprintf(1,'%g\n',-0)
0==0
0==-0
scatter_fit(res.c_buoy.air_t,fwyf1.nd
fwyf1 = get_station_from_station_name('fwyf1'); fwyf1 = load_all_ndbc_data(fwyf1);
scatter_fit(res.c_buoy.air_t,fwyf1.ndbc_air_t)
scatter_fit_ts(res.c_buoy.air_t,fwyf1.ndbc_air_t)
scatter_fit_ts(res.ne_buoy.air_t,fwyf1.ndbc_air_t)
fwyf1=[]; clear fwyf1
res
clear ans axs badix begdt c0 c5 c10 c15 cbeg cend center_ix cix cstnm cstnms dep deps derr dix e5 e10 e15 e20 e30 enddt fid fname hiix loix ne0 ne5 ne10 ne15 ne20 ne30 ne40 ne45 oterflds se26 se35 se98 stnm sw0 sw5 sw10 sw15 ta tau tbeg tend ts xp
clear ans axs badix begdt c0 c5 c10 c15 cbeg cend center_ix cix corispath cstnm cstnms dep deps derr dix doPlots e5 e10 e15 e20 e30 enddt fid fname hiix loix ne0 ne5 ne10 ne15 ne20 ne30 ne40 ne45 otherflds se26 se35 se98 stnm sw0 sw5 sw10 sw15 ta tau tbeg tend ts xp
dir *fw*m
type load_fwc_fdep.m
type plot_fwc_event.m
datestr(begdt),datestr(enddt)
dir *fw*m
type plot_all_fwri.m
dir *fw*m
type plot_fwc_fdep_bathymetry.m
type plot_fwc_fdep_stations.m
load_fwc_fdep
fwc=[]; fwyf1=[]; lkwf1=[]; clear ans fwc fwyf1 lkwf1
load_fwc_fdep
fmg; plot_ts(fwc.UPDB.sea_t);
fmg; plot_ts(fwc.UPDB.sea_t,fwc.TUCA.sea_t);
fwc
help plot_hires_bathymetry
if ( ~isfield(fwc,'bath') || ~isfield(fwc.bath,'ngdc_hires_bathy') )
fwc.bath = [];
fwc = rmfield(fwc,'bath');
fwc.bath.lon = -80.08;
fwc.bath.lat = 26.62;
end;
fwc.bath = plot_hires_bathymetry(fwc.bath,-[2:2:60],[20e3,70e3],true,@contour);
colorbar;
daspect([1,cosd(fwc.bath.lat),1]);
fwc.bath
fwc.bath = plot_hires_bathymetry(fwc.bath,-[2:2:60],[20e3,70e3],true,@contour);
fwc.bath = plot_hires_bathymetry(fwc.bath,-[2:2:60],[20e3,70e3],true);
fwc.bath = plot_hires_bathymetry(fwc.bath,-[2:2:60],[20e3,70e3]);
fwc.bath = plot_hires_bathymetry(fwc.bath,-[2:2:60],[10e3,30e3]);
dbstop plot_hires_bathymetry
fwc.bath = plot_hires_bathymetry(fwc.bath,-[2:2:60],[20e3,70e3],true,@contour);
bbox
dbstop plot_hires_coastline 32
fwc.bath = plot_hires_bathymetry(fwc.bath,-[2:2:60],[20e3,70e3],true,@contour);
help bboxint
help isinside
useix = find(bbox(1)<=sofla_coast(:,1)&sofla_coast(:,1)<=bbox(2) & bbox(3)<=sofla_coast(:,2)&sofla_coast(:,2)<=bbox(4));
ch=fill(sofla_coast(useix,1), sofla_coast(useix,2), [0.0 0.0 0.0], 'LineWidth',2);
ch=fill(sofla_coast(:,1), sofla_coast(:,2), [0.0 0.0 0.0], 'LineWidth',2);
help bbox2rect
axes(bbox)
axis(bbox)
dbquit
fwc.bath = plot_hires_bathymetry(fwc.bath,-[2:2:60],[20e3,70e3],true,@contour);
dbclear all
fwc.bath = plot_hires_bathymetry(fwc.bath,-[2:2:60],[20e3,70e3],true,@contour);
dbclear all
fwc.bath = plot_hires_bathymetry(fwc.bath,-[2:2:60],[20e3,70e3],true,@contour);
load_fwc_fdep
fwc=[]; fwyf1=[]; lkwf1=[]; clear ans fwc fwyf1 lkwf1
close all
load_fwc_fdep
fwc
res
res.sw_buoy
fwc
fwc.Rodeo
fwc.UPDB
res.sw_buoy
res.pier_cc
res.sw_buoy
fwc
fwc = rmfield(fwc,'bath')
plot_fwc_fdep_bathymetry
fwc = rmfield(fwc,'bath')
plot_fwc_fdep_bathymetry
fwc = rmfield(fwc,'bath')
plot_fwc_fdep_bathymetry
find_date_ranges(res.nw_w_btm.adcp_seatemp)
find_date_ranges(res.nw_w_btm.adcp_seatemp.date)
find_date_ranges(fwc.UPDB.date)
find_date_ranges(fwc.UPDB.sea_t.date)
fmg; plot_ts(res.nw_w_btm,fwc.Noula,fwc.Rodeo,fwc.TUCA,fwc.UPDB); legend('NW/W','Noula','Rodeo','TUCA','UPDB');
fmg; plot_ts(res.nw_w_btm.adcp_seatemp,fwc.Noula.sea_t,fwc.Rodeo.sea_t,fwc.TUCA.sea_t,fwc.UPDB.sea_t); legend('NW/W','Noula','Rodeo','TUCA','UPDB');
datetick3;
fwc
fmg; plot_ts(res.nw_w_btm.adcp_seatemp,fwc.Rodeo.sea_t,fwc.Alpha.sea_t,fwc.Noula.sea_t,fwc.TUCA.sea_t,fwc.UPDB.sea_t); legend('NW/W','Rodeo','Alpha','Noula','TUCA','UPDB');
datetick3;
plot_ts(lkwf1.ndbc_air_t,'k-')
plot_fwc_fdep_bathymetry
load_all_sefcri
scatter_fit_ts(fwc.UPDB.sea_t,updb.t)
fwc.UPDB.sea_t,updb.t
pela
station_dist(pela,fwc.UPDB)
station_dist(evca,fwc.UPDB)
station_dist(slib,fwc.UPDB)
station_dist(updb,fwc.UPDB)
clear ans
clear bc1 bc3 dc1 dc3 evca mc2 pb1 pb2 pela slib updb
fwc=[]; fwyf1=[]; lkwf1=[]; clear ans fwc fwyf1 lkwf1
plot_all_sefcri_analyses
fmg; plot_ts(res.nw_w_btm.adcp_seatemp,updb.t); legend('NW/W','UPDB');
datetick3;
f
fmg; plot_ts(res.nw_w_btm.adcp_seatemp,pb2.t,updb.t); legend('NW/W','PB2','UPDB');
datetick3;
plot_ts(lkwf1.ndbc_air_t,'k-')
plot_ts(fwyf1.ndbc_air_t,'k-')
plot_ts(lkwf1.ndbc_air_t,'r-')
plot_ts(lkwf1.ndbc_air_t,'-','color',[.5,.5,.5])
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
analyze_sefcri
type plot_sefcri_upwellings.m
fh=fmg;
plot_ts(ftpf1.ndbc_sea_t, cnnf1.ndbc_sea_t, canf1.ndbc_sea_t, sauf1.ndbc_sea_t, ...
ftpf1.ndbc_sigwavehgt,'m.-', ...
sauf1.ndbc_sigwavehgt,'m-', ...
lkwf1.ndbc_wind1_speed_3_h_median,'.-','Color',[0.75,0.75,0], ...
sauf1.ndbc_wind1_speed_3_h_median,'-','Color',[0.75,0.75,0], ...
sefcri.updb.t,'g-', sefcri.pb2.t,'y-','Color',[0.75,0.75,0], sefcri.dc3.t,'k-');
grid minor;
legend('FtP','CnN','Can','Sau','FtP Hs','Sau Hs','LkW U_3_h','Sau U_3_h',...
'UPDB','PB2','DC3', 'Location','EastOutside');
xlim(datenum([2008,2013],[5,9],1)); datetick3;
for dt = datenum(2010,[[5:9]+(0*12),[5:9]+(1*12),[5:9]+(2*12),[5:9]+(3*12)],1);
xlim([dt,dt+32]);
datetick3;
disp('Hit "Enter" to move to next time-frame');
pause;
if ( ishandle(fh) )
figure(fh);
else
disp('Quit');
break;
end;
end;
find_date_ranges(sefcri.updb.hourly_t.date,10)
find_date_ranges(sefcri.pb2.hourly_t.date,10)
find_date_ranges(sefcri.dc3.hourly_t.date,10)
fh=fmg;
plot_ts(ftpf1.ndbc_sea_t, cnnf1.ndbc_sea_t, canf1.ndbc_sea_t, sauf1.ndbc_sea_t, ...
ftpf1.ndbc_sigwavehgt,'m.-', ...
sauf1.ndbc_sigwavehgt,'m-', ...
lkwf1.ndbc_wind1_speed_3_h_median,'.-','Color',[0.75,0.75,0], ...
sauf1.ndbc_wind1_speed_3_h_median,'-','Color',[0.75,0.75,0], ...
sefcri.updb.t,'g-', sefcri.pb2.t,'y-',sefcri.bc3.t,'k.-',sefcri.dc3.t,'k-');
grid minor;
legend('FtP','CnN','Can','Sau','FtP Hs','Sau Hs','LkW U_3_h','Sau U_3_h',...
'UPDB','PB2','BC3','DC3', 'Location','EastOutside');
xlim(datenum([2008,2013],[5,9],1)); datetick3;
% % May-September, 2010-2013
% for dt = datenum(2010,[[5:9]+(0*12),[5:9]+(1*12),[5:9]+(2*12),[5:9]+(3*12)],1);
% July-September, 2010-2013
for dt = datenum(2010,[[7:9]+(0*12),[7:9]+(1*12),[7:9]+(2*12),[7:9]+(3*12)],1);
xlim([dt,dt+32]);
datetick3;
disp('Hit "Enter" to move to next time-frame');
pause;
if ( ishandle(fh) )
figure(fh);
else
disp('Quit');
break;
end;
end;
fh=fmg;
plot_ts(ftpf1.ndbc_sea_t, cnnf1.ndbc_sea_t, canf1.ndbc_sea_t, sauf1.ndbc_sea_t, ...
ftpf1.ndbc_sigwavehgt,'m.-', sauf1.ndbc_sigwavehgt,'m-', ...
lkwf1.ndbc_wind1_speed_3_h_median,'.-','Color',[0.75,0.75,0], ...
sauf1.ndbc_wind1_speed_3_h_median,'-','Color',[0.75,0.75,0], ...
lkwf1.ndbc_air_t,'k-', ...
sefcri.updb.t,'g-', sefcri.pb2.t,'b-',sefcri.bc3.t,'g.-',sefcri.dc3.t,'b.-');
grid minor;
legend('FtP','CnN','Can','Sau','FtP Hs','Sau Hs','LkW U_3_h','Sau U_3_h',...
'LkW T_a','UPDB','PB2','BC3','DC3', 'Location','EastOutside');
xlim(datenum([2008,2013],[5,9],1)); datetick3;
% % May-September, 2010-2013
% for dt = datenum(2010,[[5:9]+(0*12),[5:9]+(1*12),[5:9]+(2*12),[5:9]+(3*12)],1);
% July-September, 2010-2013
for dt = datenum(2010,[[7:9]+(0*12),[7:9]+(1*12),[7:9]+(2*12),[7:9]+(3*12)],1);
xlim([dt,dt+32]);
datetick3;
disp('Hit "Enter" to move to next time-frame');
pause;
if ( ishandle(fh) )
figure(fh);
else
disp('Quit');
break;
end;
end;
fh=fmg;
plot_ts(ftpf1.ndbc_sea_t, cnnf1.ndbc_sea_t, canf1.ndbc_sea_t, sauf1.ndbc_sea_t, ...
ftpf1.ndbc_sigwavehgt,'m.-', sauf1.ndbc_sigwavehgt,'m-', ...
lkwf1.ndbc_wind1_speed_3_h_median,'.-','Color',[0.75,0.75,0], ...
sauf1.ndbc_wind1_speed_3_h_median,'-','Color',[0.75,0.75,0], ...
lkwf1.ndbc_air_t,'k-', ...
sefcri.updb.t,'g-', sefcri.pb2.t,'g.-',sefcri.bc3.t,'b-',sefcri.dc3.t,'b.-');
grid minor;
legend('FtP','CnN','Can','Sau','FtP Hs','Sau Hs','LkW U_3_h','Sau U_3_h',...
'LkW T_a','UPDB','PB2','BC3','DC3', 'Location','EastOutside');
xlim(datenum([2008,2013],[5,9],1)); datetick3;
% % May-September, 2010-2013
% for dt = datenum(2010,[[5:9]+(0*12),[5:9]+(1*12),[5:9]+(2*12),[5:9]+(3*12)],1);
% July-September, 2010-2013
for dt = datenum(2010,[[7:9]+(0*12),[7:9]+(1*12),[7:9]+(2*12),[7:9]+(3*12)],1);
xlim([dt,dt+32]);
datetick3;
disp('Hit "Enter" to move to next time-frame');
pause;
if ( ishandle(fh) )
figure(fh);
else
disp('Quit');
break;
end;
end;
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
analyze_sefcri
daspect([1,cosd(25),1]);
daspect
daspect([1,1,1]);
axis default
help axis
axis auto
sefcri
clear ans b1 b3 d1 d3 ix locs m p1 p2 s u
sefcri.slib = plot_hires_bathymetry(sefcri.slib,-[0:2:150],[30e3,15e3]);
plot_all_fwri;
sefcri.slib = rmfield(sefcri.slib,'ngdc_hires_bathy')
sefcri.slib = plot_hires_bathymetry(sefcri.slib,-[0:2:150],[30e3,15e3]);
plot_all_fwri;
sefcri.slib = rmfield(sefcri.slib,'ngdc_hires_bathy')
sefcri.slib = plot_hires_bathymetry(sefcri.slib,-[0:2:150],[15e3,30e3]);
plot_all_fwri;
sefcri.slib = rmfield(sefcri.slib,'ngdc_hires_bathy')
sefcri.slib = plot_hires_bathymetry(sefcri.slib,-[0:2:150],[15e3,35e3]);
plot_all_fwri;
sefcri.slib = rmfield(sefcri.slib,'ngdc_hires_bathy')
sefcri.slib = plot_hires_bathymetry(sefcri.slib,-[0:2:150],[15e3,40e3]);
plot_all_fwri;
sefcri.updb = rmfield(sefcri.updb,'ngdc_hires_bathy')
sefcri.updb = plot_hires_bathymetry(sefcri.updb,-[0:5:50,100:50:300],[40e3,160e3],[],@contourf);
plot_all_fwri;
axis square
axis tight
daspect([1,cosd(25),1]);
clear ans b1 b3 d1 d3 ix locs m p1 p2 s u
sefcri
fieldnames(sefcri)
sefcri.dc3
sefcri.dc1
sefcri.bc1
sefcri.bc3
sefcri.pb1
sefcri.pb2
sefcri.mc2
sefcri.updb
sefcri.slib
sefcri.evca
sefcri
sefcri.pela
sefcri.updb
station_dist(sefcri.updb,sauf1)
station_dist(sefcri.updb,mlrf1)
station_dist(sefcri.updb,ftpf1)
station_dist(sefcri.updb,canf1)
station_dist(sefcri.updb,cnnf1)
dir *fw*m
help grep
help vrep
help grepfields
edit log4j.properties
edit javaclasspath.m
edit javaaddpath.m
edit javacomponent.m
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
extract_sfomc_etc
q
close all
clear ans axs badix begdt c0 c5 c10 c15 cbeg cend center_ix cix corispath cstnm cstnms dep deps derr dix doPlots e5 e10 e15 e20 e30 enddt fid fname hiix loix ne0 ne5 ne10 ne15 ne20 ne30 ne40 ne45 otherflds se26 se35 se98 stnm sw0 sw5 sw10 sw15 ta tau tbeg tend ts xp
analyze_sefcri
reviewanim([],0,0)
reviewanim([],0,0,0)
close all
clear ans b1 b3 d1 d3 dt fh ix locs m matfname p1 p2 s u
help pause
help jumpax
type jumpax
plot_sefcri_upwellings
help any
size(any(res.e_buoy.seatemp.prof<17))
size(any(res.e_buoy.seatemp.prof<17,2))
numel(find(any(res.e_buoy.seatemp.prof<17,2)))
numel(find(any(res.ne_buoy.seatemp.prof<17,2)))
res.e_buoy.seatemp
numel(find(any(res.se_buoy.seatemp.prof<17,2)))
numel(find(any(res.sw_buoy.seatemp.prof<17,2)))
fmg; contour(res.e_buoy.seatemp.date,res.e_buoy.seatemp.depths,res.e_buoy.seatemp.prof); datetick3;
fmg;
fmg; contourf(res.e_buoy.seatemp.date,res.e_buoy.seatemp.depths,res.e_buoy.seatemp.prof'); datetick3; colorbar;
datetick3;
help contour_field
fmg; contour(res.e_buoy.seatemp);
dbstop if error
fmg; contour(res.e_buoy.seatemp);
dbup
dbquit
fmg; contour_field(res.e_buoy.seatemp);
dbquit
whos
clear ans
fmg; contour_field(res.e_buoy.seatemp,'prof');
dbquit
dbclear all
help contourf
help contour
fmg; contour_field(res.e_buoy.seatemp,'prof');
fmg; contour_profile(res.e_buoy.seatemp,'prof');
datetick3
fmg; contour_profile(res.e_buoy.seatemp,'prof');
datetick3;
fmg; contour_profile(res.ne_buoy.seatemp,'prof');
shading interp
colorbar
shading facet
shading faceted
shading flat
shading interp
help contourf
help contour
help contourf
help contour
help contourf
set(gca)
set(gca,'linew',0)
set(gca,'linew',0.1)
set(gca,'color',0)
set(gca,'color',[0,0,0])
set(gca,'color',[1,1,1])
fmg; plot_ts(res.ne_buoy.at_45m.sbe_seatemp);
fmg; plot_ts(res.e_buoy.at_45m.sbe_seatemp);
fmg; plot_ts(res.e_buoy.at_30m.sbe_seatemp);
fmg; plot_ts(res.ne_buoy.at_30m.sbe_seatemp);
fmg; plot_ts(res.ne_buoy.at_30m.sbe_seatemp); plot_ts(res.e_buoy.at_30m.sbe_seatemp); plot_ts(res.se_buoy.at_35m.sbe_seatemp); plot_ts(res.sw_buoy.at_15m.sbe_seatemp);
fmg; plot_ts(res.ne_buoy.at_30m.sbe_seatemp,res.e_buoy.at_30m.sbe_seatemp,res.se_buoy.at_35m.sbe_seatemp,res.sw_buoy.at_15m.sbe_seatemp);
res.sw_buoy.at_
fmg; plot_ts(res.se_buoy.at_98m,res.ne_buoy.at_45m);
fmg; plot_ts(res.se_buoy.at_98m.sbe_seatemp,res.ne_buoy.at_45m.sbe_seatemp);
fmg; plot_ts(res.se_buoy.at_35m.sbe_seatemp,res.ne_buoy.at_45m.sbe_seatemp);
scatter_fit_ts_seasons(res.ne_buoy.at_45m.sbe_seatemp,res.c_buoy.at_15m.sbe_seatemp)
dbstop if error
scatter_fit_ts_seasons(res.ne_buoy.at_45m.sbe_seatemp,res.c_buoy.at_15m.sbe_seatemp)
( min(diff(ix1)) < 0 || min(diff(ix2)) < 0 )
( min(diff(ix1)) < 0),( min(diff(ix2)) < 0 )
size(ix1)
size(ix2)
dbquit
clear axLms fh fix1 fix2 fld1 fld2 ix1 ix2 lag legendArgs statToPlot tol xlbl ylbl
clear axLms dt fh fix1 fix2 fld1 fld2 inp ix1 ix2 lag legendArgs statToPlot tol xlbl ylbl
clear ans axLms dt fh fix1 fix2 fld1 fld2 inp ix1 ix2 lag legendArgs statToPlot tol xlbl ylbl
fmg; plot_ts(res.se_buoy.at_98m,res.ne_buoy.at_45m);
dbquit
fmg; plot_ts(res.se_buoy.at_98m.sbe_seatemp,res.ne_buoy.at_45m.sbe_seatemp);
fmg; plot_ts(res.se_buoy.at_98m.sbe_seatemp,res.se_buoy.at_35m.sbe_seatemp,res.ne_buoy.at_45m.sbe_seatemp);
fmg; plot_ts(res.se_buoy.at_98m.sbe_seatemp,res.se_buoy.at_35m.sbe_seatemp,res.ne_buoy.at_45m.sbe_seatemp); legend('98m SE','35m SE','45m NE');
fmg; plot_ts(res.se_buoy.at_98m.sbe_seatemp,res.se_buoy.at_35m.sbe_seatemp,res.ne_buoy.at_45m.sbe_seatemp,res.ne_buoy.at_30m.sbe_seatemp); legend('98m SE','35m SE','45m NE','30m NE');
datetick3;
find_date_ranges(res.se_buoy.at_98m.sbe_seatemp.date,1)
find_date_ranges(res.se_buoy.at_35m.sbe_seatemp.date,1)
find_date_ranges(res.se_buoy.at_26m.sbe_seatemp.date,1)
res.se_buoy
res.sw_buoy
res.se_buoy
min(res.ne_buoy.seatemp.prof(:))
min(res.e_buoy.seatemp.prof(:))
find_date_ranges(res.e_buoy.at_30m.sbe_seatemp.date,1)
find_date_ranges(res.ne_buoy.at_30m.sbe_seatemp.date,1)
res.e_buoy
res.ne_buoy
station_dist(res.ne_buoy,res.e_buoy)
280+165+50
0.54*340
ans+50
280+165+235
280*3
840+165+235
ans-280-280
sefcri
for cstnm=fieldnames(sefcri); stnm=cstnm{:}; numel(find(sefcri.(stnm).t.data<17)), end;
for cstnm=fieldnames(sefcri)'; stnm=cstnm{:}; numel(find(sefcri.(stnm).t.data<17)), end;
timenow
for cstnm=fieldnames(sefcri)'; stnm=cstnm{:}; disp({stnm,numel(find(sefcri.(stnm).t.data<17))}); end;
for cstnm=fieldnames(sefcri)'; stnm=cstnm{:}; disp({stnm,numel(find(sefcri.(stnm).t.data<17&get_season(sefcri.(stnm).t.date)==3))}); end;
for cstnm=fieldnames(sefcri)'; stnm=cstnm{:}; disp({stnm,numel(find(sefcri.(stnm).t.data<17&ismember(get_month(sefcri.(stnm).t.date),[6:9])))}); end;
type get_season
for cstnm=fieldnames(sefcri)'; stnm=cstnm{:}; disp({stnm,numel(find(sefcri.(stnm).t.data<17&ismember(get_month(sefcri.(stnm).t.date),[6:9])))}); end;
for cstnm=fieldnames(sefcri)'; stnm=cstnm{:}; disp({stnm,numel(find(sefcri.(stnm).t.data<17&ismember(get_month(sefcri.(stnm).t.date),[5:9])))}); end;
for cstnm=fieldnames(sefcri)'; stnm=cstnm{:}; disp({stnm,numel(find(sefcri.(stnm).t.data<17&ismember(get_month(sefcri.(stnm).t.date),[5:10])))}); end;
for cstnm=fieldnames(sefcri)'; stnm=cstnm{:}; disp({stnm,numel(unique(floor(sefcri.(stnm).t.date(sefcri.(stnm).t.data<17&ismember(get_month(sefcri.(stnm).t.date),[5:10])))))}); end;
for cstnm=fieldnames(sefcri)'; stnm=cstnm{:}; disp({stnm,numel(unique(floor(sefcri.(stnm).t.date(sefcri.(stnm).t.data<17&ismember(get_month(sefcri.(stnm).t.date),[5:10]))))) / numel(unique(floor(sefcri.(stnm).t.date))) }); end;
for cstnm=fieldnames(sefcri)'; stnm=cstnm{:}; disp({stnm,numel(unique(floor(sefcri.(stnm).t.date(sefcri.(stnm).t.data<17&ismember(get_month(sefcri.(stnm).t.date),[1:10]))))) / numel(unique(floor(sefcri.(stnm).t.date))) }); end;
for cstnm=fieldnames(sefcri)'; stnm=cstnm{:}; disp({stnm,numel(unique(floor(sefcri.(stnm).t.date(sefcri.(stnm).t.data<17&ismember(get_month(sefcri.(stnm).t.date),[5:10])))))}); end;
for cstnm=fieldnames(sefcri)'; stnm=cstnm{:}; disp({stnm,numel(unique(floor(sefcri.(stnm).t.date(sefcri.(stnm).t.data<17&ismember(get_month(sefcri.(stnm).t.date),[1:10])))))}); end;
for cstnm=fieldnames(sefcri)'; stnm=cstnm{:}; disp({stnm,numel(unique(floor(sefcri.(stnm).t.date(sefcri.(stnm).t.data<17&ismember(get_month(sefcri.(stnm).t.date),[1:12])))))}); end;
for cstnm=fieldnames(sefcri)'; stnm=cstnm{:}; disp({stnm,numel(unique(floor(sefcri.(stnm).t.date(sefcri.(stnm).t.data<17&ismember(get_month(sefcri.(stnm).t.date),[5:12])))))}); end;
for cstnm=fieldnames(sefcri)'; stnm=cstnm{:}; disp({stnm,numel(unique(floor(sefcri.(stnm).t.date(sefcri.(stnm).t.data<17&ismember(get_month(sefcri.(stnm).t.date),[5:11])))))}); end;
sefcri.pb1
sefcri.pb2
for cstnm=fieldnames(sefcri)'; stnm=cstnm{:}; disp({stnm,numel(unique(floor(sefcri.(stnm).t.date(sefcri.(stnm).t.data<17&ismember(get_month(sefcri.(stnm).t.date),[5:10])))))}); end;
for cstnm=fieldnames(sefcri)'; stnm=cstnm{:}; disp({stnm,numel(unique(floor(sefcri.(stnm).t.date(sefcri.(stnm).t.data<17&ismember(get_month(sefcri.(stnm).t.date),[5:12])))))}); end;
for cstnm=fieldnames(sefcri)'; stnm=cstnm{:}; disp({stnm,numel(unique(floor(sefcri.(stnm).t.date(sefcri.(stnm).t.data<17&ismember(get_month(sefcri.(stnm).t.date),[1:10])))))}); end;
for cstnm=fieldnames(sefcri)'; stnm=cstnm{:}; disp({stnm,numel(unique(floor(sefcri.(stnm).t.date(sefcri.(stnm).t.data<17&ismember(get_month(sefcri.(stnm).t.date),[5:12])))))}); end;
for cstnm=fieldnames(sefcri)'; stnm=cstnm{:}; disp({stnm,numel(unique(floor(sefcri.(stnm).t.date(sefcri.(stnm).t.data<17&ismember(get_month(sefcri.(stnm).t.date),[1:4])))))}); end;
for cstnm=fieldnames(sefcri)'; stnm=cstnm{:}; disp({stnm,numel(unique(floor(sefcri.(stnm).t.date(sefcri.(stnm).t.data<17&ismember(get_month(sefcri.(stnm).t.date),[1:3])))))}); end;
for cstnm=fieldnames(sefcri)'; stnm=cstnm{:}; disp({stnm,numel(unique(floor(sefcri.(stnm).t.date(sefcri.(stnm).t.data<17&ismember(get_month(sefcri.(stnm).t.date),[4])))))}); end;
fmg; plot_ts(sefcri.pb1.t,sefcri.pb2.t);
datetick3;
plot_ts(res.nw_w_btm.adcp_seatemp,'k-')
find_date_ranges(res.nw_w_btm.adcp_seatemp.date)
find_date_ranges(res.nw_w_btm.adcp_seatemp.date,10)
find_date_ranges(mc2.t.date,10)
find_date_ranges(sefcri.mc2.t.date,10)
find_date_ranges(res.nw_w_btm.sbe_seatemp.date,10)
plot_ts(sefcri.mc2.t,'k-')
find_date_ranges(lkwf1.ndbc_air_t.date,10)
find_date_ranges(lkwf1.ndbc_air_t.date,1)
plot_ts(lkwf1.ndbc_air_t,'k-')
find_date_ranges(fwyf1.ndbc_air_t.date,1)
plot_ts(fwyf1.ndbc_air_t,'k-')
station_dist(sefcri.pb2,ftpf1)
station_dist(sefcri.pb2,lkwf1)
station_dist(sefcri.pb2,fwyf1)
station_dist(sefcri.pb2,sauf1)
station_dist(sefcri.pb2,canf1)
station_dist(sefcri.pb2,cnnf1)
find_date_ranges(ftpf1.ndbc_air_t.date,1)
ftpf1
find_date_ranges(ftpf1.ndbc_sigwavehgt.date,1)
station_dist(sefcri.pb2,ftpf1)
plot_ts(fwyf1.ndbc_air_t,'k-')
station_dist(sefcri.pb2,cnnf1)
station_dist(sefcri.pb2,canf1)
station_dist(sefcri.pb2,sauf1)
for cstnm=fieldnames(sefcri)'; stnm=cstnm{:}; disp({stnm,numel(unique(floor(sefcri.(stnm).t.date(sefcri.(stnm).t.data<17&ismember(get_month(sefcri.(stnm).t.date),[4])))))}); end; clear cstnm stnm
for cstnm=fieldnames(sefcri)'; stnm=cstnm{:}; disp({stnm,numel(unique(floor(sefcri.(stnm).t.date(sefcri.(stnm).t.data<17&ismember(get_month(sefcri.(stnm).t.date),[4])))))}); end; clear ans cstnm stnm
plot_ts(ftpf1.ndbc_sigwavehgt,'m:')
plot_ts(ts_op(ftpf1.ndbc_sigwavehgt,28,'+'),'m:')
plot_ts(ts_op(ftpf1.ndbc_sigwavehgt,27,'+'),'m:')
dbstop verify_variable
pd ../../FACE
which face
face
analyze_face_nutrients
dbquit
clear ans dtendoff endoff
xl=[]; clear xl
face=[]; clear face
analyze_face_nutrients
face
face=[]; clear face
xl=[]; clear xl
clear ans dtendoff endoff
clear ans coastpath dix dtendoff endoff hlat hlon ix laterr latix lonerr lonix
clear ans coastpath dix dtendoff endoff hlat hlon ix laterr latix lonerr lonix non_boil_stnix
dbquit
analyze_face_nutrients
scatter_fit(face.t,face.n,'T','N0_x');
face
datestr(face.dt)
datestr(face.dt(end))
face
%-- 2/27/2016 6:02 PM --%
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
analyze_face_nutrients\
analyze_face_nutrients
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
analyze_face_nutrients
datestr(face.dt)
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
analyze_face_nutrients
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
quit
%-- 2/29/2016 1:29 PM --%
analyze_face_nutrients
nansummary(face.lon)
nansummary(face.lat)
nansummary(face.lat(ix))
nansummary(face.lon(ix))
bath.lon = mean(face.lon(ix)); bath.lat = mean(face.lat(ix));
bath = plot_hires_bathymetry(bath,[10e3,60e3],-[0:5:200]);
bath
help plot_hires_bathymetry
bath = plot_hires_bathymetry(bath,-[0:5:200],[10e3,60e3]);
plot(face.lon(ix),face.lat(ix),'.',face.lon(non_boil_stnix),face.lat(non_boil_stnix),'r.');
bath=[]; clear bath
bath.lon = mean(face.lon(ix)); bath.lat = mean(face.lat(ix));
bath = plot_hires_bathymetry(bath,[20e3,50e3],-[0:5:80,100:20:300]);
plot(face.lon(ix),face.lat(ix),'.',face.lon(non_boil_stnix),face.lat(non_boil_stnix),'r.');
%non_boil_stnix = union(strmatch('PE',face.stnm),strmatch('B',face.stnm));
non_boil_stnix = union(strmatch('PE',face.stnm),strmatch('B',face.stnm));
bath.lon = mean(face.lon(ix)); bath.lat = mean(face.lat(ix));
bath = plot_hires_bathymetry(bath,-[0:5:80,100:20:300],[20e3,50e3]);
plot(face.lon(ix),face.lat(ix),'.',face.lon(non_boil_stnix),face.lat(non_boil_stnix),'r.');
bath.lon = mean(face.lon(ix)); bath.lat = mean(face.lat(ix));
bath = plot_hires_bathymetry(bath,-[0:5:80,100:20:300,400:100:800],[20e3,60e3]);
plot(face.lon(ix),face.lat(ix),'.',face.lon(non_boil_stnix),face.lat(non_boil_stnix),'r.');
bath=[]; clear bath
clear ans coastpath dix dtendoff endoff hlat hlon ix laterr latix lonerr lonix non_boil_stnix
clear MYDOCS
xl=[]; clear xl
face=[]; clear face
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
analyze_face_nutrients
plot_hires_coastline([mean(face.lon),mean(face.lat)],[60e3,150e3]);
[mean(face.lon),mean(face.lat)],
plot_hires_coastline([nanmean(face.lon),nanmean(face.lat)],[60e3,150e3]);
[nanmean(face.lon),nanmean(face.lat)],
help plot_hires_coastline
lobath.lon=nanmean(bath.lon); lobath.lat=nanmean(bath.lat);
lobath = plot_hires_bathymetry(lobath,[60e3,150e3]);
lobath = plot_hires_bathymetry(lobath,[150e3,60e3]);
lobath=[]; clear lobath
lobath.lon=nanmean(bath.lon); lobath.lat=nanmean(bath.lat);
lobath = plot_hires_bathymetry(lobath,[150e3,60e3]);
lobath=[]; clear lobath
lobath.lon=nanmean(bath.lon); lobath.lat=nanmean(bath.lat);
lobath = plot_hires_bathymetry(lobath,[150e3,150e3]);
lobath=[]; clear lobath
clear ans coastpath dix dtendoff endoff hlat hlon ix laterr latix lonerr lonix non_boil_stnix
xl=[]; clear xl
lobath.lon=nanmean(bath.lon); lobath.lat=nanmean(bath.lat);
lobath = plot_hires_bathymetry(lobath,[150e3,150e3],-[0:50:1000]);
lobath = plot_hires_bathymetry(lobath,-[0:50:1000],[150e3,150e3]);
lobath=[]; clear lobath
pack
lobath.lon=nanmean(bath.lon); lobath.lat=nanmean(bath.lat);
lobath = plot_hires_bathymetry(lobath,-[0:50:1000],[150e3,60e3]);
lobath=[]; clear lobath
pack
lobath.lon=nanmean(bath.lon); lobath.lat=nanmean(bath.lat);
lobath = plot_hires_bathymetry(lobath,-[0:50:1000],[60e3,150e3]);
lobath=[]; clear lobath
pack
lobath.lon=nanmean(face.lon); lobath.lat=nanmean(face.lat);
lobath = plot_hires_bathymetry(lobath,-[0:50:1000],[60e3,150e3]);
face.z(face.t<25)
face.z(face.t<20)
nansummary(face.z(face.t<25))
nansummary(face.z(face.t<26))
nansummary(face.z(face.t<27))
nansummary(face.z(face.t<28))
zix=find(face.t(ix)<25);
ix = find(face.lat>25);
zix=find(face.t(ix)<25);
nansummary(face.z(zix))
nansummary(face.z(ix(zix)))
face.z(1:10)
face.h(1:10)
scatter_fit_ts(face.t,face.z)
scatter_fit(face.t,face.z)
scatter_fit(face.t(ix),face.z(ix))
non_boil_stnix = union(strmatch('PE',face.stnm),strmatch('B',face.stnm));
scatter_fit(face.t(non_boil_stnix),face.z(non_boil_stnix))
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
analyze_face_nutrients
reviewanim([],0,0,0)
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
analyze_face_nutrients
xl
xl.textdata
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
analyze_face_nutrients
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
analyze_face_nutrients
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
%xl = importdata(fullfile(get_relative_path,'Master_Data_Sheet_21aug09.xlsx'));
xl = importdata(fullfile(get_relative_path,'Master_Data_Sheet_01Feb10.xls'));
sht = 'Sheet3';
endoff=0; dtendoff=endoff+21;
face.dt=datenum(xl.textdata.(sht)(4:end-dtendoff,1))+xl.data.(sht)(1:end-endoff,1);
face.stnm=xl.textdata.(sht)(4:end-dtendoff,3);
face.lat=xl.data.(sht)(1:end-endoff,3)+(xl.data.(sht)(1:end-endoff,4)./60);
face.lon=-( xl.data.(sht)(1:end-endoff,5)+(xl.data.(sht)(1:end-endoff,6)./60) );
face.z=xl.data.(sht)(1:end-endoff,7);
face.t=xl.data.(sht)(1:end-endoff,8);
face.s=xl.data.(sht)(1:end-endoff,9);
face.n=xl.data.(sht)(1:end-endoff,13);
face.no2=xl.data.(sht)(1:end-endoff,14);
face.no3=face.n-face.no2;
face.nh4=xl.data.(sht)(1:end-endoff,15);
face.p=xl.data.(sht)(1:end-endoff,16);
face.si=xl.data.(sht)(1:end-endoff,17);
face.n(1:10)
face.no2(1:10)
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
analyze_face_nutrients
face.lat
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
analyze_face_nutrients
bath=[]; face=[]; hlat=[]; hlon=[];
clear all
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
analyze_face_nutrients
help plot_hires_bathymetry
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
analyze_face_nutrients
help scatter_curve_fit
scatter_curve_fit(face.t(ix),face.n(ix),'exp1')
scatter_fit(face.t,face.s)
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
analyze_face_nutrients
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
quit
%-- 3/2/2016 4:41 PM --%
analyze_face_nutrients
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
analyze_face_nutrients
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
analyze_face_nutrients
face.lon
face.lat
sampleCruise
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
analyze_face_nutrients
scatter_fit(face.t,face.den)
non_boil_stnix
scatter_fit(face.t,face.den)
scatter_curve_fit(face.t,face.den,'exp1')
scatter_curve_fit(face.t,face.den,'lo1')
scatter_curve_fit(face.t,face.den,'log1')
help fittype
scatter_curve_fit(face.t,face.den,'log')
scatter_curve_fit(face.t,face.den,'logit')
lookfor logit
lookfor logistic
help loglc
figure(13)
figure(15)
figure(14)
cruise
help fittype
doc fittype
scatter_curve_fit(face.t,face.den,'exp1')
scatter_curve_fit(face.t,face.den,'exp1','T','\rho')
scatter_curve_fit(face.t,face.den,'power2','T','\rho')
scatter_curve_fit(face.t,face.den,'power1','T','\rho')
scatter_curve_fit(face.t,face.den,'exp2','T','\rho')
scatter_curve_fit(face.t,face.den,'gaussian','T','\rho')
scatter_curve_fit(face.t,face.den,'gauss1','T','\rho')
scatter_curve_fit(face.t,face.den,'gauss2','T','\rho')
scatter_curve_fit(face.t,face.den,'gauss3','T','\rho')
face
face.stnm
face
timenow
reviewanim([],0,0,0)
dbstop if error
reviewanim([],0,0,0)
fhs
which isfinite
fhs(1)
fhs.Name
fhs
fhs(3)
help isvalid
isvalid(fhs)
isvalid(nan)
dbquit
whos
face
oldface=face;
clear oldface
nansummary(face.n)
fmg; hist(face.n,100)
fmg; hist(face.n,100); titlename('N');
close all
fmg; hist(face.n,100); titlename('N');
fmg; hist(face.p,100); titlename('P');
bath=[]; face=[]; hlat=[]; hlon=[];
clear all
pack
analyze_face_nutrients
which reviewanim
reviewanim([],0,0,0)
dbstop if error
reviewanim([],0,0,0)
timenow
reviewanim([],0,0)
dbquit
reviewanim([],0,0)
dbquit
reviewanim([],0,0)
reviewanim([],0,0,0)
fhs = get(0,'child')
ishandle(fhs(1))
ishandle(fhs(3))
fhs
ishandle(fhs(3))
ishandle(fhs(2))
reviewanim([],0,0,0)
face
fmg; hist(face.n,100); titlename('N');
fmg; hist(face.p,100); titlename('P');
reviewanim([],0,0,0)
reviewanim([3:17],0,0,0)
reviewanim([3:19],0,0,0)
reviewanim([3:16],0,0,0)
close(3:16)
fmg; hist(face.n,100); titlename('N');
fmg; hist(face.p,100); titlename('P');
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
xl = importdata(fullfile(get_relative_path,'Master_Data_Sheet_21aug09.xlsx'));
%xl = importdata(fullfile(get_relative_path,'Master_Data_Sheet_01Feb10.xls'));
sht = 'EddyXmptE1';
%sht = 'EddyXmptE2';
endoff=0; dtendoff=endoff+5;
xl.textdata.EddyXmptE1
xl.data.EddyXmptE1
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
xl = importdata(fullfile(get_relative_path,'Master_Data_Sheet_21aug09.xlsx'));
%xl = importdata(fullfile(get_relative_path,'Master_Data_Sheet_01Feb10.xls'));
for sht={'EddyXmptE1','EddyXmptE2'};
disp(sht);
%endoff=0; dtendoff=endoff+5;
endoff=0; dtendoff=endoff;
face.dt=datenum(xl.textdata.(sht)(3:end-dtendoff,1));
face.stnm=xl.textdata.(sht)(3:end-dtendoff,2);
face.lat=xl.data.(sht)(3:end-endoff,1);
face.lon=xl.data.(sht)(3:end-endoff,2);
face.z=xl.data.(sht)(3:end-endoff,3);
face.t=xl.data.(sht)(3:end-endoff,4);
face.s=xl.data.(sht)(3:end-endoff,5);
face.tss=xl.data.(sht)(3:end-endoff,6);
face.chl=xl.data.(sht)(3:end-endoff,7);
face.pha=xl.data.(sht)(3:end-endoff,8);
face.n=xl.data.(sht)(3:end-endoff,9);
face.p=xl.data.(sht)(3:end-endoff,10);
face.si=xl.data.(sht)(3:end-endoff,11);
end;
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
help struct
xl = importdata(fullfile(get_relative_path,'Master_Data_Sheet_21aug09.xlsx'));
%xl = importdata(fullfile(get_relative_path,'Master_Data_Sheet_01Feb10.xls'));
face = struct('dt',[],'stnm',[],'lat',[],'lon',[],'z',[],'t',[],'s',[],...
'tss',[],'chl',[],'pha',[],'n',[],'p',[],'si',[]);
for csht={'EddyXmptE1','EddyXmptE2'};
sht = csht{:};
disp(sht);
%endoff=0; dtendoff=endoff+5;
endoff=0; dtendoff=endoff;
N = length(xl.data.(sht)(3:end-endoff,1));
face.dt(end+1:end+N,1)   = datenum(xl.textdata.(sht)(3:end-dtendoff,1));
face.stnm(end+1:end+N,1) = xl.textdata.(sht)(3:end-dtendoff,2);
face.lat(end+1:end+N,1)  = xl.data.(sht)(3:end-endoff,1);
face.lon(end+1:end+N,1)  = xl.data.(sht)(3:end-endoff,2);
face.z(end+1:end+N,1)    = xl.data.(sht)(3:end-endoff,3);
face.t(end+1:end+N,1)    = xl.data.(sht)(3:end-endoff,4);
face.s(end+1:end+N,1)    = xl.data.(sht)(3:end-endoff,5);
face.tss(end+1:end+N,1)  = xl.data.(sht)(3:end-endoff,6);
face.chl(end+1:end+N,1)  = xl.data.(sht)(3:end-endoff,7);
face.pha(end+1:end+N,1)  = xl.data.(sht)(3:end-endoff,8);
face.n(end+1:end+N,1)    = xl.data.(sht)(3:end-endoff,9);
face.p(end+1:end+N,1)    = xl.data.(sht)(3:end-endoff,10);
face.si(end+1:end+N,1)   = xl.data.(sht)(3:end-endoff,11);
end;
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
xl = importdata(fullfile(get_relative_path,'Master_Data_Sheet_21aug09.xlsx'));
%xl = importdata(fullfile(get_relative_path,'Master_Data_Sheet_01Feb10.xls'));
face = struct('dt',[],'stnm',[],'lat',[],'lon',[],'z',[],'t',[],'s',[],...
'tss',[],'chl',[],'pha',[],'n',[],'p',[],'si',[]);
for csht={'EddyXmptE1','EddyXmptE2'};
sht = csht{:};
disp(sht);
%endoff=0; dtendoff=endoff+5;
endoff=5; dtendoff=endoff;
N = length(xl.data.(sht)(3:end-endoff,1));
face.dt(end+1:end+N,1)   = datenum(xl.textdata.(sht)(3:end-dtendoff,1));
face.stnm(end+1:end+N,1) = xl.textdata.(sht)(3:end-dtendoff,2);
face.lat(end+1:end+N,1)  = xl.data.(sht)(3:end-endoff,1);
face.lon(end+1:end+N,1)  = xl.data.(sht)(3:end-endoff,2);
face.z(end+1:end+N,1)    = xl.data.(sht)(3:end-endoff,3);
face.t(end+1:end+N,1)    = xl.data.(sht)(3:end-endoff,4);
face.s(end+1:end+N,1)    = xl.data.(sht)(3:end-endoff,5);
face.tss(end+1:end+N,1)  = xl.data.(sht)(3:end-endoff,6);
face.chl(end+1:end+N,1)  = xl.data.(sht)(3:end-endoff,7);
face.pha(end+1:end+N,1)  = xl.data.(sht)(3:end-endoff,8);
face.n(end+1:end+N,1)    = xl.data.(sht)(3:end-endoff,9);
face.p(end+1:end+N,1)    = xl.data.(sht)(3:end-endoff,10);
face.si(end+1:end+N,1)   = xl.data.(sht)(3:end-endoff,11);
face,
end;
face
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
xl = importdata(fullfile(get_relative_path,'Master_Data_Sheet_21aug09.xlsx'));
%xl = importdata(fullfile(get_relative_path,'Master_Data_Sheet_01Feb10.xls'));
face = struct('dt',[],'stnm',{},'lat',[],'lon',[],'z',[],'t',[],'s',[],...
'tss',[],'chl',[],'pha',[],'n',[],'p',[],'si',[]);
for csht={'EddyXmptE1','EddyXmptE2'};
sht = csht{:};
disp(sht);
%endoff=0; dtendoff=endoff+5;
endoff=5; dtendoff=endoff;
N = length(xl.data.(sht)(3:end-endoff,1));
face.dt(end+1:end+N,1)   = datenum(xl.textdata.(sht)(3:end-dtendoff,1));
face.stnm(end+1:end+N,1) = xl.textdata.(sht)(3:end-dtendoff,2);
face.lat(end+1:end+N,1)  = xl.data.(sht)(3:end-endoff,1);
face.lon(end+1:end+N,1)  = xl.data.(sht)(3:end-endoff,2);
face.z(end+1:end+N,1)    = xl.data.(sht)(3:end-endoff,3);
face.t(end+1:end+N,1)    = xl.data.(sht)(3:end-endoff,4);
face.s(end+1:end+N,1)    = xl.data.(sht)(3:end-endoff,5);
face.tss(end+1:end+N,1)  = xl.data.(sht)(3:end-endoff,6);
face.chl(end+1:end+N,1)  = xl.data.(sht)(3:end-endoff,7);
face.pha(end+1:end+N,1)  = xl.data.(sht)(3:end-endoff,8);
face.n(end+1:end+N,1)    = xl.data.(sht)(3:end-endoff,9);
face.p(end+1:end+N,1)    = xl.data.(sht)(3:end-endoff,10);
face.si(end+1:end+N,1)   = xl.data.(sht)(3:end-endoff,11);
face,
end;
face
face(1)
face.dt
N
face.dt
face = struct('dt',[],'stnm',[],'lat',[],'lon',[],'z',[],'t',[],'s',[],...
'tss',[],'chl',[],'pha',[],'n',[],'p',[],'si',[]);
face.stnm = {};
for csht={'EddyXmptE1','EddyXmptE2'};
sht = csht{:};
disp(sht);
%endoff=0; dtendoff=endoff+5;
endoff=5; dtendoff=endoff;
N = length(xl.data.(sht)(3:end-endoff,1));
face.dt(end+1:end+N,1)   = datenum(xl.textdata.(sht)(3:end-dtendoff,1));
face.stnm(end+1:end+N,1) = xl.textdata.(sht)(3:end-dtendoff,2);
face.lat(end+1:end+N,1)  = xl.data.(sht)(3:end-endoff,1);
face.lon(end+1:end+N,1)  = xl.data.(sht)(3:end-endoff,2);
face.z(end+1:end+N,1)    = xl.data.(sht)(3:end-endoff,3);
face.t(end+1:end+N,1)    = xl.data.(sht)(3:end-endoff,4);
face.s(end+1:end+N,1)    = xl.data.(sht)(3:end-endoff,5);
face.tss(end+1:end+N,1)  = xl.data.(sht)(3:end-endoff,6);
face.chl(end+1:end+N,1)  = xl.data.(sht)(3:end-endoff,7);
face.pha(end+1:end+N,1)  = xl.data.(sht)(3:end-endoff,8);
face.n(end+1:end+N,1)    = xl.data.(sht)(3:end-endoff,9);
face.p(end+1:end+N,1)    = xl.data.(sht)(3:end-endoff,10);
face.si(end+1:end+N,1)   = xl.data.(sht)(3:end-endoff,11);
face,
end;
below_25_ix = find(face.t<25);
scatter_fit(face.t,face.n,[sampleCruise,' T'],'N0_x \mumol'); legend('Location','NorthEast');
scatter_fit(face.t,face.p,[sampleCruise,' T'],'P \mumol'); legend('Location','NorthEast');
scatter_fit(face.t,face.si,[sampleCruise,' T'],'Si \mumol'); legend('Location','NorthEast');
scatter_fit(face.t,face.n,[cruiseName,' T'],'N0_x \mumol'); legend('Location','NorthEast');
scatter_fit(face.t,face.p,[cruiseName,' T'],'P \mumol'); legend('Location','NorthEast');
scatter_fit(face.t,face.si,[cruiseName,' T'],'Si \mumol'); legend('Location','NorthEast');
cruiseName = 'NF08';
scatter_fit(face.t,face.n,[cruiseName,' T'],'N0_x \mumol'); legend('Location','NorthEast');
scatter_fit(face.t,face.p,[cruiseName,' T'],'P \mumol'); legend('Location','NorthEast');
scatter_fit(face.t,face.si,[cruiseName,' T'],'Si \mumol'); legend('Location','NorthEast');
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
cruiseName = 'NF08';
if ( strcmpi(cruiseName,'NF09') )
%xl = importdata(fullfile(get_relative_path,'Master_Data_Sheet_21aug09.xlsx'));
xl = importdata(fullfile(get_relative_path,'Master_Data_Sheet_01Feb10.xls'));
sht = 'Sheet3';
endoff=0; dtendoff=endoff+21;
face.dt=datenum(xl.textdata.(sht)(4:end-dtendoff,1))+xl.data.(sht)(1:end-endoff,1);
face.stnm=xl.textdata.(sht)(4:end-dtendoff,3);
face.lat=xl.data.(sht)(1:end-endoff,3)+(xl.data.(sht)(1:end-endoff,4)./60);
face.lon=-( xl.data.(sht)(1:end-endoff,5)+(xl.data.(sht)(1:end-endoff,6)./60) );
face.z=xl.data.(sht)(1:end-endoff,7);
face.t=xl.data.(sht)(1:end-endoff,8);
face.s=xl.data.(sht)(1:end-endoff,9);
face.chl=xl.data.(sht)(1:end-endoff,11);
face.pha=xl.data.(sht)(1:end-endoff,12);
face.n=xl.data.(sht)(1:end-endoff,13);
face.no2=xl.data.(sht)(1:end-endoff,14);
face.no3=face.n-face.no2;
face.nh4=xl.data.(sht)(1:end-endoff,15);
face.p=xl.data.(sht)(1:end-endoff,16);
face.si=xl.data.(sht)(1:end-endoff,17);
elseif ( strcmpi(cruiseName,'NF08') )
xl = importdata(fullfile(get_relative_path,'Master_Data_Sheet_21aug09.xlsx'));
%xl = importdata(fullfile(get_relative_path,'Master_Data_Sheet_01Feb10.xls'));
face = struct('dt',[],'stnm',[],'lat',[],'lon',[],'z',[],'t',[],'s',[],...
'tss',[],'chl',[],'pha',[],'n',[],'p',[],'si',[]);
face.stnm = {};
for csht={'EddyXmptE1','EddyXmptE2'};
sht = csht{:};
disp(sht);
endoff=5;
N = length(xl.data.(sht)(3:end-endoff,1));
face.dt(end+1:end+N,1)   = datenum(xl.textdata.(sht)(3:end-endoff,1));
face.stnm(end+1:end+N,1) = xl.textdata.(sht)(3:end-endoff,2);
face.lat(end+1:end+N,1)  = xl.data.(sht)(3:end-endoff,1);
face.lon(end+1:end+N,1)  = xl.data.(sht)(3:end-endoff,2);
face.z(end+1:end+N,1)    = xl.data.(sht)(3:end-endoff,3);
face.t(end+1:end+N,1)    = xl.data.(sht)(3:end-endoff,4);
face.s(end+1:end+N,1)    = xl.data.(sht)(3:end-endoff,5);
face.tss(end+1:end+N,1)  = xl.data.(sht)(3:end-endoff,6);
face.chl(end+1:end+N,1)  = xl.data.(sht)(3:end-endoff,7);
face.pha(end+1:end+N,1)  = xl.data.(sht)(3:end-endoff,8);
face.n(end+1:end+N,1)    = xl.data.(sht)(3:end-endoff,9);
face.p(end+1:end+N,1)    = xl.data.(sht)(3:end-endoff,10);
face.si(end+1:end+N,1)   = xl.data.(sht)(3:end-endoff,11);
face,
end;
else
error('Which cruise??');
end;
face.den=sw_dens(face.s,face.t,face.z);
xl=[]; clear xl sht endoff dtendoff
%if ( ~exist('h','var') )
if ( ~isfield(face,'h') )
coastpath = get_ecoforecasts_path('coast');
nc = mDataset(fullfile(coastpath,'fl_east_gom_crm_v1.nc'));
if ( ~isempty(nc) )
hlon=cast(nc{'x'}(:),'double');
hlat=cast(nc{'y'}(:),'double');
%xix=find(min(face.lon)-0.1<=hlon & hlon<=max(face.lon)+0.1);
%yix=find(min(face.lat)-0.1<=hlat & hlat<=max(face.lat)+0.1);
%hz = cast(nc{'z'}(yix(1):yix(end),xix(1):xix(end)),'double');
zvar = nc{'z'};
face.h = repmat(nan,size(face.lon));
for dix=1:numel(face.lon)
[lonerr,lonix]=min(abs(hlon-face.lon(dix)));
[laterr,latix]=min(abs(hlat-face.lat(dix)));
if ( lonerr > (190/111e3) || laterr > (190/111e3) )
clear zvar; close(nc); clear nc;
error('Could not locate data point %d in NGDC 3" data: %g,%g',hix,face.lon,face.lat);
end;
face.h(dix) = cast(zvar(latix,lonix),'double');
end;
clear zvar; close(nc);
end;
clear nc;
end;
scatter_fit(face.t,face.n,[cruiseName,' T'],'N0_x \mumol'); legend('Location','NorthEast');
scatter_fit(face.t,face.p,[cruiseName,' T'],'P \mumol'); legend('Location','NorthEast');
scatter_fit(face.t,face.si,[cruiseName,' T'],'Si \mumol'); legend('Location','NorthEast');
face=[]; hlat=[]; hlon=[]; clear face hlat hlon
clear all
analyze_face_nutrients
dbquit
clear all
analyze_face_nutrients
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
dbquit
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
analyze_face_nutrients
dbquit
clear all
analyze_face_nutrients
dbquit
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
analyze_face_nutrients
dbup
dbdown
dbcont
min(face.no3)
max(face.no3)
max(face.no2)
min(face.z)
max(face.z)
cruiseName
face = struct('dt',[],'stnm',[],'lat',[],'lon',[],'z',[],'t',[],'s',[],...
'tss',[],'chl',[],'pha',[],'n',[],'p',[],'si',[]);
face.fname = fullfile(get_relative_path,'Master_Data_Sheet_21aug09.xlsx');
face
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
face = struct('fname',[],'dt',[],'stnm',[],'lat',[],'lon',[],'z',[],'t',[],'s',[],...
'tss',[],'chl',[],'pha',[],'n',[],'p',[],'si',[]);
face.fname = fullfile(get_relative_path,'Master_Data_Sheet_21aug09.xlsx');
face
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
analyze_face_nutrients
dbcont
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
analyze_face_nutrients
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
analyze_face_nutrients
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
analyze_face_nutrients
max(face.p)
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
analyze_face_nutrients
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
analyze_face_nutrients
face
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
dbstop analyze_face_nutrients
analyze_face_nutrients
dbquit
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
analyze_face_nutrients
max(face.p)
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
max(face.si)
analyze_face_nutrients
reviewanim([],0,0,0)
res
which res
help res
help resid
which resid
help ident
help toolbox ident
pd
type extract_sfomc_etc
doPlots=false; extract_sfomc_etc; clear doPlots
help dmin
help dmax
clear dmin dmax
help dmax
help dmin
nansummary(face.si)
face
fmg; hist(face.si,100);
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
analyze_face_nutrients
pd
analyze_face_nutrients
pd
dbstop if error
doPlots=false; extract_sfomc_etc; clear doPlots
reviewanim([],0,0,0)
dbup
res.c_buoy.at_0_6m.sbe_seatemp_diff_1_h_sum
fmg; plot_ts(res.c_buoy.at_0_6m.sbe_seatemp_diff_1_h_sum);
dbdown
size(s),size(t),size(d)
dbquit
min(res.c_buoy.at_15m.sbe_seatemp.data)
min(res.e_buoy.at_30m.sbe_seatemp.data)
min(res.ne_buoy.at_30m.sbe_seatemp.data)
min(res.ne_buoy.at_45m.sbe_seatemp.data)
min(res.nw_w_btm.sbe_seatemp.data)
min(res.nw_w_btm.adcp_seatemp.data)
min(res.se_buoy.at_98m.adcp_seatemp.data)
min(res.se_buoy.at_98m.sbe_seatemp.data)
min(res.sw_buoy.at_15m.sbe_seatemp.data)
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
dir *sefcri*.m
type analyze_sefcri.m
load_all_sefcri
sefcri
sefcri.updb
min(sefcri.updb.t.data)
min(sefcri.dc1.t.data)
min(sefcri.bc1.t.data)
min(sefcri.pb1.t.data)
min(sefcri.mc2.t.data)
min(sefcri.pela.t.data)
min(sefcri.evca.t.data)
min(sefcri.slib.t.data)
min(sefcri.evca.t.data)
fmg; plot_ts(sefcri.evca.t);
sefcri.evca
sefcri.updb
station_dist(sefcri.updb,sefcri.evca)
help station_dist
[d,degT] = station_dist(sefcri.updb,sefcri.evca)
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
dbstop window_func 177
doPlots=false; extract_sfomc_etc; clear doPlots
size(tdts),size(tsrc)
size(dts(1):(1.0/24.0):dts(end))
size([dts(1):(1.0/24.0):dts(end)]')
size(tdts)
size(tdts(:))
dbquit
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
analyze_face_nutrients
pd
analyze_face_nutrients
nansummary(face.si)
doPlots=false; extract_sfomc_etc; clear doPlots
pd
doPlots=false; extract_sfomc_etc; clear doPlots
fmg; plot_ts(res.nw_w_btm.sbe_implied_heatflux);
fmg; plot_ts(res.nw_w_btm.sbe_seatemp,res.nw_w_btm.sbe_salin);
fwyf1 = get_station_from_station_name('fwyf1'); fwyf1 = load_all_ndbc_data(fwyf1);
fmg; plot_ts(res.nw_w_btm.sbe_seatemp,res.nw_w_btm.sbe_salin,fwyf1.ndbc_air_t);
fmg; plot_ts(res.nw_w_btm.sbe_seatemp,res.nw_w_btm.sbe_salin,fwyf1.ndbc_air_t); xlim(res.nw_w_btm.sbe_seatemp.date([1,end])); datetick3;
fmg; plot_ts(res.nw_w_btm.sbe_seatemp,res.nw_w_btm.sbe_salin,fwyf1.ndbc_air_t,fwyf1.ndbc_wind1_speed); xlim(res.nw_w_btm.sbe_seatemp.date([1,end])); datetick3;
fmg; plot_ts(fwyf1.ndbc_air_t,fwyf1.ndbc_wind1_speed,res.nw_w_btm.sbe_seatemp,res.nw_w_btm.sbe_salin); xlim(res.nw_w_btm.sbe_seatemp.date([1,end])); datetick3;
fmg; plot_ts(fwyf1.ndbc_air_t,fwyf1.ndbc_wind1_speed,res.nw_w_btm.sbe_seatemp,res.nw_w_btm.adcp_seatemp,res.nw_w_btm.sbe_salin); xlim(res.nw_w_btm.sbe_seatemp.date([1,end])); datetick3;
fmg; plot_ts(fwyf1.ndbc_air_t,fwyf1.ndbc_wind1_speed,res.nw_w_btm.sbe_seatemp,res.nw_w_btm.adcp_seatemp,res.nw_w_btm.sbe_salin); xlim(datenum(1999,[6,8],1)); datetick3;
fmg; plot_ts(fwyf1.ndbc_air_t,fwyf1.ndbc_wind1_speed,res.nw_w_btm.sbe_seatemp,res.nw_w_btm.adcp_seatemp,res.nw_w_btm.sbe_salin); xlim(datenum(1999,[7,9],1)); datetick3;
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
quit
%-- 3/4/2016 9:16 PM --%
lonf1 = get_station_from_station_name('lonf1'); lonf1 = load_all_ndbc_data(lonf1);
plot(lonf1.ndbc_sea_t.data(get_month(lonf1.ndbc_sea_t.date)==3));
plot(get_dom(lonf1.ndbc_sea_t.date),lonf1.ndbc_sea_t.data(get_month(lonf1.ndbc_sea_t.date)==3),'.');
plot(lonf1.ndbc_sea_t.data(get_month(lonf1.ndbc_sea_t.date)==4));
plot(lonf1.ndbc_sea_t.data(get_month(lonf1.ndbc_sea_t.date
plot(lonf1.ndbc_sea_t.data(get_month(lonf1.ndbc_sea_t.date)))
ix=find(get_month(lonf1.ndbc_sea_t.date),[34]);
mo=get_month(lonf1.ndbc_sea_t.date);
ix=find(mo,[34]);
mo=get_month(lonf1.ndbc_air_t.date);
ix=find(mo,[34]);
fmg; plot(mo(ix),lonf1.ndbc_air_t.data(ix),'.');
ix=find(ismember(get_month(lonf1.ndbc_sea_t.date),[34]));
fmg; plot(get_dom(lonf1.ndbc_air_t.date(ix)),lonf1.ndbc_air_t.data(ix),'.');
ix=find(ismember(get_month(lonf1.ndbc_air_t.date),[34]));
fmg; plot(get_monthday(lonf1.ndbc_air_t.date(ix)),lonf1.ndbc_air_t.data(ix),'.');
ix=find(ismember(get_month(lonf1.ndbc_air_t.date),[34]));
lonf1.ndbc_air_t
unique(get_month(lonf1.ndbc_air_t.date))
help ismember
more on
help ismember
ix=find(ismember(get_month(lonf1.ndbc_air_t.date),[34]));
get_month(lonf1.ndbc_air_t.date)
ismember(get_month(lonf1.ndbc_air_t.date),[34])
ismember(get_month(lonf1.ndbc_air_t.date),[3 4])
ix=find(ismember(get_month(lonf1.ndbc_air_t.date),[3 4]));
fmg; plot(get_monthday(lonf1.ndbc_air_t.date(ix)),lonf1.ndbc_air_t.data(ix),'.');
lookfor celsius
lookfor convert
lookfor unit
help convtemp
convtemp(10,'k','f')
convtemp(10,'K','F')
convtemp(10,'c','f')
convtemp(10:12,'c','f')
fmg; plot(get_monthday(lonf1.ndbc_air_t.date(ix)),convtemp(lonf1.ndbc_air_t.data(ix),'c','f'),'.');
fmg; plot(get_monthday(lonf1.ndbc_air_t.date(ix)),convtemp(lonf1.ndbc_air_t.data(ix),'c','f'),'.'); datetick3('x','mm/dd');
ix=find(ismember(get_month(lonf1.ndbc_air_t.date),[1 2 3 4]));
fmg; plot(get_monthday(lonf1.ndbc_air_t.date(ix)),convtemp(lonf1.ndbc_air_t.data(ix),'c','f'),'.'); datetick3('x','mm/dd');
ix=find(ismember(get_month(lonf1.ndbc_air_t.date),[1 2 3 4])&ismember(get_hour(lonf1.ndbc_air_t.date),[13:21]));
fmg; plot(get_monthday(lonf1.ndbc_air_t.date(ix)),convtemp(lonf1.ndbc_air_t.data(ix),'c','f'),'.'); datetick3('x','mm/dd');
@
fmg; plot(get_hour(lonf1.ndbc_air_t.date(ix)),convtemp(lonf1.ndbc_air_t.data(ix),'c','f'),'.');
ix=find(ismember(get_month(lonf1.ndbc_air_t.date),[1 2 3 4]));
fmg; plot(get_hour(lonf1.ndbc_air_t.date(ix)),convtemp(lonf1.ndbc_air_t.data(ix),'c','f'),'.');
ix=find(ismember(get_month(lonf1.ndbc_air_t.date),[1 2 3 4])&ismember(get_hour(lonf1.ndbc_air_t.date),[17:23]));
fmg; plot(get_monthday(lonf1.ndbc_air_t.date(ix)),convtemp(lonf1.ndbc_air_t.data(ix),'c','f'),'.'); datetick3('x','mm/dd');
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
help pm_value
help pm_cunit
help pm_unit
which pm_unit
help physmod
help units
help pm_units
help pm_cunit
help pm_value
help pm_units
help pm_unit
help pm_commensurate
help pm_fundamentalunit
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
pd
analyze_face_nutrients
pd
load_all_sefcri
doPlots=false; extract_sfomc_etc; clear doPlots
reviewanim([],0,0,0)
nansummary(face.t)
scatter_curve_fit(face.t,face.n,'exp1',[face.cruiseName,' T'],'N0_x \mumol');
axis([tmin,tmax,nmin,nmax]); legend('Location','NorthEast');
scatter_curve_fit(face.t,face.p,'exp1',[face.cruiseName,' T'],'P \mumol');
axis([tmin,tmax,pmin,pmax]); legend('Location','NorthEast');
scatter_curve_fit(face.t,face.si,'exp1',[face.cruiseName,' T'],'Si \mumol');
axis([tmin,tmax,simin,simax]); legend('Location','NorthEast');
exp(0)
scatter_curve_fit(face.t,face.n,'exp2',[face.cruiseName,' T'],'N0_x \mumol');
scatter_curve_fit(face.t,face.n,'exp2',[face.cruiseName,' T'],'N0_x \mumol'); axis([tmin,tmax,nmin,nmax]); legend('Location','NorthEast');
help scatterhist
help scatterhist_ts
scatter_fit(face.t,face.n,'exp2',[face.cruiseName,' T'],'N0_x \mumol'); axis([tmin,tmax,nmin,nmax]); legend('Location','NorthEast');
scatter_fit(face.t,face.n,[face.cruiseName,' T'],'N0_x \mumol'); axis([tmin,tmax,nmin,nmax]); legend('Location','NorthEast');
fmg; scatterhist(face.t,face.n,[face.cruiseName,' T'],'N0_x \mumol'); axis([tmin,tmax,nmin,nmax]); legend('Location','NorthEast');
fmg; scatterhist(face.t,face.n); axis([tmin,tmax,nmin,nmax]); legend('Location','NorthEast');
fmg; scatterhist(face.t,face.n,'nbins',100,'); axis([tmin,tmax,nmin,nmax]); legend('Location','NorthEast');
fmg; scatterhist(face.t,face.n,'nbins',100); axis([tmin,tmax,nmin,nmax]); legend('Location','NorthEast');
numel(ix)
fmg; scatterhist(face.t(non_boil_stnix),face.n(non_boil_stnix),'nbins',100); axis([tmin,tmax,nmin,nmax]); legend('Location','NorthEast');
fmg; scatterhist(face.t(below_25_ix),face.n(below_25_ix),'nbins',100); axis([tmin,tmax,nmin,nmax]); legend('Location','NorthEast');
fmg; scatterhist(face.t(non_boil_stnix),face.n(non_boil_stnix),'nbins',100); axis([tmin,tmax,nmin,nmax]); legend('Location','NorthEast');
fmg; scatterhist(face.t(non_boil_stnix),face.n(non_boil_stnix),'nbins',100);
fmg; scatterhist(face.t(non_boil_stnix),face.n(non_boil_stnix),'nbins',100); axis([tmin,tmax,nmin,nmax]); legend('Location','NorthEast');
fmg; scatterhist(face.t(non_boil_stnix),face.n(non_boil_stnix),'nbins',100);
fmg; scatterhist(face.t(below_25_ix),face.n(below_25_ix),'nbins',100);
below_22_ix = find(face.t<22);
fmg; scatterhist(face.t(below_22_ix),face.n(below_22_ix),'nbins',100);
below_22_ix = find(face.t(non_boil_stnix)<22); below_22_ix = non_boil_stnix(below_22_ix);
fmg; scatterhist(face.t(below_22_ix),face.n(below_22_ix),'nbins',100);
fmg; scatterhist(face.t(below_25_ix),face.n(below_25_ix),'nbins',100);
below_24_ix = find(face.t(non_boil_stnix)<=24); below_24_ix = non_boil_stnix(below_24_ix);
fmg; scatterhist(face.t(below_24_ix),face.n(below_24_ix),'nbins',100);
numel(below_24_ix)
timenow
reviewanim([],0,0,0)
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
analyze_face_nutrients
pd
analyze_face_nutrients
clear all
pack
analyze_face_nutrients
reviewanim([],0,0,0)
datenum(2016,3,1)-datenum(2015,10,1)
152*.9
get_jday(datenum(2016,3,1)),get_jday(datenum(2015,10,1))
51+47
xl
sht
hlat = [];
hlon = [];
clear coastpath dix hlat hlon laterr latix lonerr lonix
clear ans csht dmax dmin doPrint ix N nmax nmin no2max no2min non_boil_stnix;
clear pmax pmin simax simin tmax tmin zmax zmin
face
reviewanim([],0,0,0)
fmg; scatterhist(face.t(below_24_ix),face.p(below_24_ix),'nbins',100); titlename('T<24 vs. P \mumol');
fmg; scatterhist(face.t(below_24_ix),face.si(below_24_ix),'nbins',100); titlename('T<24 vs. Si \mumol');
fmg; scatterhist(face.t(below_24_ix),face.p(below_24_ix),'nbins',100); titlename('T<24 vs. P \mumol');
fmg; scatterhist(face.t(below_24_ix),face.n(below_24_ix),'nbins',100); titlename('T<24 vs. NO_x \mumol');
fmg; scatterhist(face.t(below_24_ix),face.p(below_24_ix),'nbins',100); titlename('T<24 vs. P \mumol');
fmg; scatterhist(face.t(below_24_ix),face.si(below_24_ix),'nbins',100); titlename('T<24 vs. Si \mumol');
nanmax(face.n)
nanmax(face.p)
nanmax(face.so)
nanmax(face.si)
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
analyze_face_nutrients
clear all; pack
analyze_face_nutrients
fwyf1 = get_station_from_station_name('fwyf1'); fwyf1 = load_all_ndbc_data(fwyf1);
fmg; plot_ts(fwyf1.ndbc_air_t); plot(face.dt,face.t); xlim(face.dt([1,end])); datetick3;
fmg; plot_ts(fwyf1.ndbc_air_t); plot(face.dt,face.t,'o'); xlim(face.dt([1,end])); datetick3;
face
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
dir *.png
help dir
ls
help ls
ls -l
ls('-l')
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
pd ../Bleach
x = importdata('Drury_Coordinates.xlsx');
x
x.data
x.textdata
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
x = importdata('Drury_Coordinates.xlsx');
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
x = importdata('Drury_Coordinates_LJG.xlsx');
x
x.data
x.textdata
edit read_drury.m
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
read_drury
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
read_drury
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
read_drury
clear ans cstnms ix stnm x
drury
drury.Stephs_Reef
drury.Bowl
bath.station_name = 'Drury'; bath.lat=25.6; bath.lon=-80.15;
bath = plot_hires_bathymetry(bath,-[0:5:80],[10e3,60e3],[],@contour);
clear ans
for ix=1:numel(cstnms);
plot(drury.(cstnms{ix}).lon,drury.(cstnms{ix}).lon,'sk');
end;
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
read_drury
bath
bath.ngdc_hires_bathy
bath.lat([1,end])
help plot_hires_bathymetry
nansummary(bath.lat)
nansummary(bath.lon)
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
help plot_hires_bathymetry
dbstop plot_hires_bathymetry
bath.station_name = 'Drury'; bath.lat=25.6; bath.lon=-80.15;
bath = plot_hires_bathymetry(bath,-[0:5:80],[10e3,60e3],[],@contour);
stn_or_stnm_or_loc
stn
stn.ngdc_hires_bathy
stn.ngdc_hires_bathy.lon
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
read_drury
drury
bath
bath.ngdc_hires_bathy
nansummary(bath.ngdc_hires_bathy.lon)
nansummary(bath.ngdc_hires_bathy.lat)
axis tight
bath=[]; clear bath
bath.station_name = 'Drury'; bath.lat=25.6; bath.lon=-80.15;
bath = plot_hires_bathymetry(bath,-[0:5:80],[10e3,60e3],[],@contour);
for ix=1:numel(cstnms);
plot(drury.(cstnms{ix}).lon,drury.(cstnms{ix}).lon,'sk');
end;
x = importdata('Drury_Coordinates_LJG.xlsx');
cstnms = x.textdata(2:end,1);
clear ans ix stnm x
bath=[]; clear bath
bath.station_name = 'Drury'; bath.lat=25.6; bath.lon=-80.15;
bath = plot_hires_bathymetry(bath,-[0:5:80],[10e3,60e3],[],@contour);
cstnms = fieldnames(drury);
for ix=1:numel(cstnms);
plot(drury.(cstnms{ix}).lon,drury.(cstnms{ix}).lon,'sk');
end;
bath=[]; clear bath
bath.station_name = 'Drury'; bath.lat=25.6; bath.lon=-80.15;
bath = plot_hires_bathymetry(bath,-[0:5:80],[10e3,60e3],[],@contour);
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
read_drury
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
read_drury
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
read_drury
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
read_drury
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
read_drury
bath
help interp_field
pd ../../MATLAB/ecoforecasts
dir *interp*m
pd
help interp_field
cstnms = fieldnames(drury);
for ix=1:numel(cstnms);
stnm = cstnms{ix};
plot(drury.(stnm).lon,drury.(stnm).lat,'sk');
drury.(stnm).depth = interp_field(bath.ngdc_hires_bathymetry.lat,bath.ngdc_hires_bathymetry.lon,...
bath.ngdc_hires_bathymetry.field,drury.(stnm).lat,drury.(stnm).lon);
end;
bath.ngdc_hires_bathy
cstnms = fieldnames(drury);
for ix=1:numel(cstnms);
stnm = cstnms{ix};
plot(drury.(stnm).lon,drury.(stnm).lat,'sk');
drury.(stnm).depth = interp_field(bath.ngdc_hires_bathy.lat,bath.ngdc_hires_bathy.lon,...
bath.ngdc_hires_bathy.field,drury.(stnm).lat,drury.(stnm).lon);
end;
drury.Bowl
cstnms = fieldnames(drury);
for ix=1:numel(cstnms);
stnm = cstnms{ix};
plot(drury.(stnm).lon,drury.(stnm).lat,'sk');
drury.(stnm).depth = interp_field(bath.ngdc_hires_bathy.lat,bath.ngdc_hires_bathy.lon,...
bath.ngdc_hires_bathy.field,drury.(stnm).lat,drury.(stnm).lon);
disp({stnm,drury.(stnm).depth});
end;
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
pd
pd +3
pd
analyze_face_nutrients
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
analyze_face_nutrients
reviewanim([],0,0,0)
face
face.NF08
face.NF09
scatter_fit([face.NF08.t,face.NF09.t],[face.NF08.n,face.NF09.n]);
scatter_fit([face.NF08.t;face.NF09.t],[face.NF08.n;face.NF09.n]);
scatter_fit([face.NF08.t;face.NF09.t(ix)],[face.NF08.n;face.NF09.n(ix)]);
scatter_fit([face.NF08.t;face.NF09.t(fix)],[face.NF08.n;face.NF09.n(fix)]);
ix = find([face.NF08.t;face.NF09.t]>24)
ix = find([face.NF08.t;face.NF09.t]>24);
scatter_fit([face.NF08.t;face.NF09.t](ix),[face.NF08.n;face.NF09.n](ix));
t=[face.NF08.t;face.NF09.t]; n=[face.NF08.n;face.NF09.n];
scatter_fit(t(ix),n(ix);
scatter_fit(t(ix),n(ix));
ix = find([face.NF08.t;face.NF09.t]<=24);
scatter_fit(t(ix),n(ix));
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
pd
pd ../Bleach
read_drury
drury
drury.Jons_Reef
drury.Jons_Reef = plot_hires_bathymetry(drury.Jons_Reef,-[0:2:80],[10e3,10e3],[],@contour);
drury.Jons_Reef = rmfield(drury.Jons_Reef,'ngdc_hires_bathy');
drury.Jons_Reef = plot_hires_bathymetry(drury.Jons_Reef,-[0:1:80],[5e3,5e3],[],@contour);
drury.Jons_Reef = rmfield(drury.Jons_Reef,'ngdc_hires_bathy');
drury.Jons_Reef = plot_hires_bathymetry(drury.Jons_Reef,-[0:1:40],[4e3,4e3],[],@contour);
drury.Jons_Reef = rmfield(drury.Jons_Reef,'ngdc_hires_bathy');
drury.Bowl = plot_hires_bathymetry(drury.Bowl,-[0:1:40],[4e3,4e3],[],@contour);
print('-dpng','drury.Bowl.bathymetry.png');
drury.Jons_Reef = plot_hires_bathymetry(drury.Jons_Reef,-[0:1:40],[4e3,4e3],[],@contour);
print('-dpng','drury.Jons_Reef.bathymetry.png');
cstnms = fieldnames(drury);
for ix=1:numel(cstnms)
stnm = cstnms{ix};
drury.(stnm) = rmfield(drury.(stnm),'ngdc_hires_bathy');
drury.(stnm) = plot_hires_bathymetry(drury.(stnm),-[0:1:30],[4e3,4e3],[],@contour);
print('-dpng',['drury.',stnm,'.bathymetry.png']);
end;
drury
cstnms = fieldnames(drury);
for ix=1:numel(cstnms)
stnm = cstnms{ix};
if ( isfield(drury.(stnm),'ngdc_hires_bathy') )
drury.(stnm) = rmfield(drury.(stnm),'ngdc_hires_bathy');
end;
drury.(stnm) = plot_hires_bathymetry(drury.(stnm),-[0:1:30],[4e3,4e3],[],@contour);
print('-dpng',['drury.',stnm,'.bathymetry.png']);
end;
clear ans cstnms ix stnm x
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
pd
analyze_face_nutrients
badpix
ylim([0,2])
nansummary(face.NF09.chl)
nansummary(face.NF08.chl)
fmg; hist(face.NF08.chl,100);
fmg; hist(face.NF09.chl,100);
fmg; hist(face.NF09.nh4,100);
fmg; hist(face.NF09.pha,100);
cN='NF08';
scatter_fit(face.(cN).t(below_24_ix),face.(cN).n(below_24_ix),[cN,' T, non-boil T<24'],'N0_x \mumol');
below_24_ix = find(face.(cN).t<=24);
scatter_fit(face.(cN).t(below_24_ix),face.(cN).n(below_24_ix),[cN,' T, non-boil T<24'],'N0_x \mumol');
help scatter_fit
help scatter_fit_ts
scatter_fit(face.(cN).t(below_24_ix),face.(cN).n(below_24_ix),[cN,' T, non-boil T<24'],'N0_x \mumol');
scatter_fit(face.(cN).t,face.(cN).n,[cN,' T, non-boil T<24'],'N0_x \mumol');
badnix = find(face.NF08.t<15 & face.NF08.n<15);
face.z(badnix)
face.NF08.z(badnix)
bath = plot_hires_bathymetry(bath,-[0:5:80],[10e3,60e3],[],@contour);
plot(face.NF08.lon(badnix),face.NF08.lat(badnix),'.');
bath = plot_hires_bathymetry(bath,-[0:5:300],[10e3,60e3],[],@contour);
plot(face.NF08.lon(badnix),face.NF08.lat(badnix),'.');
plot(face.NF08.lon(badnix),face.NF08.lat(badnix),'r.');
badnix
{face.NF08.lon(badnix),face.NF08.lat(badnix),}
{face.NF08.lon(badnix);face.NF08.lat(badnix)}
{face.NF08.lon(badnix)',face.NF08.lat(badnix)',}
face.NF08.lon(badnix),face.NF08.lat(badnix),
badnix
face.NF08.z(55:62)
face.NF08.z(55:59)
face.NF08.n(55:59)
face.NF08.p(55:59)
face.NF08.si(55:59)
face.NF08.chl(55:59)
fmg; contourf(face.NF08.lon,face.NF08.z,face.NF08.n);
help grid2
help grid
help griddata
N = griddata(face.NF08.lon,face.NF08.z,face.NF08.n);
help griddata
unique(face.NF08.z)
size(min(face.NF08.lon):0.05:max(face.NF08.lon))
size(min(face.NF08.lon):0.01:max(face.NF08.lon))
size(min(face.NF08.lat):0.01:max(face.NF08.lat))
[L,Z] = meshgrid(min(face.NF08.lat):0.01:max(face.NF08.lat),[0:2:10,15:5:50,60:10:200]);
Z
[L,Z] = meshgrid(min(face.NF08.lat):0.01:max(face.NF08.lat),[0:2:10,15:5:50,60:10:300]);
N = griddata(face.NF08.lat,face.NF08.z,face.NF08.n,L,Z);
numel(find(isnan(face.NF08.lat)))
numel(find(isnan(face.NF08.z)))
find(isnan(face.NF08.z))
face.NF08.t(65)
face.NF08.h(65)
face.NF08.z(60:65)
face.NF08.z(60:66)
Z
numel(find(isnan(face.NF08.lon)))
numel(find(isnan(face.NF08.lat)))
numel(find(isnan(face.NF08.n)))
numel(find(isnan(face.NF08.p)))
cN
badix = find(isnan(face.(cN).z) | isnan(face.(cN).t));
face.(cN).dt(badix) = [];
face.(cN).stnm(badix) = [];
face.(cN).lat(badix) = [];
face.(cN).lon(badix) = [];
face.(cN).z(badix) = [];
face.(cN).t(badix) = [];
face.(cN).s(badix) = [];
face.(cN).tss(badix) = [];
face.(cN).chl(badix) = [];
face.(cN).pha(badix) = [];
face.(cN).n(badix) = [];
face.(cN).p(badix) = [];
face.(cN).si(badix) = [];
numel(find(isnan(face.NF08.z)))
badix
N = griddata(face.NF08.lat,face.NF08.z,face.NF08.n,L,Z);
fmg; contourf(face.NF08.lon,face.NF08.z,N);
fmg; contourf(L,Z,N);
fmg; contourf(L,-Z,N);
fmg; contourf(L,-Z,N); colorbar;
d
help distance_wgs84
pd ../../MATLAB/ecoforecasts
dir *track*
dir *sect*
help station_field_transect
help transect_wgs84
dir *wgs84*
help distance_wgs84
help reckon_wgs84
[d,degT] = distance_wgs84(face.NF08.lat(1:end-1),face.NF08.lon(1:end-1),face.NF08.lat(2:end),face.NF08.lon(2:end),
[d,degT] = distance_wgs84(face.NF08.lat(1:end-1),face.NF08.lon(1:end-1),face.NF08.lat(2:end),face.NF08.lon(2:end));
face.NF08
[d,degT] = distance_wgs84(face.NF08.lat(1:end-1),face.NF08.lon(1:end-1),face.NF08.lat(2:end),face.NF08.lon(2:end)); d = [0;d(:)]; degT = [0;degT(:)];
[D,Z] = meshgrid(d,[0:2:10,15:5:50,60:10:300]);
N = griddata(d,face.NF08.z,face.NF08.n,D,Z);
min(d)
fmg; hist(d,100);
numel(find(d==0))
fmg; contourf(D,-Z,N); colorbar;
[d,degT] = distance_wgs84(face.NF08.lat(1),face.NF08.lon(1),face.NF08.lat,face.NF08.lon);
fmg; hist(d,100);
help distance_wgs84
numel(find(d==0))
numel(find(d<eps))
numel(find(d<0.001))
numel(find(d<0.01))
numel(find(d<0.1))
numel(find(d<1))
numel(find(d<10))
face.NF08.z(1:10)
face.NF08.z(3:8)
face.NF08.z(4:8)
face.NF08.lon(4:8)
face.NF08.lat(4:8)
fmg; hist(diff(d),100);
nansummary(d)
[D,Z] = meshgrid(d,[0:2:10,15:5:50,60:10:300]);
N = griddata(d,face.NF08.z,face.NF08.n,D,Z);
fmg; contourf(D,-Z,N); colorbar;
shading interp
shading flat
fmg; contourf(D,-Z,N); colorbar;
fmg; contour(D,-Z,N); colorbar;
[D,Z] = meshgrid(d,-[0:2:10,15:5:50,60:10:300]);
fmg; contour(D,-Z,N); colorbar;
fmg; contourf(D,-Z,N); colorbar;
fmg; contourf(D,Z,N); colorbar;
N = griddata(d,face.NF08.z,face.NF08.n,D,Z);
nansummary(face.NF08.z)
[D,Z] = meshgrid(d,[0:2:10,15:5:50,60:10:300]);
N = griddata(d,face.NF08.z,face.NF08.n,D,Z);
1.03^10
1/ans
N
timenow
help griddata
[D,Z] = meshgrid(d,[0:2:10,15:5:50,60:10:300],'natural');
help meshgrid
[D,Z] = meshgrid(d,[0:2:10,15:5:50,60:10:300]);
N = griddata(d,face.NF08.z,face.NF08.n,D,Z,'natural');
fmg; contourf(D,Z,N); colorbar;
fmg; contourf(D,-Z,N); colorbar;
help griddata
fmg; contourf(D,-Z,N); colorbar;
timenow
[D,Z] = meshgrid(d,[0:2:300]);
N = griddata(d,face.NF08.z,face.NF08.n,D,Z,'natural');
fmg; contourf(D,-Z,N); colorbar;
N = griddata(d,face.NF08.z,face.NF08.n,D,Z,'cubic');
fmg; contourf(D,-Z,N); colorbar;
fmg; contourf(D,-Z,N,[0:1:20]); colorbar;
fmg; contourf(D,-Z,N,[0:1:10]); colorbar;
help comment
help punct
help sw_dens
nansummary(face.NF08.h)
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
rehash
pd
pd +2
x = importdata('Drury_Composite Temp Data.xlsx');
x
x.textdata
dts = datenum(x.textdata(2:end,1));
clear dts ans
x.textdata
dts = datenum(x.textdata(2:end,1));
drury.inshore.t.date=dts; drury.inshore.t.data=x.data(:,1);
drury.coopers.t.date=dts; drury.coopers.t.data=x.data(:,2);
drury.cvfd.t.date=dts; drury.cvfd.t.data=x.data(:,3);
drury.stephs.t.date=dts; drury.stephs.t.data=x.data(:,4);
x.textdata
x.textdata(1:2,:)
drury.grounding.t.date=dts; drury.grounding.t.data=x.data(:,5);
drury.miamibeach.t.date=dts; drury.miamibeach.t.data=x.data(:,6);
drury.jons.t.date=dts; drury.jos.t.data=x.data(:,7);
drury = rmfield(drury,'jos');
drury.jons.t.date=dts; drury.jons.t.data=x.data(:,7);
fmg; plot_ts(drury.inshore.t,drury.coopers.t,drury.cvfd.t,drury.stephs.t,drury.grounding.t,drury.miamibeach.t,drury.jons.t); legend('Inshore','Coopers','CVFD','Stephs','Grounding','MiamiBeach','Jons');
drury.inshore
drury.inshore.t
(3*26) +
(3*26) + (25) + (3*13)
x = importdata('Drury_Composite Temp Data_LJG.xlsx');
x.textdata
dts = datenum(x.textdata(2:end,1));
drury.inshore.t.date=dts; drury.inshore.t.data=x.data(:,1);
drury.coopers.t.date=dts; drury.coopers.t.data=x.data(:,2);
drury.cvfd.t.date=dts; drury.cvfd.t.data=x.data(:,3);
drury.stephs.t.date=dts; drury.stephs.t.data=x.data(:,4);
drury.grounding.t.date=dts; drury.grounding.t.data=x.data(:,5);
drury.miamibeach.t.date=dts; drury.miamibeach.t.data=x.data(:,6);
drury.jons.t.date=dts; drury.jons.t.data=x.data(:,7);
fmg; plot_ts(drury.inshore.t,drury.coopers.t,drury.cvfd.t,drury.stephs.t,drury.grounding.t,drury.miamibeach.t,drury.jons.t); legend('Inshore','Coopers','CVFD','Stephs','Grounding','MiamiBeach','Jons');
drury
drury.jons
drury.jons.t
x.textdata
x
x.textdata
x.data
x=[]; clear x
x = importdata('Drury_Composite Temp Data_LJG.xlsx');
x
dts = datenum(x.textdata(2:end,1));
drury
x=[]; clear x
x = importdata('Drury_Coordinates_LJG.xlsx');
cstnms = x.textdata(2:end,1);
cstnms
x=[]; clear x
dts=[]; clear dts
clear ans cstnms
x = importdata('Drury_Composite Temp Data_LJG.xlsx');
dts = datenum(x.textdata(2:end,1));
drury.Inshore.t.date=dts;		drury.Inshore.t.data=x.data(:,1);
drury.Coopers_Reef.t.date=dts;		drury.Coopers_Reef.t.data=x.data(:,2);
drury.CV_FD.t.date=dts;			drury.CV_FD.t.data=x.data(:,3);
drury.Stephs_Reef.t.date=dts;		drury.Stephs_Reef.t.data=x.data(:,4);
drury.Grounding.t.date=dts;		drury.Grounding.t.data=x.data(:,5);
drury.Miami_Beach.t.date=dts;		drury.Miami_Beach.t.data=x.data(:,6);
drury.Jons_Reef.t.date=dts;		drury.Jons_Reef.t.data=x.data(:,7);
dts=[]; x=[]; clear dts x
drury
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
read_drury
reviewanim([],0,0,0)
x = importdata('Drury_Composite Temp Data_LJG.xlsx');
x.textdata
x.data
x
x.textdata
x = importdata('Drury_Composite Temp Data_LJG.xlsx');
dts = datenum(x.textdata(2:end,1));
drury.Inshore.t.date=dts;		drury.Inshore.t.data=x.data(:,1);
drury.Coopers_Reef.t.date=dts;		drury.Coopers_Reef.t.data=x.data(:,2);
drury.CV_FD.t.date=dts;			drury.CV_FD.t.data=x.data(:,3);
drury.Site_454.t.date=dts;		drury.Site_454.t.data=x.data(:,4);
drury.Stephs_Reef.t.date=dts;		drury.Stephs_Reef.t.data=x.data(:,5);
drury.Grounding.t.date=dts;		drury.Grounding.t.data=x.data(:,6);
drury.Miami_Beach.t.date=dts;		drury.Miami_Beach.t.data=x.data(:,7);
drury.Jons_Reef.t.date=dts;		drury.Jons_Reef.t.data=x.data(:,8);
dts=[]; x=[]; clear dts x
if 1;
fmg;
plot_ts(drury.Inshore.t,drury.Coopers_Reef.t,drury.CV_FD.t,drury.Stephs_Reef.t,drury.Grounding.t,drury.Miami_Beach.t,drury.Jons_Reef.t); legend('Inshore','Coopers','CVFD','Stephs','Grounding','MiamiBeach','Jons');
end;
fmg;
plot_ts(drury.Inshore.t,drury.Coopers_Reef.t,drury.CV_FD.t,drury.Site_454.t,drury.Stephs_Reef.t,drury.Grounding.t,drury.Miami_Beach.t,drury.Jons_Reef.t);
legend('Inshore','Coopers','CV_FD','454','Stephs','Grounding','MiamiBeach','Jons');
fmg;
plot_ts(drury.Inshore.t,drury.Coopers_Reef.t,drury.CV_FD.t,drury.Site_454.t,drury.Stephs_Reef.t,drury.Grounding.t,drury.Miami_Beach.t,drury.Jons_Reef.t,'k-');
legend('Inshore','Coopers','CV_FD','454','Stephs','Grounding','MiamiBeach','Jons');
x = importdata('Drury_Composite Temp Data_LJG.xlsx');
dts = datenum(x.textdata(2:end,1));
drury.Inshore.t.date=dts;		drury.Inshore.t.data=x.data(:,1);
drury.Coopers_Reef.t.date=dts;		drury.Coopers_Reef.t.data=x.data(:,2);
drury.CV_FD.t.date=dts;			drury.CV_FD.t.data=x.data(:,3);
drury.Site_454.t.date=dts;		drury.Site_454.t.data=x.data(:,4);
drury.Stephs_Reef.t.date=dts;		drury.Stephs_Reef.t.data=x.data(:,5);
drury.Grounding.t.date=dts;		drury.Grounding.t.data=x.data(:,6);
drury.Miami_Beach.t.date=dts;		drury.Miami_Beach.t.data=x.data(:,7);
drury.Jons_Reef.t.date=dts;		drury.Jons_Reef.t.data=x.data(:,8);
if 1;
fmg;
plot_ts(drury.Inshore.t,drury.Coopers_Reef.t,drury.CV_FD.t,drury.Site_454.t,drury.Stephs_Reef.t,drury.Grounding.t,drury.Miami_Beach.t,drury.Jons_Reef.t,'k-');
xlim(dts([1,end])); datetick3;
legend('Inshore','Coopers','CV_FD','454','Stephs','Grounding','MiamiBeach','Jons');
end;
dts=[]; x=[]; clear dts x
clear ans cstnms ix stnm
print('-dpng','Drury_time_series.png');
drury.Stephs_Reef
drury.Inshore
drury.Inshore_1196
for cstnm = fieldnames(drury); disp({cstnm{:},drury.(cstnm{:}).depth}); end;
for cstnm = fieldnames(drury)'; disp({cstnm{:},drury.(cstnm{:}).depth}); end;
drury
drury.Inshore_1196
x = importdata('Drury_Composite Temp Data_LJG.xlsx');
x.textdata
x.textdata(1:2,:)
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
read_drury
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
read_drury
scatter_fit_ts(drury.Stephs_Reef.t,drury.Inshore_1196.t)
scatter_fit_ts(drury.Stephs_Reef.t,drury.Inshore_1196.t,[],[],'Steph','Inshore')
scatter_fit_ts(drury.Stephs_Reef.t,drury.Inshore_1196.t,[],[],'Steph','Inshore',[],[],true)
scatter_fit_ts_seasons(drury.Stephs_Reef.t,drury.Inshore_1196.t,[],[],'Steph','Inshore',[],[],true)
for cstnm = fieldnames(drury)'; disp({cstnm{:},drury.(cstnm{:}).depth}); end;
scatter_fit_ts_seasons(drury.Stephs_Reef.t,drury.Site_454.t,[],[],'Steph','454',[],[],true)
scatter_fit_ts_seasons(drury.Stephs_Reef.t,drury.Bowl.t,[],[],'Steph','Bowl',[],[],true)
drury.Bowl
x = importdata('Drury_Composite Temp Data_LJG.xlsx');
x.textdata(1:2,:)
drury.Bowl
for cstnm = fieldnames(drury)'; disp({cstnm{:},drury.(cstnm{:}).depth}); drury.(cstnm{:}), end;
for cstnm = fieldnames(drury)'; drury.(cstnm{:}), end;
scatter_fit_ts_seasons(drury.Stephs_Reef.t,drury.Coopers_Reef.t,[],[],'Steph','Coopers',[],[],true)
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
pd
analyze_face_nutrients
re
reviewanim([],0,0,0)
close all
help nearest
lookfor nearest
help dsearchn feasibl knnsearch knnr find_water find_water_xy
help dsearchn
help tsearchn
help feasibl
help knnsearch
help knnr
help find_water
help find_water_xy
which find_water.m
help tmd
help knnsearch
dx = face.NF08.lon;
[DX,Z] = meshgrid(dx,-[0:2:10,15:5:50,60:10:300]);
T = griddata(dx,face.NF08.z,face.NF08.t,DX,Z);
fmg; contourf(DX,Z,T); colorbar; titlename('Along-track T vs. depth');
N = griddata(dx,face.NF08.z,face.NF08.n,DX,Z);
fmg; contourf(DX,Z,N); colorbar; titlename('Along-track NO_x vs. depth');
P = griddata(dx,face.NF08.z,face.NF08.p,DX,Z);
fmg; contourf(DX,Z,P); colorbar; titlename('Along-track P vs. depth');
SI = griddata(dx,face.NF08.z,face.NF08.si,DX,Z);
fmg; contourf(DX,Z,SI); colorbar; titlename('Along-track Si vs. depth');
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
quit
%-- 3/12/2016 5:18 PM --%
analyze_face_nutrients
reviewanim([],0,0,0)
close all
% Longitude - proxy for distance from shore
dx = face.NF08.lon;
% % Along-track distance - proxy for distance from shore
% [dx,degT] = distance_wgs84(face.NF08.lat(1),face.NF08.lon(1),face.NF08.lat,face.NF08.lon);
% % Actually CALCULATE distance from shore
% cst = plot_hires_coastline(bath);
% [idx,dx] = knnsearch(cst,stn);
[DX,Z] = meshgrid(dx,-[0:2:10,15:5:50,60:10:300]);
%{
[DX,Z] = meshgrid(dx,-[0:2:300]);
N = griddata(d,face.NF08.z,face.NF08.n,DX,Z,'natural');
N = griddata(d,face.NF08.z,face.NF08.n,DX,Z,'cubic');
fmg; contourf(DX,Z,N); colorbar; titlename('Along-track NO_x vs. depth');
%}
T = griddata(dx,face.NF08.z,face.NF08.t,DX,Z);
fmg; contourf(DX,Z,T); colorbar; titlename('Along-track T vs. depth');
N = griddata(dx,face.NF08.z,face.NF08.n,DX,Z);
fmg; contourf(DX,Z,N); colorbar; titlename('Along-track NO_x vs. depth');
P = griddata(dx,face.NF08.z,face.NF08.p,DX,Z);
fmg; contourf(DX,Z,P); colorbar; titlename('Along-track P vs. depth');
SI = griddata(dx,face.NF08.z,face.NF08.si,DX,Z);
fmg; contourf(DX,Z,SI); colorbar; titlename('Along-track Si vs. depth');
help plot_hires_coastline
type plot_hires_coastline
more on
type plot_hires_coastline
help contour
[c,h] =
help contour
help contourf
help clabel
help contour
type contour
doc contour
help contourc
bath
bath.ngdc_hires_bathy
CST = contour(bath.ngdc_hires_bathy.lon,bath.ngdc_hires_bathy.lat,bath.ngdc_hires_bathy.field,[0,0]);
CST = contourc(bath.ngdc_hires_bathy.lon,bath.ngdc_hires_bathy.lat,bath.ngdc_hires_bathy.field,[0,0]);
nansummary(CST(1,:))
nansummary(CST(2,:))
fmg hist(CST(1,:),100)
fmg; hist(CST(1,:),100)
fmg; hist(CST(2,:),100)
CST(2,1:10)
CST(2,1:20)
nanmedian(CST(2,:))
median(CST(2,:))
iqr(CST(2,:))
CST(2,1:20)
CST(1,1:20)
help knnsearch
CST(2,1:20)
CST = contourc(bath.ngdc_hires_bathy.lon,bath.ngdc_hires_bathy.lat,bath.ngdc_hires_bathy.field,[0,0]);
CST(1,CST(2,:)==0) = []; CST(2,CST(2,:)==0) = [];
size(CST)
CST(:,CST(2,:)==0) = [];
size(CST)
CST(2,1:20)
CST(1,1:20)
CST(:,CST(1,:)==0) = [];
size(CST)
[idx,dx] = knnsearch(CST,[face.lon,face.lat]);
face
face.NF08
[idx,dx] = knnsearch(CST,[face.NF08.lon,face.NF08.lat]');
[idx,dx] = knnsearch(CST',[face.NF08.lon,face.NF08.lat]);
dx
[idx,dx] = knnsearch(CST',[face.NF08.lon,face.NF08.lat]);
[DX,Z] = meshgrid(dx,-[0:2:10,15:5:50,60:10:300]);
%{
[DX,Z] = meshgrid(dx,-[0:2:300]);
N = griddata(d,face.NF08.z,face.NF08.n,DX,Z,'natural');
N = griddata(d,face.NF08.z,face.NF08.n,DX,Z,'cubic');
fmg; contourf(DX,Z,N); colorbar; titlename('Along-track NO_x vs. depth');
%}
T = griddata(dx,face.NF08.z,face.NF08.t,DX,Z);
fmg; contourf(DX,Z,T); colorbar; titlename('Along-track T vs. depth');
fmg; contourf(DX,Z,T); colorbar; titlename('Along-track T vs. depth');
size(T)
face.NF08.stnm
face.NF09.stnm
face.NF08.stnm
face.NF08.lon
d
dx
diff(dx)
secix = find(diff(dx)<0)
dx(10:16)
diff(dx(10:16))
dx(10:16)
dx(12:14)
secix = [0,find(diff(dx)<0)+1]
secix = [0;find(diff(dx)<0)+1]
secix = [1;find(diff(dx)<0)+1]
diff(secix)
diff(secix)-1
nansummary(diff(secix)-1)
nansummary(dx)
help knnsearch
CST(1,1:20)
dx = distance_wgs84(face.NF08.lat,face.NF08.lon,CST(2,idx),CST(1,idx));
dx = distance_wgs84(face.NF08.lat,face.NF08.lon,CST(2,idx)',CST(1,idx)');
nansummary(dx)
secix = [0;find(diff(dx)<0);length(dx)]+1;
secln = diff(secix)-1;
secln
size(secln)
size(secix)
fmg; contourf(dx(secix),face.NF08.z(secix),face.NF08.t(secix)); colorbar;
clear secix secln
brkix = [0;find(diff(dx)<0);length(dx)]+1;
brkln = diff(brkix)-1;
ixix=1;
secix = brkix(ixix):brkix(ixix)+brkln(ixix);
fmg; contourf(dx(secix),face.NF08.z(secix),face.NF08.t(secix)); colorbar;
for ixix=1;
secix = brkix(ixix):brkix(ixix)+brkln(ixix);
[DX,Z] = meshgrid(dx(secix),-[0:2:10,15:5:50,60:10:300]);
T = griddata(dx(secix),face.NF08.z(secix),face.NF08.t(secix),DX,Z);
fmg; contourf(DX,Z,T); colorbar; titlename('Along-track T vs. depth');
end;
for ixix=2;
secix = brkix(ixix):brkix(ixix)+brkln(ixix);
[DX,Z] = meshgrid(dx(secix),-[0:2:10,15:5:50,60:10:300]);
T = griddata(dx(secix),face.NF08.z(secix),face.NF08.t(secix),DX,Z);
fmg; contourf(DX,Z,T); colorbar; titlename('Along-track T vs. depth');
end;
for ixix=1:length(brkln);
secix = brkix(ixix):brkix(ixix)+brkln(ixix);
[DX,Z] = meshgrid(dx(secix),-[0:2:10,15:5:50,60:10:300]);
T = griddata(dx(secix),face.NF08.z(secix),face.NF08.t(secix),DX,Z);
fmg; contourf(DX,Z,T); colorbar; titlename('Along-track T vs. depth');
titlename(['Section ',num2str(ixix)]);
DX=[]; Z=[]; T=[]; N=[]; P=[]; SI=[]; clear DX Z T N P SI
end;
ixix
nansummary(face.NF08.z)
secix
brkix
for ixix=1:length(brkln);
disp(ixix);
secix = brkix(ixix):brkix(ixix)+brkln(ixix);
%[DX,Z] = meshgrid(dx(secix),-[0:2:10,15:5:50,60:10:300]);
[DX,Z] = meshgrid(dx(secix),-[0:2:300]);
T = griddata(dx(secix),face.NF08.z(secix),face.NF08.t(secix),DX,Z);
fmg; contourf(DX,Z,T); colorbar; titlename('Along-track T vs. depth');
titlename(['Section ',num2str(ixix)]);
DX=[]; Z=[]; T=[]; N=[]; P=[]; SI=[]; clear DX Z T N P SI
end;
% Actually CALCULATE distance from shore
CST = contourc(bath.ngdc_hires_bathy.lon,bath.ngdc_hires_bathy.lat,bath.ngdc_hires_bathy.field,[0,0]);
CST(:,CST(1,:)==0) = [];
%[idx,dx] = knnsearch(CST',[face.NF08.lon,face.NF08.lat]);
idx = knnsearch(CST',[face.NF08.lon,face.NF08.lat]);
dx = distance_wgs84(face.NF08.lat,face.NF08.lon,CST(2,idx)',CST(1,idx)');
% Divide cruise data into cross-shore sections
brkix = [0;find(diff(dx)<0);length(dx)]+1;
brkln = diff(brkix)-1;
for ixix=1:length(brkln);
disp(ixix);
secix = brkix(ixix):brkix(ixix)+brkln(ixix);
%[DX,Z] = meshgrid(dx(secix),-[0:2:10,15:5:50,60:10:300]);
[DX,Z] = meshgrid(dx(secix),-[0:2:300]);
T = griddata(dx(secix),face.NF08.z(secix),face.NF08.t(secix),DX,Z);
fmg; contourf(DX,Z,T); colorbar; titlename('Along-track T vs. depth');
titlename(['Section ',num2str(ixix)]);
DX=[]; Z=[]; T=[]; N=[]; P=[]; SI=[]; clear DX Z T N P SI
end;
close all
% Actually CALCULATE distance from shore
CST = contourc(bath.ngdc_hires_bathy.lon,bath.ngdc_hires_bathy.lat,bath.ngdc_hires_bathy.field,[0,0]);
CST(:,CST(1,:)==0) = [];
%[idx,dx] = knnsearch(CST',[face.NF08.lon,face.NF08.lat]);
idx = knnsearch(CST',[face.NF08.lon,face.NF08.lat]);
dx = distance_wgs84(face.NF08.lat,face.NF08.lon,CST(2,idx)',CST(1,idx)');
% Divide cruise data into cross-shore sections
brkix = [0;find(diff(dx)<0);length(dx)]+1;
brkln = diff(brkix)-1;
for ixix=1:length(brkln);
disp(ixix);
secix = brkix(ixix):brkix(ixix)+brkln(ixix);
%[DX,Z] = meshgrid(dx(secix),-[0:2:10,15:5:50,60:10:300]);
[DX,Z] = meshgrid(dx(secix),-[0:2:300]);
T = griddata(dx(secix),face.NF08.z(secix),face.NF08.t(secix),DX,Z);
fmg; contourf(DX,Z,T); colorbar; titlename('Along-track T vs. depth');
titlename(['Section ',num2str(ixix)]);
DX=[]; Z=[]; T=[]; N=[]; P=[]; SI=[]; clear DX Z T N P SI
end;
for ixix=1:length(brkln);
disp(ixix);
secix = brkix(ixix):brkix(ixix)+brkln(ixix);
%[DX,Z] = meshgrid(dx(secix),-[0:2:10,15:5:50,60:10:300]);
[DX,Z] = meshgrid(dx(secix),-[0:2:300]);
T = griddata(dx(secix),face.NF08.z(secix),face.NF08.t(secix),DX,Z);
fmg; contourf(DX,Z,T); colorbar; titlename('Along-track T vs. depth');
titlename(['Section ',num2str(ixix)]);
DX=[]; Z=[]; T=[]; N=[]; P=[]; SI=[]; clear DX Z T N P SI
end;
reviewanim([],0,0,0)
reviewanim([1:10],0,0,0)
close all
for ixix=1:length(brkln);
disp(ixix);
secix = brkix(ixix):brkix(ixix)+brkln(ixix);
%[DX,Z] = meshgrid(dx(secix),-[0:2:10,15:5:50,60:10:300]);
%[DX,Z] = meshgrid(dx(secix),-[0:2:300]);
[DX,Z] = meshgrid(dx(secix),unique(face.z(secix)));
T = griddata(dx(secix),face.NF08.z(secix),face.NF08.t(secix),DX,Z);
fmg; contourf(DX,Z,T); colorbar; titlename('Along-track T vs. depth');
titlename(['Section ',num2str(ixix)]);
DX=[]; Z=[]; T=[]; N=[]; P=[]; SI=[]; clear DX Z T N P SI
end;
for ixix=1:length(brkln);
disp(ixix);
secix = brkix(ixix):brkix(ixix)+brkln(ixix);
%[DX,Z] = meshgrid(dx(secix),-[0:2:10,15:5:50,60:10:300]);
%[DX,Z] = meshgrid(dx(secix),-[0:2:300]);
[DX,Z] = meshgrid(dx(secix),unique(face.NF08.z(secix)));
T = griddata(dx(secix),face.NF08.z(secix),face.NF08.t(secix),DX,Z);
fmg; contourf(DX,Z,T); colorbar; titlename('Along-track T vs. depth');
titlename(['Section ',num2str(ixix),' T']);
DX=[]; Z=[]; T=[]; N=[]; P=[]; SI=[]; clear DX Z T N P SI
end;
for ixix=1:length(brkln);
disp(ixix);
secix = brkix(ixix):brkix(ixix)+brkln(ixix);
%[DX,Z] = meshgrid(dx(secix),-[0:2:10,15:5:50,60:10:300]);
%[DX,Z] = meshgrid(dx(secix),-[0:2:300]);
[DX,Z] = meshgrid(dx(secix),unique(face.NF08.z(secix)));
% T = griddata(dx(secix),face.NF08.z(secix),face.NF08.t(secix),DX,Z);
% fmg; contourf(DX,Z,T); colorbar; titlename(['Section ',num2str(ixix),' T']);
N = griddata(dx(secix),face.NF08.z(secix),face.NF08.n(secix),DX,Z);
fmg; contourf(DX,Z,N); colorbar; titlename(['Section ',num2str(ixix),' T']);
DX=[]; Z=[]; T=[]; N=[]; P=[]; SI=[]; clear DX Z T N P SI
end;
brklen
brkln
for ixix=1:length(brkln); secix = brkix(ixix):brkix(ixix)+brkln(ixix); disp({ixix,min(face.NF08.z(secix)),max(face.NF08.z(secix)),}); end;
for ixix=1:length(brkln); secix = brkix(ixix):brkix(ixix)+brkln(ixix); disp({ixix,brkix,brkln,min(face.NF08.z(secix)),max(face.NF08.z(secix)),}); end;
for ixix=1:length(brkln); secix = brkix(ixix):brkix(ixix)+brkln(ixix); disp({ixix,brkix(ixix),brkln(ixix),min(face.NF08.z(secix)),max(face.NF08.z(secix)),}); end;
for ixix=1:length(brkln); secix = brkix(ixix):brkix(ixix)+brkln(ixix); disp({ixix,brkix(ixix),brkln(ixix),prctile(face.NF08.z(secix),25),prctile(face.NF08.z(secix),75),}); end;
for ixix=1:length(brkln);
disp(ixix);
secix = brkix(ixix):brkix(ixix)+brkln(ixix);
%[DX,Z] = meshgrid(dx(secix),-[0:2:10,15:5:50,60:10:300]);
%[DX,Z] = meshgrid(dx(secix),-[0:2:300]);
[DX,Z] = meshgrid(dx(secix),unique(face.NF08.z(secix)));
T  = griddata(dx(secix),face.NF08.z(secix),face.NF08.t(secix),DX,Z);
N  = griddata(dx(secix),face.NF08.z(secix),face.NF08.n(secix),DX,Z);
P  = griddata(dx(secix),face.NF08.z(secix),face.NF08.p(secix),DX,Z);
SI = griddata(dx(secix),face.NF08.z(secix),face.NF08.si(secix),DX,Z);
% fmg; contourf(DX,Z,T); colorbar; titlename(['Section ',num2str(ixix),' T']);
fmg; contourf(DX,Z,N); colorbar; titlename(['Section ',num2str(ixix),' N']);
[c,h]=contour(DX,Z,T); clabel(c,h);
DX=[]; Z=[]; T=[]; N=[]; P=[]; SI=[]; clear DX Z T N P SI
end;
for ixix=1:length(brkln); secix = brkix(ixix):brkix(ixix)+brkln(ixix); disp({ixix,brkix(ixix),brkln(ixix),prctile(face.NF08.z(secix),25),prctile(face.NF08.z(secix),75),}); end;
close all
for ixix=1:length(brkln);
disp(ixix);
secix = brkix(ixix):brkix(ixix)+brkln(ixix);
%[DX,Z] = meshgrid(dx(secix),-[0:2:10,15:5:50,60:10:300]);
%[DX,Z] = meshgrid(dx(secix),-[0:2:300]);
[DX,Z] = meshgrid(dx(secix),unique(face.NF08.z(secix)));
D  = griddata(dx(secix),face.NF08.z(secix),face.NF08.dena(secix),DX,Z);
T  = griddata(dx(secix),face.NF08.z(secix),face.NF08.t(secix),DX,Z);
N  = griddata(dx(secix),face.NF08.z(secix),face.NF08.n(secix),DX,Z);
P  = griddata(dx(secix),face.NF08.z(secix),face.NF08.p(secix),DX,Z);
SI = griddata(dx(secix),face.NF08.z(secix),face.NF08.si(secix),DX,Z);
CH = griddata(dx(secix),face.NF08.z(secix),face.NF08.chl(secix),DX,Z);
% fmg; contourf(DX,Z,T); colorbar; titlename(['Section ',num2str(ixix),' T']);
fmg; contourf(DX,Z,N);  colorbar; titlename(['Section ',num2str(ixix),' N']);
[c,h]=contour(DX,Z,D);  clabel(c,h);
%{
fmg; contourf(DX,Z,P);  colorbar; titlename(['Section ',num2str(ixix),' P']);
[c,h]=contour(DX,Z,D);  clabel(c,h);
fmg; contourf(DX,Z,SI); colorbar; titlename(['Section ',num2str(ixix),' Si']);
[c,h]=contour(DX,Z,D);  clabel(c,h);
fmg; contourf(DX,Z,CH); colorbar; titlename(['Section ',num2str(ixix),' Chl_a']);
[c,h]=contour(DX,Z,D);  clabel(c,h);
%}
disp({ixix,brkix(ixix),brkln(ixix),prctile(face.NF08.z(secix),25),prctile(face.NF08.z(secix),75),});
DX=[]; Z=[]; D=[]; T=[]; N=[]; P=[]; SI=[]; CH=[]; clear DX Z D T N P SI CH
end;
close all
for ixix=1:length(brkln);
secix = brkix(ixix):brkix(ixix)+brkln(ixix);
disp({ixix,brkix(ixix),brkln(ixix),min(face.NF08.z(secix)),max(face.NF08.z(secix)),});
%[DX,Z] = meshgrid(dx(secix),-[0:2:10,15:5:50,60:10:300]);
%[DX,Z] = meshgrid(dx(secix),-[0:2:300]);
[DX,Z] = meshgrid(dx(secix),unique(face.NF08.z(secix)));
D  = griddata(dx(secix),face.NF08.z(secix),face.NF08.dena(secix),DX,Z);
T  = griddata(dx(secix),face.NF08.z(secix),face.NF08.t(secix),DX,Z);
N  = griddata(dx(secix),face.NF08.z(secix),face.NF08.n(secix),DX,Z);
P  = griddata(dx(secix),face.NF08.z(secix),face.NF08.p(secix),DX,Z);
SI = griddata(dx(secix),face.NF08.z(secix),face.NF08.si(secix),DX,Z);
CH = griddata(dx(secix),face.NF08.z(secix),face.NF08.chl(secix),DX,Z);
% fmg; contourf(DX,Z,T); colorbar; titlename(['Section ',num2str(ixix),' T']);
fmg; contourf(DX,Z,N);  colorbar; titlename(['Section ',num2str(ixix),' N']);
[c,h]=contour(DX,Z,D);  clabel(c,h);
%{
fmg; contourf(DX,Z,P);  colorbar; titlename(['Section ',num2str(ixix),' P']);
[c,h]=contour(DX,Z,D);  clabel(c,h);
fmg; contourf(DX,Z,SI); colorbar; titlename(['Section ',num2str(ixix),' Si']);
[c,h]=contour(DX,Z,D);  clabel(c,h);
fmg; contourf(DX,Z,CH); colorbar; titlename(['Section ',num2str(ixix),' Chl_a']);
[c,h]=contour(DX,Z,D);  clabel(c,h);
%}
DX=[]; Z=[]; D=[]; T=[]; N=[]; P=[]; SI=[]; CH=[]; clear DX Z D T N P SI CH
end;
help contour
for ixix=1:length(brkln);
secix = brkix(ixix):brkix(ixix)+brkln(ixix);
disp({ixix,brkix(ixix),brkln(ixix),min(face.NF08.z(secix)),max(face.NF08.z(secix)),});
[DX,Z] = meshgrid(dx(secix),unique(face.NF08.z(secix)));
D  = griddata(dx(secix),face.NF08.z(secix),face.NF08.dena(secix),DX,Z);
T  = griddata(dx(secix),face.NF08.z(secix),face.NF08.t(secix),DX,Z);
N  = griddata(dx(secix),face.NF08.z(secix),face.NF08.n(secix),DX,Z);
P  = griddata(dx(secix),face.NF08.z(secix),face.NF08.p(secix),DX,Z);
SI = griddata(dx(secix),face.NF08.z(secix),face.NF08.si(secix),DX,Z);
CH = griddata(dx(secix),face.NF08.z(secix),face.NF08.chl(secix),DX,Z);
%fmg; contourf(DX,Z,T); colorbar; titlename(['Section ',num2str(ixix),' T']);
fmg; contourf(DX,Z,N); colorbar; titlename(['Section ',num2str(ixix),' N']);
[c,h]=contour(DX,Z,D,'k-'); clabel(c,h); [c,h]=contour(DX,Z,T,'k--','Color',[.5,.5,.5]); clabel(c,h);
xlabel('Distance offshore [km]'); ylabel('Depth [m]');
%{
fmg; contourf(DX,Z,P); colorbar; titlename(['Section ',num2str(ixix),' P']);
[c,h]=contour(DX,Z,D); clabel(c,h);
xlabel('Distance offshore [km]'); ylabel('Depth [m]');
fmg; contourf(DX,Z,SI); colorbar; titlename(['Section ',num2str(ixix),' Si']);
[c,h]=contour(DX,Z,D); clabel(c,h);
xlabel('Distance offshore [km]'); ylabel('Depth [m]');
fmg; contourf(DX,Z,CH); colorbar; titlename(['Section ',num2str(ixix),' Chl_a']);
[c,h]=contour(DX,Z,D); clabel(c,h);
xlabel('Distance offshore [km]'); ylabel('Depth [m]');
%}
DX=[]; Z=[]; D=[]; T=[]; N=[]; P=[]; SI=[]; CH=[]; clear DX Z D T N P SI CH
end;
Z
close all
for ixix=1:length(brkln);
DX=[]; Z=[]; D=[]; T=[]; N=[]; P=[]; SI=[]; CH=[]; clear DX Z D T N P SI CH
secix = brkix(ixix):brkix(ixix)+brkln(ixix);
disp({ixix,brkix(ixix),brkln(ixix),min(face.NF08.z(secix)),max(face.NF08.z(secix)),});
[DX,Z] = meshgrid(dx(secix),unique(face.NF08.z(secix)));
D  = griddata(dx(secix),face.NF08.z(secix),face.NF08.dena(secix),DX,Z);
T  = griddata(dx(secix),face.NF08.z(secix),face.NF08.t(secix),DX,Z);
N  = griddata(dx(secix),face.NF08.z(secix),face.NF08.n(secix),DX,Z);
P  = griddata(dx(secix),face.NF08.z(secix),face.NF08.p(secix),DX,Z);
SI = griddata(dx(secix),face.NF08.z(secix),face.NF08.si(secix),DX,Z);
CH = griddata(dx(secix),face.NF08.z(secix),face.NF08.chl(secix),DX,Z);
%fmg; contourf(DX,Z,T); colorbar; titlename(['Section ',num2str(ixix),' T']);
fmg; contourf(DX,Z,N); colorbar; titlename(['Section ',num2str(ixix),' N']);
[c,h]=contour(DX,Z,D,'k-'); clabel(c,h); [c,h]=contour(DX,Z,T,'k--','Color',[.5,.5,.5]); clabel(c,h);
xlabel('Distance offshore [km]'); ylabel('Depth [m]');
%{
fmg; contourf(DX,Z,P); colorbar; titlename(['Section ',num2str(ixix),' P']);
[c,h]=contour(DX,Z,D); clabel(c,h);
xlabel('Distance offshore [km]'); ylabel('Depth [m]');
fmg; contourf(DX,Z,SI); colorbar; titlename(['Section ',num2str(ixix),' Si']);
[c,h]=contour(DX,Z,D); clabel(c,h);
xlabel('Distance offshore [km]'); ylabel('Depth [m]');
fmg; contourf(DX,Z,CH); colorbar; titlename(['Section ',num2str(ixix),' Chl_a']);
[c,h]=contour(DX,Z,D); clabel(c,h);
xlabel('Distance offshore [km]'); ylabel('Depth [m]');
%}
end;
Z
N
close all
for ixix=1:length(brkln);
DX=[]; Z=[]; D=[]; T=[]; N=[]; P=[]; SI=[]; CH=[]; clear DX Z D T N P SI CH
secix = brkix(ixix):brkix(ixix)+brkln(ixix);
disp({ixix,brkix(ixix),brkln(ixix),min(face.NF08.z(secix)),max(face.NF08.z(secix)),});
[DX,Z] = meshgrid(dx(secix),unique(face.NF08.z(secix)));
D  = griddata(dx(secix),face.NF08.z(secix),face.NF08.dena(secix),DX,Z);
T  = griddata(dx(secix),face.NF08.z(secix),face.NF08.t(secix),DX,Z);
N  = griddata(dx(secix),face.NF08.z(secix),face.NF08.n(secix),DX,Z);
P  = griddata(dx(secix),face.NF08.z(secix),face.NF08.p(secix),DX,Z);
SI = griddata(dx(secix),face.NF08.z(secix),face.NF08.si(secix),DX,Z);
CH = griddata(dx(secix),face.NF08.z(secix),face.NF08.chl(secix),DX,Z);
%fmg; contourf(DX,Z,T); colorbar; titlename(['Section ',num2str(ixix),' T']);
fmg; contourf(DX,Z,N); colorbar; titlename(['Section ',num2str(ixix),' N']);
[c,h]=contour(DX,Z,D,'k-'); clabel(c,h); [c,h]=contour(DX,Z,T,'k--','Color',[.5,.5,.5]); clabel(c,h);
xlabel('Distance offshore [km]'); ylabel('Depth [m]');
%{
fmg; contourf(DX,Z,P); colorbar; titlename(['Section ',num2str(ixix),' P']);
[c,h]=contour(DX,Z,D); clabel(c,h);
xlabel('Distance offshore [km]'); ylabel('Depth [m]');
fmg; contourf(DX,Z,SI); colorbar; titlename(['Section ',num2str(ixix),' Si']);
[c,h]=contour(DX,Z,D); clabel(c,h);
xlabel('Distance offshore [km]'); ylabel('Depth [m]');
fmg; contourf(DX,Z,CH); colorbar; titlename(['Section ',num2str(ixix),' Chl_a']);
[c,h]=contour(DX,Z,D); clabel(c,h);
xlabel('Distance offshore [km]'); ylabel('Depth [m]');
%}
end;
for ixix=[1,2,3,4,5,9,10]; disp(ixix); disp('...'); end;
close all
%for ixix=1:length(brkln);
% Only plot cross-shore sections
for ixix=[1,2,3,4,5,9,10];
DX=[]; Z=[]; D=[]; T=[]; N=[]; P=[]; SI=[]; CH=[]; clear DX Z D T N P SI CH
secix = brkix(ixix):brkix(ixix)+brkln(ixix);
disp({ixix,brkix(ixix),brkln(ixix),min(face.NF08.z(secix)),max(face.NF08.z(secix)),});
[DX,Z] = meshgrid(dx(secix),unique(face.NF08.z(secix)));
D  = griddata(dx(secix),face.NF08.z(secix),face.NF08.dena(secix),DX,Z);
T  = griddata(dx(secix),face.NF08.z(secix),face.NF08.t(secix),DX,Z);
N  = griddata(dx(secix),face.NF08.z(secix),face.NF08.n(secix),DX,Z);
P  = griddata(dx(secix),face.NF08.z(secix),face.NF08.p(secix),DX,Z);
SI = griddata(dx(secix),face.NF08.z(secix),face.NF08.si(secix),DX,Z);
CH = griddata(dx(secix),face.NF08.z(secix),face.NF08.chl(secix),DX,Z);
%fmg; contourf(DX,Z,T); colorbar; titlename(['Section ',num2str(ixix),' T']);
fmg; contourf(DX,Z,N); colorbar; titlename(['Section ',num2str(ixix),' N']);
[c,h]=contour(DX,Z,D,'k-'); clabel(c,h); [c,h]=contour(DX,Z,T,'k--','Color',[.5,.5,.5]); clabel(c,h);
xlabel('Distance offshore [km]'); ylabel('Depth [m]');
%{
fmg; contourf(DX,Z,P); colorbar; titlename(['Section ',num2str(ixix),' P']);
[c,h]=contour(DX,Z,D); clabel(c,h);
xlabel('Distance offshore [km]'); ylabel('Depth [m]');
fmg; contourf(DX,Z,SI); colorbar; titlename(['Section ',num2str(ixix),' Si']);
[c,h]=contour(DX,Z,D); clabel(c,h);
xlabel('Distance offshore [km]'); ylabel('Depth [m]');
fmg; contourf(DX,Z,CH); colorbar; titlename(['Section ',num2str(ixix),' Chl_a']);
[c,h]=contour(DX,Z,D); clabel(c,h);
xlabel('Distance offshore [km]'); ylabel('Depth [m]');
%}
end;
close all
xsecix=[1,2,3,4,5,9];
for ixix=xsecix; disp(ixix); disp('...'); end;
timenow
%for ccN = {'NF08','NF09'};
for ccN = {'NF09'};
cN = ccN{:};
% % Longitude - proxy for distance from shore
% face.(cN).dx = face.(cN).lon;
% % Along-track distance - proxy for distance from shore
% [face.(cN).dx,degT] = distance_wgs84(face.(cN).lat(1),face.(cN).lon(1),face.(cN).lat,face.(cN).lon);
%[idx,face.(cN).dx] = knnsearch(CST',[face.(cN).lon,face.(cN).lat]);
idx = knnsearch(CST',[face.(cN).lon,face.(cN).lat]);
face.(cN).dx = distance_wgs84(face.(cN).lat,face.(cN).lon,CST(2,idx)',CST(1,idx)');
% Divide cruise data into cross-shore sections
brkix = [0;find(diff(face.(cN).dx)<0);length(face.(cN).dx)]+1;
brkln = diff(brkix)-1;
xsecix=1:numel(brkln);
% Only plot cross-shore sections
if ( strcmpi(cN,'NF08') )
xsecix=[1,2,3,4,5,9];
elseif ( strcmpi(cN,'NF09') )
end;
for ixix=xsecix;
DX=[]; Z=[]; D=[]; T=[]; N=[]; P=[]; SI=[]; CH=[]; clear DX Z D T N P SI CH
secix = brkix(ixix):brkix(ixix)+brkln(ixix);
disp({ixix,brkix(ixix),brkln(ixix),min(face.(cN).z(secix)),max(face.(cN).z(secix)),});
[DX,Z] = meshgrid(face.(cN).dx(secix),unique(face.(cN).z(secix)));
D  = griddata(face.(cN).dx(secix),face.(cN).z(secix),face.(cN).dena(secix),DX,Z);
T  = griddata(face.(cN).dx(secix),face.(cN).z(secix),face.(cN).t(secix),DX,Z);
N  = griddata(face.(cN).dx(secix),face.(cN).z(secix),face.(cN).n(secix),DX,Z);
P  = griddata(face.(cN).dx(secix),face.(cN).z(secix),face.(cN).p(secix),DX,Z);
SI = griddata(face.(cN).dx(secix),face.(cN).z(secix),face.(cN).si(secix),DX,Z);
CH = griddata(face.(cN).dx(secix),face.(cN).z(secix),face.(cN).chl(secix),DX,Z);
%fmg; contourf(DX,Z,T); colorbar; titlename([cN,' Section ',num2str(ixix),' T']);
fmg; contourf(DX,Z,N); colorbar; titlename([cN,' Section ',num2str(ixix),' N']);
[c,h]=contour(DX,Z,D,'k-'); clabel(c,h); [c,h]=contour(DX,Z,T,'k--','Color',[.5,.5,.5]); clabel(c,h);
xlabel('Distance offshore [km]'); ylabel('Depth [m]');
fmg; contourf(DX,Z,P); colorbar; titlename([cN,' Section ',num2str(ixix),' P']);
[c,h]=contour(DX,Z,D,'k-'); clabel(c,h); [c,h]=contour(DX,Z,T,'k--','Color',[.5,.5,.5]); clabel(c,h);
xlabel('Distance offshore [km]'); ylabel('Depth [m]');
fmg; contourf(DX,Z,SI); colorbar; titlename([cN,' Section ',num2str(ixix),' Si']);
[c,h]=contour(DX,Z,D,'k-'); clabel(c,h); [c,h]=contour(DX,Z,T,'k--','Color',[.5,.5,.5]); clabel(c,h);
xlabel('Distance offshore [km]'); ylabel('Depth [m]');
fmg; contourf(DX,Z,CH); colorbar; titlename([cN,' Section ',num2str(ixix),' Chl_a']);
[c,h]=contour(DX,Z,D,'k-'); clabel(c,h); [c,h]=contour(DX,Z,T,'k--','Color',[.5,.5,.5]); clabel(c,h);
xlabel('Distance offshore [km]'); ylabel('Depth [m]');
end; %for ixix=xsecix;
end; %for ccN = {'NF08','NF09'};
cN
xsecix
ccN = {'NF09'};
cN = ccN{:};
% % Longitude - proxy for distance from shore
% face.(cN).dx = face.(cN).lon;
% % Along-track distance - proxy for distance from shore
% [face.(cN).dx,degT] = distance_wgs84(face.(cN).lat(1),face.(cN).lon(1),face.(cN).lat,face.(cN).lon);
%[idx,face.(cN).dx] = knnsearch(CST',[face.(cN).lon,face.(cN).lat]);
idx = knnsearch(CST',[face.(cN).lon,face.(cN).lat]);
face.(cN).dx = distance_wgs84(face.(cN).lat,face.(cN).lon,CST(2,idx)',CST(1,idx)');
% Divide cruise data into cross-shore sections
brkix = [0;find(diff(face.(cN).dx)<0);length(face.(cN).dx)]+1;
brkln = diff(brkix)-1;
xsecix=1:numel(brkln);
% Only plot cross-shore sections
if ( strcmpi(cN,'NF08') )
xsecix=[1,2,3,4,5,9];
elseif ( strcmpi(cN,'NF09') )
end;
xsecix
dx
xsecix
brkix
brkln
brkln(4:10)
brkln(11:end)
brkln(16:end)
% Only plot cross-shore sections
if ( strcmpi(cN,'NF08') )
xsecix=[1,2,3,4,5,9];
elseif ( strcmpi(cN,'NF09') )
xsecix=[6,7,8,9,11,12,16,18,19,23];
end;
for ixix=xsecix;
DX=[]; Z=[]; D=[]; T=[]; N=[]; P=[]; SI=[]; CH=[]; clear DX Z D T N P SI CH
secix = brkix(ixix):brkix(ixix)+brkln(ixix);
disp({ixix,brkix(ixix),brkln(ixix),min(face.(cN).z(secix)),max(face.(cN).z(secix)),});
[DX,Z] = meshgrid(face.(cN).dx(secix),unique(face.(cN).z(secix)));
D  = griddata(face.(cN).dx(secix),face.(cN).z(secix),face.(cN).dena(secix),DX,Z);
T  = griddata(face.(cN).dx(secix),face.(cN).z(secix),face.(cN).t(secix),DX,Z);
N  = griddata(face.(cN).dx(secix),face.(cN).z(secix),face.(cN).n(secix),DX,Z);
P  = griddata(face.(cN).dx(secix),face.(cN).z(secix),face.(cN).p(secix),DX,Z);
SI = griddata(face.(cN).dx(secix),face.(cN).z(secix),face.(cN).si(secix),DX,Z);
CH = griddata(face.(cN).dx(secix),face.(cN).z(secix),face.(cN).chl(secix),DX,Z);
%fmg; contourf(DX,Z,T); colorbar; titlename([cN,' Section ',num2str(ixix),' T']);
fmg; contourf(DX,Z,N); colorbar; titlename([cN,' Section ',num2str(ixix),' N']);
[c,h]=contour(DX,Z,D,'k-'); clabel(c,h); [c,h]=contour(DX,Z,T,'k--','Color',[.5,.5,.5]); clabel(c,h);
xlabel('Distance offshore [km]'); ylabel('Depth [m]');
fmg; contourf(DX,Z,P); colorbar; titlename([cN,' Section ',num2str(ixix),' P']);
[c,h]=contour(DX,Z,D,'k-'); clabel(c,h); [c,h]=contour(DX,Z,T,'k--','Color',[.5,.5,.5]); clabel(c,h);
xlabel('Distance offshore [km]'); ylabel('Depth [m]');
fmg; contourf(DX,Z,SI); colorbar; titlename([cN,' Section ',num2str(ixix),' Si']);
[c,h]=contour(DX,Z,D,'k-'); clabel(c,h); [c,h]=contour(DX,Z,T,'k--','Color',[.5,.5,.5]); clabel(c,h);
xlabel('Distance offshore [km]'); ylabel('Depth [m]');
fmg; contourf(DX,Z,CH); colorbar; titlename([cN,' Section ',num2str(ixix),' Chl_a']);
[c,h]=contour(DX,Z,D,'k-'); clabel(c,h); [c,h]=contour(DX,Z,T,'k--','Color',[.5,.5,.5]); clabel(c,h);
xlabel('Distance offshore [km]'); ylabel('Depth [m]');
end; %for ixix=xsecix;
for ixix=xsecix;
DX=[]; Z=[]; D=[]; T=[]; N=[]; P=[]; SI=[]; CH=[]; clear DX Z D T N P SI CH
secix = brkix(ixix):brkix(ixix)+brkln(ixix);
disp({ixix,brkix(ixix),brkln(ixix),min(face.(cN).z(secix)),max(face.(cN).z(secix)),});
[DX,Z] = meshgrid(face.(cN).dx(secix),unique(face.(cN).z(secix)));
D  = griddata(face.(cN).dx(secix),face.(cN).z(secix),face.(cN).dena(secix),DX,Z);
T  = griddata(face.(cN).dx(secix),face.(cN).z(secix),face.(cN).t(secix),DX,Z);
N  = griddata(face.(cN).dx(secix),face.(cN).z(secix),face.(cN).n(secix),DX,Z);
P  = griddata(face.(cN).dx(secix),face.(cN).z(secix),face.(cN).p(secix),DX,Z);
SI = griddata(face.(cN).dx(secix),face.(cN).z(secix),face.(cN).si(secix),DX,Z);
CH = griddata(face.(cN).dx(secix),face.(cN).z(secix),face.(cN).chl(secix),DX,Z);
%fmg; contourf(DX,Z,T); colorbar; titlename([cN,' Section ',num2str(ixix),' T']);
fmg; contourf(DX,Z,N); colorbar; titlename([cN,' Section ',num2str(ixix),' N']);
[c,h]=contour(DX,Z,D,'k-'); clabel(c,h); [c,h]=contour(DX,Z,T,'k--','Color',[.5,.5,.5]); clabel(c,h);
xlabel('Distance offshore [km]'); ylabel('Depth [m]');
%{
fmg; contourf(DX,Z,P); colorbar; titlename([cN,' Section ',num2str(ixix),' P']);
[c,h]=contour(DX,Z,D,'k-'); clabel(c,h); [c,h]=contour(DX,Z,T,'k--','Color',[.5,.5,.5]); clabel(c,h);
xlabel('Distance offshore [km]'); ylabel('Depth [m]');
fmg; contourf(DX,Z,SI); colorbar; titlename([cN,' Section ',num2str(ixix),' Si']);
[c,h]=contour(DX,Z,D,'k-'); clabel(c,h); [c,h]=contour(DX,Z,T,'k--','Color',[.5,.5,.5]); clabel(c,h);
xlabel('Distance offshore [km]'); ylabel('Depth [m]');
fmg; contourf(DX,Z,CH); colorbar; titlename([cN,' Section ',num2str(ixix),' Chl_a']);
[c,h]=contour(DX,Z,D,'k-'); clabel(c,h); [c,h]=contour(DX,Z,T,'k--','Color',[.5,.5,.5]); clabel(c,h);
xlabel('Distance offshore [km]'); ylabel('Depth [m]');
%}
end; %for ixix=xsecix;
shading interp
help contourf
help surf
help surfc
help surfl
help patch
fmg; patch(DX,Z,N); colorbar; titlename([cN,' Section ',num2str(ixix),' N']);
[DX,Z] = meshgrid(face.(cN).dx(secix),unique(face.(cN).z(secix)));
D  = griddata(face.(cN).dx(secix),face.(cN).z(secix),face.(cN).dena(secix),DX,Z);
T  = griddata(face.(cN).dx(secix),face.(cN).z(secix),face.(cN).t(secix),DX,Z);
N  = griddata(face.(cN).dx(secix),face.(cN).z(secix),face.(cN).n(secix),DX,Z);
fmg; patch(DX,Z,N); colorbar; titlename([cN,' Section ',num2str(ixix),' N']);
N  = griddata(face.(cN).dx(secix),face.(cN).z(secix),face.(cN).n(secix),DX,Z);
numel(find(isnan(DX)))
numel(find(isnan(Z)))
numel(find(isnan(face.NF09.z)))
size(DX)
numel(find(isnan(face.NF09.t)))
numel(find(isnan(face.NF09.n)))
numel(find(isnan(face.NF09.t) & isnan(face.NF09.z)))
numel(find(isnan(face.NF09.t) & isnan(face.NF09.n)))
numel(find(isnan(face.NF09.t) & isnan(face.NF09.p)))
nansummary(face.NF09.z)
nansummary(face.NF08.z)
ix
nansummary(face.NF09.z(ix))
max(ix)
face.NF09
cN
ix = find(face.(cN).lat>25);
nansummary(face.NF09.z(ix))
nansummary(face.NF09.z()
nansummary(face.NF09.z)
nansummary(face.NF08.z)
bath
ix = find(t<=24);
bath = plot_hires_bathymetry(bath,-[0:5:80,100:20:300,400:100:800],[20e3,60e3],true,@contour);
plot(lon(ix),lat(ix),'.');
bath = plot_hires_bathymetry(bath,-[0:5:80,100:20:300,400:100:800],[20e3,60e3],true,@contour);
ix = find(lon>25);
plot(lon(ix),lat(ix),'.');
plot(lon,lat,'.');
bath = plot_hires_bathymetry(bath,-[0:5:80,100:20:300,400:100:800],[20e3,60e3],true,@contour);
plot(lon(ix),lat(ix),'.');
ix
bath = plot_hires_bathymetry(bath,-[0:5:80,100:20:300,400:100:800],[20e3,60e3],true,@contour);
ix = find(lat>25);
plot(lon(ix),lat(ix),'.');
text(lon(ix),lat(ix),stnm(ix));
stnm={face.NF08.stnm;face.NF09.stnm};
text(lon(ix),lat(ix),stnm(ix));
size(stnm)
stnm={face.NF08.stnm{:};face.NF09.stnm{:}};
stnm={face.NF08.stnm{:},face.NF09.stnm{:}};
size(stnm)
size(lon)
stnm={face.NF08.stnm{:}';face.NF09.stnm{:}'};
stnm={face.NF08.stnm{:},face.NF09.stnm{:}}';
size(stnm)
text(lon(ix),lat(ix),stnm(ix));
brkix
bath = plot_hires_bathymetry(bath,-[0:5:80,100:20:300,400:100:800],[20e3,60e3],true,@contour);
ix = find(lat>25);
below_24_ix = find(t(ix)<=24);
below_24_ix = ix(below_24_ix);
plot(lon(ix),lat(ix),'r.');
plot(lon(below_24_ix),lat(below_24_ix),'b.');
badnix = find(face.NF08.t<15 & face.NF08.n<15);
plot(lon(below_24_ix),lat(below_24_ix),'g.');
text(lon(ix),lat(ix),stnm(ix));
% Plot bathymetric map showing our station locations
bath.station_name = 'NF';
bath.lon = mean(lon(ix)); bath.lat = mean(lat(ix));
bath = plot_hires_bathymetry(bath,-[0:5:80,100:20:300,400:100:800],[20e3,60e3],true,@contour);
ix = find(lat>25);
below_24_ix = find(t(ix)<=24);
below_24_ix = ix(below_24_ix);
plot(lon(ix),lat(ix),'r.');
plot(lon(below_24_ix),lat(below_24_ix),'b.');
badnix = find(t<15 & n<15);
plot(lon(badnix),lat(badnix),'g.');
% Too Much Information
%text(lon(ix),lat(ix),stnm(ix));
bath = plot_hires_bathymetry(bath,-[0:5:80,100:20:300,400:100:800],[20e3,60e3],true,@contour);
% Plot bathymetric map showing our station locations
bath.station_name = 'NF';
bath.lon = mean(lon(ix)); bath.lat = mean(lat(ix));
bath = plot_hires_bathymetry(bath,-[0:5:80,100:20:300,400:100:800],[20e3,60e3],true,@contour);
ix = find(lat>25);
plot(lon(ix),lat(ix),'r.');
below_24_ix = find(t(ix)<=24);
below_24_ix = ix(below_24_ix);
plot(lon(below_24_ix),lat(below_24_ix),'b.');
badnix = find(t(ix)<15 & n(ix)<15);
badnix = ix(badnix);
plot(lon(badnix),lat(badnix),'g.');
bath=[]; clear bath
below_24_ix = find(t<=24);
scatter_fit(t(below_24_ix),n(below_24_ix),'T<24','Both NO_x \mumol');
scatter_fit(t(below_24_ix),p(below_24_ix),'T<24','Both P \mumol');
scatter_fit(t(below_24_ix),si(below_24_ix),'T<24','Both Si \mumol');
close all
ix = find(lat>25);
below_24_ix = find(t(ix)<=24);
below_24_ix = ix(below_24_ix);
weird_n_ix = find(t(ix)<15 & n(ix)<15);
weird_n_ix = ix(weird_n_ix);
scatter_fit(t(below_24_ix),n(below_24_ix),'T<24','Both NO_x \mumol');
scatter_fit(t(below_24_ix),p(below_24_ix),'T<24','Both P \mumol');
scatter_fit(t(below_24_ix),si(below_24_ix),'T<24','Both Si \mumol');
scatter_fit(t(ix),n(ix),'T','Both NO_x \mumol');
scatter_fit(t(ix),p(ix),'T','Both P \mumol');
scatter_fit(t(ix),si(ix),'T','Both Si \mumol');
close all
for below_T=[25,24,23,22,21,20,19,18,17];
below_T_ix = find(t(ix)<=below_T);
below_T_ix = ix(below_T_ix);
scatter_fit(t(below_T_ix),n(below_T_ix),['T<',num2str(below_T)],'Both NO_x \mumol');
end;
close all
pack
for below_T=[27,26,25,24,23,22];
below_T_ix = find(t(ix)<=below_T);
below_T_ix = ix(below_T_ix);
scatter_fit(t(below_T_ix),n(below_T_ix),['T<',num2str(below_T)],'Both NO_x \mumol');
end;
close all
for below_T=[26,25,24,23,22];
below_T_ix = find(t(ix)<=below_T);
below_T_ix = ix(below_T_ix);
scatter_fit(t(below_T_ix),n(below_T_ix),['T<',num2str(below_T)],'Both NO_x \mumol');
scatter_fit(t(below_T_ix),p(below_T_ix),['T<',num2str(below_T)],'Both P \mumol');
scatter_fit(t(below_T_ix),si(below_T_ix),['T<',num2str(below_T)],'Both Si \mumol');
end;
clear below_T below_T_ix
close all
clear ccN cN
clear c h
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
analyze_face_nutrients
cN
face.(cN)
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
analyze_face_nutrients
cN
face.(cN)
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
analyze_face_nutrients
size(DX),size(Z),size(N)
cN
reviewanim([1:10],0,0,0)
reviewanim(],0,0,0)
reviewanim([],0,0,0)
reviewanim([1:37],0,0,0)
brkln
brkln(ixix)
ixix
reviewanim([1:37],0,0,0)
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
analyze_face_nutrients
reviewanim([1:37],0,0,0)
reviewanim([],0,0,0)
get(0,'child')
figure(1)
timenow
reviewanim([1:21],0,0,0)
reviewanim([1:25],0,0,0)
reviewanim([],0,0,0)
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
analyze_face_nutrients
close all
for ccN = {'NF08','NF09'};
cN = ccN{:};
disp(cN);
% Divide cruise data into cross-shore sections
brkix = [0;find(diff(face.(cN).dx)<0);length(face.(cN).dx)]+1;
brkln = diff(brkix)-1;
xsecix=1:numel(brkln);
% Only plot cross-shore sections
if ( strcmpi(cN,'NF08') )
xsecix=[1,2,3,4,5,9];
elseif ( strcmpi(cN,'NF09') )
xsecix=[6,7,8,11,9,12,16,18,19,23];
end;
for ixix=xsecix;
if ( brkln(ixix) < 3 )
disp(['SKIPPING short section ',num2str(ixix)]);
%%%% EARLY CONTINUE
continue;
end;
secix = brkix(ixix):brkix(ixix)+brkln(ixix);
disp({ixix,brkix(ixix),brkln(ixix),min(face.(cN).z(secix)),max(face.(cN).z(secix)),});
[DX,Z] = meshgrid(face.(cN).dx(secix),unique(face.(cN).z(secix)));
D  = griddata(face.(cN).dx(secix),face.(cN).z(secix),face.(cN).dena(secix),DX,Z);
T  = griddata(face.(cN).dx(secix),face.(cN).z(secix),face.(cN).t(secix),DX,Z);
N  = griddata(face.(cN).dx(secix),face.(cN).z(secix),face.(cN).n(secix),DX,Z);
P  = griddata(face.(cN).dx(secix),face.(cN).z(secix),face.(cN).p(secix),DX,Z);
SI = griddata(face.(cN).dx(secix),face.(cN).z(secix),face.(cN).si(secix),DX,Z);
CH = griddata(face.(cN).dx(secix),face.(cN).z(secix),face.(cN).chl(secix),DX,Z);
%fmg; contourf(DX,Z,T); colorbar; titlename([cN,' Section ',num2str(ixix),' T']);
fmg; contourf(DX,Z,N); colorbar; titlename([cN,' Section ',num2str(ixix),' N']);
[c,h]=contour(DX,Z,D,'k-'); clabel(c,h); [c,h]=contour(DX,Z,T,'k--','Color',[.5,.5,.5]); clabel(c,h);
xlabel('Distance offshore [km]'); ylabel('Depth [m]');
%{
fmg; contourf(DX,Z,P); colorbar; titlename([cN,' Section ',num2str(ixix),' P']);
[c,h]=contour(DX,Z,D,'k-'); clabel(c,h); [c,h]=contour(DX,Z,T,'k--','Color',[.5,.5,.5]); clabel(c,h);
xlabel('Distance offshore [km]'); ylabel('Depth [m]');
fmg; contourf(DX,Z,SI); colorbar; titlename([cN,' Section ',num2str(ixix),' Si']);
[c,h]=contour(DX,Z,D,'k-'); clabel(c,h); [c,h]=contour(DX,Z,T,'k--','Color',[.5,.5,.5]); clabel(c,h);
xlabel('Distance offshore [km]'); ylabel('Depth [m]');
fmg; contourf(DX,Z,CH); colorbar; titlename([cN,' Section ',num2str(ixix),' Chl_a']);
[c,h]=contour(DX,Z,D,'k-'); clabel(c,h); [c,h]=contour(DX,Z,T,'k--','Color',[.5,.5,.5]); clabel(c,h);
xlabel('Distance offshore [km]'); ylabel('Depth [m]');
%}
clear c h
DX=[]; Z=[]; D=[]; T=[]; N=[]; P=[]; SI=[]; CH=[]; clear DX Z D T N P SI CH
end; %for ixix=xsecix;
end; %for ccN = {'NF08','NF09'};
close all
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
analyze_face_nutrients
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
analyze_face_nutrients
face.NF08
nansummary(face.NF08.dx)
nansummary(face.NF09.dx)
nansummary(face.NF09.dx(above_25N_ix))
dx=[face.NF08.dx;face.NF09.dx];
nansummary(dx(above_25N_ix))
brkix
face.(cN).lat(brkix(10))
face.NF09.lat(brkix(10))
face.NF09.lat(brkix(5))
face.NF09.lat(brkix(7))
face.NF09.lat(brkix(8))
face.NF09.lat(brkix(9))
unique(face.NF09.lat(brkix(8):brkix(9))
unique(face.NF09.lat(brkix(8):brkix(9)))
help clim
close all
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
numel(find(~isnan([])))
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
analyze_face_nutrients
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
tic,
analyze_face_nutrients
toc,
bath = plot_hires_bathymetry(bath,-[0:5:80,100:20:300,400:100:800],[20e3,60e3],true,@contour);
for ix=1:numel(lon)-1; plot([lon(ix),lon(ix+1)],[lat(ix),lat(ix+1)]); pause(0.5); end;
for ix=above_25N_ix(1):numel(lon); plot([lon(ix),lon(ix+1)],[lat(ix),lat(ix+1)]); pause(0.5); end;
bath = plot_hires_bathymetry(bath,-[0:5:80,100:20:300,400:100:800],[20e3,60e3],true,@contour);
for ix=1:length(lon(above_25N_ix)); plot([lon(above_25N_ix(ix)),lon(above_25N_ix(ix)+1)],[lat(above_25N_ix(ix)),lat(above_25N_ix(ix)+1)]); pause(0.5); end;
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
qyut
quitquit
quit
%-- 3/16/2016 11:16 AM --%
timenow
analyze_face_nutrients
f
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
analyze_face_nutrients
scatter_fit(z,t)
close all
for ccN = {'NF09'};
%for ccN = {'NF08','NF09'};
cN = ccN{:};
disp(cN);
% Divide cruise data into cross-shore sections
brkix = [0;find(diff(face.(cN).dx)<0);length(face.(cN).dx)]+1;
brkln = diff(brkix)-1;
xsecix=1:numel(brkln);
% Only plot cross-shore sections
if ( strcmpi(cN,'NF08') )
xsecix=[1,2,3,4,5,9];
elseif ( strcmpi(cN,'NF09') )
% %xsecix=[6,7,8,11,9,12,16,18,19,23];
% xsecix=[6,7,8,9,12,13,18,19,20,23];
% There's something amiss with "Section #6"
% Sections 1-8 are in the Keys anyway
%    xsecix=[9,12,13,19,20,23];
xsecix=[9:23];
end;
for ixix=xsecix;
if ( brkln(ixix) < 3 )
disp(['SKIPPING short (<3 station) section ',num2str(ixix)]);
%%%% EARLY CONTINUE
continue;
end;
secix = brkix(ixix):brkix(ixix)+brkln(ixix);
disp({ixix,brkix(ixix),brkln(ixix),min(face.(cN).z(secix)),max(face.(cN).z(secix)),});
[DX,Z] = meshgrid(face.(cN).dx(secix),unique(face.(cN).z(secix)));
D  = griddata(face.(cN).dx(secix),face.(cN).z(secix),face.(cN).dena(secix),DX,Z);
T  = griddata(face.(cN).dx(secix),face.(cN).z(secix),face.(cN).t(secix),DX,Z);
N  = griddata(face.(cN).dx(secix),face.(cN).z(secix),face.(cN).n(secix),DX,Z);
P  = griddata(face.(cN).dx(secix),face.(cN).z(secix),face.(cN).p(secix),DX,Z);
SI = griddata(face.(cN).dx(secix),face.(cN).z(secix),face.(cN).si(secix),DX,Z);
CH = griddata(face.(cN).dx(secix),face.(cN).z(secix),face.(cN).chl(secix),DX,Z);
%fmg; contourf(DX,Z,T); colorbar; titlename([cN,' Section ',num2str(ixix),' T']);
if ( numel(find(~isnan(N))) == 0 )
disp(['SKIPPING section with no valid concentrations ',num2str(ixix)]);
else
fmg; contourf(DX,Z,N); colorbar; titlename([cN,' Section ',num2str(ixix),' N']);
[c,h]=contour(DX,Z,D,'k--','Color',[.5,.5,.5]); clabel(c,h); [c,h]=contour(DX,Z,T,'k-'); clabel(c,h);
xlim([0,23]); ylim([-360,0]); xlabel('Distance offshore [km]'); ylabel('Depth [m]');
%{
fmg; contourf(DX,Z,P); colorbar; titlename([cN,' Section ',num2str(ixix),' P']);
[c,h]=contour(DX,Z,D,'k--','Color',[.5,.5,.5]); clabel(c,h); [c,h]=contour(DX,Z,T,'k-'); clabel(c,h);
xlim([0,23]); ylim([-360,0]); xlabel('Distance offshore [km]'); ylabel('Depth [m]');
fmg; contourf(DX,Z,SI); colorbar; titlename([cN,' Section ',num2str(ixix),' Si']);
[c,h]=contour(DX,Z,D,'k--','Color',[.5,.5,.5]); clabel(c,h); [c,h]=contour(DX,Z,T,'k-'); clabel(c,h);
xlim([0,23]); ylim([-360,0]); xlabel('Distance offshore [km]'); ylabel('Depth [m]');
fmg; contourf(DX,Z,CH); colorbar; titlename([cN,' Section ',num2str(ixix),' Chl_a']);
[c,h]=contour(DX,Z,D,'k--','Color',[.5,.5,.5]); clabel(c,h); [c,h]=contour(DX,Z,T,'k-'); clabel(c,h);
xlim([0,23]); ylim([-360,0]); xlabel('Distance offshore [km]'); ylabel('Depth [m]');
%}
clear c h
end;
DX=[]; Z=[]; D=[]; T=[]; N=[]; P=[]; SI=[]; CH=[]; clear DX Z D T N P SI CH
end; %for ixix=xsecix;
end; %for ccN = {'NF08','NF09'};
clear ccN cN
clear ccN cN
close all
% Plot property sections (distance vs. depth vs. concentration) of each cruise
for ccN = {'NF09'};
%for ccN = {'NF08','NF09'};
cN = ccN{:};
disp(cN);
% Divide cruise data into cross-shore sections
brkix = [0;find(diff(face.(cN).dx)<0);length(face.(cN).dx)]+1;
brkln = diff(brkix)-1;
xsecix=1:numel(brkln);
% Only plot cross-shore sections
if ( strcmpi(cN,'NF08') )
xsecix=[1,2,3,4,5,9];
elseif ( strcmpi(cN,'NF09') )
% %xsecix=[6,7,8,11,9,12,16,18,19,23];
% xsecix=[6,7,8,9,12,13,18,19,20,23];
% There's something amiss with "Section #6"
% Sections 1-8 are in the Keys anyway
%    xsecix=[9,12,13,19,20,23];
xsecix=[9:23];
end;
for ixix=xsecix;
% if ( brkln(ixix) < 3 )
%   disp(['SKIPPING short (<3 station) section ',num2str(ixix)]);
%   %%%% EARLY CONTINUE
%   continue;
% end;
secix = brkix(ixix):brkix(ixix)+brkln(ixix);
disp({ixix,brkix(ixix),brkln(ixix),min(face.(cN).z(secix)),max(face.(cN).z(secix)),});
[DX,Z] = meshgrid(face.(cN).dx(secix),unique(face.(cN).z(secix)));
D  = griddata(face.(cN).dx(secix),face.(cN).z(secix),face.(cN).dena(secix),DX,Z);
T  = griddata(face.(cN).dx(secix),face.(cN).z(secix),face.(cN).t(secix),DX,Z);
N  = griddata(face.(cN).dx(secix),face.(cN).z(secix),face.(cN).n(secix),DX,Z);
P  = griddata(face.(cN).dx(secix),face.(cN).z(secix),face.(cN).p(secix),DX,Z);
SI = griddata(face.(cN).dx(secix),face.(cN).z(secix),face.(cN).si(secix),DX,Z);
CH = griddata(face.(cN).dx(secix),face.(cN).z(secix),face.(cN).chl(secix),DX,Z);
%fmg; contourf(DX,Z,T); colorbar; titlename([cN,' Section ',num2str(ixix),' T']);
% if ( numel(find(~isnan(N))) == 0 )
%   disp(['SKIPPING section with no valid concentrations ',num2str(ixix)]);
% else
fmg; contourf(DX,Z,N); colorbar; titlename([cN,' Section ',num2str(ixix),' N']);
[c,h]=contour(DX,Z,D,'k--','Color',[.5,.5,.5]); clabel(c,h); [c,h]=contour(DX,Z,T,'k-'); clabel(c,h);
xlim([0,23]); ylim([-360,0]); xlabel('Distance offshore [km]'); ylabel('Depth [m]');
%{
fmg; contourf(DX,Z,P); colorbar; titlename([cN,' Section ',num2str(ixix),' P']);
[c,h]=contour(DX,Z,D,'k--','Color',[.5,.5,.5]); clabel(c,h); [c,h]=contour(DX,Z,T,'k-'); clabel(c,h);
xlim([0,23]); ylim([-360,0]); xlabel('Distance offshore [km]'); ylabel('Depth [m]');
fmg; contourf(DX,Z,SI); colorbar; titlename([cN,' Section ',num2str(ixix),' Si']);
[c,h]=contour(DX,Z,D,'k--','Color',[.5,.5,.5]); clabel(c,h); [c,h]=contour(DX,Z,T,'k-'); clabel(c,h);
xlim([0,23]); ylim([-360,0]); xlabel('Distance offshore [km]'); ylabel('Depth [m]');
fmg; contourf(DX,Z,CH); colorbar; titlename([cN,' Section ',num2str(ixix),' Chl_a']);
[c,h]=contour(DX,Z,D,'k--','Color',[.5,.5,.5]); clabel(c,h); [c,h]=contour(DX,Z,T,'k-'); clabel(c,h);
xlim([0,23]); ylim([-360,0]); xlabel('Distance offshore [km]'); ylabel('Depth [m]');
%}
clear c h
end;
DX=[]; Z=[]; D=[]; T=[]; N=[]; P=[]; SI=[]; CH=[]; clear DX Z D T N P SI CH
end; %for ixix=xsecix;
end; %for ccN = {'NF08','NF09'};
clear ccN cN
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
quit
%-- 3/16/2016 1:57 PM --%
234/7
`
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
analyze_face_nutrients
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
analyze_face_nutrients
clear ccN cN
% Plot property sections (distance vs. depth vs. concentration) of each cruise
for ccN = {'NF09'};
%for ccN = {'NF08','NF09'};
cN = ccN{:};
disp(cN);
% Divide cruise data into cross-shore sections
brkix = [0;find(diff(face.(cN).dx)<0);length(face.(cN).dx)]+1;
brkln = diff(brkix)-1;
xsecix=1:numel(brkln);
% Only plot cross-shore sections
if ( strcmpi(cN,'NF08') )
xsecix=[1,2,3,4,5,9];
elseif ( strcmpi(cN,'NF09') )
% %xsecix=[6,7,8,11,9,12,16,18,19,23];
% xsecix=[6,7,8,9,12,13,18,19,20,23];
% There's something amiss with "Section #6"
% Sections 1-8 are in the Keys anyway
%    xsecix=[9,12,13,19,20,23];
end;
for ixix=xsecix;
if ( brkln(ixix) < 3 )
disp(['SKIPPING short (<3 station) section ',num2str(ixix)]);
%%%% EARLY CONTINUE
continue;
end;
secix = brkix(ixix):brkix(ixix)+brkln(ixix);
disp({ixix,brkix(ixix),brkln(ixix),min(face.(cN).z(secix)),max(face.(cN).z(secix)),});
[DX,Z] = meshgrid(face.(cN).dx(secix),unique(face.(cN).z(secix)));
D  = griddata(face.(cN).dx(secix),face.(cN).z(secix),face.(cN).dena(secix),DX,Z);
T  = griddata(face.(cN).dx(secix),face.(cN).z(secix),face.(cN).t(secix),DX,Z);
N  = griddata(face.(cN).dx(secix),face.(cN).z(secix),face.(cN).n(secix),DX,Z);
P  = griddata(face.(cN).dx(secix),face.(cN).z(secix),face.(cN).p(secix),DX,Z);
SI = griddata(face.(cN).dx(secix),face.(cN).z(secix),face.(cN).si(secix),DX,Z);
CH = griddata(face.(cN).dx(secix),face.(cN).z(secix),face.(cN).chl(secix),DX,Z);
%fmg; contourf(DX,Z,T); colorbar; titlename([cN,' Section ',num2str(ixix),' T']);
if ( numel(find(~isnan(N))) == 0 )
disp(['SKIPPING section with no valid concentrations ',num2str(ixix)]);
else
fmg; contourf(DX,Z,N); colorbar; titlename([cN,' Section ',num2str(ixix),' N']);
[c,h]=contour(DX,Z,D,'k--','Color',[.5,.5,.5]); clabel(c,h); [c,h]=contour(DX,Z,T,'k-'); clabel(c,h);
xlim([0,23]); ylim([-360,0]); xlabel('Distance offshore [km]'); ylabel('Depth [m]');
%{
fmg; contourf(DX,Z,P); colorbar; titlename([cN,' Section ',num2str(ixix),' P']);
[c,h]=contour(DX,Z,D,'k--','Color',[.5,.5,.5]); clabel(c,h); [c,h]=contour(DX,Z,T,'k-'); clabel(c,h);
xlim([0,23]); ylim([-360,0]); xlabel('Distance offshore [km]'); ylabel('Depth [m]');
fmg; contourf(DX,Z,SI); colorbar; titlename([cN,' Section ',num2str(ixix),' Si']);
[c,h]=contour(DX,Z,D,'k--','Color',[.5,.5,.5]); clabel(c,h); [c,h]=contour(DX,Z,T,'k-'); clabel(c,h);
xlim([0,23]); ylim([-360,0]); xlabel('Distance offshore [km]'); ylabel('Depth [m]');
fmg; contourf(DX,Z,CH); colorbar; titlename([cN,' Section ',num2str(ixix),' Chl_a']);
[c,h]=contour(DX,Z,D,'k--','Color',[.5,.5,.5]); clabel(c,h); [c,h]=contour(DX,Z,T,'k-'); clabel(c,h);
xlim([0,23]); ylim([-360,0]); xlabel('Distance offshore [km]'); ylabel('Depth [m]');
%}
clear c h
end;
DX=[]; Z=[]; D=[]; T=[]; N=[]; P=[]; SI=[]; CH=[]; clear DX Z D T N P SI CH
end; %for ixix=xsecix;
end; %for ccN = {'NF08','NF09'};
clear ccN cN
% Plot property sections (distance vs. depth vs. concentration) of each cruise
for ccN = {'NF09'};
%for ccN = {'NF08','NF09'};
cN = ccN{:};
disp(cN);
% Divide cruise data into cross-shore sections
brkix = [0;find(diff(face.(cN).dx)<0);length(face.(cN).dx)]+1;
brkln = diff(brkix)-1;
xsecix=1:numel(brkln);
% Only plot cross-shore sections
if ( strcmpi(cN,'NF08') )
xsecix=[1,2,3,4,5,9];
elseif ( strcmpi(cN,'NF09') )
xsecix=[9,12,13,19,20];
end;
bath = plot_hires_bathymetry(bath,-[0:5:80,100:20:300,400:100:800],[20e3,60e3],true,@contour);
titlename(['Sections used from ',cN]);
bfig = gcf;
for ixix=xsecix;
if ( brkln(ixix) < 3 )
disp(['SKIPPING short (<3 station) section ',num2str(ixix)]);
%%%% EARLY CONTINUE
continue;
end;
secix = brkix(ixix):brkix(ixix)+brkln(ixix);
disp({ixix,brkix(ixix),brkln(ixix),min(face.(cN).z(secix)),max(face.(cN).z(secix)),});
figure(bfig);
for secixix=1:length(secix)-1
plot(face.(cN).lon(secix(secixix)),face.(cN).lat(secix(secixix)));
text(face.(cN).lon(secix(secixix)),face.(cN).lat(secix(secixix)),num2str(ixix));
end;
[DX,Z] = meshgrid(face.(cN).dx(secix),unique(face.(cN).z(secix)));
D  = griddata(face.(cN).dx(secix),face.(cN).z(secix),face.(cN).dena(secix),DX,Z);
T  = griddata(face.(cN).dx(secix),face.(cN).z(secix),face.(cN).t(secix),DX,Z);
N  = griddata(face.(cN).dx(secix),face.(cN).z(secix),face.(cN).n(secix),DX,Z);
P  = griddata(face.(cN).dx(secix),face.(cN).z(secix),face.(cN).p(secix),DX,Z);
SI = griddata(face.(cN).dx(secix),face.(cN).z(secix),face.(cN).si(secix),DX,Z);
CH = griddata(face.(cN).dx(secix),face.(cN).z(secix),face.(cN).chl(secix),DX,Z);
%fmg; contourf(DX,Z,T); colorbar; titlename([cN,' Section ',num2str(ixix),' T']);
if ( numel(find(~isnan(N))) == 0 )
disp(['SKIPPING section with no valid concentrations ',num2str(ixix)]);
else
fmg; contourf(DX,Z,N); colorbar; titlename([cN,' Section ',num2str(ixix),' N']);
[c,h]=contour(DX,Z,D,'k--','Color',[.5,.5,.5]); clabel(c,h); [c,h]=contour(DX,Z,T,'k-'); clabel(c,h);
xlim([0,23]); ylim([-360,0]); xlabel('Distance offshore [km]'); ylabel('Depth [m]');
%{
fmg; contourf(DX,Z,P); colorbar; titlename([cN,' Section ',num2str(ixix),' P']);
[c,h]=contour(DX,Z,D,'k--','Color',[.5,.5,.5]); clabel(c,h); [c,h]=contour(DX,Z,T,'k-'); clabel(c,h);
xlim([0,23]); ylim([-360,0]); xlabel('Distance offshore [km]'); ylabel('Depth [m]');
fmg; contourf(DX,Z,SI); colorbar; titlename([cN,' Section ',num2str(ixix),' Si']);
[c,h]=contour(DX,Z,D,'k--','Color',[.5,.5,.5]); clabel(c,h); [c,h]=contour(DX,Z,T,'k-'); clabel(c,h);
xlim([0,23]); ylim([-360,0]); xlabel('Distance offshore [km]'); ylabel('Depth [m]');
fmg; contourf(DX,Z,CH); colorbar; titlename([cN,' Section ',num2str(ixix),' Chl_a']);
[c,h]=contour(DX,Z,D,'k--','Color',[.5,.5,.5]); clabel(c,h); [c,h]=contour(DX,Z,T,'k-'); clabel(c,h);
xlim([0,23]); ylim([-360,0]); xlabel('Distance offshore [km]'); ylabel('Depth [m]');
%}
clear c h
end;
DX=[]; Z=[]; D=[]; T=[]; N=[]; P=[]; SI=[]; CH=[]; clear DX Z D T N P SI CH
end; %for ixix=xsecix;
end; %for ccN = {'NF08','NF09'};
clear ccN cN
close all
% Plot property sections (distance vs. depth vs. concentration) of each cruise
for ccN = {'NF09'};
%for ccN = {'NF08','NF09'};
cN = ccN{:};
disp(cN);
% Divide cruise data into cross-shore sections
brkix = [0;find(diff(face.(cN).dx)<0);length(face.(cN).dx)]+1;
brkln = diff(brkix)-1;
xsecix=1:numel(brkln);
% Only plot cross-shore sections
if ( strcmpi(cN,'NF08') )
xsecix=[1,2,3,4,5,9];
elseif ( strcmpi(cN,'NF09') )
xsecix=[9,12,13,19,20];
end;
bath = plot_hires_bathymetry(bath,-[0:5:80,100:20:300,400:100:800],[20e3,60e3],true,@contour);
titlename(['Sections used from ',cN]);
bfig = gcf;
for ixix=xsecix;
if ( brkln(ixix) < 3 )
disp(['SKIPPING short (<3 station) section ',num2str(ixix)]);
%%%% EARLY CONTINUE
continue;
end;
secix = brkix(ixix):brkix(ixix)+brkln(ixix);
disp({ixix,brkix(ixix),brkln(ixix),min(face.(cN).z(secix)),max(face.(cN).z(secix)),});
figure(bfig);
for secixix=1:length(secix)-1
plot(face.(cN).lon(secix(secixix)),face.(cN).lat(secix(secixix)));
text(face.(cN).lon(secix(secixix)),face.(cN).lat(secix(secixix)),num2str(ixix));
end;
[DX,Z] = meshgrid(face.(cN).dx(secix),unique(face.(cN).z(secix)));
D  = griddata(face.(cN).dx(secix),face.(cN).z(secix),face.(cN).dena(secix),DX,Z);
T  = griddata(face.(cN).dx(secix),face.(cN).z(secix),face.(cN).t(secix),DX,Z);
N  = griddata(face.(cN).dx(secix),face.(cN).z(secix),face.(cN).n(secix),DX,Z);
P  = griddata(face.(cN).dx(secix),face.(cN).z(secix),face.(cN).p(secix),DX,Z);
SI = griddata(face.(cN).dx(secix),face.(cN).z(secix),face.(cN).si(secix),DX,Z);
CH = griddata(face.(cN).dx(secix),face.(cN).z(secix),face.(cN).chl(secix),DX,Z);
%fmg; contourf(DX,Z,T); colorbar; titlename([cN,' Section ',num2str(ixix),' T']);
if ( numel(find(~isnan(N))) == 0 )
disp(['SKIPPING section with no valid concentrations ',num2str(ixix)]);
else
fmg; contourf(DX,Z,N); colorbar; titlename([cN,' Section ',num2str(ixix),' N']);
[c,h]=contour(DX,Z,D,'k--','Color',[.5,.5,.5]); clabel(c,h); [c,h]=contour(DX,Z,T,'k-'); clabel(c,h);
xlim([0,23]); ylim([-360,0]); xlabel('Distance offshore [km]'); ylabel('Depth [m]');
%{
fmg; contourf(DX,Z,P); colorbar; titlename([cN,' Section ',num2str(ixix),' P']);
[c,h]=contour(DX,Z,D,'k--','Color',[.5,.5,.5]); clabel(c,h); [c,h]=contour(DX,Z,T,'k-'); clabel(c,h);
xlim([0,23]); ylim([-360,0]); xlabel('Distance offshore [km]'); ylabel('Depth [m]');
fmg; contourf(DX,Z,SI); colorbar; titlename([cN,' Section ',num2str(ixix),' Si']);
[c,h]=contour(DX,Z,D,'k--','Color',[.5,.5,.5]); clabel(c,h); [c,h]=contour(DX,Z,T,'k-'); clabel(c,h);
xlim([0,23]); ylim([-360,0]); xlabel('Distance offshore [km]'); ylabel('Depth [m]');
fmg; contourf(DX,Z,CH); colorbar; titlename([cN,' Section ',num2str(ixix),' Chl_a']);
[c,h]=contour(DX,Z,D,'k--','Color',[.5,.5,.5]); clabel(c,h); [c,h]=contour(DX,Z,T,'k-'); clabel(c,h);
xlim([0,23]); ylim([-360,0]); xlabel('Distance offshore [km]'); ylabel('Depth [m]');
%}
clear c h
end;
DX=[]; Z=[]; D=[]; T=[]; N=[]; P=[]; SI=[]; CH=[]; clear DX Z D T N P SI CH
end; %for ixix=xsecix;
end; %for ccN = {'NF08','NF09'};
clear ccN cN
clear ccN cN
close all
% Plot property sections (distance vs. depth vs. concentration) of each cruise
for ccN = {'NF09'};
%for ccN = {'NF08','NF09'};
cN = ccN{:};
disp(cN);
% Divide cruise data into cross-shore sections
brkix = [0;find(diff(face.(cN).dx)<0);length(face.(cN).dx)]+1;
brkln = diff(brkix)-1;
xsecix=1:numel(brkln);
% Only plot cross-shore sections
if ( strcmpi(cN,'NF08') )
xsecix=[1,2,3,4,5,9];
elseif ( strcmpi(cN,'NF09') )
xsecix=[9,12,13,19,20];
end;
bath = plot_hires_bathymetry(bath,-[0:5:80,100:20:300,400:100:800],[20e3,60e3],true,@contour);
titlename(['Sections used from ',cN]);
bfig = gcf;
for ixix=xsecix;
if ( brkln(ixix) < 3 )
disp(['SKIPPING short (<3 station) section ',num2str(ixix)]);
%%%% EARLY CONTINUE
continue;
end;
secix = brkix(ixix):brkix(ixix)+brkln(ixix);
disp({ixix,brkix(ixix),brkln(ixix),min(face.(cN).z(secix)),max(face.(cN).z(secix)),});
figure(bfig);
for secixix=1:length(secix)-1
plot(face.(cN).lon(secix(secixix:secixix+1)),face.(cN).lat(secix(secixix:secixix+1)));
text(mean(face.(cN).lon(secix(secixix:secixix+1))),mean(face.(cN).lat(secix(secixix:secixix+1))),num2str(ixix));
end;
[DX,Z] = meshgrid(face.(cN).dx(secix),unique(face.(cN).z(secix)));
D  = griddata(face.(cN).dx(secix),face.(cN).z(secix),face.(cN).dena(secix),DX,Z);
T  = griddata(face.(cN).dx(secix),face.(cN).z(secix),face.(cN).t(secix),DX,Z);
N  = griddata(face.(cN).dx(secix),face.(cN).z(secix),face.(cN).n(secix),DX,Z);
P  = griddata(face.(cN).dx(secix),face.(cN).z(secix),face.(cN).p(secix),DX,Z);
SI = griddata(face.(cN).dx(secix),face.(cN).z(secix),face.(cN).si(secix),DX,Z);
CH = griddata(face.(cN).dx(secix),face.(cN).z(secix),face.(cN).chl(secix),DX,Z);
%fmg; contourf(DX,Z,T); colorbar; titlename([cN,' Section ',num2str(ixix),' T']);
if ( numel(find(~isnan(N))) == 0 )
disp(['SKIPPING section with no valid concentrations ',num2str(ixix)]);
else
fmg; contourf(DX,Z,N); colorbar; titlename([cN,' Section ',num2str(ixix),' N']);
[c,h]=contour(DX,Z,D,'k--','Color',[.5,.5,.5]); clabel(c,h); [c,h]=contour(DX,Z,T,'k-'); clabel(c,h);
xlim([0,23]); ylim([-360,0]); xlabel('Distance offshore [km]'); ylabel('Depth [m]');
%{
fmg; contourf(DX,Z,P); colorbar; titlename([cN,' Section ',num2str(ixix),' P']);
[c,h]=contour(DX,Z,D,'k--','Color',[.5,.5,.5]); clabel(c,h); [c,h]=contour(DX,Z,T,'k-'); clabel(c,h);
xlim([0,23]); ylim([-360,0]); xlabel('Distance offshore [km]'); ylabel('Depth [m]');
fmg; contourf(DX,Z,SI); colorbar; titlename([cN,' Section ',num2str(ixix),' Si']);
[c,h]=contour(DX,Z,D,'k--','Color',[.5,.5,.5]); clabel(c,h); [c,h]=contour(DX,Z,T,'k-'); clabel(c,h);
xlim([0,23]); ylim([-360,0]); xlabel('Distance offshore [km]'); ylabel('Depth [m]');
fmg; contourf(DX,Z,CH); colorbar; titlename([cN,' Section ',num2str(ixix),' Chl_a']);
[c,h]=contour(DX,Z,D,'k--','Color',[.5,.5,.5]); clabel(c,h); [c,h]=contour(DX,Z,T,'k-'); clabel(c,h);
xlim([0,23]); ylim([-360,0]); xlabel('Distance offshore [km]'); ylabel('Depth [m]');
%}
clear c h
end;
DX=[]; Z=[]; D=[]; T=[]; N=[]; P=[]; SI=[]; CH=[]; clear DX Z D T N P SI CH
end; %for ixix=xsecix;
end; %for ccN = {'NF08','NF09'};
clear ccN cN
timenow
,
timenow
,
;
close all
% Plot property sections (distance vs. depth vs. concentration) of each cruise
for ccN = {'NF09'};
%for ccN = {'NF08','NF09'};
cN = ccN{:};
disp(cN);
% Divide cruise data into cross-shore sections
brkix = [0;find(diff(face.(cN).dx)<0);length(face.(cN).dx)]+1;
brkln = diff(brkix)-1;
xsecix=1:numel(brkln);
% Only plot cross-shore sections
if ( strcmpi(cN,'NF08') )
xsecix=[1,2,3,4,5,9];
elseif ( strcmpi(cN,'NF09') )
xsecix=[9,12,13,19,20];
end;
bath = plot_hires_bathymetry(bath,-[0:5:80,100:20:300,400:100:800],[20e3,60e3],true,@contour);
titlename(['Sections used from ',cN]);
bfig = gcf;
for ixix=xsecix;
if ( brkln(ixix) < 3 )
disp(['SKIPPING short (<3 station) section ',num2str(ixix)]);
%%%% EARLY CONTINUE
continue;
end;
secix = brkix(ixix):brkix(ixix)+brkln(ixix);
disp({ixix,brkix(ixix),brkln(ixix),min(face.(cN).z(secix)),max(face.(cN).z(secix)),});
figure(bfig);
for secixix=1:length(secix)-1
plot(face.(cN).lon(secix(secixix:secixix+1)),face.(cN).lat(secix(secixix:secixix+1)));
text(mean(face.(cN).lon(secix(secixix:secixix+1))),mean(face.(cN).lat(secix(secixix:secixix+1))),num2str(ixix));
end;
[DX,Z] = meshgrid(face.(cN).dx(secix),unique(face.(cN).z(secix)));
D  = griddata(face.(cN).dx(secix),face.(cN).z(secix),face.(cN).dena(secix),DX,Z);
T  = griddata(face.(cN).dx(secix),face.(cN).z(secix),face.(cN).t(secix),DX,Z);
N  = griddata(face.(cN).dx(secix),face.(cN).z(secix),face.(cN).n(secix),DX,Z);
P  = griddata(face.(cN).dx(secix),face.(cN).z(secix),face.(cN).p(secix),DX,Z);
SI = griddata(face.(cN).dx(secix),face.(cN).z(secix),face.(cN).si(secix),DX,Z);
CH = griddata(face.(cN).dx(secix),face.(cN).z(secix),face.(cN).chl(secix),DX,Z);
%fmg; contourf(DX,Z,T); colorbar; titlename([cN,' Section ',num2str(ixix),' T']);
if ( numel(find(~isnan(N))) == 0 )
disp(['SKIPPING section with no valid concentrations ',num2str(ixix)]);
else
fmg; contourf(DX,Z,N); colorbar; titlename([cN,' Section ',num2str(ixix),' N']);
[c,h]=contour(DX,Z,D,'k--','Color',[.5,.5,.5]); clabel(c,h); [c,h]=contour(DX,Z,T,'k-'); clabel(c,h);
xlim([0,23]); ylim([-360,0]); xlabel('Distance offshore [km]'); ylabel('Depth [m]');
for secixix=1:length(secix)-1
plot(face.(cN).lon(secix(secixix:secixix+1)),face.(cN).z(secix(secixix:secixix+1)));
end;
%{
fmg; contourf(DX,Z,P); colorbar; titlename([cN,' Section ',num2str(ixix),' P']);
[c,h]=contour(DX,Z,D,'k--','Color',[.5,.5,.5]); clabel(c,h); [c,h]=contour(DX,Z,T,'k-'); clabel(c,h);
xlim([0,23]); ylim([-360,0]); xlabel('Distance offshore [km]'); ylabel('Depth [m]');
fmg; contourf(DX,Z,SI); colorbar; titlename([cN,' Section ',num2str(ixix),' Si']);
[c,h]=contour(DX,Z,D,'k--','Color',[.5,.5,.5]); clabel(c,h); [c,h]=contour(DX,Z,T,'k-'); clabel(c,h);
xlim([0,23]); ylim([-360,0]); xlabel('Distance offshore [km]'); ylabel('Depth [m]');
fmg; contourf(DX,Z,CH); colorbar; titlename([cN,' Section ',num2str(ixix),' Chl_a']);
[c,h]=contour(DX,Z,D,'k--','Color',[.5,.5,.5]); clabel(c,h); [c,h]=contour(DX,Z,T,'k-'); clabel(c,h);
xlim([0,23]); ylim([-360,0]); xlabel('Distance offshore [km]'); ylabel('Depth [m]');
%}
clear c h
end;
DX=[]; Z=[]; D=[]; T=[]; N=[]; P=[]; SI=[]; CH=[]; clear DX Z D T N P SI CH
end; %for ixix=xsecix;
end; %for ccN = {'NF08','NF09'};
clear ccN cN
close all
% Plot property sections (distance vs. depth vs. concentration) of each cruise
for ccN = {'NF09'};
%for ccN = {'NF08','NF09'};
cN = ccN{:};
disp(cN);
% Divide cruise data into cross-shore sections
brkix = [0;find(diff(face.(cN).dx)<0);length(face.(cN).dx)]+1;
brkln = diff(brkix)-1;
xsecix=1:numel(brkln);
% Only plot cross-shore sections
if ( strcmpi(cN,'NF08') )
xsecix=[1,2,3,4,5,9];
elseif ( strcmpi(cN,'NF09') )
xsecix=[9,12,13,19,20];
end;
bath = plot_hires_bathymetry(bath,-[0:5:80,100:20:300,400:100:800],[20e3,60e3],true,@contour);
titlename(['Sections used from ',cN]);
bfig = gcf;
for ixix=xsecix;
if ( brkln(ixix) < 3 )
disp(['SKIPPING short (<3 station) section ',num2str(ixix)]);
%%%% EARLY CONTINUE
continue;
end;
secix = brkix(ixix):brkix(ixix)+brkln(ixix);
disp({ixix,brkix(ixix),brkln(ixix),min(face.(cN).z(secix)),max(face.(cN).z(secix)),});
figure(bfig);
for secixix=1:length(secix)-1
plot(face.(cN).lon(secix(secixix:secixix+1)),face.(cN).lat(secix(secixix:secixix+1)));
text(mean(face.(cN).lon(secix(secixix:secixix+1))),mean(face.(cN).lat(secix(secixix:secixix+1))),num2str(ixix));
end;
[DX,Z] = meshgrid(face.(cN).dx(secix),unique(face.(cN).z(secix)));
D  = griddata(face.(cN).dx(secix),face.(cN).z(secix),face.(cN).dena(secix),DX,Z);
T  = griddata(face.(cN).dx(secix),face.(cN).z(secix),face.(cN).t(secix),DX,Z);
N  = griddata(face.(cN).dx(secix),face.(cN).z(secix),face.(cN).n(secix),DX,Z);
P  = griddata(face.(cN).dx(secix),face.(cN).z(secix),face.(cN).p(secix),DX,Z);
SI = griddata(face.(cN).dx(secix),face.(cN).z(secix),face.(cN).si(secix),DX,Z);
CH = griddata(face.(cN).dx(secix),face.(cN).z(secix),face.(cN).chl(secix),DX,Z);
%fmg; contourf(DX,Z,T); colorbar; titlename([cN,' Section ',num2str(ixix),' T']);
if ( numel(find(~isnan(N))) == 0 )
disp(['SKIPPING section with no valid concentrations ',num2str(ixix)]);
else
fmg; contourf(DX,Z,N); colorbar; titlename([cN,' Section ',num2str(ixix),' N']);
[c,h]=contour(DX,Z,D,'k--','Color',[.5,.5,.5]); clabel(c,h); [c,h]=contour(DX,Z,T,'k-'); clabel(c,h);
xlim([0,23]); ylim([-360,0]); xlabel('Distance offshore [km]'); ylabel('Depth [m]');
for secixix=1:length(secix)-1
plot(face.(cN).lon(secix(secixix:secixix+1)),face.(cN).h(secix(secixix:secixix+1)));
end;
%{
fmg; contourf(DX,Z,P); colorbar; titlename([cN,' Section ',num2str(ixix),' P']);
[c,h]=contour(DX,Z,D,'k--','Color',[.5,.5,.5]); clabel(c,h); [c,h]=contour(DX,Z,T,'k-'); clabel(c,h);
xlim([0,23]); ylim([-360,0]); xlabel('Distance offshore [km]'); ylabel('Depth [m]');
fmg; contourf(DX,Z,SI); colorbar; titlename([cN,' Section ',num2str(ixix),' Si']);
[c,h]=contour(DX,Z,D,'k--','Color',[.5,.5,.5]); clabel(c,h); [c,h]=contour(DX,Z,T,'k-'); clabel(c,h);
xlim([0,23]); ylim([-360,0]); xlabel('Distance offshore [km]'); ylabel('Depth [m]');
fmg; contourf(DX,Z,CH); colorbar; titlename([cN,' Section ',num2str(ixix),' Chl_a']);
[c,h]=contour(DX,Z,D,'k--','Color',[.5,.5,.5]); clabel(c,h); [c,h]=contour(DX,Z,T,'k-'); clabel(c,h);
xlim([0,23]); ylim([-360,0]); xlabel('Distance offshore [km]'); ylabel('Depth [m]');
%}
clear c h
end;
DX=[]; Z=[]; D=[]; T=[]; N=[]; P=[]; SI=[]; CH=[]; clear DX Z D T N P SI CH
end; %for ixix=xsecix;
end; %for ccN = {'NF08','NF09'};
clear ccN cN
close all
% Plot property sections (distance vs. depth vs. concentration) of each cruise
for ccN = {'NF09'};
%for ccN = {'NF08','NF09'};
cN = ccN{:};
disp(cN);
% Divide cruise data into cross-shore sections
brkix = [0;find(diff(face.(cN).dx)<0);length(face.(cN).dx)]+1;
brkln = diff(brkix)-1;
xsecix=1:numel(brkln);
% Only plot cross-shore sections
if ( strcmpi(cN,'NF08') )
xsecix=[1,2,3,4,5,9];
elseif ( strcmpi(cN,'NF09') )
xsecix=[9,12,13,19,20];
end;
bath = plot_hires_bathymetry(bath,-[0:5:80,100:20:300,400:100:800],[20e3,60e3],true,@contour);
titlename(['Sections used from ',cN]);
bfig = gcf;
for ixix=xsecix;
if ( brkln(ixix) < 3 )
disp(['SKIPPING short (<3 station) section ',num2str(ixix)]);
%%%% EARLY CONTINUE
continue;
end;
secix = brkix(ixix):brkix(ixix)+brkln(ixix);
disp({ixix,brkix(ixix),brkln(ixix),min(face.(cN).z(secix)),max(face.(cN).z(secix)),});
figure(bfig);
for secixix=1:length(secix)-1
plot(face.(cN).lon(secix(secixix:secixix+1)),face.(cN).lat(secix(secixix:secixix+1)));
text(mean(face.(cN).lon(secix(secixix:secixix+1))),mean(face.(cN).lat(secix(secixix:secixix+1))),num2str(ixix));
end;
[DX,Z] = meshgrid(face.(cN).dx(secix),unique(face.(cN).z(secix)));
D  = griddata(face.(cN).dx(secix),face.(cN).z(secix),face.(cN).dena(secix),DX,Z);
T  = griddata(face.(cN).dx(secix),face.(cN).z(secix),face.(cN).t(secix),DX,Z);
N  = griddata(face.(cN).dx(secix),face.(cN).z(secix),face.(cN).n(secix),DX,Z);
P  = griddata(face.(cN).dx(secix),face.(cN).z(secix),face.(cN).p(secix),DX,Z);
SI = griddata(face.(cN).dx(secix),face.(cN).z(secix),face.(cN).si(secix),DX,Z);
CH = griddata(face.(cN).dx(secix),face.(cN).z(secix),face.(cN).chl(secix),DX,Z);
%fmg; contourf(DX,Z,T); colorbar; titlename([cN,' Section ',num2str(ixix),' T']);
if ( numel(find(~isnan(N))) == 0 )
disp(['SKIPPING section with no valid concentrations ',num2str(ixix)]);
else
fmg; contourf(DX,Z,N); colorbar; titlename([cN,' Section ',num2str(ixix),' N']);
[c,h]=contour(DX,Z,D,'k--','Color',[.5,.5,.5]); clabel(c,h); [c,h]=contour(DX,Z,T,'k-'); clabel(c,h);
xlim([0,23]); ylim([-360,0]); xlabel('Distance offshore [km]'); ylabel('Depth [m]');
for secixix=1:length(secix)-1
plot(face.(cN).lon(secix(secixix:secixix+1)),face.(cN).h(secix(secixix:secixix+1)),'LineWidth',2);
end;
%{
fmg; contourf(DX,Z,P); colorbar; titlename([cN,' Section ',num2str(ixix),' P']);
[c,h]=contour(DX,Z,D,'k--','Color',[.5,.5,.5]); clabel(c,h); [c,h]=contour(DX,Z,T,'k-'); clabel(c,h);
xlim([0,23]); ylim([-360,0]); xlabel('Distance offshore [km]'); ylabel('Depth [m]');
fmg; contourf(DX,Z,SI); colorbar; titlename([cN,' Section ',num2str(ixix),' Si']);
[c,h]=contour(DX,Z,D,'k--','Color',[.5,.5,.5]); clabel(c,h); [c,h]=contour(DX,Z,T,'k-'); clabel(c,h);
xlim([0,23]); ylim([-360,0]); xlabel('Distance offshore [km]'); ylabel('Depth [m]');
fmg; contourf(DX,Z,CH); colorbar; titlename([cN,' Section ',num2str(ixix),' Chl_a']);
[c,h]=contour(DX,Z,D,'k--','Color',[.5,.5,.5]); clabel(c,h); [c,h]=contour(DX,Z,T,'k-'); clabel(c,h);
xlim([0,23]); ylim([-360,0]); xlabel('Distance offshore [km]'); ylabel('Depth [m]');
%}
clear c h
end;
DX=[]; Z=[]; D=[]; T=[]; N=[]; P=[]; SI=[]; CH=[]; clear DX Z D T N P SI CH
end; %for ixix=xsecix;
end; %for ccN = {'NF08','NF09'};
clear ccN cN
close all
% Plot property sections (distance vs. depth vs. concentration) of each cruise
for ccN = {'NF09'};
%for ccN = {'NF08','NF09'};
cN = ccN{:};
disp(cN);
% Divide cruise data into cross-shore sections
brkix = [0;find(diff(face.(cN).dx)<0);length(face.(cN).dx)]+1;
brkln = diff(brkix)-1;
xsecix=1:numel(brkln);
% Only plot cross-shore sections
if ( strcmpi(cN,'NF08') )
xsecix=[1,2,3,4,5,9];
elseif ( strcmpi(cN,'NF09') )
xsecix=[9,12,13,19,20];
end;
bath = plot_hires_bathymetry(bath,-[0:5:80,100:20:300,400:100:800],[20e3,60e3],true,@contour);
titlename(['Sections used from ',cN]);
bfig = gcf;
for ixix=xsecix;
if ( brkln(ixix) < 3 )
disp(['SKIPPING short (<3 station) section ',num2str(ixix)]);
%%%% EARLY CONTINUE
continue;
end;
secix = brkix(ixix):brkix(ixix)+brkln(ixix);
disp({ixix,brkix(ixix),brkln(ixix),min(face.(cN).z(secix)),max(face.(cN).z(secix)),});
figure(bfig);
for secixix=1:length(secix)-1
plot(face.(cN).lon(secix(secixix:secixix+1)),face.(cN).lat(secix(secixix:secixix+1)));
text(mean(face.(cN).lon(secix(secixix:secixix+1))),mean(face.(cN).lat(secix(secixix:secixix+1))),num2str(ixix));
end;
[DX,Z] = meshgrid(face.(cN).dx(secix),unique(face.(cN).z(secix)));
D  = griddata(face.(cN).dx(secix),face.(cN).z(secix),face.(cN).dena(secix),DX,Z);
T  = griddata(face.(cN).dx(secix),face.(cN).z(secix),face.(cN).t(secix),DX,Z);
N  = griddata(face.(cN).dx(secix),face.(cN).z(secix),face.(cN).n(secix),DX,Z);
P  = griddata(face.(cN).dx(secix),face.(cN).z(secix),face.(cN).p(secix),DX,Z);
SI = griddata(face.(cN).dx(secix),face.(cN).z(secix),face.(cN).si(secix),DX,Z);
CH = griddata(face.(cN).dx(secix),face.(cN).z(secix),face.(cN).chl(secix),DX,Z);
%fmg; contourf(DX,Z,T); colorbar; titlename([cN,' Section ',num2str(ixix),' T']);
if ( numel(find(~isnan(N))) == 0 )
disp(['SKIPPING section with no valid concentrations ',num2str(ixix)]);
else
fmg; contourf(DX,Z,N); colorbar; titlename([cN,' Section ',num2str(ixix),' N']);
[c,h]=contour(DX,Z,D,'k--','Color',[.5,.5,.5]); clabel(c,h); [c,h]=contour(DX,Z,T,'k-'); clabel(c,h);
xlim([0,23]); ylim([-360,0]); xlabel('Distance offshore [km]'); ylabel('Depth [m]');
for secixix=1:length(secix)-1
plot(face.(cN).dx(secix(secixix:secixix+1)),face.(cN).h(secix(secixix:secixix+1)),'LineWidth',2);
end;
%{
fmg; contourf(DX,Z,P); colorbar; titlename([cN,' Section ',num2str(ixix),' P']);
[c,h]=contour(DX,Z,D,'k--','Color',[.5,.5,.5]); clabel(c,h); [c,h]=contour(DX,Z,T,'k-'); clabel(c,h);
xlim([0,23]); ylim([-360,0]); xlabel('Distance offshore [km]'); ylabel('Depth [m]');
fmg; contourf(DX,Z,SI); colorbar; titlename([cN,' Section ',num2str(ixix),' Si']);
[c,h]=contour(DX,Z,D,'k--','Color',[.5,.5,.5]); clabel(c,h); [c,h]=contour(DX,Z,T,'k-'); clabel(c,h);
xlim([0,23]); ylim([-360,0]); xlabel('Distance offshore [km]'); ylabel('Depth [m]');
fmg; contourf(DX,Z,CH); colorbar; titlename([cN,' Section ',num2str(ixix),' Chl_a']);
[c,h]=contour(DX,Z,D,'k--','Color',[.5,.5,.5]); clabel(c,h); [c,h]=contour(DX,Z,T,'k-'); clabel(c,h);
xlim([0,23]); ylim([-360,0]); xlabel('Distance offshore [km]'); ylabel('Depth [m]');
%}
clear c h
end;
DX=[]; Z=[]; D=[]; T=[]; N=[]; P=[]; SI=[]; CH=[]; clear DX Z D T N P SI CH
end; %for ixix=xsecix;
end; %for ccN = {'NF08','NF09'};
clear ccN cN
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
analyze_face_nutrients
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
quit
%-- 3/16/2016 9:48 PM --%
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
analyze_face_nutrients
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
analyze_face_nutrients
get(gca,'clim')
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
analyze_face_nutrients
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
analyze_face_nutrients
reviewanim([1:39],0,0,0)
1/16
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
help clabel
analyze_face_nutrients
get(0,'child')
`1.5/16
1.5/16
1/16
figsnamed
get(0,'child')
get(ans(1))
figs=get(0,'child')
get(figs,'Name')
class(get(figs,'Name'))
help regexp
help cellfun
figsnamed('Chl_a')
reviewanim(figsnamed('Chl_a'),[],[],[])
which reviewanmi
which reviewanim
ishandle(42)
ishandle(43)
ishandle(99)
reviewanim(figsnamed('Chl_a'),[],[],[])
reviewanim(figsnamed('Chl_a'),0,0,0)
reviewanim(figsnamed('Section 9 N'),0,0,0)
figsnamed('Section 9 N')
reviewanim(figsnamed('Section 9 N'),0,0,0)
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
help get
figure/get
help figure/get
help handle/get
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
% Plot property sections (distance vs. depth vs. concentration) of each cruise
%for ccN = {'NF09'};
for ccN = {'NF08','NF09'};
cN = ccN{:};
disp(cN);
% Divide cruise data into cross-shore sections
brkix = [0;find(diff(face.(cN).dx)<0);length(face.(cN).dx)]+1;
brkln = diff(brkix)-1;
xsecix=1:numel(brkln);
% Only plot cross-shore sections
if ( strcmpi(cN,'NF08') )
xsecix=[1,2,3,4,5,9];
elseif ( strcmpi(cN,'NF09') )
xsecix=[9,12,13,19,20];
end;
bath = plot_hires_bathymetry(bath,-[0:5:80,100:20:300,400:100:800],[20e3,60e3],true,@contour);
titlename(['Sections used from ',cN]);
bfig = gcf;
for ixix=xsecix;
if ( brkln(ixix) < 3 )
disp(['SKIPPING short (<3 station) section ',num2str(ixix)]);
%%%% EARLY CONTINUE
continue;
end;
secix = brkix(ixix):brkix(ixix)+brkln(ixix);
disp({ixix,brkix(ixix),brkln(ixix),min(face.(cN).z(secix)),max(face.(cN).z(secix)),});
figure(bfig);
for secixix=1:length(secix)-1
plot(face.(cN).lon(secix(secixix:secixix+1)),face.(cN).lat(secix(secixix:secixix+1)));
text(mean(face.(cN).lon(secix(secixix:secixix+1))),mean(face.(cN).lat(secix(secixix:secixix+1))),num2str(ixix));
end;
nd;
end;
end;
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
analyze_face_nutrients
reviewanim(figsnamed('NO_'),0,0,0)
reviewanim(figsnamed('[0-9] N'),0,0,0)
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
analyze_face_nutrients
nansummary(face.NF08.s)
nansummary(face.NF09.s)
size(s)
size(S)
s = [face.NF08.s;face.NF09.s];
fmg; hist(s,100);
numel(find(s<35))
numel(find(s<34))
numel(find(s<34.5))
numel(find(s<34.9))
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
analyze_face_nutrients
\
figsnamed('Section 9 N')
figsnamed('NF09.*N')
reviewanim(ans,0,0,0)
nansummary(h)
nansummary(face.NF08.h)
nansummary(face.NF09.h)
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
quit
%-- 3/17/2016 2:26 PM --%
load_all_sefcri
pd
pd +#
pd +3
pd
load_all_sefcri
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
load_all_sefcri
sefcri
sefcri.dc1
find_date_ranges(sefcri.dc1.t.date,10)
find_date_ranges(sefcri.bc1.t.date,10)
find_date_ranges(sefcri.updb.t.date,10)
sefcri
foreach cs=fieldnames(sefcri); sefcri.(cs{:}).depth, end;
foreach cs=fieldnames(sefcri); disp({cs{:},sefcri.(cs{:}).depth,}); end;
for cs=fieldnames(sefcri); disp({cs{:},sefcri.(cs{:}).depth,}); end;
for cs=fieldnames(sefcri)'; disp({cs{:},sefcri.(cs{:}).depth,}); end;
timenow
fmg; for cs=fieldnames(sefcri)'; disp({cs{:},sefcri.(cs{:}).depth,}); plot(sefcri.(cs{:}).t.date,sefcri.(cs{:}).t.data); end; datetick3;
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
analyze_face_nutrients
pd
analyze_face_nutrients
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
quit
%-- 3/17/2016 4:14 PM --%
extract_sfomc_etc
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
pd
extract_sfomc_etc
q
reviewanim(ans,0,0,0)
ans=[]
reviewanim(ans,0,0,0)
which ans
clear ans
which ans
type ans
sumix=find(get_season(res.ne_buoy.at_45m.sbe_seatemp.date)==3);
res.ne_buoy.at_45m.sbe_seatemp
min(diff(res.ne_buoy.at_45m.sbe_seatemp.date)
min(diff(res.ne_buoy.at_45m.sbe_seatemp.date))
min(diff(res.ne_buoy.at_45m.sbe_seatemp.date))*24
min(diff(res.ne_buoy.at_45m.sbe_seatemp.date))*24*60
median(diff(res.ne_buoy.at_45m.sbe_seatemp.date))*24*60
sumix=find(get_season(res.ne_buoy.at_45m.sbe_seatemp.date)==3);
numel(find(res.ne_buoy.at_45m.sbe_seatemp.data(sumix)<20))
numel(find(res.ne_buoy.at_45m.sbe_seatemp.data(sumix)<23))
oldsumix=sumix;
sumix=find(get_season(res.ne_buoy.at_40m.sbe_seatemp.date)==3);
all(sumix==oldsumix)
[at45,at40] = intersect_tses(res.ne_buoy.at_45m.sbe_seatemp,res.ne_buoy.at_40m.sbe_seatemp);
sumix=find(get_season(at40.date)==3);
clear oldsumix
sumix=find(get_season(at45.date)==3);
sumix=find(get_season(at40.date)==3);
numel(find(at40.data<23))
numel(find(at40.data<23))*5/60
numel(find(at40.data<23))*5/60/24
5450*5
5450*5/60/24
numel(find(at40.data<20))*5/60/24
[at30,at45,at40] = intersect_tses(res.ne_buoy.at_30m.sbe_seatemp,res.ne_buoy.at_45m.sbe_seatemp,res.ne_buoy.at_40m.sbe_seatemp);
sumix=find(get_season(at30.date)==3);
numel(find(at30.data<20))*5/60/24
numel(find(at30.data(sumix)<20))*5/60/24
numel(find(at30.data(sumix)<23))*5/60/24
[at20,at30,at45,at40] = intersect_tses(res.ne_buoy.at_20m.sbe_seatemp,res.ne_buoy.at_30m.sbe_seatemp,res.ne_buoy.at_45m.sbe_seatemp,res.ne_buoy.at_40m.sbe_seatemp);
sumix=find(get_season(at20.date)==3);
numel(find(at20.data(sumix)<23))*5/60/24
numel(find(at20.data(sumix)<24))*5/60/24
numel(find(at20.data(sumix)<24))*5/60
numel(unique(floor(at20.date(sumix))))
unique(get_year(at20.date(sumix)))
numel(find(at30.data(sumix)<23))*5/60/24
sumix=find(get_season(res.e_buoy.at_30m.sbe_seatemp.date)==3);
sumix=find(get_season(at20.date)==3);
esumix=find(get_season(res.e_buoy.at_30m.sbe_seatemp.date)==3);
numel(find(res.e_buoy.at_30m.sbe_seatemp.data(esumix)<23))*5/60/24
numel(find(res.e_buoy.at_30m.sbe_seatemp.data(esumix)<23))
median(diff(res.e_buoy.at_30m.sbe_seatemp.date))*24*60
find_date_ranges(res.e_buoy.at_30m.sbe_seatemp.date,10)
res.se_buoy
find_date_ranges(res.se_buoy.at_35m.sbe_seatemp.date,10)
res
find_date_ranges(res.sw_buoy.at_15m.sbe_seatemp.date,10)
fmg; plot_ts(res.sw_buoy.at_15m.sbe_seatemp);
fmg; plot_ts(res.sw_buoy.at_10m.sbe_seatemp);
fmg; plot_ts(res.c_buoy.at_15m.sbe_seatemp);
fmg; boxplot_ts(res.nw_w_btm.adcp_seatemp);
fmg; plot_ts(res.nw_w_btm.adcp_seatemp);
fmg; plot_ts(res.ne_buoy.at_45m.sbe_seatemp);
numel(find(at40.data(sumix)<23))*5/60/24
numel(find(at45.data(sumix)<23))*5/60/24
numel(find((24-at45.data(sumix))>5))
numel(find((24-at45.data(sumix))>5))*5/60/24
elevix=find((24-at45.data(sumix))>1);
elevix=find((24-at45.data)>1);
nt = (24-at45.data(sumix))>1;
nt = (24-at45.data(sumix));
timenow
elevix=find(nt>1);
nt(nt<0)=0;
elevix=find(nt>1);
fmg; plot(nt);
fmg; plot(nt); datetick3;
fmg; plot(at45.date(sumix),nt); datetick3;
sum(nt)
sum(nt)/1000
sum(nt)*20/1000
help convtemp
help convacc
help convmass
help convlength
sum(nt)*20*14
sum(nt)*20*14/1e3
ans/16
1000/16
sum(nt)*20*14/1e3
station_dist('fwyf1','lkwf1')
station_dist('fwyf1','lkwf1')*1e3
station_dist('fwyf1','lkwf1')*1e3*20
sum(nt)*20*14 % day-degrees below 24C times 20 m high water times N atomic mass
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
pd ../Sediment
dir *.m
ttpe ef_eval_ool_250m
type ef_eval_ool_250m
stn = ef_eval_ool_250m('pomf1');
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
more off
stn = ef_eval_ool_250m('pomf1');
reviewanim([],0,0,0)
[d,m] = lastwarn
oversample_ww3_to_ool
ef_eval_ool_250m_do_figures
sit
dbstop catchwarn
ef_eval_ool_250m_do_figures
msg
e
dbquit
close all
dbstatus
ef_eval_ool_250m_do_figures
e
e.stack
e.stack(:)
e.stack.name
dbup
fld1
fld2
numel(find(~isnan(fld1.data)))
numel(find(~isnan(fld2.data)))
dbup
ts
numel(find(~isnan(ts.data)))
stn.(dataset).hs
stn.(dataset)
dbcont
e.identifier
e.message
warnstr
e
warnstr
e.identifier
help warning
dbcont
close all
clear ans cbh d dataset edgeix lonix m seas seasix ts
close all
ef_eval_ool_250m_do_figures
dbstop catchwarn
close all
clear ans cbh seas seasix ts
dbstop catchwarn
ef_eval_ool_250m_do_figures
warnstr
warning('Foo c:\P');
warning('foo:bar','Foo c:\P');
warning('foo:bar','Foo c:\\P');
dbquit
ef_eval_ool_250m_do_figures
warning(warnstr)
e.identifier
warnstr
warning(e.identifier,warnstr)
warning(e.identifier,'%s',warnstr)
e.identifier
warning(e.identifier,warnstr)
warnstr
warning(e.identifier,'%s',warnstr)
warning(e.identifier,'W %s',warnstr)
help warning
more on
help warning
help sprintf
warning(e.identifier,sprintf(warnstr))
warning(e.identifier,sprintf('%s',warnstr))
warnstr
help sprintf
sprintf('%s',warnstr)
sprintf(warnstr)
warning(e.identifier,sprintf('%s',warnstr))
warnstr = sprintf('%s\n%s\n%s\n',warnstr,e.identifier,e.message);
warning(e.identifier,sprintf('%s',warnstr))
warning(e.identifier,warnstr)
warnstr
e.identifier
warnstr = ['CAUGHT ERROR: ',msg];
warnstr = sprintf('%s\n%s\n%s\n',warnstr,e.identifier,e.message);
warning(e.identifier,warnstr)
warnstr
warning(warnstr)
warning(e.identifier)
warning(e.identifier,warnstr)
e
help warning
s = warning('query',e.identifier)
dbquit
clear ans cbh seas seasix ts
ef_eval_ool_250m_do_figures
s
s = warning('query','foobie:bletch')
dbquit
close all
clear ans cbh seas seasix ts
ef_eval_ool_250m_do_figures
help set_more
close all
ef_eval_ool_250m_do_figures
close all
ef_eval_ool_250m_do_figures
close all
ef_eval_ool_250m_do_figures
dbstop catchwarn
close all
ef_eval_ool_250m_do_figures
warnstr
e.identifier
e.identifier{:}
e.identifier(:)
e.identifier
warnstr
warning(e.identifier,warnstr)
[d,m] = lastwarn
dbquit
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
quit
%-- 3/18/2016 2:09 PM --%
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
extract_sfomc_etc
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
stn = ef_eval_ool_250m('pomf1');
dbstop catchwarn
ef_eval_ool_250m_do_figures
e
warnstr
help findstr
findstr(warnstr,'\P')
[d,m] = lastwarn
type warning
findstr(warnstr,'\\P')
sprintf('\P')
[d,m] = lastwarn
dbquit
close all
help warning
more on
help warning
close all
clear ans cbh seas seasix ts
ef_eval_ool_250m_do_figures
close all
ef_eval_ool_250m_do_figures
dbclear all
dbstop scatter_fit_ts_seasons 196
close all
ef_eval_ool_250m_do_figures
axLms
mn
mx
isvector([])
isvector([1])
isvector(1)
ismatrix(1)
ismatrix([1])
ismatrix([1,2;3,4])
dbup
dbquit
close all
clear ans cbh seas seasix ts
help nanmin
diff([1,3,2,nan,6,5])
ef_eval_ool_250m_do_figures
which catchwarn
close all
more on
ts.date = stn.(dataset).date; ts.data = stn.(dataset).hs(:,stn.(dataset).latix,stn.(dataset).lonix);
numel(find(isnan(ts.data)))
numel(find(~isnan(ts.data)))
help interp_field
ts=[]; clear ts
ts.data = interp_field(stn.(dataset).lon,stn.(dataset).lat,stn.(dataset).hs,stn.lon,stn.lat);
stn.(dataset)
help interp_field
ts.data = interp_field(stn.(dataset).lat,stn.(dataset).lon,stn.(dataset).hs,stn.lat,stn.lon);
ts
ts.date = stn.(dataset).date; ts.data = interp_field(stn.(dataset).lat,stn.(dataset).lon,stn.(dataset).hs,stn.lat,stn.lon);
numel(find(~isnan(ts.data)))
ts.date = stn.(dataset).date; ts.data = interp_field(stn.(dataset).lat,stn.(dataset).lon,stn.(dataset).hs,stn.lat,stn.lon,@nanmean);
numel(find(~isnan(ts.data)))
fmg; contour_field(stn.(dataset),'hs',1);
fmg; contour_field(stn.(dataset),'hs',@nanmean);
help contour_field
ishandle(@nanmean)
isfunction(@nanmean)
help isa
isa(@nanmean,'function_handle')
fmg; contour_field(stn.(dataset),'hs',@nanmean);
fmg; contour_field(stn.(dataset),'hs',@nanmean); plot(stn.lon,stn.lat,'kp',stn.lon,stn.lat,'wd');
fmg; contour_field(stn.(dataset),'hs',{@nanmean,3});
fmg; contour_field(stn.(dataset),'hs',@nanmean); plot(stn.lon,stn.lat,'kp',stn.lon,stn.lat,'wd');
fmg; contour_field(stn.ww3,[],@nanmean); plot(stn.lon,stn.lat,'kp',stn.lon,stn.lat,'wd');
stn.ww3
fmg; contour_field(stn.ww3,'hs',@nanmean); plot(stn.lon,stn.lat,'kp',stn.lon,stn.lat,'wd');
stn.ww3.hs
fmg; contour_field(stn.ww3,stn.ww3.hs.field,@nanmean); plot(stn.lon,stn.lat,'kp',stn.lon,stn.lat,'wd');
fmg; contour_field(stn.ww3,stn.ww3.hs.field,@nanmean); plot(stn.lon,stn.lat,'kp',stn.lon,stn.lat,'wd'); axis(min(stn.(dataset).lon),max(stn.(dataset).lon),min(stn.(dataset).lat),max(stn.(dataset).lat));
fmg; contour_field(stn.ww3,stn.ww3.hs.field,@nanmean); plot(stn.lon,stn.lat,'kp',stn.lon,stn.lat,'wd'); axis([min(stn.(dataset).lon),max(stn.(dataset).lon),min(stn.(dataset).lat),max(stn.(dataset).lat)]);
newfld = oversample_attenuate_field(stn.ww3,5,stn.(dataset),stn.(dataset).depth)
newfld = oversample_attenuate_field({stn.ww3.lon,stn.ww3.lat,stn.ww3.hs.date,stn.ww3.hs.date,},5,stn.(dataset),stn.(dataset).depth)
dbstop if error
newfld = oversample_attenuate_field({stn.ww3.lon,stn.ww3.lat,stn.ww3.hs.date,stn.ww3.hs.date,},5,stn.(dataset),stn.(dataset).depth)
lostr
histr
[flat,fdts,flon] = meshgrid(fldlat,histr.date,fldlon);
dbquit
newfld = oversample_attenuate_field({stn.ww3.lon,stn.ww3.lat,stn.ww3.hs.date,stn.ww3.hs.date,},5,stn.(dataset),stn.(dataset).depth)
dbquit
newfld = oversample_attenuate_field({stn.ww3.lon,stn.ww3.lat,stn.ww3.hs.date,stn.ww3.hs.date,},5,stn.(dataset),stn.(dataset).depth)
dbup
lostr
dbquit
newfld = oversample_attenuate_field({stn.ww3.lon,stn.ww3.lat,stn.ww3.hs.date,stn.ww3.hs.field},5,stn.(dataset),stn.(dataset).depth)
dep
dbquit
stn.(dataset)
help ismatrix
newfld = oversample_attenuate_field({stn.ww3.lon,stn.ww3.lat,stn.ww3.hs.date,stn.ww3.hs.field},5,stn.(dataset),stn.(dataset).depth)
ans=[]; clear ans
clear ts
size(newfld)
tic,
newfld = oversample_attenuate_field({stn.ww3.lon,stn.ww3.lat,stn.ww3.hs.date,stn.ww3.hs.field},5,stn.(dataset),stn.(dataset).depth);
toc,
newfld=[]; clear newfld
clear ans cbh seas seasix ts
ts.date = stn.(dataset).date; ts.data = interp_field(stn.(dataset).lat,stn.(dataset).lon,stn.(dataset).hs,stn.lat,stn.lon,@nanmean);
numel(find(~isnan(ts.data)))
fmg; contour_field(stn.(dataset),'hs',@nanmean); plot(stn.lon,stn.lat,'kp',stn.lon,stn.lat,'wd');
stn.(dataset).hs = oversample_attenuate_field(stn.ww3,5,stn.(dataset),stn.(dataset).depth);
dbquit
stn.(dataset).hs = oversample_attenuate_field({stn.ww3.lon,stn.ww3.lat,stn.ww3.hs.date,stn.ww3.hs.field},5,stn.(dataset),stn.(dataset).depth);
stn.(dataset)
stn.(dataset).hs = oversample_attenuate_field({stn.ww3.lon,stn.ww3.lat,stn.ww3.hs.date,stn.ww3.hs.field},5,stn.(dataset),stn.(dataset).depth);
ans=[]; clear ans
ts.date = stn.(dataset).date; ts.data = interp_field(stn.(dataset).lat,stn.(dataset).lon,stn.(dataset).hs,stn.lat,stn.lon,@nanmean);
fmg; contour_field(stn.(dataset),'hs',@nanmean); plot(stn.lon,stn.lat,'kp',stn.lon,stn.lat,'wd');
oversample_ww3_to_ool
fmg; contour_field(stn.(dataset),'hs',@nanmean); plot(stn.lon,stn.lat,'kp',stn.lon,stn.lat,'wd');
stn.(dataset).hs = oversample_attenuate_field({stn.ww3.lon,stn.ww3.lat,stn.ww3.hs.date,stn.ww3.hs.field},6,stn.(dataset),stn.(dataset).depth);
fmg; contour_field(stn.(dataset),'hs',@nanmean); plot(stn.lon,stn.lat,'kp',stn.lon,stn.lat,'wd');
stn.(dataset).hs = oversample_attenuate_field({stn.ww3.lon,stn.ww3.lat,stn.ww3.hs.date,stn.ww3.hs.field},10,stn.(dataset),stn.(dataset).depth);
fmg; contour_field(stn.(dataset),'hs',@nanmean); plot(stn.lon,stn.lat,'kp',stn.lon,stn.lat,'wd');
stn.(dataset).hs = oversample_attenuate_field({stn.ww3.lon,stn.ww3.lat,stn.ww3.hs.date,stn.ww3.hs.field},5,stn.(dataset),stn.(dataset).depth);
fmg; contour_field(stn.(dataset),'hs',@nanmean); plot(stn.lon,stn.lat,'kp',stn.lon,stn.lat,'wd'); appendtitlename(' 5x');
ts.date = stn.(dataset).date; ts.data = interp_field(stn.(dataset).lat,stn.(dataset).lon,stn.(dataset).hs,stn.lat,stn.lon,@nanmean);
numel(find(~isnan(ts.data)))
stn.(dataset).hs = oversample_attenuate_field({stn.ww3.lon,stn.ww3.lat,stn.ww3.hs.date,stn.ww3.hs.field},10,stn.(dataset),stn.(dataset).depth);
ts.date = stn.(dataset).date; ts.data = interp_field(stn.(dataset).lat,stn.(dataset).lon,stn.(dataset).hs,stn.lat,stn.lon,@nanmean);
numel(find(~isnan(ts.data)))
tic,
stn.(dataset).hs = oversample_attenuate_field({stn.ww3.lon,stn.ww3.lat,stn.ww3.hs.date,stn.ww3.hs.field},10,stn.(dataset),stn.(dataset).depth);
toc,
tic,
stn.(dataset).hs = oversample_attenuate_field({stn.ww3.lon,stn.ww3.lat,stn.ww3.hs.date,stn.ww3.hs.field},5,stn.(dataset),stn.(dataset).depth);
toc,
ts.date = stn.(dataset).date; ts.data = interp_field(stn.(dataset).lat,stn.(dataset).lon,stn.(dataset).hs,stn.lat,stn.lon,{@nanmean,2});
numel(find(~isnan(ts.data)))
ts.date = stn.(dataset).date; ts.data = interp_field(stn.(dataset).lat,stn.(dataset).lon,stn.(dataset).hs,stn.lat,stn.lon,{@nanmean,1});
numel(find(~isnan(ts.data)))
ts.date = stn.(dataset).date; ts.data = interp_field(stn.(dataset).lat,stn.(dataset).lon,stn.(dataset).hs,stn.lat,stn.lon,{@nanmean,1,2});
numel(find(~isnan(ts.data)))
ts.date = stn.(dataset).date; ts.data = interp_field(stn.(dataset).lat,stn.(dataset).lon,stn.(dataset).hs,stn.lat,stn.lon,{@nanmean,2,1});
numel(find(~isnan(ts.data)))
help interp_field
ts.date = stn.(dataset).date; ts.data = interp_field(stn.(dataset).lat,stn.(dataset).lon,stn.(dataset).hs,stn.lat,stn.lon,{@nanmean,2,1},[],true);
numel(find(~isnan(ts.data)))
ts.date = stn.(dataset).date; ts.data = interp_field(stn.(dataset).lat,stn.(dataset).lon,stn.(dataset).hs,stn.lat,stn.lon,{@nanmean,1,2},[],true);
numel(find(~isnan(ts.data)))
ts.date = stn.(dataset).date; ts.data = interp_field(stn.(dataset).lat,stn.(dataset).lon,stn.(dataset).hs,stn.lat,stn.lon,{@nanmean,2},[],true);
numel(find(~isnan(ts.data)))
ts.date = stn.(dataset).date; ts.data = interp_field(stn.(dataset).lat,stn.(dataset).lon,stn.(dataset).hs,stn.lat,stn.lon,{@nanmean,2,1});
numel(find(~isnan(ts.data)))
dbstop oversample_attenuate_field 109
clear ans cbh seas seasix ts
clear edgeix lonix
close all
ef_eval_ool_250m_do_figures
dbstatus
dbclear all
dbstop oversample_attenuate_field 109
dbstop oversample_attenuate_field
dbclear all
dbstop oversample_attenuate_field
dbclear all
dbstop oversample_attenuate_field 109
ef_eval_ool_250m_do_figures
55*4
263.20
263.20 + 734.46 + 313.90 + 220 + 100
roundn( 263.20 + 734.46 + 313.90 + 220 + 100 , -2 )
roundn( 263.20 + 734.46 + 313.90 + 220 + 100 , -1 )
roundn( 263.20 + 734.46 + 313.90 + 220 + 100 , 0 )
roundn( 263.20 + 750 + 313.90 + 220 + 100 , 0 )
roundn( 263.20 + 750 + 313.90 + 220 + 103 , 0 )
roundn( 263.20 + 750 + 313.90 + 220 + 103 , 2 )
roundn( 263.20 + 750 + 313.90 + 220 + 103 , -2 )
roundn( 264 + 750 + 314 + 220 + 103 , -2 )
stn
res.CI
stn.CI
ans=[]; clear ans
stn.CI
close all
ef_eval_ool_250m_do_figures
0.15^2
0.3^2
type ef_eval_ool_250m.m
type ef_ool_250m.m
type ef_eval_ool_250m_do_figures.m
type ef_eval_ool_250m.m
ef_ool_250m
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
quit
%-- 3/19/2016 6:12 PM --%
stn = ef_eval_ool_250m('pomf1');
ef_eval_ool_250m_do_figures
which oversample_attenuate_field
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
stn = ef_eval_ool_250m('pomf1');
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
dbstop oversample_attenuate_field 107
stn = ef_eval_ool_250m('pomf1');
size(newfld)
dbquit
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
stn = ef_eval_ool_250m('pomf1');
ef_eval_ool_250m_do_figures
ef_eval_ool_250m_run_efs
close all
size(stn.ww3.ef_value)
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
stn = ef_eval_ool_250m('pvgf1');
ef_eval_ool_250m_do_figures
fmg; contour_field(stn.(dataset),'hs',@nanmean); plot(stn.lon,stn.lat,'kp',stn.lon,stn.lat,'wd');
fmg; contour_field(stn.ww3,stn.ww3.hs.field,@nanmean); plot(stn.lon,stn.lat,'kp',stn.lon,stn.lat,'wd'); axis([min(stn.(dataset).lon),max(stn.(dataset).lon),min(stn.(dataset).lat),max(stn.(dataset).lat)]);
set(gca,'clim',[0,.8])
fmg; contour_field(stn.ww3,stn.ww3.hs.field,@nanmean); plot(stn.lon,stn.lat,'kp',stn.lon,stn.lat,'wd'); axis([min(stn.(dataset).lon),max(stn.(dataset).lon),min(stn.(dataset).lat),max(stn.(dataset).lat)]); set(gca,'clim',[0,1])
fmg; contour_field(stn.(dataset),'hs',@nanmean); plot(stn.lon,stn.lat,'kp',stn.lon,stn.lat,'wd'); set(gca,'clim',[0,1])
fmg; contour_field(stn.(dataset),'hs',@nanmean); plot(stn.lon,stn.lat,'kp',stn.lon,stn.lat,'wd'); set(gca,'clim',[0,1]); titlename('Attenuated');
fmg; contour_field(stn.ww3,stn.ww3.hs.field,@nanmean); plot(stn.lon,stn.lat,'kp',stn.lon,stn.lat,'wd'); axis([min(stn.(dataset).lon),max(stn.(dataset).lon),min(stn.(dataset).lat),max(stn.(dataset).lat)]); set(gca,'clim',[0,1]) titlename('Model');
fmg; contour_field(stn.ww3,stn.ww3.hs.field,@nanmean); plot(stn.lon,stn.lat,'kp',stn.lon,stn.lat,'wd'); axis([min(stn.(dataset).lon),max(stn.(dataset).lon),min(stn.(dataset).lat),max(stn.(dataset).lat)]); set(gca,'clim',[0,1]); titlename('Model');
fmg; contour_field(stn.ww3,stn.ww3.hs.field,@nanmean); plot(stn.lon,stn.lat,'kp',stn.lon,stn.lat,'wd'); axis([min(stn.(dataset).lon),max(stn.(dataset).lon),min(stn.(dataset).lat),max(stn.(dataset).lat)]); set(gca,'clim',[0,1]) titlename('Model');
fmg; contour_field(stn.ww3,stn.ww3.hs.field,@nanmean); plot(stn.lon,stn.lat,'kp',stn.lon,stn.lat,'wd'); axis([min(stn.(dataset).lon),max(stn.(dataset).lon),min(stn.(dataset).lat),max(stn.(dataset).lat)]); set(gca,'clim',[0,1]); titlename('Model');
fmg; contour_field(stn.(dataset),'hs',@nanmean); plot(stn.lon,stn.lat,'kp',stn.lon,stn.lat,'wd'); set(gca,'clim',[0,1]); titlename('Attenuated');
ef_eval_ool_250m_run_efs
size(stn.ww3.ef_value)
size(stn.(Dataset).ef_value)
size(stn.(dataset).ef_value)
stn.(dataset)
size(stn.(dataset).ef_value)
ef_eval_ool_250m_run_efs
size(stn.(dataset).ef_value)
size(stn.ww3.ef_value)
find_date_ranges(stn.(dataset).hs.date,10)
find_date_ranges(stn.(dataset).date,10)
stn.(dataset)
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
pd ../../FAV
quit
%-- 3/20/2016 11:33 AM --%
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
pbbe1 = get_station_from_station_name('pbbe1');
pbbe1 = plot_hires_bathymetry(pbbe1,-[0:2:30,50:50:300],[20e3,20e3]);
pbbe1 = plot_hires_bathymetry(pbbe1,-[0:2:80],[20e3,20e3]);
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
pbbe1 = get_station_from_station_name('pbbe1');
extract_ww3_multi_region([],'berm',[]);
pbbe1 = get_ww3_multi_region([],'berm',[],[]);
pbbe1
old_pbbe1=pbbe1
pbbe1=[]; clear pbbe1
berm = get_ww3_multi_region([],'berm',[],[]);
berm
help save
berm = get_ww3_multi_region([],'berm',[],[]);
%-- 3/20/2016 3:08 PM --%
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
berm = get_ww3_multi_region([],'berm',[],[]);
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
berm = get_ww3_multi_region([],'berm',[],[]);
pbbe1 = get_station_from_station_name('pbbe1');
pbbe1.ww3 = get_ww3_multi_region([],'berm',[],[]);
berm=[]; clear berm
pbbe1 = plot_hires_bathymetry(pbbe1,-[0:2:80],[20e3,20e3]);
print('-dpng','pbbe1-bermuda-bathymetry.png');
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
pbbe1 = plot_hires_bathymetry(pbbe1,-[0:2:80],[30e3,30e3]);
pbbe1 = get_station_from_station_name('pbbe1');
pbbe1 = plot_hires_bathymetry(pbbe1,-[0:2:80],[30e3,30e3]);
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
swbe1 = get_station_from_station_name('swbe1');
swbe1 = plot_hires_bathymetry(swbe1,-[0:2:80],[30e3,30e3]);
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
swbe1 = get_station_from_station_name('swbe1');
swbe1 = plot_hires_bathymetry(swbe1,-[2:2:80],[30e3,30e3]);
print('-dpng','swbe1-bermuda-bathymetry.png');
swbe1.ww3 = seasonalize_ww3_region('berm',false,false,false);
swbe1.ww3 = seasonalize_ww3_region('ber',false,false,false);
swbe1.ww
swbe1.ww3
swbe1
ans=[]; clear ans
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
swbe1 = get_station_from_station_name('swbe1');
swbe1 = plot_hires_bathymetry(swbe1,-[2:2:80],[30e3,30e3]);
swbe1.ww3 = seasonalize_ww3_region('ber',true,true,false);
swbe1.ww3
swbe1
swbe1.ww3 = seasonalize_ww3_region('ber',false,false,false);
swbe1.ww3 = seasonalize_ww3_region('ber',true,true,false);
swbe1.ww3 = seasonalize_ww3_region(swbe1,true,true,false);
dbstop if error
swbe1.ww3 = seasonalize_ww3_region(swbe1,true,true,false);
dbup
clm
clm.hs
numel(find(~isnan(clm.hs)))
numel(find(~isnan(clm.hs.field)))
dbup
numel(find(~isnan(clm.hs.field)))
dbquit
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
swbe1 = get_station_from_station_name('swbe1');
swbe1.ww3 = seasonalize_ww3_region('ber',false,false,false);
swbe1.ww3
swbe1.ww3.hs
numelswbe1.ww3.hs.field))
numel(find(~isnan(swbe1.ww3.hs.field)))
pd ../../Postdoc/data
dir *2015*at_4*
dir *at_4*2015*
numel(find(~isnan(swbe1.ww3.dp.field)))
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
load ww3_at_4m_berm_-67.000_-62.000_30.000_34.000_200502010000_201507312100.mat
fld
fmg; contourf(fld.lon,fld.lat,squeeze(nanmean(fld.hs.field))); colorbar;
dir ww3/*201507*
numel(find(~isnan(fld.dp.field)))
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
nc = mDataset('multi_1.at_4m.hs.201507.grb2');
nj_info(nc)
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
nc = mDataset('ww3/multi_1.at_4m.hs.201507.grb2');
nj_info(nc)
lon = nc{'lon'}(:);
lat = nc{'lat'}(:);
min(lon),max(lon),min(lat),max(lat)
min(360-lon),max(360-lon),min(lat),max(lat)
min(lon-360),max(lon-360),min(lat),max(lat)
nj_info(nc)
hs = cast('double',pi)
hs = cast(pi,'double')
hs = cast(nc{'Significant_height_of_combined_wind_waves_and_swell'}(1,:,:),'double');
fmg; contourf(lon,lat,hs); colorbar;
dirs
pd +2
dir *.png
dir *BERMUDA*.png
close(nc); clear nc
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
dir *BERMUDA*.png
x = load('A20160801735.1KM.BERMUDA.PASS.L3D_RRC.CI.png
x = load('A20160801735.1KM.BERMUDA.PASS.L3D_RRC.CI.png');
x = imread('A20160801735.1KM.BERMUDA.PASS.L3D_RRC.CI.png');
size(x)
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
swbe1 = get_station_from_station_name('swbe1');
swbe1 = analyze_ool_250m(swbe1,[0.30,0.30],'BERMUDA');
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
swbe1 = get_station_from_station_name('swbe1');
swbe1 = analyze_ool_250m(swbe1,[0.30,0.30],'BERMUDA','SST');
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
swbe1 = get_station_from_station_name('swbe1');
swbe1 = analyze_ool_250m(swbe1,[0.30,0.30],'BERMUDA','CHL');
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
swbe1 = get_station_from_station_name('swbe1');
swbe1 = analyze_ool_250m(swbe1,[0.30,0.30],'BERMUDA','AFAI');
swbe1 = analyze_ool_250m(swbe1,[0.30,0.30],'BERMUDA','AFAI');
fmg; contourf(swbe1.AFAI.lon,swbe1.AFAI.lat,squeeze(nanmean(swbe1.AFAI.field))); colorbar;
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
swbe1 = get_station_from_station_name('swbe1'); swbe1 = analyze_ool_250m(swbe1,[0.30,0.30],'BERMUDA','CI');
fmg; contourf(swbe1.AFAI.lon,swbe1.AFAI.lat,squeeze(nanmean(swbe1.AFAI.field))); colorbar;
fmg; contourf(swbe1.CI.lon,swbe1.CI.lat,squeeze(nanmean(swbe1.CI.field))); colorbar;
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
swbe1 = get_station_from_station_name('swbe1'); swbe1 = analyze_ool_250m(swbe1,[0.30,0.30],'BERMUDA','SST');
help contour_field
fmg; contour_field(swbe1.SST,[],@nanmean);
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
type ef_eval_ool_250m
swbe1 = ef_eval_ool_250m('swbe1','SST',true,false);
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
swbe1 = ef_eval_ool_250m('swbe1','SST',true,false);
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
swbe1 = ef_eval_ool_250m('swbe1','SST',true,false);
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
more on
swbe1 = ef_eval_ool_250m('swbe1','SST',true,false);
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
swbe1 = ef_eval_ool_250m('swbe1','SST',true,false);
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
swbe1 = ef_eval_ool_250m('swbe1','SST',true,true);
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
quit
%-- 3/25/2016 4:34 PM --%
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
swbe1=ext
pd ../../Postdoc
pd ../../../Postdoc
swbe1 = get_station_from_station_name('swbe1');
help extract_ww3_multi_region
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
extract_ww3_multi_region([],'berm',[]);
swbe1 = get_station_from_station_name('swbe1');
help get_ww3_multi_region
pd
swbe1.ww3 = seasonalize_ww3_region(swbe1,true,false,false);
numel(find(~isnan(swbe1.ww3.dp.field)))
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
type get_all_station_metadata.m
pbbe1 = get_station_from_station_name('pbbe1');
pbbe1.ww3 = seasonalize_ww3_region(pbbe1,true,false,false);
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
stn = ef_eval_ool_250m('swbe1');
ef_eval_ool_250m_do_figures
axis([min(stn.(dataset).lon),max(stn.(dataset).lon),min(stn.(dataset).lat),max(stn.(dataset).lat)]);
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
stn = ef_eval_ool_250m('swbe1');
ef_eval_ool_250m_run_efs
stn.ww3
stn.ww3.hs
find_date_ranges(stn.ww3.hs.date,10)
find_date_ranges(stn.CI.date,10)
dbquit
ans=[];
clear ans dataset fld ndts sz2d sz3d
stn.ww3
stn.ww3.hs
stn.CI
stn.ww3
ef_eval_ool_250m_run_efs
doPrint
close all
ef_eval_ool_250m_run_efs
close all
ef_eval_ool_250m_run_efs
colorbar;
roi
ef_roi
axis(ef_roi)
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
stn = ef_eval_ool_250m('pvgf1');
ef_eval_ool_250m_do_figures
close all
ef_eval_ool_250m_run_efs
stn
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
dbstop ef_eval_ool_250m     plot(stn.lon,stn.lat,'rp');
dbstop ef_eval_ool_250m 73
help catchwarn
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
dbstop ef_eval_ool_250m 73
stn = ef_eval_ool_250m('pvgf1');
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
dbstop ef_eval_ool_250m 73
dbstop if error
stn = ef_eval_ool_250m('pvgf1');
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
dbstop get_ww3_multi_region
stn = ef_eval_ool_250m('pvgf1');
dbquit
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
help extract_ww3_multi_region
extract_ww3_multi_region([],'fknms',[]);
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
extract_ww3_multi_region([],'amsam',[]);
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
extract_ww3_multi_region([],'cnmi',[]);
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
extract_ww3_multi_region([],'hi',[]);
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
extract_ww3_multi_region([],'nwhi',[]);
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
extract_ww3_multi_region([],'pr_vi',[]);
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
stn = ef_eval_ool_250m('pvgf1');
ef_eval_ool_250m_run_efs
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
stn = get_station_from_station_name('pvgf1');
stn.ww3 = seasonalize_ww3_region(stn,doFigs,false,doPrint);
stn.ww3 = seasonalize_ww3_region(stn,false,false,false);
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
stn = get_station_from_station_name('pomf1');
stn.ww3 = seasonalize_ww3_region(stn,false,false,false);
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
stn = get_station_from_station_name('pgsp6');
stn.ww3 = seasonalize_ww3_region(stn,false,false,false);
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
vesp6 = get_station_from_station_name('vesp6');
stn.ww3 = seasonalize_ww3_region(stn,false,false,false);
stn = get_station_from_station_name('vesp6');
stn.ww3 = seasonalize_ww3_region(stn,false,false,false);
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
stn = get_station_from_station_name('vwsp6');
stn.ww3 = seasonalize_ww3_region(stn,false,false,false);
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
stn = get_station_from_station_name('ofsp6');
stn.ww3 = seasonalize_ww3_region(stn,false,false,false);
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
stn = get_station_from_station_name('llbp7');
stn.ww3 = seasonalize_ww3_region(stn,false,false,false);
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
stn = get_station_from_station_name('manm1');
stn.ww3 = seasonalize_ww3_region(stn,false,false,false);
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
stn = get_station_from_station_name('sasp7');
stn.ww3 = seasonalize_ww3_region(stn,false,false,false);
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
stn = get_station_from_station_name('tasp7');
stn.ww3 = seasonalize_ww3_region(stn,false,false,false);
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
stn = get_station_from_station_name('tinp7');
stn.ww3 = seasonalize_ww3_region(stn,false,false,false);
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
stn = get_station_from_station_name('lppr1');
stn.ww3 = seasonalize_ww3_region(stn,false,false,false);
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
% SE_FL
stn = ef_eval_ool_250m('pvgf1',[],true,true); ef_eval_ool_250m_run_efs
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
stn = ef_eval_ool_250m('pomf1',[],true,true); ef_eval_ool_250m_run_efs
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
% SAMOA
stn = ef_eval_ool_250m('pgsp6',[],true,true); ef_eval_ool_250m_run_efs
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
stn = ef_eval_ool_250m('pvgf1',[],true,true); ef_eval_ool_250m_run_efs
stn.ww3
stn.ww3.hs
stn.(dataset).date(end) >= stn.ww3.hs.date(end)
stn.(dataset).date(end), stn.ww3.hs.date(end)
datestr([stn.(dataset).date(end), stn.ww3.hs.date(end)])
datestr([stn.(dataset).date(end), stn.(dataset).ef_date(end)])
dbquit
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
dbstop if error\
dbstop if error
stn = ef_eval_ool_250m('pomf1',[],true,true); ef_eval_ool_250m_run_efs
[stn.station_name,'-',dataset,'-EF-events-N.png']
doPrint=true;
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
stn = ef_eval_ool_250m('pomf1',[],true,true,true);
reviewanim([],0,0,0)
stn
stn.CI
stn.CI.ef_events
doPrint=true;
dbstop ef_eval_ool_250m_run_efs
ef_eval_ool_250m_run_efs
newdat
newdat.CI
numel(find(~isnan(newdat.CI.hs)))
stn.ww3.ef_value
dbquit
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
dbstop ef_eval_ool_250m_run_efs 123
stn = ef_eval_ool_250m('pomf1',[],true,true,true);
stn.ww3.ef_criterion
numel(find(~isnan(stn.ww3.ef_criterion)))
numel(find(~isnan(stn.(dataset).ef_criterion)))
stn.(dataset)
150*53
numel(find(~isnan(stn.ww3.ef_criterion)))
stn.ww3
numel(find(~isnan(stn.ww3.ef_value)))
dbquit
close all
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
stn = ef_eval_ool_250m('pomf1',[],false,false,true);
stn.(dataset)
dataset='CI';
stn.(dataset)
stn.(dataset).hs = oversample_attenuate_field({stn.ww3.lon,stn.ww3.lat,stn.ww3.hs.date,stn.ww3.hs.field},5,stn.(dataset),stn.(dataset).depth);
numel(find(~isnan(stn.(dataset).hs)))
fmg; contourf(stn.(dataset).hs);
fmg; contour_field(stn.(dataset),'hs',@nanmean);
fmg; contour_field(stn.ww3,[],@nanmean);
fmg; contour_field(stn.ww3.hs,[],@nanmean);
stn.ww3
fmg; contour_field(stn.ww3,stn.ww3.hs.field,@nanmean);
fmg; contour_field(stn.ww3,stn.ww3.hs.field,@nanmean); plot_hires_coastline(stn);
fmg; contour_field(stn.ww3,stn.ww3.hs.field,@nanmean); plot_hires_coastline(stn.ngdc_hires_bathy);
dbstop oversample_attenuate_field
stn.(dataset).hs = oversample_attenuate_field({stn.ww3.lon,stn.ww3.lat,stn.ww3.hs.date,stn.ww3.hs.field},5,stn.(dataset),stn.(dataset).depth);
fmg; contourf(squeeze(nanmean(fldfld)));
fmg; contourf(flon,flat,squeeze(nanmean(fldfld)));
dbstatus
size(flat)
fmg; contourf(fldlon,fldlat,squeeze(nanmean(fldfld)));
fldfld = interp3(lostr.lat,lostr.date,lostr.lon,lostr.field,flat,fdts,flon,'spline');
fmg; contourf(fldlon,fldlat,squeeze(nanmean(fldfld)));
fldfld = interp3(lostr.lat,lostr.date,lostr.lon,lostr.field,flat,fdts,flon,'linear');
fmg; contourf(fldlon,fldlat,squeeze(nanmean(fldfld)));
help interp3
fldfld = interp3(lostr.lat,lostr.date,lostr.lon,lostr.field,flat,fdts,flon,'cubic');
fmg; contourf(fldlon,fldlat,squeeze(nanmean(fldfld)));
dbquit
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
quit
%-- 3/29/2016 10:11 PM --%
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
stn = ef_eval_ool_250m('pomf1',[],false,false,true);
fmg; contour_field(stn.(dataset),'hs',@nanmean); plot_hires_coastline(stn.ngdc_hires_bathy);
fmg; contour_field(stn.(dataset),'hs',@nanmean);
stn.(dataset)
dataset='CI';
fmg; contour_field(stn.(dataset),'hs',@nanmean); plot_hires_coastline(stn.ngdc_hires_bathy);
stn.(dataset).hs = oversample_attenuate_field({stn.ww3.lon,stn.ww3.lat,stn.ww3.hs.date,stn.ww3.hs.field},stn.(dataset),stn.(dataset).depth);
fmg; contour_field(stn.(dataset),'hs',@nanmean); plot_hires_coastline(stn.ngdc_hires_bathy);
which bbox2points
help bbox2points
help bbox2rect
fmg; contour_field(stn.(dataset),'hs',@nanmean); plot_hires_coastline(stn.ngdc_hires_bathy); plot(stn.lon,stn.lat,'kp',stn.lon,stn.lat,'wd'); axis(fieldbbox(stn.(dataset)));
fmg; contour_field(stn.(dataset),'hs',@nanmean); plot_hires_coastline(stn.ngdc_hires_bathy); plot(stn.lon,stn.lat,'kp',stn.lon,stn.lat,'wd'); axis(field_bbox(stn.(dataset)));
if ( 1 ); disp yeah; end;
if ( [1 0] ); disp yeah; end;
if ( [1 1] ); disp yeah; end;
if ( [1 1 ; 1 0] ); disp yeah; end;
if ( [1 1 ; 1 1] ); disp yeah; end;
if ( [1 1 ; 1 1] ~= 0 ); disp yeah; end;
if ( [1 1 ; 1 0] ~= 0 ); disp yeah; end;
help isscalar
size(pi)
lgutter = gutter(1);
rgutter = gutter(1);
bgutter = gutter(2);
tgutter = gutter(2);
fmg; contour_field(stn.(dataset),'hs',@nanmean); plot_hires_coastline(stn.ngdc_hires_bathy); plot(stn.lon,stn.lat,'kp',stn.lon,stn.lat,'wd'); axis(field_bbox(stn.(dataset)));
tic,
stn.(dataset).hs = oversample_attenuate_field({stn.ww3.lon,stn.ww3.lat,stn.ww3.hs.date,stn.ww3.hs.field},10,stn.(dataset),stn.(dataset).depth);
toc,
fmg; contour_field(stn.(dataset),'hs',@nanmean); plot_hires_coastline(stn.ngdc_hires_bathy); plot(stn.lon,stn.lat,'kp',stn.lon,stn.lat,'wd'); axis(field_bbox(stn.(dataset))); titlename('10x');
tic,
stn.(dataset).hs = oversample_attenuate_field({stn.ww3.lon,stn.ww3.lat,stn.ww3.hs.date,stn.ww3.hs.field},stn.(dataset),stn.(dataset).depth);
toc,
fmg; contour_field(stn.(dataset),'hs',@nanmean); plot_hires_coastline(stn.ngdc_hires_bathy); plot(stn.lon,stn.lat,'kp',stn.lon,stn.lat,'wd'); axis(field_bbox(stn.(dataset))); titlename('Full');
fmg; contour_field(stn.ww3,stn.ww3.hs,@nanmean); plot_hires_coastline(stn.ngdc_hires_bathy); plot(stn.lon,stn.lat,'kp',stn.lon,stn.lat,'wd'); axis(field_bbox(stn.(dataset))); titlename('Orig');
fmg; contour_field(stn.ww3,stn.ww3.hs.field,@nanmean); plot_hires_coastline(stn.ngdc_hires_bathy); plot(stn.lon,stn.lat,'kp',stn.lon,stn.lat,'wd'); axis(field_bbox(stn.(dataset))); titlename('Orig');
field_bbox(stn.(dataset))
field_bbox(stn.(dataset),'gutter')
field_bbox(stn.(dataset),1)
field_bbox(stn.(dataset),[1,2[)
field_bbox(stn.(dataset),[1,2])
field_bbox(stn.(dataset),[1,2,3,4])
dbstop oversample_attenuate_field 116
stn.(dataset).hs = oversample_attenuate_field({stn.ww3.lon,stn.ww3.lat,stn.ww3.hs.date,stn.ww3.hs.field},stn.(dataset),stn.(dataset).depth);
oldfld
oldfld = newfld;
newfld = interp3(flat,fdts,flon,fldfld,...
depstr.lat,histr.date,depstr.lon);
tic,
for lonix=1:numel(depstr.lon)
newfld(:,:,lonix) = interp3(flat,fdts,flon,fldfld,...
depstr.lat,histr.date,depstr.lon(lonix));
end;
toc,
newfld = interp3(flat,fdts,flon,fldfld,...
depstr.lat,histr.date,depstr.lon);
numel(find(newfld==oldfld))
numel(find(newfld~=oldfld))
numel(find(newfld(oldfld==0)==oldfld(oldfld==0)))
numel(find(newfld(oldfld==0)~=oldfld(oldfld==0)))
size(newfld), size(fldfld)
fmg; contourf(fldlon,fldlat,squeeze(nanmean(newfld)));
fmg; contourf(fldlon,fldlat,squeeze(nanmean(oldfld)));
newfld = interp3(flat,fdts,flon,fldfld,depstr.lat,histr.date,depstr.lon);
newfld = interp3(flat,fdts,flon,fldfld,depstr.lat,histr.date,depstr.lon,'spline');
newfld=[]; clear newfld
newfld = interp3(flat,fdts,flon,fldfld,depstr.lat,histr.date,depstr.lon,'spline');
newfld = interp3(flat,fdts,flon,fldfld,depstr.lat,histr.date,depstr.lon,'cubic');
help interp3
newfld = interp3(flat,fdts,flon,fldfld,depstr.lat,histr.date,depstr.lon);
tic,
newfld = interp3(flat,fdts,flon,fldfld,depstr.lat,histr.date,depstr.lon);
toc,
help memory
size(pi)
help whos
whos pi
x
x=pi;
whos x
clear x
help size
help sizem
which sizem
help sizedat
which sizedat
help idutils
help arx
doc argx
help whos
help memory
x = memory
x = memory; x.MaxPossibleArrayBytes
x = memory; x.MaxPossibleArrayBytes/8
dbclear all
dbcont
clear ans
help byte2any
help prod
stn.(dataset).hs = oversample_attenuate_field({stn.ww3.lon,stn.ww3.lat,stn.ww3.hs.date,stn.ww3.hs.field},stn.(dataset),stn.(dataset).depth);
fmg; contour_field(stn.(dataset),'hs',@nanmean); plot_hires_coastline(stn.ngdc_hires_bathy); plot(stn.lon,stn.lat,'kp',stn.lon,stn.lat,'wd'); axis(field_bbox(stn.(dataset))); titlename('Full 2');
stn.(dataset).hs = oversample_attenuate_field({stn.ww3.lon,stn.ww3.lat,stn.ww3.hs.date,stn.ww3.hs.field},stn.(dataset),stn.(dataset).depth);
fmg; contour_field(stn.(dataset),'hs',@nanmean); plot_hires_coastline(stn.ngdc_hires_bathy); plot(stn.lon,stn.lat,'kp',stn.lon,stn.lat,'wd'); axis(field_bbox(stn.(dataset))); titlename('Nearest');
titlename('Full 1');
stn.(dataset).hs = oversample_attenuate_field({stn.ww3.lon,stn.ww3.lat,stn.ww3.hs.date,stn.ww3.hs.field},stn.(dataset),stn.(dataset).depth);
fmg; contour_field(stn.(dataset),'hs',@nanmean); plot_hires_coastline(stn.ngdc_hires_bathy); plot(stn.lon,stn.lat,'kp',stn.lon,stn.lat,'wd'); axis(field_bbox(stn.(dataset))); titlename('Full');
stn.(dataset).hs = oversample_attenuate_field({stn.ww3.lon,stn.ww3.lat,stn.ww3.hs.date,stn.ww3.hs.field},stn.(dataset),stn.(dataset).depth);
tic,
stn.(dataset).hs = oversample_attenuate_field({stn.ww3.lon,stn.ww3.lat,stn.ww3.hs.date,stn.ww3.hs.field},stn.(dataset),stn.(dataset).depth);
toc,
tic,
stn.(dataset).hs = oversample_attenuate_field({stn.ww3.lon,stn.ww3.lat,stn.ww3.hs.date,stn.ww3.hs.field},stn.(dataset),stn.(dataset).depth);
toc,
dbstop oversample_attenuate_field 124
stn.(dataset).hs = oversample_attenuate_field({stn.ww3.lon,stn.ww3.lat,stn.ww3.hs.date,stn.ww3.hs.field},stn.(dataset),stn.(dataset).depth);
fmg; contourf(fldlon,fldlat,squeeze(nanmean(newfld)));
fmg; contourf(fldlon,fldlat,squeeze(nanmean(newfld))); daspect([1,cosd(fldlat(1)),1]);
fmg; contourf(fldlon,fldlat,squeeze(nanmean(fldfld))); daspect([1,cosd(fldlat(1)),1]);
fmg; contourf(lostr.lon,lostr.lon,squeeze(nanmean(lostr.field))); daspect([1,cosd(fldlat(1)),1]);
fmg; contourf(fldlon,fldlat,squeeze(nanmean(fldfld))); daspect([1,cosd(fldlat(1)),1]);
fmg; contourf(lostr.lon,lostr.lat,squeeze(nanmean(lostr.field))); daspect([1,cosd(fldlat(1)),1]);
fmg; contourf(fldlon,fldlat,squeeze(nanmean(fldfld))); daspect([1,cosd(fldlat(1)),1]);
fmg; contourf(lostr.lon,lostr.lat,squeeze(nanmean(lostr.field))); axis([min(fldlon),max(fldlon),min(fldlat),max(fldlat)]); daspect([1,cosd(fldlat(1)),1]);
help interp3
help griddata
help griddatan
help ndgrid
help slice
help griddatan
help griddata
tic,
newfld = griddata(flat,fdts,flon,fldfld,depstr.lat,histr.date,depstr.lon);
toc,
%-- 3/30/2016 10:03 PM --%
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
stn = ef_eval_ool_250m('pomf1',[],false,false,true);
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
stn = ef_eval_ool_250m('pomf1',[],false,false,true);
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
stn = ef_eval_ool_250m('pomf1',[],false,false,true);
fmg; contour_field(stn.(dataset),'hs',@nanmean); plot_hires_coastline(stn.ngdc_hires_bathy); plot(stn.lon,stn.lat,'kp',stn.lon,stn.lat,'wd'); axis(field_bbox(stn.(dataset))); titlename('Full');
dataset='CI';
fmg; contour_field(stn.(dataset),'hs',@nanmean); plot_hires_coastline(stn.ngdc_hires_bathy); plot(stn.lon,stn.lat,'kp',stn.lon,stn.lat,'wd'); axis(field_bbox(stn.(dataset))); titlename('Full');
stn.(dataset).hs = oversample_attenuate_field({stn.ww3.lon,stn.ww3.lat,stn.ww3.hs.date,stn.ww3.hs.field},stn.(dataset),stn.(dataset).depth,[],true);
fmg; contour_field(stn.(dataset),'hs',@nanmean); plot_hires_coastline(stn.ngdc_hires_bathy); plot(stn.lon,stn.lat,'kp',stn.lon,stn.lat,'wd'); axis(field_bbox(stn.(dataset))); titlename('No NaNs');
stn.(dataset).hs = oversample_attenuate_field({stn.ww3.lon,stn.ww3.lat,stn.ww3.hs.date,stn.ww3.hs.field},stn.(dataset),stn.(dataset).depth,[],false);
fmg; contour_field(stn.(dataset),'hs',@nanmean); plot_hires_coastline(stn.ngdc_hires_bathy); plot(stn.lon,stn.lat,'kp',stn.lon,stn.lat,'wd'); axis(field_bbox(stn.(dataset))); titlename('Full');
help contour_field
fmg; contour_field(stn.(dataset),'hs',@nanmean,0:0.05:1); plot_hires_coastline(stn.ngdc_hires_bathy); plot(stn.lon,stn.lat,'kp',stn.lon,stn.lat,'wd'); axis(field_bbox(stn.(dataset))); titlename('Full');
stn.(dataset).hs = oversample_attenuate_field({stn.ww3.lon,stn.ww3.lat,stn.ww3.hs.date,stn.ww3.hs.field},stn.(dataset),stn.(dataset).depth,[],true);
fmg; contour_field(stn.(dataset),'hs',@nanmean,0:0.05:1); plot_hires_coastline(stn.ngdc_hires_bathy); plot(stn.lon,stn.lat,'kp',stn.lon,stn.lat,'wd'); axis(field_bbox(stn.(dataset))); titlename('Full');
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
quit
%-- 4/1/2016 2:42 PM --%
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
stn = ef_eval_ool_250m('pomf1',[],false,false,true);
fmg; contour_field(stn.(dataset),'hs',@nanmean,0:0.05:1); plot_hires_coastline(stn.ngdc_hires_bathy); plot(stn.lon,stn.lat,'kp',stn.lon,stn.lat,'wd'); axis(field_bbox(stn.(dataset))); titlename('Full');
dataset='CI';
fmg; contour_field(stn.(dataset),'hs',@nanmean,0:0.05:1); plot_hires_coastline(stn.ngdc_hires_bathy); plot(stn.lon,stn.lat,'kp',stn.lon,stn.lat,'wd'); axis(field_bbox(stn.(dataset))); titlename('Full');
numel(find(~isnan(stn.(dataset).hs)))
numel(find(isnan(stn.(dataset).hs)))
numel(find(~isnan(stn.(dataset).hs)))
numel(find(~isnan(stn.(dataset).ef_value)))
numel(find(~isnan(stn.ww3.ef_value)))
numel(find(~isnan(stn.ww3.hs.field)))
numel(find(isnan(stn.ww3.hs.field)))
numel(find(~isnan(stn.ww3.ef_value)))
dbstop ef_eval_ool_250m_run_efs 28
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
dbstop ef_eval_ool_250m_run_efs 28
stn = ef_eval_ool_250m('pomf1',[],false,false,true);
fmg; contour_field(stn.(dataset),'hs',@nanmean,0:0.05:1); plot_hires_coastline(stn.ngdc_hires_bathy); plot(stn.lon,stn.lat,'kp',stn.lon,stn.lat,'wd'); axis(field_bbox(stn.(dataset))); titlename('Full');
fmg; contour_field(newdat.(dataset),'hs',@nanmean,0:0.05:1); plot_hires_coastline(stn.ngdc_hires_bathy); plot(stn.lon,stn.lat,'kp',stn.lon,stn.lat,'wd'); axis(field_bbox(stn.(dataset))); titlename('Full');
fmg; contour_field(stn.(dataset),'hs',@nanmean,0:0.05:1); plot_hires_coastline(stn.ngdc_hires_bathy); plot(stn.lon,stn.lat,'kp',stn.lon,stn.lat,'wd'); axis(field_bbox(stn.(dataset))); titlename('Full');
numel(find(~isnan(stn.(dataset).hs)))
numel(find(~isnan(stn.ww3.ef_value)))
stn.ww3
dbquit
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
stn = ef_eval_ool_250m('pomf1',[],false,false,true);
reviewanim([],0,0,0)
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
stn = ef_eval_ool_250m('pomf1',[],false,false,true);
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
quit
%-- 4/1/2016 9:33 PM --%
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
stn = ef_eval_ool_250m('pomf1',[],false,false,true);
ef_eval_ool_250m_run_efs
close(3:7)
ef_eval_ool_250m_run_efs
close(5:9)
ef_eval_ool_250m_run_efs
close(7:11)
stn=[];
clear all
stn = ef_eval_ool_250m('pomf1',[],false,false,true);
set(gca,'clim',[-1,4]);
set(gca,'clim',[-1,2]);
set(gca,'clim',[-1,1]);
set(gca,'clim',[-.8,.8]);
set(gca,'clim',[-.1,.1]);
set(gca,'clim',[-1,0]);
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
stn = ef_eval_ool_250m('pomf1',[],false,false,true);
axis tight
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
stn = ef_eval_ool_250m('pvgf1',[],false,false,true);
set(gca,'clim',[200,405]);
set(gca,'clim',[300,420]);
set(gca,'clim',[300,405]);
set(gca,'clim',[300,401]);
set(gca,'clim',[300,400]);
set(gca,'clim',[0,400]);
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
stn
stn = ef_eval_ool_250m('pvgf1',[],false,false,true);
stn.(dataset)
dataset='CI';
stn.(dataset)
ans=[]; clear ans
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
stn = ef_eval_ool_250m('pomf1',[],trye,false,true);
stn = ef_eval_ool_250m('pomf1',[],try,false,true);
stn = ef_eval_ool_250m('pomf1',[],trye,false,true);
stn = ef_eval_ool_250m('pomf1',[],try,false,true);
stn = ef_eval_ool_250m('pomf1',[],true,false,true);
set(gca,'clim',[0,1]);
set(gca,'clim',[0,0.1]);
set(gca,'clim',[0,0.05]);
set(gca,'clim',[0,0.02]);
fmg; boxplot_ts(stn.ww3.hs);
fmg; boxplot_ts(stn.ww3.tp);
nanmean(stn.(dataset).ef_value(:))
dataset='CI';
nanmean(stn.(dataset).ef_value(:))
nanmean(stn.(dataset).ef_events.n(:))
stn.(dataset).ef_events
nanmean(stn.(dataset).ef_events.field_n(:))
nanmean(stn.(dataset).data_density(:))
stn.(dataset)
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
stn = ef_eval_ool_250m('vwsp1',[],true,false,true);
stn = ef_eval_ool_250m('vwsp6',[],true,false,true);
reviewanim([],0,0,0)
fmg; boxplot_ts(stn.ww3.hs);
fmg; boxplot_ts(stn.ww3.tp);
dataset='CI';
stn.(dataset).hs = oversample_attenuate_field({stn.ww3.lon,stn.ww3.lat,stn.ww3.hs.date,stn.ww3.hs.field},stn.(dataset),stn.(dataset).depth,[],true);
fmg; contour_field(stn.(dataset),'hs',@nanmean); plot_hires_coastline(stn.ngdc_hires_bathy); axis([min(stn.(dataset).lon),max(stn.(dataset).lon),min(stn.(dataset).lat),max(stn.(dataset).lat)]); set(gca,'clim',[0,1]); titlename('Attenuated');
stn.(dataset).hs = oversample_attenuate_field({stn.ww3.lon,stn.ww3.lat,stn.ww3.hs.date,stn.ww3.hs.field},stn.(dataset),stn.(dataset).depth,[],false);
fmg; contour_field(stn.(dataset),'hs',@nanmean); plot_hires_coastline(stn.ngdc_hires_bathy); axis([min(stn.(dataset).lon),max(stn.(dataset).lon),min(stn.(dataset).lat),max(stn.(dataset).lat)]); set(gca,'clim',[0,1]); titlename('Attenuated (old)');
stn.(dataset).hs = oversample_attenuate_field({stn.ww3.lon,stn.ww3.lat,stn.ww3.hs.date,stn.ww3.hs.field},stn.(dataset),5,stn.(dataset).depth,[],false);
stn.(dataset).hs = oversample_attenuate_field({stn.ww3.lon,stn.ww3.lat,stn.ww3.hs.date,stn.ww3.hs.field},stn.(dataset),stn.(dataset).depth,5,[],false);
stn.(dataset).hs = oversample_attenuate_field({stn.ww3.lon,stn.ww3.lat,stn.ww3.hs.date,stn.ww3.hs.field},5,stn.(dataset),stn.(dataset).depth,[],false);
help oversample_attenuate_field
stn.(dataset).hs = oversample_attenuate_field({stn.ww3.lon,stn.ww3.lat,stn.ww3.hs.date,stn.ww3.hs.field},5,stn.(dataset),stn.(dataset).depth,[],false);
fmg; contour_field(stn.(dataset),'hs',@nanmean); plot_hires_coastline(stn.ngdc_hires_bathy); axis([min(stn.(dataset).lon),max(stn.(dataset).lon),min(stn.(dataset).lat),max(stn.(dataset).lat)]); set(gca,'clim',[0,1]); titlename('Attenuated (5x)');
stn.(dataset).hs = oversample_attenuate_field({stn.ww3.lon,stn.ww3.lat,stn.ww3.hs.date,stn.ww3.hs.field},10,stn.(dataset),stn.(dataset).depth,[],false);
fmg; contour_field(stn.(dataset),'hs',@nanmean); plot_hires_coastline(stn.ngdc_hires_bathy); axis([min(stn.(dataset).lon),max(stn.(dataset).lon),min(stn.(dataset).lat),max(stn.(dataset).lat)]); set(gca,'clim',[0,1]); titlename('Attenuated (10x)');
nanmax(stn.ww3.hs)
nanmax(stn.ww3.hs.field(:))
fmg; hist(stn.ww3.hs.field(:),100);
nansummary(stn.ww3.hs.field(:));
fmg; hist(stn.ww3.hs.field(:),100);
prctile(stn.ww3.hs.field(:),93);
prctile(stn.ww3.hs.field(:),93)
prctile(stn.ww3.hs.field(:),99)
set(gca,'CLim',[0,prctile(stn.ww3.hs.field(:),99)]);
set(gca,'CLim',[0,prctile(stn.ww3.hs.field(:),93)]);
set(gca,'CLim',[0,prctile(stn.ww3.hs.field(:),87)]);
set(gca,'CLim',prctile(stn.ww3.hs.field(:),[7,97]));
set(gca,'CLim',prctile(stn.ww3.hs.field(:),[7,93]));
stn.(dataset).hs = oversample_attenuate_field({stn.ww3.lon,stn.ww3.lat,stn.ww3.hs.date,stn.ww3.hs.field},10,stn.(dataset),stn.(dataset).depth,[],false);
fmg; contour_field(stn.(dataset),'hs',@nanmean); plot_hires_coastline(stn.ngdc_hires_bathy); axis([min(stn.(dataset).lon),max(stn.(dataset).lon),min(stn.(dataset).lat),max(stn.(dataset).lat)]); set(gca,'clim',[0,1]); titlename('Attenuated (10x)');
set(gca,'CLim',prctile(stn.ww3.hs.field(:),[7,93]));
close all
fmg; contour_field(stn.(dataset),'hs',@nanmean,[0:0.2:2]); plot_hires_coastline(stn.ngdc_hires_bathy); axis([min(stn.(dataset).lon),max(stn.(dataset).lon),min(stn.(dataset).lat),max(stn.(dataset).lat)]); set(gca,'clim',[0,1]); titlename('Attenuated (10x)');
stn.(dataset).hs = oversample_attenuate_field({stn.ww3.lon,stn.ww3.lat,stn.ww3.hs.date,stn.ww3.hs.field},5,stn.(dataset),stn.(dataset).depth,[],false);
fmg; contour_field(stn.(dataset),'hs',@nanmean,[0:0.2:2]); plot_hires_coastline(stn.ngdc_hires_bathy); axis([min(stn.(dataset).lon),max(stn.(dataset).lon),min(stn.(dataset).lat),max(stn.(dataset).lat)]); set(gca,'clim',[0,1]); titlename('Attenuated (5x)');
stn.(dataset).hs = oversample_attenuate_field({stn.ww3.lon,stn.ww3.lat,stn.ww3.hs.date,stn.ww3.hs.field},stn.(dataset),stn.(dataset).depth,[],false);
fmg; contour_field(stn.(dataset),'hs',@nanmean,[0:0.2:2]); plot_hires_coastline(stn.ngdc_hires_bathy); axis([min(stn.(dataset).lon),max(stn.(dataset).lon),min(stn.(dataset).lat),max(stn.(dataset).lat)]); set(gca,'clim',[0,1]); titlename('Attenuated (OLD)');
stn.(dataset).hs = oversample_attenuate_field({stn.ww3.lon,stn.ww3.lat,stn.ww3.hs.date,stn.ww3.hs.field},stn.(dataset),stn.(dataset).depth,[],true);
fmg; contour_field(stn.(dataset),'hs',@nanmean,[0:0.2:2]); plot_hires_coastline(stn.ngdc_hires_bathy); axis([min(stn.(dataset).lon),max(stn.(dataset).lon),min(stn.(dataset).lat),max(stn.(dataset).lat)]); set(gca,'clim',[0,1]); titlename('Attenuated');
stn.(dataset).hs = oversample_attenuate_field({stn.ww3.lon,stn.ww3.lat,stn.ww3.hs.date,stn.ww3.hs.field},stn.(dataset),stn.(dataset).depth);
fmg; contour_field(stn.(dataset),'hs',@nanmean,[0:0.2:2]); plot_hires_coastline(stn.ngdc_hires_bathy); axis([min(stn.(dataset).lon),max(stn.(dataset).lon),min(stn.(dataset).lat),max(stn.(dataset).lat)]); set(gca,'clim',[0,1]); titlename('Attenuated');
help strfind
regexprep('depth','[,]*depth[,]*','')
regexprep('depth,,','[,]*depth[,]*','')
regexprep('depth,*linear','[,]*depth[,]*','')
regexprep('*linear,depth','[,]*depth[,]*','')
clear ans
stn.(dataset).hs = oversample_attenuate_field({stn.ww3.lon,stn.ww3.lat,stn.ww3.hs.date,stn.ww3.hs.field},stn.(dataset),stn.(dataset).depth,'depth',true);
dbstop if error
stn.(dataset).hs = oversample_attenuate_field({stn.ww3.lon,stn.ww3.lat,stn.ww3.hs.date,stn.ww3.hs.field},stn.(dataset),stn.(dataset).depth,'depth',true);
size(fldfld),size(depscl)
depscl = flddep ./ 20;
fmg; hist(depscl,100);
depscl = -flddep ./ 20;
depscl(depscl > 1) = 1; depscl(depscl < 0) = 0;
fmg; hist(depscl,100);
depscl = repmat(depscl,[size(fldfld,1),size(depscl)]);
size(fldfld,1)
size(depscl)
size(fldfld,1)*88*88
size(fldfld,2)
size(fldfld,3)
depscl = repmat(depscl,size(fldfld));
size(fldfld)
depscl = repmat(depscl,size(fldfld,1));
size(fldfld,1)
help repmat
depscl(1,:,:) = depscl;
depscl(1,:,:) = depscl(:,:);
depscl(2,:,:) = 0;
size(depscl)
size(fldfld,1)*88*88
size(fldfld,1)*88*88/1e9
size(fldfld,1)*88*88/1e6
prod(size(fldfld)/1e6
prod(size(fldfld))/1e6
prod(size(fldfld))/(2^20)
prod(size(fldfld)*1e3)/(2^30)
prod(size(fldfld)*1e3)/1e9
prod(size(fldfld))*1e3/(2^30)
prod(size(fldfld))*1e3/1e9
help repmat
size(depscl)
depscl(2,:,:) = 0;
size(depscl)
depscl(2,:,:) = depscl;
depscl(2,size(depscl,1),size(depscl,2)) = depscl;
depscl(2,1:size(depscl,1),1:size(depscl,2)) = depscl;
size(depscl)
depscl = -flddep ./ 20;
depscl(depscl < 0) = 0;
depscl(depscl > 1) = 1;
size(depscl)
help repmat
help bsxfun
size(depscl)
size(bsxfun(@times,depscl,fldfld))
help times
help mtimes
help repmat
size( repmat(depscl,[1,1,size(fldfld,1)]) )
size( repmat(depscl,[1,1,size(fldfld,1)])' )
size( permute(repmat(depscl,[1,1,size(fldfld,1)]),[3,1,2]) )
depscl(:,end+1) = depscl(:,end);
size(depscl)
size( permute(repmat(depscl,[1,1,size(fldfld,1)]),[3,1,2]) )
depscl = -flddep ./ 20;
depscl(depscl < 0) = 0;
depscl(depscl > 1) = 1;
depscl = permute(repmat(depscl,[1,1,size(fldfld,1)]),[3,1,2]);
newfld = fldfld .* depscl;
depscl=[]; clear depscl;
size(newfld)
dbquit
stn.(dataset).hs = oversample_attenuate_field({stn.ww3.lon,stn.ww3.lat,stn.ww3.hs.date,stn.ww3.hs.field},stn.(dataset),stn.(dataset).depth,'depth',true);
fmg; contour_field(stn.(dataset),'hs',@nanmean,[0:0.2:2]); plot_hires_coastline(stn.ngdc_hires_bathy); axis([min(stn.(dataset).lon),max(stn.(dataset).lon),min(stn.(dataset).lat),max(stn.(dataset).lat)]); set(gca,'clim',[0,1]); titlename('Attenuated (depth)');
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
stn = ef_eval_ool_250m('pvgf1',[],false,false,true);
fmg; contour_field(stn.(dataset),'hs',@nanmean,[0:0.2:2]); plot_hires_coastline(stn.ngdc_hires_bathy); axis([min(stn.(dataset).lon),max(stn.(dataset).lon),min(stn.(dataset).lat),max(stn.(dataset).lat)]); set(gca,'clim',[0,1]); titlename('Attenuated');
dataset='CI';
fmg; contour_field(stn.(dataset),'hs',@nanmean,[0:0.2:2]); plot_hires_coastline(stn.ngdc_hires_bathy); axis([min(stn.(dataset).lon),max(stn.(dataset).lon),min(stn.(dataset).lat),max(stn.(dataset).lat)]); set(gca,'clim',[0,1]); titlename('Attenuated');
stn.(dataset).hs = oversample_attenuate_field({stn.ww3.lon,stn.ww3.lat,stn.ww3.hs.date,stn.ww3.hs.field},stn.(dataset),stn.(dataset).depth,[],false);
fmg; contour_field(stn.(dataset),'hs',@nanmean,[0:0.2:2]); plot_hires_coastline(stn.ngdc_hires_bathy); axis([min(stn.(dataset).lon),max(stn.(dataset).lon),min(stn.(dataset).lat),max(stn.(dataset).lat)]); set(gca,'clim',[0,1]); titlename('Attenuated (OLD)');
stn.(dataset).hs = oversample_attenuate_field({stn.ww3.lon,stn.ww3.lat,stn.ww3.hs.date,stn.ww3.hs.field},stn.(dataset),stn.(dataset).depth,'depth',true);
fmg; contour_field(stn.(dataset),'hs',@nanmean,[0:0.2:2]); plot_hires_coastline(stn.ngdc_hires_bathy); axis([min(stn.(dataset).lon),max(stn.(dataset).lon),min(stn.(dataset).lat),max(stn.(dataset).lat)]); set(gca,'clim',[0,1]); titlename('Attenuated (depth)');
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
stn = ef_eval_ool_250m('pomf1',[],false,false,true);
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
stn = ef_eval_ool_250m('pomf1',[],false,false,true);
fmg; contour_field(stn.(dataset),'hs',@nanmean,[0:0.2:2]); plot_hires_coastline(stn.ngdc_hires_bathy); axis([min(stn.(dataset).lon),max(stn.(dataset).lon),min(stn.(dataset).lat),max(stn.(dataset).lat)]); set(gca,'clim',[0,1]); titlename('Attenuated');
dataset='CI';
fmg; contour_field(stn.(dataset),'hs',@nanmean,[0:0.2:2]); plot_hires_coastline(stn.ngdc_hires_bathy); axis([min(stn.(dataset).lon),max(stn.(dataset).lon),min(stn.(dataset).lat),max(stn.(dataset).lat)]); set(gca,'clim',[0,1]); titlename('Attenuated');
stn.(dataset).hs = oversample_attenuate_field({stn.ww3.lon,stn.ww3.lat,stn.ww3.hs.date,stn.ww3.hs.field},stn.(dataset),stn.(dataset).depth,[],false);
fmg; contour_field(stn.(dataset),'hs',@nanmean,[0:0.2:2]); plot_hires_coastline(stn.ngdc_hires_bathy); axis([min(stn.(dataset).lon),max(stn.(dataset).lon),min(stn.(dataset).lat),max(stn.(dataset).lat)]); set(gca,'clim',[0,1]); titlename('Attenuated (OLD)');
stn.(dataset).hs = oversample_attenuate_field({stn.ww3.lon,stn.ww3.lat,stn.ww3.hs.date,stn.ww3.hs.field},stn.(dataset),stn.(dataset).depth,'depth',true);
fmg; contour_field(stn.(dataset),'hs',@nanmean,[0:0.2:2]); plot_hires_coastline(stn.ngdc_hires_bathy); axis([min(stn.(dataset).lon),max(stn.(dataset).lon),min(stn.(dataset).lat),max(stn.(dataset).lat)]); set(gca,'clim',[0,1]); titlename('Attenuated (depth)');
stn.(dataset).hs = oversample_attenuate_field({stn.ww3.lon,stn.ww3.lat,stn.ww3.hs.date,stn.ww3.hs.field},stn.(dataset),stn.(dataset).depth,'depth',false);
fmg; contour_field(stn.(dataset),'hs',@nanmean,[0:0.2:2]); plot_hires_coastline(stn.ngdc_hires_bathy); axis([min(stn.(dataset).lon),max(stn.(dataset).lon),min(stn.(dataset).lat),max(stn.(dataset).lat)]); set(gca,'clim',[0,1]); titlename('Attenuated (depth)');
stn.(dataset).hs = oversample_attenuate_field({stn.ww3.lon,stn.ww3.lat,stn.ww3.hs.date,stn.ww3.hs.field},stn.(dataset),stn.(dataset).depth,'depth',true);
fmg; contour_field(stn.(dataset),'hs',@nanmean,[0:0.2:2]); plot_hires_coastline(stn.ngdc_hires_bathy); axis([min(stn.(dataset).lon),max(stn.(dataset).lon),min(stn.(dataset).lat),max(stn.(dataset).lat)]); set(gca,'clim',[0,1]); titlename('Attenuated (depth)');
stn.(dataset).hs = oversample_attenuate_field({stn.ww3.lon,stn.ww3.lat,stn.ww3.hs.date,stn.ww3.hs.field},stn.(dataset),stn.(dataset).depth,'depth',true);
fmg; contour_field(stn.(dataset),'hs',@nanmean,[0:0.2:2]); plot_hires_coastline(stn.ngdc_hires_bathy); axis([min(stn.(dataset).lon),max(stn.(dataset).lon),min(stn.(dataset).lat),max(stn.(dataset).lat)]); set(gca,'clim',[0,1]); titlename('Attenuated (depth)');
fmg; contour_field(stn.(dataset),'hs',@nanmean,[0:0.05:2]); plot_hires_coastline(stn.ngdc_hires_bathy); axis([min(stn.(dataset).lon),max(stn.(dataset).lon),min(stn.(dataset).lat),max(stn.(dataset).lat)]); set(gca,'clim',[0,1]); titlename('Attenuated (depth)');
stn.(dataset).hs = oversample_attenuate_field({stn.ww3.lon,stn.ww3.lat,stn.ww3.hs.date,stn.ww3.hs.field},stn.(dataset),stn.(dataset).depth,'depth',true);
fmg; contour_field(stn.(dataset),'hs',@nanmean,[0:0.05:2]); plot_hires_coastline(stn.ngdc_hires_bathy); axis([min(stn.(dataset).lon),max(stn.(dataset).lon),min(stn.(dataset).lat),max(stn.(dataset).lat)]); set(gca,'clim',[0,1]); titlename('Attenuated (depth)');
x=0:0.05:1; fmg; plot(x,exp(-x/20))
x=0:0.05:20; fmg; plot(x,exp(-x/20))
x=0:0.05:20; fmg; plot(x,ln(x/20))
x=0:0.05:20; fmg; plot(x,log(x/20))
x=0:0.05:20; fmg; plot(x,1/log(x/20))
x=0:0.05:20; fmg; plot(x,1./log(x/20))
x=0:0.05:20; fmg; plot(x,log(-x/20))
x=0:0.05:20; fmg; plot(x,exp(-x/20))
x=0:0.05:20; fmg; plot(x,exp(x/20))
x=0:0.05:20; fmg; plot(x,exp(x/20)-1)
x=0:0.05:20; fmg; plot(x,exp(-x/20))
x=0:0.05:20; fmg; plot(x,exp((20-x)/20))
x=0:0.05:20; fmg; plot(x,exp((x-20)/20))
x=0:0.05:20; fmg; plot(x,exp(1)/exp((x-20)/20))
x=0:0.05:20; fmg; plot(x,exp(1)./exp((x-20)/20))
x=0:0.05:20; fmg; plot(x,exp((x-20)/20)./exp(1))
x=0:0.05:20; fmg; plot(x,exp((x-20)/20))
x=0:0.05:20; fmg; plot(x,cosd((x-20)/pi))
x=0:0.05:20; fmg; plot(x,cosd((20-x)/pi))
x=0:0.05:20; fmg; plot(x,cosd(2*(20-x)/pi))
x=0:0.05:20; fmg; plot(x,cosd((20-x)))
x=0:0.05:20; fmg; plot(x,cosd((20-x)/pi/20))
x=0:0.05:20; fmg; plot(x,cosd(20*(20-x)/pi))
x=0:0.05:20; fmg; plot(x,cosd(2*(20-x)/pi))
x=0:0.05:20; fmg; plot(x,cosd((20-x)/(20*pi)))
x=0:0.05:20; fmg; plot(x,cosd((20-x)/(20/pi)))
x=0:0.05:20; fmg; plot(x,cosd((20-x)/(pi/20)))
x=0:0.05:20; fmg; plot(x,cosd((20-x)/pi))
x=0:0.05:20; fmg; plot(x, cosd( ((20-x)/pi) + (pi/2) ) );
x=0:0.05:20; fmg; plot(x, cosd( ((20-x)/pi) ));
stn.(dataset).hs = oversample_attenuate_field({stn.ww3.lon,stn.ww3.lat,stn.ww3.hs.date,stn.ww3.hs.field},stn.(dataset),stn.(dataset).depth,'depth',true);
fmg; contour_field(stn.(dataset),'hs',@nanmean,[0:0.05:2]); plot_hires_coastline(stn.ngdc_hires_bathy); axis([min(stn.(dataset).lon),max(stn.(dataset).lon),min(stn.(dataset).lat),max(stn.(dataset).lat)]); set(gca,'clim',[0,1]); titlename('Attenuated (depth)');
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
stn = ef_eval_ool_250m('pagp6',[],true,true,true);
stn = ef_eval_ool_250m('ppgp6',[],true,true,true);
stn = ef_eval_ool_250m('pgop6',[],true,true,true);
type get_all_station_metadata
stn = ef_eval_ool_250m('vwsp6',[],true,true,true);
reviewanim([],0,0,0)
fmg; contourf(nanmean(stn.(dataset).field)); colorbar;
dataset='CI';
fmg; contourf(nanmean(stn.(dataset).field)); colorbar;
fmg; contourf(squeeze(nanmean(stn.(dataset).field))); colorbar;
fmg; contourf(squeeze(iqr(stn.(dataset).field))); colorbar;
fmg; contourf(squeeze(nanmean(stn.(dataset).norm_field))); colorbar;
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
quit
%-- 4/4/2016 3:23 PM --%
timenow
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
stn = ef_eval_ool_250m('vwsp6',[],false,false,true);
reviewanim([],0,0,0)
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
stn = ef_eval_ool_250m('pgsp6',[],false,false,true);
fmg; contour_field(stn.(dataset),'hs',@nanmean,[0:0.05:2]); plot_hires_coastline(stn.ngdc_hires_bathy); axis([min(stn.(dataset).lon),max(stn.(dataset).lon),min(stn.(dataset).lat),max(stn.(dataset).lat)]); set(gca,'clim',[0,1]); titlename('Attenuated (depth)');
dataset='CI';
fmg; contour_field(stn.(dataset),'hs',@nanmean,[0:0.05:2]); plot_hires_coastline(stn.ngdc_hires_bathy); axis([min(stn.(dataset).lon),max(stn.(dataset).lon),min(stn.(dataset).lat),max(stn.(dataset).lat)]); set(gca,'clim',[0,1]); titlename('Attenuated (depth)');
fmg; contourf(squeeze(iqr(stn.(dataset).field))); colorbar;
nansummary(stn.(dataset).hs)
nansummary(stn.(dataset).field)
nansummary(stn.(dataset).norm_field)
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
stn = ef_eval_ool_250m('llbp7',[],true,false,true);
reviewanim([],0,0,0)
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
stn = ef_eval_ool_250m('tasp7',[],true,false,true);
reviewanim([],0,0,0)
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
stn = ef_eval_ool_250m('tasp7',[],true,false,true,true);
reviewanim([],0,0,0)
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
stn = ef_eval_ool_250m('ofsp6',[],false,false,true);
reviewanim([],0,0,0)
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
stn = ef_eval_ool_250m('ofsp6',[],false,false,true);
reviewanim([],0,0,0)
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
stn = ef_eval_ool_250m('sasp7',[],true,false,true,true);
reviewanim([],0,0,0)
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
stn = ef_eval_ool_250m('llbp7',[],false,false,true);
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
stn = ef_eval_ool_250m('vesp6',[],false,false,true);
reviewanim([],0,0,0)
set(gca,'CLim',[0,prctile(stn.(dataset).ef_pct(:),93)]);
dataset='CI';
set(gca,'CLim',[0,prctile(stn.(dataset).ef_events.pct(:),93)]);
set(gca,'CLim',[0,nanmax(stn.(dataset).ef_events.pct(:)]);
set(gca,'CLim',[0,nanmax(stn.(dataset).ef_events.pct(:))]);
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
quit
quit
%-- 4/8/2016 6:00 PM --%
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
lonf1 = load_all_ndbc_data('lonf1');
lonf1 = get_station_from_station_name('lonf1'); lonf1 = load_all_ndbc_data(lonf1);
cd tin_020509\
dir
[lon,lat,dat] = asc2mat('tinlidballmos','.','.',[]);
nansummary(lon)
nansummary(dat)
help history
help utmgeoid
help utmzoneui
zn = utmzone
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
dir *.mat
delete tinlidballmos.mat
[lon,lat,dat] = asc2mat('tinlidballmos','.','.',[],'55P');
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
quit
%-- 4/8/2016 7:41 PM --%
po
lonf1 = get_station_from_station_name('lonf1'); lonf1 = load_all_ndbc_data(lonf1);
dir ../*slope*
type find_station_ngdc_offshore_slope.m
more on
type find_station_ngdc_offshore_slope.m
dir ../*slope*
dir ../../../RSMAS/Coastal/thesis/src/*slope*
type station_ngdc_offshore_slope.m
lonf1 = get_station_from_station_name('lonf1'); lonf1 = load_all_ndbc_data(lonf1); lonf1 = station_ngdc_offshore_slope(lonf1);
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
edit ../../../Postdoc/plot_beta_vs_t_season.m
fwyf1 = get_station_from_station_name('fwyf1'); fwyf1 = load_all_ndbc_data(fwyf1); fwyf1 = station_ngdc_offshore_slope(fwyf1);
lonf1 = get_station_from_station_name('lonf1'); lonf1 = load_all_ndbc_data(lonf1); lonf1 = station_ngdc_offshore_slope(lonf1);
mlrf1 = get_station_from_station_name('mlrf1'); mlrf1 = load_all_ndbc_data(mlrf1); mlrf1 = station_ngdc_offshore_slope(mlrf1);
smkf1 = get_station_from_station_name('smkf1'); smkf1 = load_all_ndbc_data(smkf1); smkf1 = station_ngdc_offshore_slope(smkf1);
help boxplot_ts
type boxplot_ts
fwyf1
bets = [fwyf1.ngdc_offshore.slope,lonf1.ngdc_offshore.slope, ...
mlrf1.ngdc_offshore.slope,smkf1.ngdc_offshore.slope,];
for jd=1:365;
jdix = find(get_jday(fwyf1.ndbc_sea_t.date)==jd);
mn(jd,1) = min(fwyf1.ndbc_sea_t.data(jdix));
mn(jd,2) = min(fwyf1.ndbc_sea_t.data(jdix));
mn(jd,3) = min(fwyf1.ndbc_sea_t.data(jdix));
mn(jd,4) = min(fwyf1.ndbc_sea_t.data(jdix));
mx(jd,1) = max(fwyf1.ndbc_sea_t.data(jdix));
mx(jd,2) = max(fwyf1.ndbc_sea_t.data(jdix));
mx(jd,3) = max(fwyf1.ndbc_sea_t.data(jdix));
mx(jd,4) = max(fwyf1.ndbc_sea_t.data(jdix));
end;
[ig,betix] = sort(bets);
fmg;
for ix=betix(:)';
plot(mn(:,ix),bets(ix)+(jd/1e3),'k.');
plot(mx(:,ix),bets(ix)+(jd/1e3),'k.');
end;
bets = [fwyf1.ngdc_offshore_slope,lonf1.ngdc_offshore_slope, ...
mlrf1.ngdc_offshore_slope,smkf1.ngdc_offshore_slope,];
bets
for jd=1:365;
jdix = find(get_jday(fwyf1.ndbc_sea_t.date)==jd);
mn(jd,1) = min(fwyf1.ndbc_sea_t.data(jdix));
mn(jd,2) = min(fwyf1.ndbc_sea_t.data(jdix));
mn(jd,3) = min(fwyf1.ndbc_sea_t.data(jdix));
mn(jd,4) = min(fwyf1.ndbc_sea_t.data(jdix));
mx(jd,1) = max(fwyf1.ndbc_sea_t.data(jdix));
mx(jd,2) = max(fwyf1.ndbc_sea_t.data(jdix));
mx(jd,3) = max(fwyf1.ndbc_sea_t.data(jdix));
mx(jd,4) = max(fwyf1.ndbc_sea_t.data(jdix));
end;
[ig,betix] = sort(bets);
fmg;
for ix=betix(:)';
plot(mn(:,ix),bets(ix)+(jd/1e3),'k.');
plot(mx(:,ix),bets(ix)+(jd/1e3),'k.');
end;
fmg; for ix=betix(:)'; plot(mn(:,ix),bets(ix)+(jd/1e2),'k.'); plot(mx(:,ix),bets(ix)+(jd/1e2),'k.'); end;
jds = 1:365;
for jd=jds(:)';
jdix = find(get_jday(fwyf1.ndbc_sea_t.date)==jd);
mn(jd,1) = min(fwyf1.ndbc_sea_t.data(jdix));
mn(jd,2) = min(fwyf1.ndbc_sea_t.data(jdix));
mn(jd,3) = min(fwyf1.ndbc_sea_t.data(jdix));
mn(jd,4) = min(fwyf1.ndbc_sea_t.data(jdix));
mx(jd,1) = max(fwyf1.ndbc_sea_t.data(jdix));
mx(jd,2) = max(fwyf1.ndbc_sea_t.data(jdix));
mx(jd,3) = max(fwyf1.ndbc_sea_t.data(jdix));
mx(jd,4) = max(fwyf1.ndbc_sea_t.data(jdix));
end;
[ig,betix] = sort(bets);
fmg; for ix=betix(:)'; plot(mn(:,ix),bets(ix)+(jd/1e4),'k.'); plot(mx(:,ix),bets(ix)+(jd/1er),'k.'); end;
fmg; for ix=betix(:)'; plot(mn(:,ix),bets(ix)+(jds./1e4),'k.'); plot(mx(:,ix),bets(ix)+(jds./1e4),'k.'); end;
fmg; for ix=betix(:)'; plot(mn(:,ix),bets(ix)+(jds./5e5),'k.'); plot(mx(:,ix),bets(ix)+(jds./5e5),'k.'); end;
fmg; for ix=betix(:)'; plot(mn(:,ix),bets(ix)+(jds./8e5),'k.'); plot(mx(:,ix),bets(ix)+(jds./8e5),'k.'); end;
fmg; for ix=betix(:)'; plot(mn(:,ix),bets(ix)+(jds./10e5),'k.'); plot(mx(:,ix),bets(ix)+(jds./10e5),'k.'); end;
fmg; for ix=betix(:)'; plot(mn(:,ix),bets(ix)+(jds./5e4),'k.'); plot(mx(:,ix),bets(ix)+(jds./5e4),'k.'); end;
bets = [fwyf1.ngdc_offshore_slope,lonf1.ngdc_offshore_slope, ...
mlrf1.ngdc_offshore_slope,smkf1.ngdc_offshore_slope,];
jds = 1:365;
for jd=jds(:)';
jdix = find(get_jday(fwyf1.ndbc_sea_t.date)==jd);
mn(jd,1) = min(fwyf1.ndbc_sea_t.data(jdix));
mn(jd,2) = min(fwyf1.ndbc_sea_t.data(jdix));
mn(jd,3) = min(fwyf1.ndbc_sea_t.data(jdix));
mn(jd,4) = min(fwyf1.ndbc_sea_t.data(jdix));
mx(jd,1) = max(fwyf1.ndbc_sea_t.data(jdix));
mx(jd,2) = max(fwyf1.ndbc_sea_t.data(jdix));
mx(jd,3) = max(fwyf1.ndbc_sea_t.data(jdix));
mx(jd,4) = max(fwyf1.ndbc_sea_t.data(jdix));
mnd(jd,1) = mean(fwyf1.ndbc_sea_t.data(jdix)) - min(fwyf1.ndbc_sea_t.data(jdix));
mnd(jd,2) = mean(fwyf1.ndbc_sea_t.data(jdix)) - min(fwyf1.ndbc_sea_t.data(jdix));
mnd(jd,3) = mean(fwyf1.ndbc_sea_t.data(jdix)) - min(fwyf1.ndbc_sea_t.data(jdix));
mnd(jd,4) = mean(fwyf1.ndbc_sea_t.data(jdix)) - min(fwyf1.ndbc_sea_t.data(jdix));
mxd(jd,1) = max(fwyf1.ndbc_sea_t.data(jdix)) - mean(fwyf1.ndbc_sea_t.data(jdix));
mxd(jd,2) = max(fwyf1.ndbc_sea_t.data(jdix)) - mean(fwyf1.ndbc_sea_t.data(jdix));
mxd(jd,3) = max(fwyf1.ndbc_sea_t.data(jdix)) - mean(fwyf1.ndbc_sea_t.data(jdix));
mxd(jd,4) = max(fwyf1.ndbc_sea_t.data(jdix)) - mean(fwyf1.ndbc_sea_t.data(jdix));
end;
[ig,betix] = sort(bets);
fmg; for ix=betix(:)'; plot(mnd(:,ix),bets(ix)+(jds./5e4),'k.'); plot(mxd(:,ix),bets(ix)+(jds./5e4),'k.'); end;
bets = [fwyf1.ngdc_offshore_slope,lonf1.ngdc_offshore_slope, ...
mlrf1.ngdc_offshore_slope,smkf1.ngdc_offshore_slope,];
jds = 1:365;
for jd=jds(:)';
jdix = find(get_jday(fwyf1.ndbc_sea_t.date)==jd);
mn(jd,1) = min(fwyf1.ndbc_sea_t.data(jdix));
mn(jd,2) = min(fwyf1.ndbc_sea_t.data(jdix));
mn(jd,3) = min(fwyf1.ndbc_sea_t.data(jdix));
mn(jd,4) = min(fwyf1.ndbc_sea_t.data(jdix));
mx(jd,1) = max(fwyf1.ndbc_sea_t.data(jdix));
mx(jd,2) = max(fwyf1.ndbc_sea_t.data(jdix));
mx(jd,3) = max(fwyf1.ndbc_sea_t.data(jdix));
mx(jd,4) = max(fwyf1.ndbc_sea_t.data(jdix));
mnd(jd,1) = min(fwyf1.ndbc_sea_t.data(jdix)) - mean(fwyf1.ndbc_sea_t.data(jdix));
mnd(jd,2) = min(fwyf1.ndbc_sea_t.data(jdix)) - mean(fwyf1.ndbc_sea_t.data(jdix));
mnd(jd,3) = min(fwyf1.ndbc_sea_t.data(jdix)) - mean(fwyf1.ndbc_sea_t.data(jdix));
mnd(jd,4) = min(fwyf1.ndbc_sea_t.data(jdix)) - mean(fwyf1.ndbc_sea_t.data(jdix));
mxd(jd,1) = max(fwyf1.ndbc_sea_t.data(jdix)) - mean(fwyf1.ndbc_sea_t.data(jdix));
mxd(jd,2) = max(fwyf1.ndbc_sea_t.data(jdix)) - mean(fwyf1.ndbc_sea_t.data(jdix));
mxd(jd,3) = max(fwyf1.ndbc_sea_t.data(jdix)) - mean(fwyf1.ndbc_sea_t.data(jdix));
mxd(jd,4) = max(fwyf1.ndbc_sea_t.data(jdix)) - mean(fwyf1.ndbc_sea_t.data(jdix));
end;
[ig,betix] = sort(bets);
fmg; for ix=betix(:)'; plot(mnd(:,ix),bets(ix)+(jds./5e4),'k.'); plot(mxd(:,ix),bets(ix)+(jds./5e4),'k.'); end;
bets = [fwyf1.ngdc_offshore_slope,lonf1.ngdc_offshore_slope, ...
mlrf1.ngdc_offshore_slope,smkf1.ngdc_offshore_slope,];
jds = 1:365;
for jd=jds(:)';
jdix = find(get_jday(fwyf1.ndbc_sea_t.date)==jd);
mn(jd,1) = min(fwyf1.ndbc_sea_t.data(jdix));
mn(jd,2) = min(fwyf1.ndbc_sea_t.data(jdix));
mn(jd,3) = min(fwyf1.ndbc_sea_t.data(jdix));
mn(jd,4) = min(fwyf1.ndbc_sea_t.data(jdix));
mx(jd,1) = max(fwyf1.ndbc_sea_t.data(jdix));
mx(jd,2) = max(fwyf1.ndbc_sea_t.data(jdix));
mx(jd,3) = max(fwyf1.ndbc_sea_t.data(jdix));
mx(jd,4) = max(fwyf1.ndbc_sea_t.data(jdix));
mnd(jd,1) = min(fwyf1.ndbc_sea_t.data(jdix)) - mean(fwyf1.ndbc_sea_t.data(jdix));
mnd(jd,2) = min(fwyf1.ndbc_sea_t.data(jdix)) - mean(fwyf1.ndbc_sea_t.data(jdix));
mnd(jd,3) = min(fwyf1.ndbc_sea_t.data(jdix)) - mean(fwyf1.ndbc_sea_t.data(jdix));
mnd(jd,4) = min(fwyf1.ndbc_sea_t.data(jdix)) - mean(fwyf1.ndbc_sea_t.data(jdix));
mxd(jd,1) = max(fwyf1.ndbc_sea_t.data(jdix)) - mean(fwyf1.ndbc_sea_t.data(jdix));
mxd(jd,2) = max(fwyf1.ndbc_sea_t.data(jdix)) - mean(fwyf1.ndbc_sea_t.data(jdix));
mxd(jd,3) = max(fwyf1.ndbc_sea_t.data(jdix)) - mean(fwyf1.ndbc_sea_t.data(jdix));
mxd(jd,4) = max(fwyf1.ndbc_sea_t.data(jdix)) - mean(fwyf1.ndbc_sea_t.data(jdix));
p7d(jd,1)  = prctile(fwyf1.ndbc_sea_t.data(jdix), 7) - mean(fwyf1.ndbc_sea_t.data(jdix));
p7d(jd,2)  = prctile(fwyf1.ndbc_sea_t.data(jdix), 7) - mean(fwyf1.ndbc_sea_t.data(jdix));
p7d(jd,3)  = prctile(fwyf1.ndbc_sea_t.data(jdix), 7) - mean(fwyf1.ndbc_sea_t.data(jdix));
p7d(jd,4)  = prctile(fwyf1.ndbc_sea_t.data(jdix), 7) - mean(fwyf1.ndbc_sea_t.data(jdix));
p93d(jd,1)  = prctile(fwyf1.ndbc_sea_t.data(jdix),93) - mean(fwyf1.ndbc_sea_t.data(jdix));
p93d(jd,2)  = prctile(fwyf1.ndbc_sea_t.data(jdix),93) - mean(fwyf1.ndbc_sea_t.data(jdix));
p93d(jd,3)  = prctile(fwyf1.ndbc_sea_t.data(jdix),93) - mean(fwyf1.ndbc_sea_t.data(jdix));
p93d(jd,4)  = prctile(fwyf1.ndbc_sea_t.data(jdix),03) - mean(fwyf1.ndbc_sea_t.data(jdix));
end;
[ig,betix] = sort(bets);
%fmg; for ix=betix(:)'; plot(mnd(:,ix),bets(ix)+(jds./5e4),'k.'); plot(mxd(:,ix),bets(ix)+(jds./5e4),'k.'); end;
fmg; for ix=betix(:)'; plot(p7d(:,ix),bets(ix)+(jds./5e4),'k.'); plot(p93d(:,ix),bets(ix)+(jds./5e4),'k.'); end;
nansummary(p93d(:,1))
nansummary(p93d(:,2))
nansummary(p93d(:,3))
nansummary(p93d(:,4))
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
min([])
plot_beta_vs_t_season
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
plot_beta_vs_t_season
stix
help struct
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
plot_beta_vs_t_season
stns
stn
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
plot_beta_vs_t_season
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
plot_beta_vs_t_season
ix
size(p7d)
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
plot_beta_vs_t_season
stns.station_name
stns(:).station_name
stns{:}.station_name
stns
fmg;
plot(p7d,bets+(jds./5e4),'k.');
for ix=betix(:)';
% plot(p7d(:,ix),bets(ix)+(jds./5e4),'k.');
% plot(p93d(:,ix),bets(ix)+(jds./5e4),'k.');
legs(ix) = stns{stix}.station_name;
end;
legend(legs);
size(p7d),size(bets+(jds./5e4))
help colormap
fmg;
clr = colormap;
for ix=betix(:)';
plot(p7d(:,ix),bets(ix)+(jds./5e4),clr(ix,:));
plot(p93d(:,ix),bets(ix)+(jds./5e4),clr{ix,:});
legs(ix) = stns{stix}.station_name;
end;
legend(legs);
fmg;
clr = colormap;
for ix=betix(:)';
plot(p7d(:,ix),bets(ix)+(jds./5e4),clr(ix,:));
plot(p93d(:,ix),bets(ix)+(jds./5e4),clr{ix,:});
legs(ix) = stns{stix}.station_name;
end;
legend(legs);
fmg;
clr = colormap;
for ix=betix(:)';
plot(p7d(:,ix),bets(ix)+(jds./5e4),'Color',clr(ix,:));
plot(p93d(:,ix),bets(ix)+(jds./5e4),'Color',clr{ix,:});
legs(ix) = stns{stix}.station_name;
end;
legend(legs);
fmg;
clr = colormap;
for ix=betix(:)';
plot(p7d(:,ix),bets(ix)+(jds./5e4),'Color',clr(ix,:));
plot(p93d(:,ix),bets(ix)+(jds./5e4),'Color',clr{ix,:});
end;
fmg;
clr = colormap;
for ix=betix(:)';
plot(p7d(:,ix),bets(ix)+(jds./5e4),'Color',clr(ix,:));
plot(p93d(:,ix),bets(ix)+(jds./5e4),'Color',clr(ix,:));
legs(ix) = stns{stix}.station_name;
end;
legend(legs);
fmg;
clr = colormap; legs = {};
for ix=betix(:)';
plot(p7d(:,ix),bets(ix)+(jds./5e4),'Color',clr(ix,:));
plot(p93d(:,ix),bets(ix)+(jds./5e4),'Color',clr(ix,:));
legs{ix} = stns{stix}.station_name;
end;
legend(legs);
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
plot_beta_vs_t_season
clr
size(clr)
fmg;
clr = colormap; legs = {}; jdoff = 1e4; %jdoff = 5e4;
for ix=betix(:)';
plot(p7d(:,ix),bets(ix)+(jds./jdoff),'Color',clr(ix,:));
plot(p93d(:,ix),bets(ix)+(jds./jdoff),'Color',clr(ix,:));
legs{end+1} = upper(stns{ix}.station_name);
end;
legend(legs);
fmg;
clr = colormap; legs = {}; jdoff = 1e4; jdoff = 5e4;
for ix=betix(:)';
plot(p7d(:,ix),bets(ix)+(jds./jdoff),'Color',clr(ix,:));
plot(p93d(:,ix),bets(ix)+(jds./jdoff),'Color',clr(ix,:));
legs{end+1} = upper(stns{ix}.station_name);
end;
legend(legs);
bets
stns
stns.station_name
for stix=1:numel(stns); disp(stns{stix}.station_name); end;
bets
fmg;
clr = colormap; legs = {}; jdoff = 1e4; jdoff = 5e4;
clr = clr(1:numel(stnms):end);
for ix=betix(:)';
plot(p7d(:,ix),bets(ix)+(jds./jdoff),'Color',clr(ix,:));
plot(p93d(:,ix),bets(ix)+(jds./jdoff),'Color',clr(ix,:));
legs{end+1} = upper(stns{ix}.station_name);
end;
legend(legs);
fmg;
clr = colormap; legs = {}; jdoff = 1e4; jdoff = 5e4;
clr = clr(1:numel(stnms)-1:end);
for ix=betix(:)';
plot(p7d(:,ix),bets(ix)+(jds./jdoff),'Color',clr(ix,:));
plot(p93d(:,ix),bets(ix)+(jds./jdoff),'Color',clr(ix,:));
legs{end+1} = upper(stns{ix}.station_name);
end;
legend(legs);
clr
fmg;
clr = colormap; legs = {}; jdoff = 1e4; jdoff = 5e4;
clr = clr(1:numel(stnms)-1:end,:);
for ix=betix(:)';
plot(p7d(:,ix),bets(ix)+(jds./jdoff),'Color',clr(ix,:));
plot(p93d(:,ix),bets(ix)+(jds./jdoff),'Color',clr(ix,:));
legs{end+1} = upper(stns{ix}.station_name);
end;
legend(legs);
clr
fmg;
clr = colormap; legs = {}; jdoff = 1e4; jdoff = 5e4;
clr = clr(1:numel(stnms):end,:);
for ix=betix(:)';
plot(p7d(:,ix),bets(ix)+(jds./jdoff),'Color',clr(ix,:));
plot(p93d(:,ix),bets(ix)+(jds./jdoff),'Color',clr(ix,:));
legs{end+1} = upper(stns{ix}.station_name);
end;
legend(legs);
clr
fmg;
clr = colormap; legs = {}; jdoff = 3e4; %jdoff = 5e4;
clr = clr(1:numel(stnms):end,:);
for ix=betix(:)';
plot(p7d(:,ix),bets(ix)+(jds./jdoff),'Color',clr(ix,:));
plot(p93d(:,ix),bets(ix)+(jds./jdoff),'Color',clr(ix,:));
legs{end+1} = upper(stns{ix}.station_name);
end;
legend(legs);
fmg;
clr = colormap; legs = {}; jdoff = 2e4; %jdoff = 5e4;
clr = clr(1:numel(stnms):end,:);
for ix=betix(:)';
plot(p7d(:,ix),bets(ix)+(jds./jdoff),'Color',clr(ix,:));
plot(p93d(:,ix),bets(ix)+(jds./jdoff),'Color',clr(ix,:));
legs{end+1} = upper(stns{ix}.station_name);
end;
legend(legs);
[ig,betix] = sort(bets);
fmg;
clr = colormap; legs = {}; jdoff = 2e4; %jdoff = 5e4;
clr = clr(1:numel(stnms):end,:);
for ix=betix(:)';
plot(p7d(:,ix),bets(ix)+(jds./jdoff),'.','Color',clr(ix,:));
plot(p93d(:,ix),bets(ix)+(jds./jdoff),'.','Color',clr(ix,:));
legs{end+1} = upper(stns{ix}.station_name);
end;
legend(legs);
fmg;
clr = colormap; legs = {}; jdoff = 5e5; %jdoff = 5e4;
clr = clr(1:numel(stnms):end,:);
for ix=betix(:)';
plot(p7d(:,ix),bets(ix)+(jds./jdoff),'.','Color',clr(ix,:));
plot(p93d(:,ix),bets(ix)+(jds./jdoff),'.','Color',clr(ix,:));
legs{end+1} = upper(stns{ix}.station_name);
end;
legend(legs);
fmg;
clr = colormap; legs = {}; jdoff = 10e5; %jdoff = 5e4;
clr = clr(1:numel(stnms):end,:);
for ix=betix(:)';
plot(p7d(:,ix),bets(ix)+(jds./jdoff),'.','Color',clr(ix,:));
plot(p93d(:,ix),bets(ix)+(jds./jdoff),'.','Color',clr(ix,:));
legs{end+1} = upper(stns{ix}.station_name);
end;
legend(legs);
fmg;
clr = colormap; legs = {}; jdoff = 3e4; %jdoff = 5e4;
clr = clr(1:numel(stnms):end,:);
for ix=betix(:)';
plot(p7d(:,ix),bets(ix)+(jds./jdoff),'.','Color',clr(ix,:));
plot(p93d(:,ix),bets(ix)+(jds./jdoff),'.','Color',clr(ix,:));
legs{end+1} = upper(stns{ix}.station_name);
end;
legend(legs);
fmg;
clr = colormap; legs = {}; jdoff = 8e4; %jdoff = 5e4;
clr = clr(1:numel(stnms):end,:);
for ix=betix(:)';
plot(p7d(:,ix),bets(ix)+(jds./jdoff),'.','Color',clr(ix,:));
plot(p93d(:,ix),bets(ix)+(jds./jdoff),'.','Color',clr(ix,:));
legs{end+1} = upper(stns{ix}.station_name);
end;
legend(legs);
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
plot_beta_vs_t_season
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
help sort
plot_beta_vs_t_season
[ig,betix] = sort(bets,1,'descend');
betix
[ig,betix] = sort(bets,2,'descend');
betix
[ig,betix] = sort(bets,2,'descend');
fmg;
legs = {}; jdoff = 8e4; %jdoff = 5e4;
clr = {'k','r',[0.0,0.5,0.0],'b','m'};
for ix=betix(:)';
plot(p7d(:,ix),bets(ix)+(jds./jdoff),'.','Color',clr{ix});
plot(p93d(:,ix),bets(ix)+(jds./jdoff),'.','Color',clr{ix});
legs{end+1} = upper(stns{ix}.station_name);
end;
legend(legs);
for stix=1:numel(stns); disp(stns{stix}.station_name); end;
stns{2}
x = intersect_tses(stns{1}.ndbc_sea_t,stns{2}.ndbc_sea_t,stns{3}.ndbc_sea_t,stns{4}.ndbc_sea_t);
x
x{1}
x{2}
x = intersect_tses(stns{1}.ndbc_sea_t,stns{2}.ndbc_sea_t,stns{3}.ndbc_sea_t,stns{4}.ndbc_sea_t,stns{5}.ndbc_sea_t);
x{1}
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
x
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
help intersect_tses
help intersect_ts_data
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
plot_beta_vs_t_season
help upper
titlename(upper(strrep(fld,'_','\_')));
clear all
pack
plot_beta_vs_t_season
clear all; pack
plot_beta_vs_t_season
%-- 4/9/2016 2:09 PM --%
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
plot_beta_vs_t_season
fmg;
legs = {}; jdoff = 10e4; %jdoff = 5e4;
clr = {'k','r',[0.0,0.5,0.0],'b','m'};
for ix=betix(:)';
lhs(end+1) = plot(p7d(:,ix),bets(ix)+(jds./jdoff),'.','Color',clr{ix});
plot(p93d(:,ix),bets(ix)+(jds./jdoff),'.','Color',clr{ix});
annotline([],bets(ix)+(  1./jdoff),'JFM');
annotline([],bets(ix)+( 91./jdoff),'AMJ');
annotline([],bets(ix)+(183./jdoff),'JAS');
annotline([],bets(ix)+(274./jdoff),'OND');
legs{end+1} = upper(stns{ix}.station_name);
end;
legend(lhs,legs,'Location','West');
titlename(upper(strrep(fld,'_','\_')));
fmg;
lhs = []; legs = {};
jdoff = 10e4; %jdoff = 5e4;
clr = {'k','r',[0.0,0.5,0.0],'b','m'};
for ix=betix(:)';
lhs(end+1) = plot(p7d(:,ix),bets(ix)+(jds./jdoff),'.','Color',clr{ix});
plot(p93d(:,ix),bets(ix)+(jds./jdoff),'.','Color',clr{ix});
annotline([],bets(ix)+(  1./jdoff),'JFM');
annotline([],bets(ix)+( 91./jdoff),'AMJ');
annotline([],bets(ix)+(183./jdoff),'JAS');
annotline([],bets(ix)+(274./jdoff),'OND');
legs{end+1} = upper(stns{ix}.station_name);
end;
legend(lhs,legs,'Location','West');
titlename(upper(strrep(fld,'_','\_')));
betix
[ig,betix] = sort(bets,2,'descend');
fmg;
lhs = []; legs = {};
jdoff = 10e4; %jdoff = 5e4;
clr = {'k','r',[0.0,0.5,0.0],'b','m'};
for ix=betix(:)';
lhs(end+1) = plot(p7d(:,ix),bets(ix)+(jds./jdoff),'.','Color',clr{ix});
plot(p93d(:,ix),bets(ix)+(jds./jdoff),'.','Color',clr{ix});
annotline([],bets(ix)+(  1./jdoff),'JFM');
annotline([],bets(ix)+( 91./jdoff),'AMJ');
annotline([],bets(ix)+(183./jdoff),'JAS');
annotline([],bets(ix)+(274./jdoff),'OND');
legs{end+1} = upper(stns{ix}.station_name);
end;
legend(lhs,legs,'Location','West');
titlename(upper(strrep(fld,'_','\_')));
[ig,betix] = sort(bets,2,'descend');
fmg;
lhs = []; legs = {};
jdoff = 10e4; %jdoff = 5e4;
clr = {'k','r',[0.0,0.5,0.0],'b','m'};
for ix=betix(:)';
lhs(end+1) = plot(p7d(:,ix),bets(ix)+(jds./jdoff),'.','Color',clr{ix});
plot(p93d(:,ix),bets(ix)+(jds./jdoff),'.','Color',clr{ix});
% annotline([],bets(ix)+(  1./jdoff),'JFM');
% annotline([],bets(ix)+( 91./jdoff),'AMJ');
% annotline([],bets(ix)+(183./jdoff),'JAS');
% annotline([],bets(ix)+(274./jdoff),'OND');
legs{end+1} = upper(stns{ix}.station_name);
end;
legend(lhs,legs,'Location','West');
titlename(upper(strrep(fld,'_','\_')));
clear all; pack
plot_beta_vs_t_season
dir *ht*m
pd ../../../Postdoc
dir *ht*.m
dir *hc*.m
plot_hc_vs_q0_vs_beta
station_dist('smkf1','sanf1')
clear all; pack
plot_beta_vs_t_season
clear all; pack
plot_beta_vs_t_season
clear all; pack
stns
stnms = { 'fwyf1', 'lonf1', 'mlrf1', 'smkf1', 'sanf1' };
stns = {}; bets = [];
for stix=1:numel(stnms)
stns{stix} = get_station_from_station_name(stnms{stix});
stns{stix} = station_ngdc_offshore_slope(stns{stix});
stns{stix} = load_all_ndbc_data(stns{stix});
stns{stix} = station_heat_flux(stns{stix},'ndbc_wind1_speed','ndbc_air_t',0.10,'ndbc_barom','ndbc_sea_t',[],[],'ndbc');
bets(stix) = stns{stix}.ngdc_offshore_slope;
end;
stns
stns(1)
stns{1}
clear all; pack
plot_beta_vs_t_season
stns{1}
nansummary(stns{1}.ndbc_relumid.data)
nansummary(stns{1}.ndbc_relhumid.data)
dbquit
clear all; pack
help hfbulktc
which hfbulktc.m
help dewp_to_relhumid
help relhumid_to_spechumid
help strrep
which strrep
plot_beta_vs_t_season
stns{1}
stns{1}.ndbc_relhumid
stns{1}.ndbc_spechumid
nansummary(stns{1}.ndbc_spechumid.data)
fmg; hist(stns{1}.ndbc_spechumid.data)
fmg; hist(stns{1}.ndbc_relhumid.data)
nansummary(stns{1}.ndbc_relhumid.data)
dbquit
clear all; pack
dir ../MATLAB/ecoforecasts/*erai*
dir ../MATLAB/ecoforecasts/*ncep*
dir ../MATLAB/ecoforecasts/*era*
dir ../RSMAS/Coastal/thesis/src/*erai*
help get_erai_station
plot_beta_vs_t_season
clear all; pack
plot_beta_vs_t_season
clear all; pack
plot_beta_vs_t_season
stns{1}
grepstruct(stns{1},'net')
clear all; pack
plot_beta_vs_t_season
clear all; pack
plot_beta_vs_t_season
size(jds)
jds = 1:365;
for stix=1:numel(stnms)
for jd=jds(:)';
jdix = find(get_jday(tses{stix}.date)==jd);
mn(jd,stix) = min(tses{stix}.data(jdix));
mx(jd,stix) = max(tses{stix}.data(jdix));
mnd(jd,stix) = min(tses{stix}.data(jdix)) - mean(tses{stix}.data(jdix));
mxd(jd,stix) = max(tses{stix}.data(jdix)) - mean(tses{stix}.data(jdix));
p7d(jd,stix)  = prctile(tses{stix}.data(jdix), 7) - mean(tses{stix}.data(jdix));
p93d(jd,stix)  = prctile(tses{stix}.data(jdix),93) - mean(tses{stix}.data(jdix));
end;
end;
clear all; pack
plot_beta_vs_t_season
tses
tses{1}
stix
tses{stix}
size(jdix)
find_date_ranges(stns{1}.ndbc_net_flux.date)
find_date_ranges(stns{2}.ndbc_net_flux.date)
find_date_ranges(stns{3}.ndbc_net_flux.date)
find_date_ranges(stns{4}.ndbc_net_flux.date)
find_date_ranges(stns{5}.ndbc_net_flux.date)
for stix=1:numel(stnms); disp(stns{stix}.station_name); end;
clear all; pack
plot_beta_vs_t_season
clear all; pack
plot_beta_vs_t_season
xlim([-600,400])
find_date_ranges(tses{1}.date)
find_date_ranges(tses{3}.date)
find_date_ranges(stns{2}.ndbc_net_flux.date)
for stix=1:numel(stnms); disp(stns{stix}.station_name); end;
find_date_ranges(stns{2}.ndbc_net_flux.date)
find_date_ranges(stns{4}.ndbc_net_flux.date)
clear all; pack
plot_beta_vs_t_season
xlim([-600,400])
for stix=1:numel(stnms); disp(stns{stix}.station_name); end;
find_date_ranges(tses{3}.date)
clear all; pack
plot_beta_vs_t_season
fld
print('-dpng',['beta_vs_','ndbc_net_flux','_season.png']);
print('-dpng',['beta_vs_','ndbc_air_t','_jd.png']);
delete(['beta_vs_','ndbc_net_flux','_season.png']);
print('-dpng',['beta_vs_','ndbc_net_flux','_jd.png']);
print('-dpng',['beta_vs_','ndbc_sea_t','_jd.png']);
print('-dpng',['beta_vs_','ndbc_wind1_speed','_jd.png']);
print('-dpng',['beta_vs_','ndbc_air_t','_jd.png']);
print('-dpng',['beta_vs_','ndbc_sea_t','_jd.png']);
print('-dpng',['beta_vs_','ndbc_wind1_speed','_jd.png']);
print('-dpng',['beta_vs_','ndbc_net_flux','_jd.png']);
print('-dpng',['beta_vs_','ndbc_air_t','_jd.png']);
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
quit
%-- 4/9/2016 7:26 PM --%
asc2mat('sai_dball',[],[],[],'55P');
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
quit
%-- 4/9/2016 7:43 PM --%
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
asc2mat('sai_mb_li_db',[],[],[],'55P');
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
quit
%-- 4/9/2016 8:36 PM --%
asc2mat('rot_dbmb_5m',[],[],[],'55P');
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
quit
%-- 4/9/2016 9:00 PM --%
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
asc2mat('tinmblidbmos',[],[],[],'55P');
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
quit
%-- 4/10/2016 1:17 PM --%
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
[lon,lat,dat] = asc2mat('tinmblidbmos');
fmg; contour(lon,lat,dat,[0 0]); colorbar;
fmg; contour(lon,lat,dat,-[0:5,30,50:50:300]); colorbar;
nansummary(dat)
fmg; hist(dat)
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
[lon,lat,dat] = asc2mat('tinlidballmos');
fmg; hist(dat)
nansummary(dat)
fmg; contour(lon,lat,dat,-[0:25]); colorbar;
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
[lon,lat,dat] = asc2mat('Rot_5m_bathymetry');
[lon,lat,dat] = asc2mat('Rota_5m_bathymetry');
nansummary(dat)
fmg; hist(dat)
dp.lon=lon; dp.lat=lat; dp.dat=dat; lon=[]; lat=[]; dat=[]; clear lon lat dat
[sh.lon,sh.lat,sh.dat] = asc2mat('rot_dbmb_5m');
min(sh.lon),min(dp.lon)
max(sh.lon),max(dp.lon)
fmg; hist(dp.dat(:)); fmg; hist(sh.dat(:));
fmg; contour(sh.lon,sh.lat,sh.dat,-[0:30,50]); colorbar; ax=axis;
fmg; contour(dp.lon,dp.lat,dp.dat,-[0:30,50]); colorbar; ax=axis;
axis(ax)
ax=axis
axis(ax)
fmg; contour(sh.lon,sh.lat,sh.dat,-[0:30,50]); colorbar; titlename('rot_dbmb_5m');
axis(ax)
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
[sh.lon,sh.lat,sh.dat] = asc2mat('rot_dbmb_5m');
[dp.lon,do.lat,do.dat] = asc2mat('Rota_5m_bathymetry');
fmg; contour(sh.lon,sh.lat,sh.dat,-[0:10:300]); colorbar; titlename('rot dbmb 5m');
fmg; contour(dp.lon,dp.lat,dp.dat,-[0:10:300]); colorbar; titlename('Rota 5m bathymetry');
dp=[]; do=[]; clear dp do
dp=[]; do=[]; clear dp do ans
[dp.lon,dp.lat,dp.dat] = asc2mat('Rota_5m_bathymetry');
fmg; contour(dp.lon,dp.lat,dp.dat,-[0:10:300]); colorbar; titlename('Rota 5m bathymetry');
ax=axis
axis(ax)
[wi.lon,wi.lat,wi.dat] = asc2mat('Rota_60m');
fmg; contour(wi.lon,wi.lat,wi.dat,-[0:10:300]); colorbar; titlename('Rota 60m');
axis(ax)
fmg; contour(wi.lon,wi.lat,wi.dat,-[0:50:2000]); colorbar; titlename('Rota 60m');
fmg; contour(dp.lon,dp.lat,dp.dat,-[0:50:2000]); colorbar; titlename('Rota 5m bathymetry');
fmg; contour(sh.lon,sh.lat,sh.dat,-[0:50:2000]); colorbar; titlename('rot dbmb 5m');
ax=axis
axis(ax)
nansummary(sh.dat),nansummary(dp.dat),nansummary(wi.dat)
fmg; contour(sh.lon,sh.lat,sh.dat,-[-5:5:80]); colorbar; titlename('rot dbmb 5m');
axis(ax)
field_bbox
field_bbox(sh)
numel(find(sh.dat(:)>0))
numel(find(sh.dat(:)>=0))
numel(find(sh.dat(:)>=-0.8))
111e3/5
360/ans
ans/60
360/ans
111e3/5
360/111e3/5
360/(111e3/5)
60*360/(111e3/5)
111e3/5
360/(111e3/5)
60/(111e3/5)
60(60/(111e3/5))
60/(60/(111e3/5))
60/(111e3/5)
60*60/(111e3/5)
1/3
60*60/(111e3/5)
field_bbox(sh)
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
quit
%-- 4/10/2016 4:15 PM --%
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
quit
%-- 4/11/2016 12:38 PM --%
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
[dp.lon,do.lat,do.dat] = asc2mat('saipan_5m');
fmg; contour(dp.lon,dp.lat,dp.dat,-[0:30,50]); colorbar; titlename('saipan_5m');
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
which do
do
[dp.lon,dp.lat,dp.dat] = asc2mat('saipan_5m');
fmg; contour(dp.lon,dp.lat,dp.dat,-[0:30,50]); colorbar; titlename('saipan_5m');
nansummary(sh.dat),nansummary(dp.dat),nansummary(wi.dat)
[sh.lon,sh.lat,sh.dat] = asc2mat('sai_mb_li_db');
fmg; contour(sh.lon,sh.lat,sh.dat,-[0:30,50]); colorbar; titlename('sai mb li db');
[wi.lon,wi.lat,wi.dat] = asc2mat('cnmi-saipan_c3-ssl-REVERSE');
fmg; contour(wi.lon,wi.lat,wi.dat,-[0:30,50]); colorbar; titlename('cnmi-saipan_c3-ssl-REVERSE');
ax=axis
axis(ax)
fmg; hist(wi.dat(:),100); titlename('cnmi-saipan...'); fmg; hist(dp.dat(:),100); titlename('saipan 5m'); fmg; hist(sh.dat(:),100); titlename('sai mb li db');
min(diff(unique(wi.lon)))
min(diff(unique(dp.lon)))
min(diff(unique(sh.lon)))
fmg; contour(sh.lon,sh.lat,sh.dat,-[0:50:2000]); colorbar; titlename('sai mb li db');
fmg; contour(dp.lon,dp.lat,dp.dat,-[0:50:2000]); colorbar; titlename('saipan 5m');
fmg; contour(wi.lon,wi.lat,wi.dat,-[0:50:2000]); colorbar; titlename('cnmi-saipan...');
ax=axis
axis(ax)
field_bbox(sh)
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
quit
%-- 4/11/2016 2:38 PM --%
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
[sh.lon,sh.lat,sh.dat] = asc2mat('tinmblibdbmos');
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
pd ../MATLAB/ecoforecasts/coast
[sh.lon,sh.lat,sh.dat] = asc2mat('tinmblibdbmos');
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
[sh.lon,sh.lat,sh.dat] = asc2mat('tinmblidbmos');
[dp.lon,dp.lat,dp.dat] = asc2mat('tinlidballmos.mat');
[dp.lon,dp.lat,dp.dat] = asc2mat('tinlidballmos');
fmg; contour(dp.lon,dp.lat,dp.dat,-[0:50:2000]); colorbar; titlename('tinlidballmos');
close all
fmg; contour(dp.lon,dp.lat,dp.dat,-[0:30,50]); colorbar; titlename('tinlidballmos');
fmg; contour(sh.lon,sh.lat,sh.dat,-[0:30,50]); colorbar; titlename('tinmblidbmos');
ax=axis
axis(ax)
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
quit
%-- 4/11/2016 2:54 PM --%
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
asc2mat('tinian_5m',[],[],[],'55P');
fcl
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
quit
%-- 4/11/2016 4:32 PM --%
[sh.lon,sh.lat,sh.dat] = asc2mat('tinmblidbmos');
[dp.lon,dp.lat,dp.dat] = asc2mat('tinlidballmos');
[wi.lon,wi.lat,wi.dat] = asc2mat('tinian_5m');
fmg; contour(sh.lon,sh.lat,sh.dat,-[0:30,50]); colorbar; titlename('tinmblidbmos');
fmg; contour(dp.lon,dp.lat,dp.dat,-[0:30,50]); colorbar; titlename('tinlidballmos');
fmg; contour(wi.lon,wi.lat,wi.dat,-[0:50:2000]); colorbar; titlename('tinian 5m');
ax=axis
axis(ax)
fmg; hist(wi.dat(:),100); titlename('tinian 5m'); fmg; hist(dp.dat(:),100); titlename('tinlidballmos'); fmg; hist(sh.dat(:),100); titlename('tinmblibdbmos');
fmg; hist(sh.dat(sh.dat>-25),100);
fmg; contour(dp.lon,dp.lat,dp.dat,-[0:30,50]); colorbar; titlename('tinlidballmos');
axis(ax)
fmg; contour(sh.lon,sh.lat,sh.dat,-[0:30,50]); colorbar; titlename('tinmblidbmos');
fmg; contour(wi.lon,wi.lat,wi.dat,-[0:30,50]); colorbar; titlename('tinian 5m');
axis(ax)
cl = get(gca,'clim')
set(gca,'CLim',cl);
set(gca,'CLim',[-50 0]);
set(gca,'CLim',[-50 0]);
set(gca,'CLim',[-50 0]);
field_bbox(sh)
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
quit
%-- 4/11/2016 5:02 PM --%
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
asc2mat('Tutuila_5m',[],[],[],'2L');
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
quit
asc2mat('tut_dball',[],[],[],'2L');
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
[lon,lat,dat] = asc2mat('tut_dball',[],[],[],'2L');
fmg; contour(lon,lat,dat,-[0:5:30,50]);
fmg; contour(lon,lat,dat,-[0:5:30,50]); colorbar;
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
quit
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
asc2mat('tut_dbmb',[],[],[],'2L');
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
asc2mat('oo_dbmb_mos3',[],[],[],'2L');
asc2mat('oo_dbmb_mos4',[],[],[],'2L');
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
asc2mat('tau_db_mos_5m',[],[],[],'2L');
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
asc2mat('tau_dbmb_mos',[],[],[],'2L');
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
asc2mat('ofuolo_5m',[],[],[],'2L');
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
asc2mat('Tau_5m',[],[],[],'2L');
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
asc2mat('tut_dbmb',[],[],[],'2L');
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
quit
%-- 4/12/2016 11:07 AM --%
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
pd
plot_beta_vs_t_season
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
plot_beta_vs_t_season
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
lookfor capitalize
lookfor capitalise
which capitalize
rehash
which capitalize
which capitalise
which capital
which cap
which caps
help caps
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
lookfor capital
help interspace
which interspace
help strfind
strfind('foo_bar bif',' _')
capitalize('foo_bar bif')
union([[],[]])
help union
union([],[])
capitalize('foo_bar bif')
capitalize('')
capitalize([])
capitalize(pi)
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
ischar('')
ischar({''})
iscellstr({''})
isempty({})
isempty({''})
iscellstr({'',['ab';'cd']})
iscellstr({'',['ab';'cd'],pi})
numel([1,3])
numel([1,3;4,5])
numel('ab')
ndims('ab')
ismatrix('ab')
isvector('ab')
isvector(['ab'])
isvector(['ab';'cd'])
ischar(['ab';'cd'])
capitalize(pi)
capitalize([])
capitalize('')
capitalize('foo_bar bif')
capitalize({'foo_bar bif'})
capitalize({'foo_bar bif','',pi})
upper({'foo_bar bif','',pi})
capitalize({'foo_bar bif',['ab';'cd'],pi})
ans{2}
help lower
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
plot_beta_vs_t_season
cal(2016,6)
calendar(2016,6)
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
plot_beta_vs_t_season
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
plot_beta_vs_t_season
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
plot_beta_vs_t_season
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
plot_beta_vs_t_season
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
plot_beta_vs_t_season
fmg;
lhs = []; legs = {};
peroff = 10e4; %peroff = 5e4;
peroff = 5e3;
peroff = numel(pers)*3e2;
clr = {'k','r',[0.0,0.5,0.0],'b','m'};
for ix=prmix(:)';
% Left-hand (7th percentile) cluster
lhs(end+1) = plot(p07d(:,ix),prms(ix)+([1:numel(pers)]./peroff),'.','Color',clr{ix});
% Right-hand (93rd percentile) cluster
plot(p93d(:,ix),prms(ix)+([1:numel(pers)]./peroff),'.','Color',clr{ix});
%annotline([],prms(ix)+(  1./peroff),[stracc,'=',num2str(pers(1))]);
legs{end+1} = upper(stns{ix}.station_name);
end;
legend(lhs,legs,'Location','West');
titlename([upper(strrep(fld,'_','\_')),' ',capitalize(strsub),' ',capitalize(stracc)]);
xlabel(strrep(capitalize(fld),'_','\_'));
ylabel(strrep(capitalize(ctlprm),'_','\_'));
%print('-dpng',[ctlprm,'_vs_',fld,'_',strsub,'_',stracc,'.png']);
for ix=prmix(:)';
annotline([],prms(ix)+(  1./peroff),[stracc,'=',num2str(pers(1))]);
end;
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
plot(p93d(:,ix),prms(ix)+([1:numel(pers)]./peroff),'.','Color',clr{ix});
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
plot_beta_vs_t_season
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
plot_beta_vs_t_season
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
plot_beta_vs_t_season
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
pd
asc2mat('tut_dbmb',[],[],[],'2L');
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
[sh.lon,sh.lat,sh.dat] = asc2mat('tut_dbmb');
[dp.lon,dp.lat,dp.dat] = asc2mat('tut_dball');
[wi.lon,wi.lat,wi.dat] = asc2mat('tinian_5m');
wi=[]; clear wi
pack
[wi.lon,wi.lat,wi.dat] = asc2mat('Tutuila_5m');
fmg; contour(sh.lon,sh.lat,sh.dat,-[0:30,50]); colorbar; titlename('tut dbmb (SH)');
fmg; contour(dp.lon,dp.lat,dp.dat,-[0:30,50]); colorbar; titlename('tut dball (DP)');
fmg; contour(wi.lon,wi.lat,wi.dat,-[0:30,50]); colorbar; titlename('Tutuila 5m (WI)');
fmg; hist(wi.dat(:),100); titlename('Tutuila 5m (WI)'); fmg; hist(dp.dat(:),100); titlename('tut dball (DP)'); fmg; hist(sh.dat(:),100); titlename('tut dbmb (SH)');
ax=axis
axis(ax)
set(gca,'CLim',[-50 0]);
nc = mDataset('Tutuila_5m.grd');
nj_info(nc),
gr
gr.lon=cast('double',nc{'x'}(:));
gr.lon=cast(nc{'x'}(:),'double');
gr.lat=cast(nc{'y'}(:),'double');
gr.dat=cast(nc{'z'}(:),'double');
gr
gr.dat=cast(nc{'z'}(:,:),'double');
nj_info(nc),
gr.dat=cast(nc{'z'}(1000:3000,3000:7000),'double');
gr.lat=cast(nc{'y'}(1000:3000),'double');
gr.lon=cast(nc{'x'}(3000:7000),'double');
fmg; contour(gr.lon,gr.lat,gr.dat); colorbar; titlename('GRD');
grd=[]; clear grd
gr.dat=cast(nc{'z'}(500:4000,3000:7000),'double');
gr=[]; clear gr
pack
gr.dat=cast(nc{'z'}(500:4000,3000:7000),'double');
nj_info(nc),
close(nc); clear nc
nc = mDataset('Tutuila_5m.grd');
nj_info(nc),
gr.dat=cast(nc{'z'}(500:4000,3000:7000),'double');
gr
gr.dat=cast(nc{'z'}(500:4000,1000:8000),'double');
gr=[]; clear gr
gr.dat=cast(nc{'z'}(500:4000,1000:8000),'double');
gr.dat=cast(nc{'z'}(:,:),'double');
nj_info(nc),
gr=[]; clear gr
close(nc); clear nc
pack
nc = mDataset('Tutuila_5m.grd');
gr.dat=cast(nc{'z'}(:,:),'double');
gr.dat=cast(nc{'z'}(500:end-500,500:end-500),'double');
gr.lon=cast(nc{'x'}(500:end-500),'double');
gr.lat=cast(nc{'y'}(500:end-500),'double');
close(nc); clear nc
fmg; contour(gr.lon,gr.lat,gr.dat); colorbar; titlename('GRD');
fmg; contour(gr.lon,gr.lat,gr.dat,-[0:30,50]); colorbar; titlename('Tutuila 5m (GRD)');
gr=rmfield(gr,'lon');
gr=rmfield(gr,'lat');
gr.x=cast(nc{'x'}(500:end-500),'double');
nc = mDataset('Tutuila_5m.grd');
gr.x=cast(nc{'x'}(500:end-500),'double');
gr.y=cast(nc{'y'}(500:end-500),'double');
close(nc); clear nc
[gr.lat,gr.lon]=utm2latlon(gr.x,gr.y,'2L');
sh
gr
gr.LON=gr.lon; gr.LAT=gr.lat; gr=rmfield(gr,'lon'); gr=rmfield(gr,'lat');
gr.lon = linspace(min(gr.LON(:)),max(gr.LON(:)),numel(unique(gr.x)));
gr.lat = linspace(min(gr.LAT(:)),max(gr.LAT(:)),numel(unique(gr.y)));
gr.dat = griddata(gr.LON,gr.LAT,gr.dat,gr.lon,gr.lat');
fmg; contour(gr.lon,gr.lat,gr.dat,-[0:30,50]); colorbar; titlename('Tutuila 5m (GRD)');
%-- 4/12/2016 1:39 PM --%
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
[sh.lon,sh.lat,sh.dat] = asc2mat('tut_dbmb');
field_bbox(sh)
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
quit
%-- 4/12/2016 2:33 PM --%
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
[wi.lon,wi.lat,wi.dat] = asc2mat('Tutuila_5m');
nc = mDataset('Tutuila_5m.grd');
gr.dat=cast(nc{'z'}(:,:),'double');
gr
gr.dat=cast(nc{'z'}(500:end-500,500:end-500),'double');
gr.lat=cast(nc{'y'}(500:end-500),'double');
gr.lon=cast(nc{'x'}(500:end-500),'double');
fmg; contour(gr.lon,gr.lat,gr.dat,-[0:30,50]); colorbar; titlename('Tutuila 5m (GRD)');
close(nc); clear nc
fmg; contour(wi.lon,wi.lat,wi.dat,-[0:30,50]); colorbar; titlename('Tutuila 5m (WI)');
axis(wi.lon([500,end-500]),wi.lat([500,end-500]))
axis([wi.lon([500,end-500]),wi.lat([500,end-500])])
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
[wi.lon,wi.lat,wi.dat] = asc2mat('ofuolo_5m');
[sh.lon,sh.lat,sh.dat] = asc2mat('oo_dbmb_mos4');
[dp.lon,dp.lat,dp.dat] = asc2mat('oo_dbmb_mos3');
fmg; contour(wi.lon,wi.lat,wi.dat,-[0:30,50]); colorbar; titlename('ofuolo 5m (WI)');
fmg; contour(sh.lon,sh.lat,sh.dat,-[0:30,50]); colorbar; titlename('oo dbmb mos4 (SH)');
fmg; contour(dp.lon,dp.lat,dp.dat,-[0:30,50]); colorbar; titlename('oo dbmb mos3 (DP)');
set(gca,'CLim',[-50 0]);
ax=axis
axis(ax)
fmg; hist(wi.dat(:),100); titlename('Tutuila 5m (WI)'); fmg; hist(dp.dat(:),100); titlename('tut dball (DP)'); fmg; hist(sh.dat(:),100); titlename('tut dbmb (SH)');
field_bbox(sh)
numel(find(sh.dat>=0))
numel(find(sh.dat>=eps))
numel(find(sh.dat>=-1))
numel(find(sh.dat>=-0.05))
numel(find(sh.dat>=-0.1))
numel(find(sh.dat>=-0.5))
fmg; contour(sh.lon,sh.lat,sh.dat,[0,0]);
fmg; [ig,ch]=contour(sh.lon,sh.lat,sh.dat,[0,0]);
ch
fmg; [ig,ch]=contour(sh.lon,sh.lat,sh.dat,[-1,-1]);
fmg; [ig,ch]=contour(sh.lon,sh.lat,sh.dat,[-0.5,-0.5]);
fmg; [ig,ch]=contour(sh.lon,sh.lat,sh.dat,[-1,-1]);
fmg; [ig,ch]=contour(sh.lon,sh.lat,sh.dat,[-2,-2]);
fmg; [ig,ch]=contour(sh.lon,sh.lat,sh.dat,[-1,-1]);
ch
help contour
fmg; [ig,ch]=contour(sh.lon,sh.lat,sh.dat,[-.5,-.5]);
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
[sh.lon,sh.lat,sh.dat] = asc2mat('tau_dbmb_mos');
[dp.lon,dp.lat,dp.dat] = asc2mat('tau_db_mos_5m');
[wi.lon,wi.lat,wi.dat] = asc2mat('Tau_5m');
fmg; contour(sh.lon,sh.lat,sh.dat,-[0:30,50]); colorbar; titlename('tau dbmb mos (SH)');
fmg; contour(dp.lon,dp.lat,dp.dat,-[0:30,50]); colorbar; titlename('tau db mos 5m (DP)');
fmg; contour(wi.lon,wi.lat,wi.dat,-[0:30,50]); colorbar; titlename('Tau 5m (WI)');
ax=axis
axis(ax)
set(gca,'CLim',[-50 0]);
field_bbox(sh)
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
help ef_eval_ool_250m
stn = ef_eval_ool_250m('tasp7',[],true,false,true,false);
dbstop read_hires_bathymetry
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
dbstop read_hires_bathymetry
stn = ef_eval_ool_250m('tasp7',[],true,false,true,false);
bathfiles
stn_covered
bathfiles
stn_covered
bathfiles
stn_covered
dbquit
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
stn = ef_eval_ool_250m('tasp7',[],true,false,true,false);
ax=axis
axis(ax)
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
stn = ef_eval_ool_250m('tasp7',[],true,true,true,false);
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
po
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
dirs
pd +2
quit
%-- 4/12/2016 9:58 PM --%
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
stn = get_station_from_station_name('tasp7')
stn = get_station_from_station_name('tasp7');
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
quit
%-- 4/13/2016 6:40 PM --%
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
stn = ef_eval_ool_250m('tasp7',[],true,true,true,false);
reviewanim([],0,0,0)
close([2:5,10:18])
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
x=1;
x{1}
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
stn = ef_eval_ool_250m('tasp7',[],false,false,false, false);
stn
x = get_isobath(stn.ngdc_hires_bathy,0);
x(1,1:10)
x(2,1:10)
x = get_isobath(stn.ngdc_hires_bathy,0);
x(1,1:10)
fmg; plot(x);
fmg; plot(x','.');
fmg; plot(x,'.');
nansummary(x(1,:))
nansummary(x(2,:))
fmg; plot(x(1,:),x(2,:),'.');
x = get_isobath(stn.ngdc_hires_bathy,-3);
fmg; plot(x(1,:),x(2,:),'.');
nansummary(x(1,:)),nansummary(x(2,:))
x = get_isobath(stn.ngdc_hires_bathy,-3);
dep
x(:,1:10)
coord(:,1:10)
dbquit
for ix=1:10; ix=ix+1; disp(ix); end;
x=[]; x(:,1)=[];
x = get_isobath(stn.ngdc_hires_bathy,-3);
nansummary(x(1,:)),nansummary(x(2,:))
fmg; plot(x(1,:),x(2,:),'.');
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
quit
%-- 4/14/2016 12:52 PM --%
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
x = get_isobath(stn.ngdc_hires_bathy,0);
stn = ef_eval_ool_250m('tasp7',[],false,false,false, false);
x = get_isobath(stn.ngdc_hires_bathy,0);
fmg; plot(x(1,:),x(2,:),'.');
x0 = get_isobath(stn.ngdc_hires_bathy,0);
clear x x0
x = get_isobath(stn.ngdc_hires_bathy,[0,-1,-2]);
x{2}
isscalar([])
isnumeric([])
y = [1,2,3]; {y(:)}
which mat2cell
mat2cell(y)
help mat2cell
help num2cell
which num2cell
num2cell(y)
ans{3}
num2cell(y)
ans{1}(1)
num2cell(y)
ans{1}
num2cell(y)
class(ans{1})
num2cell(y)
class(ans(1))
clear x y ans
num2cell(pi)
ans{1}
clear x y ans
x = get_isobath(stn.ngdc_hires_bathy,[0,-1,-2]);
x{1}
x{2}
x = get_isobath(stn.ngdc_hires_bathy,0,-2);
clear x y ans
x = get_isobath(stn.ngdc_hires_bathy,0,-2);
fmg; plot(x(1,:),x(2,:),'.');
clear x y ans
x = get_isobath(stn.ngdc_hires_bathy,0,-2);
dbstop if error
x = get_isobath(stn.ngdc_hires_bathy,0,-2);
mindep
deps
dbquit
clear x y ans
isnumeric({1})
isnumeric([1])
isnumeric([1,2])
isnumeric({1})
isnumeric({[1]})
clear x y ans
size(stn.ngdc_hires_bathy.lon([1:10]))
size(stn.ngdc_hires_bathy.lat([1:10]))
stn.ngdc_hires_bathy
x=[1:10];
x(1:3)
ix=1:3;
x(ix)
x(ix(:))
x=[1:10]';
x(ix(:))
x(ix)
x(ix(:))
x(:)
x=[1:10];
x(:)
x(ix(:))'
x=[1:10]';
x(ix(:))'
x(ix(:))
x(ix(:))'
clear x y ans
clear ix x y ans
x = get_isobath(stn.ngdc_hires_bathy,0,-2);
deps
mindep
isscalar(mindep),isnumeric(mindep),isscalar(deps),isnumeric(deps),mindep < deps)
isscalar(mindep),isnumeric(mindep),isscalar(deps),isnumeric(deps),(mindep < deps)
dbquit
x = get_isobath(stn.ngdc_hires_bathy,0,-2);
ix = find(mindep <= bath.field(:) && bath.field(:) <= deps);
ix = find(mindep <= bath.field(:) & bath.field(:) <= deps);
dbquit
x = get_isobath(stn.ngdc_hires_bathy,0,-2);
size(ix)
[ix,jx] = find(mindep <= bath.field & bath.field <= deps);
whos
ix
jx
[ix,jx] = find(mindep <= bath.field & bath.field <= deps);
bath.field(ix(1:5),jx(1:5))
bath.field(jx(1:5),ix(1:5))
bath.field(ix(1:5),jx(1:5))
[ix,jx] = find(-3 <= bath.field & bath.field <= 0);
bath.field(ix(1:5),jx(1:5))
bath.field(jx(1:5),ix(1:5))
[ix,jx] = find(-3 <= bath.field & bath.field <= 0);
bath.field(jx(1:5),ix(1:5))
bath.field(ix(1:5),jx(1:5))
bath.field( sub2ind(ix(1:5),jx(1:5)) )
help sub2ind
bath.field( sub2ind(size(bath.field),ix(1:5),jx(1:5)) )
[ix,jx] = find(-3 <= bath.field & bath.field <= -2);
bath.field( sub2ind(size(bath.field),ix(1:5),jx(1:5)) )
nansummary( bath.field( sub2ind(size(bath.field),ix(1:5),jx(1:5)) ) )
ix = find(-3 <= bath.field & bath.field <= -2);
nansummary( bath.field( ix(1:5) ) )
nansummary( bath.field( ix ) )
fmg; contour(bath.lon,bath.lat,bath.field,-[0:1:5]);
plot(bath.lon(ix),bath.lat(ix),'ro');
[ix,jx] = find(-3 <= bath.field & bath.field <= -2);
clear ix jx
whos
dbup
dbdown
ix = find(-3 <= bath.field & bath.field <= -2);
whos
[ix,jx] = find(-3 <= bath.field & bath.field <= -2);
whos
nansummary(ix)
plot(bath.lon(ix),bath.lat(jx),'ro');
plot(bath.lon(jx),bath.lat(ix),'ro');
fmg; contour(bath.lon,bath.lat,bath.field,-[0:5:20]);
plot(bath.lon(jx),bath.lat(ix),'ro');
[ix,jx] = find(-4 <= bath.field & bath.field <= -1);
whos
plot(bath.lon(jx),bath.lat(ix),'ro');
dbquit
which locus
locus
help rlocus
isnumeric(stn)
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
stn = ef_eval_ool_250m('pvgf1',[],true,true,true,false);
print('-dpng','pvgf1-event-days.png');
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
stn = ef_eval_ool_250m('tasp7',[],false,false,false, false);
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
dirs
pd +2
pd coast
nc = mDataset('usgs_dem_10m_saipan.dods');
help netcdf
ncid = netcdf.open('usgs_dem_10m_saipan.dods');
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
stn = ef_eval_ool_250m('tasp7',[],false,false,false, false);
x = get_isobath(stn.ngdc_hires_bathy,2,-2);
x = get_isobath(stn.ngdc_hires_bathy,0,-2);
x = get_isobath(stn.ngdc_hires_bathy,0,-20);
x = get_isobath(stn.ngdc_hires_bathy,-2,-20);
x = get_isobath(stn.ngdc_hires_bathy,-4,-20);
dx = get_contour_distance(cst,x);
cst = get_isobath(stn.ngdc_hires_bathy,0,-4);
dx = get_contour_distance(cst,x);
dbstop if error
dx = get_contour_distance(cst,x);
locus_or_bath
size(locus_or_bath)
dbquit
size(cst)
[1 2 ; 3 4]
ans(:)'
ans(:)
ans(:)'
cst = get_isobath(stn.ngdc_hires_bathy,0,-4);
dx = get_contour_distance(cst,x);
nDims
nDims2
dbup
size(locus')
size(fromCoords)
dbup
size(fromCoords)
dbup
size(x)
size(cst)
dbquit
x = get_isobath(stn.ngdc_hires_bathy,-4,-20);
size(x)
size(cst)
dx = get_contour_distance(cst,x);
size(fromCoords,1)
dbquit
dx = get_contour_distance(cst,x);
nDims2
nDims
dbup
help knnsearch
dbup
dbdown
idx = knnsearch(locus',fromCoords');
whos
nansummary(idx)
min(diff(unique(idx)))
max(diff(unique(idx)))
dbquit
dx = get_contour_distance(cst,x);
dx = distance_wgs84(lat,lon,nearestCoords(2,:)',nearestCoords(1,:)');
size(lat),size(lon),size(nearestCoords(2,:))
dbup
size(lat),size(lon),size(nearestCoords(2,:))
size(lat),size(lon),size(nearestCoords(1,:))
dbdown
ndims(lat1)
ndims(lat2)
ndims(lon1)
ndims(lat1)
size(lon1),size(lon2)
size(lat1),size(lat2)
dbup
dbdown
dbup
dx = distance_wgs84(lat,lon,nearestCoords(2,:),nearestCoords(1,:));
whos
dbquit
whos
x = get_isobath(stn.ngdc_hires_bathy,-4,-20);
cst = get_isobath(stn.ngdc_hires_bathy,0,-4);
dx = get_contour_distance(cst,x);
fmg; hist(dx(:));
fmg; hist(dx(:),100);
fmg; contour_field(stn.ngdc_hires_bathy);
whos
dx
6+131+131
plot(x(:,dx>=1),'r.');
x(:,1:10)
fmg; contour_field(stn.ngdc_hires_bathy);
plot(x(:,dx>=1)','r.');
fmg; contour_field(stn.ngdc_hires_bathy);
plot(x(1,dx>=1),x(2,dx>=1),'r.');
fmg; hist(dx(:),100);
help contour_field
fmg; contour_field(stn.ngdc_hires_bathy,[],[],@contour);
dbup
contourf
dbup
dbquit
whos
clear ix jx ans
cst = get_isobath(stn.ngdc_hires_bathy,3,-3);
fmg; plot(cst,'.');
fmg; plot(cst','.');
fmg; plot(cst(1,:),cst(2,:),'.');
cst = get_isobath(stn.ngdc_hires_bathy,0,-2);
fmg; plot(cst(1,:),cst(2,:),'.');
cst = get_isobath(stn.ngdc_hires_bathy,10,-4);
fmg; plot(cst(1,:),cst(2,:),'.');
cst = get_isobath(stn.ngdc_hires_bathy,10,-1);
fmg; plot(cst(1,:),cst(2,:),'.');
cst = get_isobath(stn.ngdc_hires_bathy,1,-4);
fmg; plot(cst(1,:),cst(2,:),'.');
cst = get_isobath(stn.ngdc_hires_bathy,0,-4);
fmg; plot(cst(1,:),cst(2,:),'.');
x = get_isobath(stn.ngdc_hires_bathy,-4.01,-30);
dx = get_contour_distance(cst,x);
close all;
fmg; hist(dx(:),100);
fmg; contour_field(stn.ngdc_hires_bathy,[],[],@contour);
dbquit
help isa
ismatrix('abc')
ismatrix(1)
ismatrix([1,1])
ismatrix([1,1;1,1])
isnumeric([])
ismatrix([])
ischar([])
ischar('')
isempty('')
fmg; contour_field(stn.ngdc_hires_bathy,@contour);
fmg; contour_field(stn.ngdc_hires_bathy,@contour,[],[],[0 0]);
fmg; contour_field(stn.ngdc_hires_bathy,@contour,[],[],-[0:50:450]);
fmg; contour_field(stn.ngdc_hires_bathy,@contour);
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
6+131+131
ans*0.55
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
nc = mDataset('http://oos.soest.hawaii.edu/thredds/dodsC/usgs_dem_10m_saipan');
nj_info(nc),
lat = cast(nc{'lat'}(:),'double');
lon = cast(nc{'lon'}(:),'double');
dat = cast(nc{'elev'}(:,:),'double');
close(nc); clear nc
fmg; contour(lon,lat,dat,-[0:50:450]); colorbar;
nansummary(dat)
fmg; contour(lon,lat,dat,[0:50:450]); colorbar;
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
nc = mDataset('http://oos.soest.hawaii.edu/thredds/dodsC/ngdc_bathy_10m_tutuila');
nj_info(nc),
lat = cast(nc{'lat'}(:),'double');
lon = cast(nc{'lon'}(:),'double');
dat = cast(nc{'elev'}(:,:),'double');
close(nc); clear nc
fmg; contour(lon,lat,dat,-[0:50:450]); colorbar;
newdat=round(dat);
fmg; contour(lon,lat,newdat,-[0:50:450]); colorbar;
fmg; contour(lon,lat,newdat,[0 0]); colorbar;
fmg; contour(lon,lat,newdat,[-1 -1]); colorbar;
fmg; contour(lon,lat,newdat,[-1 0]); colorbar;
round(-0.9*2)/2
round(-0.9/2)*2
round(-1/2)*2
round(-1.1/2)*2
newdat=dat; newdat(newdat>-1)=0;
nansummary(dat)
fmg; hist(dat(:),100);
fmg; hist(dat(dat>-80),100);
fmg; contour(lon,lat,dat,-[50:2:70]); colorbar;
fmg; contour(lon,lat,dat,-[30:2:50]); colorbar;
fmg; contour(lon,lat,dat,-[10:2:30]); colorbar;
fmg; contour(lon,lat,dat,-[-10:2:10]); colorbar;
fmg; hist(dat(dat>-80),100);
fmg; hist(dat(dat>-20),100);
min(diff(unique(dat)))
newdat=[]; clear ans newdat
whos('saipan_5m.mat')
whos('-mat','saipan_5m.mat')
help whos
whos('-file','saipan_5m.mat')
url = 'http://oos.soest.hawaii.edu/thredds/dodsC/ngdc_bathy_10m_tutuila';
save('ngdc_bathy_10m_tutuila.mat');
close all
save('ngdc_bathy_10m_tutuila.mat');
delete('ngdc_bathy_10m_tutuila.mat');
save('ngdc_bathy_10m_tutuila.dods.mat');
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
url = 'http://oos.soest.hawaii.edu/thredds/dodsC/pibhmc_bathy_5m_tutuila';
nc = mDataset(url);
nj_info(nc),
lon = cast(nc{'lon'}(:),'double');
lat = cast(nc{'lat'}(:),'double');
dat = cast(nc{'elev'}(:,:),'double');
close(nc); clear nc
clear ix jx ans
save('pibhmc_bathy_5m_tutuila.dods.mat');
nansummary(lon)
lat=[]; lon=[]; clear lat lon
dat=[]; lat=[]; lon=[]; clear dat lat lon ans
nc = mDataset(url);
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
quit
%-- 4/16/2016 12:16 PM --%
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
url = 'http://oos.soest.hawaii.edu/thredds/dodsC/pibhmc_bathy_5m_tutuila';
nc = mDataset(url);
dat = cast(nc{'elev'}(:,:),'double');
lon = cast(nc{'lon'}(:),'double');
lat = cast(nc{'lat'}(:),'double');
close(nc); clear nc
url
save('pibhmc_bathy_5m_tutuila.dods.mat');
fmg; contour(lon,lat,dat,-[-10:2:10]); colorbar;
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
url = 'http://oos.soest.hawaii.edu/thredds/dodsC/ngdc_bathy_90m_amsamoa';
nc = mDataset(url);
lon = cast(nc{'lon'}(:),'double');
lat = cast(nc{'lat'}(:),'double');
dat = cast(nc{'elev'}(:,:),'double');
close(nc); clear nc
save('ngdc_bathy_90m_amsamoa.dods.mat');
fmg; contour(lon,lat,dat,-[-10:2:10]); colorbar;
fmg; contour(lon,lat,dat,-[0:3]); colorbar;
fmg; contour(lon,lat,dat,[0 0]); colorbar;
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
250*6
55*7
55*8
6*10
6*10 + 1
6*10 + 10*2
6*10 + 10*7
1300+1500+300+440+150
1300+1500+300+440+160
300*6
1400+1800+300+440+160
6*10 + 10*7 + 60*0.55
60*0.55
6*10 + 10*7
1400+1800+300+440+163
250*1.15*1.20
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
x(1,1:2,1:2)=randn(2,2)
x(2,1:2,1:2)=randn(2,2)
x
ndims(x)
x(2,:,:)=[];
x
ndims(x)
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
x=randn(2,2)
ndims(x)
x(1,1:2,1:2)=randn(2,2)
x
size(x)
x(2,:,:)
x=[];
x(1,1:2,1:2)=randn(2,2)
ndims(x)
size(x)
ndims(squeeze(x))
x(2,:,:)
x(2,:,:)=[];
x(2,1:2,1:2)=randn(2,2)
size(x)
ndims(squeeze(x))
help subsindex
help subsref
ndims(pi)
pi
ndims([pi])
ndims({pi})
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
x = {1};
x(2:end)
x{2:end}
x(2:end)
ans{:}
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
ix={};
ix{:}
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
nc = mDataset('usgs_dem_10m_saipan.dods');
450*4 + 270*3 + 100*4
close(nc); clear nc
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
url = 'http://oos.soest.hawaii.edu/thredds/dodsC/usgs_dem_10m_saipan';
nc = mDataset(url);
dat = cast(nc{'elev'}(:,:),'double');
lon = cast(nc{'lon'}(:),'double');
lat = cast(nc{'lat'}(:),'double');
close(nc); clear nc
save('usgs_dem_10m_saipan.dods.mat');
fmg; contour(lon,lat,dat,[0 0]); colorbar;
132+131+2+2
ans*0.54
24008+132
23808+132
23951+131
23808+132+131+2+2
40*14
fmg; hist(dat(:),100);
nansummary(dat(:))
min(diff(unique(dat)))
newdat=round(dat);
fmg; contour(lon,lat,newdat,[0 0]); colorbar;
newdat=[]; clear ans newdat
newdat=dat; newdat(-1<newdat & newdat<1)=0;
fmg; contour(lon,lat,newdat,[0 0]); colorbar;
fmg; hist(dat(:),100);
nansummary(dat(:))
fmg; hist(dat(dat>-20),100);
fmg; hist(dat(dat>-20 & dat<20),100);
newdat=[]; clear ans newdat
dem
dem.lon=lon; dem.lat=lat; dem.dat=dat;
dat=[]; lat=[]; lon=[]; clear dat lat lon ans
clear url
url = 'http://oos.soest.hawaii.edu/thredds/dodsC/ngdc_bathy_10m_tutuila';
clear url
url = 'http://oos.soest.hawaii.edu/thredds/dodsC/usgs_dem_10m_saipan';
url = 'http://oos.soest.hawaii.edu/thredds/dodsC/pibhmc_bathy_5m_saipan';
help save
nc = mDataset(url);
dat = cast(nc{'elev'}(:,:),'double');
lon = cast(nc{'lon'}(:),'double');
lat = cast(nc{'lat'}(:),'double');
save('pibhmc_bathy_5m_saipan.plus.usgs_dem_10m_saipan.dods.mat','url','lon','lat','dat','-v7.3');
save('pibhmc_bathy_5m_saipan.dods.mat','url','lon','lat','dat','-v7.3');
bat
bar
bat.lon=lon; bat.lat=lat; bat.dat=dat;
dat=[]; lat=[]; lon=[]; clear dat lat lon ans
close(nc); clear nc
nansummary(bat.lon),nansummary(dep.lon)
nansummary(bat.lon),nansummary(dem.lon)
nansummary(bat.lat),nansummary(dem.lat)
lon = interp1(bat.lon,bat.lon,dem.lon);
dem
lat = interp1(bat.lat,bat.lat,dem.lat);
[lon,lonix] = interp1(bat.lon,bat.lon,dem.lon);
help interp1
lon = interp1(bat.lon,bat.lon,dem.lon,'nearest');
lonix = interp1(bat.lon,1:numel(bat.lon),dem.lon,'nearest');
lat = interp1(bat.lat,bat.lat,dem.lat,'nearest');
latix = interp1(bat.lat,1:numel(bat.lat),dem.lat,'nearest');
dem
bat
dat = interp2(bat.lon,bat.lat,bat.dat,dem.lon,dem.lat,'nearest');
dat = interp2(bat.lon,bat.lat,bat.dat',dem.lon,dem.lat,'nearest');
bat
dem
help interp2
dat = interp2(bat.lon,bat.lat,bat.dat,dem.lon,dem.lat,'nearest');
dat = interp2(bat.lat,bat.lon,bat.dat,dem.lat,dem.lon,'nearest');
dat = interp2(bat.lat,bat.lon,bat.dat',dem.lat,dem.lon,'nearest');
dat = interp2(bat.lon,bat.lat,bat.dat,dem.lon,dem.lat,'nearest');
[LON,LAT]=meshgrid(dem.lon,dem.lat);
dat = interp2(bat.lon,bat.lat,bat.dat,LON,LAT,'nearest');
save('pibhmc_bathy_5m_saipan.plus.usgs_dem_10m_saipan.dods.mat','url','lon','lat','dat','-v7.3');
dat(isnan(dat)) = dem.dat(isnan(dat));
save('pibhmc_bathy_5m_saipan.plus.usgs_dem_10m_saipan.dods.mat','url','lon','lat','dat','-v7.3');
fmg; contour(lon,lat,dat,[0 0]); colorbar;
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
quit
%-- 4/20/2016 4:26 PM --%
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
load('pibhmc_bathy_5m_saipan.plus.usgs_dem_10m_saipan.dods.mat');
fmg; contour(lon,lat,dat,[0 0]); colorbar;
fmg; contour(lon,lat,dat,-[0:3]); colorbar;
fmg; hist(dat(:),100);
fmg; hist(dat(abs(dat)<20),100);
fmg; hist(dat(abs(dat)<5),100);
newdat=roundn(dat,-1);
timenow
fmg; hist(dat(abs(dat)<5),100);
fmg; contour(lon,lat,dat,[0 0]); colorbar;
min(diff(unique(dat)))
min(diff(unique(mewdat)))
min(diff(unique(newdat)))
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
load('usgs_dem_10m_saipan.dods.mat')
fmg; contour(lon,lat,dat,[0 0]); colorbar;
timenow
newdat=roundn(dat,-1);
timenow
fmg; contour(lon,lat,dat,[0 0]); colorbar;
fmg; contour(lon,lat,dat,[1 1]); colorbar;
fmg; contour(lon,lat,dat,[0.1 0.1]); colorbar;
fmg; contour(lon,lat,dat,[0.2 0.2]); colorbar;
fmg; contour(lon,lat,dat,[0.3 0.3]); colorbar;
fmg; contour(lon,lat,dat,[1 1]); colorbar;
fmg; contour(lon,lat,newdat,[0.1 0.1]); colorbar;
fmg; contour(lon,lat,newdat,[0.2 0.2]); colorbar;
fmg; contour(lon,lat,newdat,[0.3 0.3]); colorbar;
fmg; contour(lon,lat,newdat,[0.5 0.5]); colorbar;
fmg; contour(lon,lat,newdat,[1 1]); colorbar;
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
dem = load('usgs_dem_10m_saipan.dods.mat');
bat = load('pibhmc_bathy_5m_saipan.dods.mat');
bat
lon = interp1(bat.lon,bat.lon,dem.lon,'nearest');
lat = interp1(bat.lat,bat.lat,dem.lat,'nearest');
lonix = interp1(bat.lon,1:numel(bat.lon),dem.lon,'nearest');
latix = interp1(bat.lat,1:numel(bat.lat),dem.lat,'nearest');
[LON,LAT]=meshgrid(dem.lon,dem.lat);
dat = interp2(bat.lon,bat.lat,bat.dat,LON,LAT,'nearest');
fmg; contour(lon,lat,dat,[0 0]); colorbar;
dat(isnan(dat)) = dem.dat(isnan(dat));
fmg; contour(lon,lat,dat,[0 0]); colorbar;
nansummary(bat.dat(:))
min(diff(unique(bat.dat(:))))
min(diff(unique(bat.dat(:,1))))
min(diff(unique(bat.dat(1,:))))
min(diff(unique(bat.dat(~isnan(bat.dat)))))
nanmin(diff(unique(bat.dat(:,1))))
nanmin(diff(unique(bat.dat(:,10))))
nanmin(diff(unique(bat.dat(:,end))))
nanmin(diff(unique(bat.dat(:,100))))
nanmin(diff(unique(bat.dat(:,1000))))
bat=[]; clear bat
LON=[]; LAT=[]; clear LON LAT
latix=[]; lonix=[]; clear latix lonix
dem
fmg; contour(bat.lon,bat.lat,bat.dat,[0 0]); colorbar;
bat = load('pibhmc_bathy_5m_saipan.dods.mat');
fmg; contour(bat.lon,bat.lat,bat.dat,[0 0]); colorbar;
fmg; contour(bat.lon,bat.lat,bat.dat,-[0.1 0.1]]); colorbar;
fmg; contour(bat.lon,bat.lat,bat.dat,-[0.1 0.1]); colorbar;
fmg; contour(bat.lon,bat.lat,bat.dat,-[1 1]); colorbar;
nansummary(bat.dat)
timenow
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
url = 'http://oos.soest.hawaii.edu/thredds/dodsC/pibhmc_bathy_5m_saipan';
nc = mDataset(url);
dat = cast(nc{'elev'}(:,:),'double');
close(nc); clear nc
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
quit
%-- 4/20/2016 6:09 PM --%
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
url = 'http://oos.soest.hawaii.edu/thredds/dodsC/pibhmc_bathy_5m_saipan';
dat = cast(nc{'elev'}(:,:),'double');
nc = mDataset(url);
dat = cast(nc{'elev'}(:,:),'double');
timenow
lon = cast(nc{'lon'}(:),'double');
lat = cast(nc{'lat'}(:),'double');
close(nc); clear nc
fmg; contour(lon,lat,dat,[0 0]); colorbar;
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
quit
%-- 4/21/2016 7:53 PM --%
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
stn = get_station_from_station_name('lonf1'); stn = load_ndbc_data(stn);
stn = get_station_from_station_name('lonf1'); stn = load_all_ndbc_data(stn);
lookfor climatol
[acl,asd] = climatologize_time_series(stn.ndbc_air_t.date,stn.ndbc_air_t.data,'Ta');
[scl,ssd] = climatologize_time_series(stn.ndbc_sea_t.date,stn.ndbc_sea_t.data,'Ts');
[acl,asd] = climatologize_time_series(stn.ndbc_air_t.date,stn.ndbc_air_t.data,'Ta');
[scl,ssd] = climatologize_time_series(stn.ndbc_sea_t.date,stn.ndbc_sea_t.data,'Ts');
fmg; plot(1:365,acl,1:365,scl);
[ta,ts] = intersect_tses(stn.ndbc_air_t,stn.ndbc_sea_t);
[acl,asd] = climatologize_time_series(ts.date,ta.data,'Ta');
[scl,ssd] = climatologize_time_series(ts.date,ts.data,'Ts');
fmg; plot(1:365,acl,1:365,scl);
fmg; plot(1:365,acl,'b-',1:365,acl-asd,'b:',1:365,acl+asd,'b:',1:365,scl,'r-',1:365,scl-ssd,'r:',1:365,scl+ssd,'r:');
scatter_fit_ts_seasons(ta,ts,[],[],'Ta','Ts');
corr2(ta.data,ts.data)
jfmix = find(get_season(ta.date)==1);
amjix = find(get_season(ta.date)==2);
jasix = find(get_season(ta.date)==3);
ondix = find(get_season(ta.date)==4);
corr2(ta.data(jfmix),ts.data(jfmix))
corr2(ta.data(amjix),ts.data(amjix))
corr2(ta.data(jasix),ts.data(jasix))
help corr2
help corrcoef
help corr2
tsa.date=ts.date; tsa.data=ts.data-scl(get_jday(ts.date));
tsa
tsa=[]; clear tsa
nansummary(get_jday(ts.date))
tsa.date=ts.date; tsa.data=ts.data-scl(get_jday_no_leap(ts.date));
size(ts.data),size(scl(get_jday_no_leap(ts.date)));
size(ts.data),size(scl(get_jday_no_leap(ts.date)))
tsa.date=ts.date; tsa.data=ts.data-scl(get_jday_no_leap(ts.date))';
taa.date=ta.date; taa.data=ta.data-scl(get_jday_no_leap(ta.date))';
corr2(taa.data,tsa.data)
scatter_fit_ts_seasons(taa,tsa,[],[],'Ta','Ts');
.74^2
edit ../corr_ts_ta.m
for mo=1:12;
c2(mo) = corr2(taa.data(get_month(taa.date)==mo),taa.data(get_month(taa.date)==mo));
end;
c2
for mo=1:12;
c2(mo) = corr2(taa.data(get_month(taa.date)==mo),tsa.data(get_month(tsa.date)==mo));
end;
c2
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
corr_ts_ta
help hist
help bar
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
corr_ts_ta
stn = verify_variable(stn,'ndbc_sea_t_1_d_max_minus_ndbc_sea_t_1_d_min');
stn = verify_variable(stn,'ndbc_sea_t_1_d_max_sub_ndbc_sea_t_1_d_min');
help corr2
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
corr_ts_ta
help months
xlabel(get_monthnames(1:12));
help get_monthnames
timenow
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
corr_all_ts_ta
stns
a = intersect_tses(stns.ndbc_sea_t);
tas
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
corr_all_ts_ta
tas
tass
tass = intersect_tses(tas);
tass = intersect_tses(tas(:));
tass = intersect_tses({tas(:)});
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
corr_all_ts_ta
help corr2
lookfor pearson
tas
tas{1}
tass = intersect_tses(tas);
dbstop if error
tass = intersect_tses(tas);
dbquit
tass = intersect_tses(tas{:});
tass
tass{1}
tass{3}
tsss = intersect_tses(tss{:});
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
corr_all_ts_ta
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
datenum(2016,5,4)-now
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
corr_all_ts_ta
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
corr_all_ts_ta
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
stn = get_station_from_station_name('sanf1'); stn = load_all_ndbc_data(stn);
find_date_ranges(stn.ndbc_sea_t.date)
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
corr_all_ts_ta
find_date_ranges(tsss{1}.date)
73583/12
90*24
30*24
5146/720
4368/720
2005-1992
find_date_ranges(tsss{1}.date,10)
find_date_ranges(stns{1}.ndbc_sea_t.date,10)
find_date_ranges(stns{2}.ndbc_sea_t.date,10)
find_date_ranges(stns{3}.ndbc_sea_t.date,10)
find_date_ranges(stns{4}.ndbc_sea_t.date,10)
find_date_ranges(stns{5}.ndbc_sea_t.date,10)
stns{3}
stns{4}
stns{6}
find_date_ranges(stns{5}.ndbc_sea_t.date,10)
find_date_ranges(stns{3}.ndbc_sea_t.date,10)
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
2005-1992
4368/720
5146/720
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
corr_all_ts_ta
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
corr_all_ts_ta
for yr=1992:2003
for cix=1:numel(cstnms)
for mo=1:12;
moix = find(get_month(taass{cix}.date)==mo & get_year(taass{cix}.date)==yr);
N(mo) = numel(moix);
R(mo) = corr2(taass{cix}.data(moix),tsass{cix}.data(moix));
end;
N,
R,
if ( min(N) < 24 )
fmg;
bar(1:12,R.^2);
axis([0.5,12.5,0,1]);
titlename([upper(cstnms{cix}),' Ta - Ts R^2 ',num2str(yr)]);
end;
end;
end;
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
corr_all_ts_ta
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
corr_all_ts_ta
yrs
corr_all_ts_ta
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
corr_all_ts_ta
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
corr_all_ts_ta
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
corr_all_ts_ta
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
help grp_ts
help ismember
help intersect
help climatologize_time_series
help intersect
help nargout
help sd
help stdev
help stddev
help nanstd
help grp_ts
stn = get_station_from_station_name('sanf1'); stn = load_all_ndbc_data(stn);
[anomts,clim,tid,sdts] = anomalize_ts(stn.ndbc_sea_t,@get_yearhour);
which grpstats
[anomts,clim,tid,sdts] = anomalize_ts(stn.ndbc_sea_t);
dbstop if error
[anomts,clim,tid,sdts] = anomalize_ts(stn.ndbc_sea_t);
size(ix)
[ig,ix,jx] = intersect(tid,cumfun(anomts.date));
whos
help intersect
help ismember
[ig,ix] = intersect(cumfun(anomts.date),tid);
who
whos
clear ig ix jx
whos
[ig,ix,jx] = intersect(cumfun(anomts.date),tid);
whos
size(cumfun(anomts.date))
[ig,ix] = ismember(cumfun(anomts.date),tid);
whos
clear ig ix jx
[ig,ix] = ismember(cumfun(anomts.date),tid);
numel(find(ix==0))
dbquit
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
dbstop if error
[anomts,clim,tid,sdts] = anomalize_ts(stn.ndbc_sea_t);
stn = get_station_from_station_name('sanf1'); stn = load_all_ndbc_data(stn);
[anomts,clim,tid,sdts] = anomalize_ts(stn.ndbc_sea_t);
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
[anomts,clim,tid,asd] = anomalize_ts(stn.ndbc_sea_t);
stn = get_station_from_station_name('sanf1'); stn = load_all_ndbc_data(stn);
[anomts,clim,tid,asd] = anomalize_ts(stn.ndbc_sea_t);
anomts
nansummary(tid)
fmg; plot(tid,cum);
fmg; plot(tid,clim);
[anomts,clim,tid,sdts] = anomalize_ts(stn.ndbc_sea_t,@get_yearhour);
tic,
[anomts,clim,tid,sdts] = anomalize_ts(stn.ndbc_sea_t,@get_yearhour);
toc,
help get_yearhour
help get_jhour
[anomts,clim,tid,sdts] = anomalize_ts(stn.ndbc_sea_t,@get_jhour);
tic,
[anomts,clim,tid,sdts] = anomalize_ts(stn.ndbc_sea_t,@get_jhour);
toc
24*365
whos
max(diff(tid))
median(diff(tid))
tid(1:10)
nansummary(clim)
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
stn = get_station_from_station_name('mlrf1'); stn = load_all_ndbc_data(stn);
[anomts,clim,tid,sdts] = anomalize_ts(stn.ndbc_sea_t,@get_jhour);
numel(unique(get_jhour(stn.ndbc_sea_t.date)))
numel(unique(get_jhour([1:8800]/24)))
find_date_ranges(stn.ndbc_sea_t.date,1)
type get_jhour
nansummary(get_jhour(stn.ndbc_sea_t.date))
24*365
numel(unique(get_jhour(stn.ndbc_sea_t.date)))
nansummary(get_jhour(stn.ndbc_sea_t.date))
min(diff(unique(get_jhour(stn.ndbc_sea_t.date))))
max(diff(unique(get_jhour(stn.ndbc_sea_t.date))))
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
quit
%-- 5/1/2016 7:37 PM --%
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
stn = get_station_from_station_name('mlrf1'); stn = load_all_ndbc_data(stn);
numel(unique(get_jhour(stn.ndbc_sea_t.date)))
nansummary(get_jhour(stn.ndbc_sea_t.date))
median(diff(unique(get_jhour(stn.ndbc_sea_t.date))))
numel(unique(get_jhour(stn.ndbc_sea_t.date)))
24*365
numel(unique(get_jhour_no_leap(stn.ndbc_sea_t.date)))
numel(get_jhour(stn.ndbc_sea_t.date))
numel(unique(get_jhour_no_leap(stn.ndbc_sea_t.date)))
numel(unique(get_jhour(stn.ndbc_sea_t.date)))
nansummary(unique(get_jhour(stn.ndbc_sea_t.date)))
nansummary(unique(get_jhour_no_leap(stn.ndbc_sea_t.date)))
[anomts,clim,tid,sdts] = anomalize_ts(stn.ndbc_sea_t,@get_jhour_no_leap);
whos
24*365
fmg; plot(tid,clim);
[anomts,clim,tid,sdts] = anomalize_ts(stn.ndbc_sea_t,@get_jday_no_leap);
fmg; plot(tid,clim);
(239*3)+(199*3)
(190.40*3)+(599)
round( (190.40*3)+(599) )
1366
(239*3)+(199*2)
timenow
fmg; plot(tid,clim);
ylim([22,31])
help grp_ts
[anomts,clim,tid,sdts] = anomalize_ts(stn.ndbc_sea_t,@get_jday_no_leap,@nanmean);
fmg; plot(tid,clim); titlename('Mean');
ylim([22,31])
[anomts,clim,tid,sdts] = anomalize_ts(stn.ndbc_sea_t,@get_jhour_no_leap,@nanmean);
fmg; plot(tid,clim); titlename('Mean');
ylim([22,31])
axis([0,366,22,31])
help datenum
help datetick
datetick('x','mmm')
datetick3('x','mmm')
axis([0,365,22,31])
axis([0,366,22,31])
datetick3('x','mmm')
axis([0,366,22,31])
datetick3('x','mmm')
help mean
help mode
help nanmode
[anomts,clim,tid,sdts] = anomalize_ts(stn.ndbc_sea_t,@get_jday_no_leap,@nanmode);
[anomts,clim,tid,sdts] = anomalize_ts(stn.ndbc_sea_t,@get_jday_no_leap,@mode);
fmg; plot(tid,clim); titlename('Mode');
axis([0,366,22,31])
datetick3('x','mmm')
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
stn = get_station_from_station_name('mlrf1'); stn = load_all_ndbc_data(stn);
[anomts,clim,tid,asd] = anomalize_ts(stn.ndbc_sea_t,@get_jhour_no_leap);
fmg; plot(tid,clim,'k-',tid,clim-asd,'k:',tid,clim+asd,'k:');
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
corr_all_ts_ta
reviewanim([],0,0,0)
axis([0,366,13,34])
axis([0,366,13,34]); datetick3('x','mmm');
close all
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
corr_all_ts_ta
reviewanim([],0,0,0)
close all
clear ans cix mo moix N R stnm
corr_all_ts_ta
reviewanim([],0,0,0)
close all
clear ans cix mo moix N R stnm
corr_all_ts_ta
clear ans cix mo moix N R stnm
close all
corr_all_ts_ta
help corrcoef
close all
clear ans cix mo moix N R stnm
corr_all_ts_ta
size(r),size(p)
r
p
R
dbquit
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
help bar
help bar3h
help histogram2
help histogram
help bar
corr_all_ts_ta
reviewanim([],0,0,0)
close all
help text
corr_all_ts_ta
num2str(P)
num2str(P')
dbquit
close all
clear ans cix mo moix N R stnm
clear ans cix mo moix n p r N P R stnm yr yrs
corr_all_ts_ta
close all
corr_all_ts_ta
close all
corr_all_ts_ta
close all
corr_all_ts_ta
help dunn
lookfor dunn
close all
help boxplot
help kruskalwallis
help multcompare
help anova1
help multcompare
help get_erai_station
help correctairspeed
which corrairspeed
which correctairspeed
help aero
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
corr_all_ts_ta
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
stn = get_station_from_station_name('mlrf1'); stn = load_all_ndbc_data(stn);
help get_ncep_narr
help get_ncep_station
stn = get_ncep_station(stn,'narr');
stn
stn = get_erai_station(stn);
stn
help strfind
help strmatch
help findstr
lookfor correct
lookfor adjust
help adjust_erai_station
stn = adjust_erai_station(stn);
stn
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
corr_all_ts_ta
clear all; pack
reviewanim([],0,0,0)
corr_all_ts_ta
clear all; pack
reviewanim([],0,0,0)
corr_all_ts_ta
reviewanim([],0,0,0)
figsnamed('LONF')
reviewanim(figsnamed('LONF'))
reviewanim(figsnamed('LONF'),0,0,0)
help figsnamed
reviewanim(figsnamed('LONF.*R^2'),0,0,0)
stns{1}
reviewanim(figsnamed('LONF1.*R^2'),0,0,0)
figsnamed('LONF1.*R^2')
figsnamed('LONF.* R')
reviewanim(figsnamed('LONF1.*R\^2'),0,0,0)
reviewanim(figsnamed('MLRF1.*R\^2'),0,0,0)
clear all; pack
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
cX = {'ERAI qa', 'erai_spechumid'};
size(cX)
cX{1,1}
cX{1,2}
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
cX = { ...
'Ta', 'ndbc_air_t' ;...
'ERAI Ta', 'erai_air_t' ;...
'NARR Ta', 'ncep_air_t' ;...
'ERAI qa', 'erai_spechumid' ;...
'NARR qa', 'ncep_spechumid' ;...
'ERAI QswI', 'erai_dsrf' ;...
'NARR QswI', 'ncep_dsrf' ;...
};
cX{1,2}
corr_all_ts_ta
help multcompare
reviewanim(figsnamed('MLRF1.*R\^2'),0,0,0)
reviewanim(figsnamed('NARR.*R\^2'),0,0,0)
reviewanim(figsnamed('R\^2'))
stns{1}
grepstruct(stns{1},'flux')
nansummary(stns{1}.ncep_latent_heat_flux.data+stns{1}.ncep_sensible_heat_flux.data-stns{1}.ncep_net_heat_flux.data);
nansummary(stns{1}.ncep_latent_heat_flux.data),nansummary(stns{1}.ncep_sensible_heat_flux.data),nansummary(stns{1}.ncep_net_heat_flux.data),
sw_dens0(36,25)
sw_dens0(36,25)/sw_dens0(36,15)
sw_dens0(36,30)/sw_dens0(36,15)
1/ans
sw_cp(36,15)
help sw_cp
sw_cp(36,15,0)
sw_cp(36,35,0)
sw_cp(36,35,0)/sw_cp(36,15,0)
1/ans
help sw_dens0
help air_dens
sw_cp(36,20)*sw_dens0(36,20)*2
sw_cp(36,20,2)*sw_dens0(36,20,2)*2
sw_cp(36,20,2)*sw_dens(36,20,2)*2
nansummary(stns{1}.ncep_latent_heat_flux.data),nansummary(stns{1}.ncep_sensible_heat_flux.data),nansummary(stns{1}.ncep_net_heat_flux.data),
830/8.18e6
3600*830/8.18e6
-1356*3600/8.18e6
830*3600/8.18e6
nansummary(stns{1}.ncep_latent_heat_flux.data),nansummary(stns{1}.ncep_sensible_heat_flux.data),nansummary(stns{1}.ncep_net_heat_flux.data),
nansummary(stns{1}.erai_latent_heat_flux.data),nansummary(stns{1}.erai_sensible_heat_flux.data),nansummary(stns{1}.erai_net_heat_flux.data),
grepstruct(stns{1},'flux')
nansummary(stns{1}.erai_latent_heat_flux.data),nansummary(stns{1}.erai_sensible_heat_flux.data),nansummary(stns{1}.erai_srf.data),nansummary(stns{1}.erai_net_heat_flux.data),
nansummary(stns{1}.ncep_latent_heat_flux.data),nansummary(stns{1}.ncep_sensible_heat_flux.data),nansummary(stns{1}.ncep_srf.data),nansummary(stns{1}.ncep_net_heat_flux.data),
scatter_fit_ts(stns{1}.ncep_net_heat_flux,stns{1}.erai_net_heat_flux)
scatter_fit_ts(stns{1}.ncep_net_heat_flux,stns{1}.raw_erai_net_heat_flux)
scatter_fit_ts(stns{1}.ncep_net_heat_flux,stns{1}.erai_net_heat_flux)
scatter_fit_ts(stns{1}.ncep_latent_heat_flux,stns{1}.erai_latent_heat_flux)
scatter_fit_ts(stns{1}.ncep_sensible_heat_flux,stns{1}.erai_sensible_heat_flux)
scatter_fit_ts(stns{1}.ncep_air_t,stns{1}.erai_air_t)
scatter_fit_ts(stns{1}.ncep_spechumid,stns{1}.erai_spechumid)
scatter_fit_ts(stns{1}.ncep_air_t,stns{1}.erai_air_t)
scatter_fit_ts_seasons(stns{1}.ncep_air_t,stns{1}.erai_air_t,[],[],[],[],[],[],'resid')
help scatter_fit_ts_seasons
scatter_fit_ts_seasons(stns{1}.ncep_air_t,stns{1}.erai_air_t,[],[],[],[],[],'resid')
scatter_fit_ts_seasons(stns{1}.ncep_air_t,stns{1}.raw_erai_air_t,[],[],[],[],[],'resid',true)
subplots_set('xlim',[5,33],'ylim',[5,33]);
scatter_fit_ts_seasons(stns{1}.ncep_air_t,stns{1}.raw_erai_air_t,[],[],'NARR','ERAI',[],'resid',true)
subplots_set('xlim',[5,33],'ylim',[5,33]);
scatter_fit_ts_seasons(stns{1}.ndbc_air_t,stns{1}.raw_erai_air_t,[],[],'NDBC','ERAI',[],'resid',true)
subplots_set('xlim',[5,33],'ylim',[5,33]);
reviewanim([28:33],0,0,0)
reviewanim([29:33],0,0,0)
reviewanim([34:39],0,0,0)
scatter_fit_ts_seasons(stns{1}.ndbc_air_t,stns{1}.ncep_air_t,[],[],'NDBC','NARR',[],'resid',true)
subplots_set('xlim',[5,33],'ylim',[5,33]);
nansummary(stns{1}.erai_latent_heat_flux.data),nansummary(stns{1}.erai_sensible_heat_flux.data),nansummary(stns{1}.erai_srf.data),nansummary(stns{1}.erai_net_heat_flux.data),
nansummary(stns{1}.ncep_latent_heat_flux.data),nansummary(stns{1}.ncep_sensible_heat_flux.data),nansummary(stns{1}.ncep_srf.data),nansummary(stns{1}.ncep_net_heat_flux.data),
scatter_fit_ts_seasons(stns{1}.ncep_dsrf,stns{1}.erai_dsrf,[],[],'NARR','ERAI',[],'resid',true)
close([29:49])
reviewanim([],0,0,0)
figsnamed('')
figsnamed('.')
figure(50)
timenow
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
corr_all_ts_ta
moix
cix
tass{cix}
taas{cix}
tsas{cix}
help intersect_tses
cix
ta{cix}
tas{cix}
ts{cix}
ta{cix}
ts{cix}
ta{2}
ts{2}
ta{3}
ts{3}
ta{4}
help intersect_tses
help intersect_all_dates
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
corr_all_ts_ta
reviewanim([],0,0,0)
close all
clear ans cix mo moix n p r N P R stnm yr yrs
clear ans cix cX cY mo moix n p r N P R stnm yr yrs
clear ans cix cX cY mo moix n p r N P R stnm xix yix X Y yr yrs
clear ans cix cX cY mo moix n p r N P R stnm xix yix X Y Xs Yx yr yrs
clear ans cix cX cY mo moix n p r N P R stnm xix yix X Y Xs Ys yr yrs
959/6
clear cstnms ans cix cX cY mo moix n p r N P R stnm xix yix X Y Xs Ys yr yrs
190*6
corr_all_ts_ta
reviewanim([],0,0,0)
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
help plotyy
corr_all_ts_ta
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
tic,
corr
toc,
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
tic,
corr_all_ts_ta
toc,
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
corr_all_ts_ta
reviewanim(figsnamed('MLRF1.*R\^2'))
reviewanim(figsnamed('MLRF1.*R\^2'),0,0,0)
reviewanim(figsnamed('LONF1.*R\^2'),0,0,0)
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
corr_all_ts_ta
reviewanim(figsnamed('LONF1.*R\^2'),0,0,0)
reviewanim([],0,0,0)
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
corr_all_ts_ta
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
quit
%-- 5/4/2016 1:54 PM --%
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
interp_dem_bathy
size(bat)
size(dat)
bat
dem
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
interp_dem_bathy
fmg; contour(lon,lat,dat,[0 0]); colorbar;
bat=[]; dem=[]; dat=[]; clear all
pack
interp_dem_bathy
fmg; contour(lon,lat,dat,[0 0]); colorbar;
help interp2
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
interp_dem_bathy
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
quit
%-- 5/6/2016 11:55 AM --%
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
interp_dem_bathy
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
interp_dem_bathy
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
interp_dem_bathy
dbquit
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
interp_dem_bathy
size(dat),bat
fmg; contour(lon,lat,dat,-20:5:10); colorbar;
fmg; contour(lon,lat,dat,0:5:10); colorbar;
%-- 5/6/2016 1:48 PM --%
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
c = load('usgs_dem_10m_saipan.dods.mat');
e
e = load('../usgs_dem_10m_saipan.dods.mat');
c
e
all(c.dat==e.dat)
all(c.dat(:)==e.dat(:))
all(c.dat(:)==e.dat(:)|isnan(c.dat(:)))
all(c.dat(:)==e.dat(:)|isnan(e.dat(:)))
all(c.lon==e.lon)
all(c.lat==e.lat)
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
interp_dem_bathy
bat=[]; dem=[]; clear bat dem
pack
numel(dat(:)==0)
numel(dat(:)==-1)
numel(dat(:)==1)
numel(find(dat(:)==0))
numel(find(dat(:)==1))
numel(find(dat(:)==-1))
numel(find( roundn(dat(:),-1) == 0 ))
numel(find( roundn(dat(:),-1) == 1 ))
numel(find( roundn(dat(:),-1) == -1 ))
numel(find( round(dat(:)) == 0 ))
numel(find( round(dat(:)) == 1 ))
numel(find( round(dat(:)) == -1 ))
for z=-1:0.1:1; disp([z, numel(find( roundn(dat(:),-1) == z ))]); end;
for z=-1:0.1:1; disp({z, numel(find( roundn(dat(:),-1) == z ))}); end;
for z=-2:0.1:2; disp({z, numel(find( roundn(dat(:),-1) == z ))}); end;
for z=-10:1:10; disp({z, numel(find( roundn(dat(:),-1) == z ))}); end;
help contourc
newdat = interp2(lon,lat,dat,lon,lat,'linear',nan);
[LON,LAT] = meshgrid(lon,lat);
newdat = interp2(lon,lat,dat,LON,LAT,'linear',nan);
LON=[]; LAT=[]; clear LON LAT
all(newdat(:)==dat(:)|isnan(newdat))
all(newdat(:)==dat(:)|isnan(newdat(:)))
newdat=[]; clear newdat
min(diff(unique(lon)))
max(diff(unique(lon)))
min(diff(unique(lat)))
max(diff(unique(lat)))
for z=-10:1:10; disp({z, numel(find( roundn(dat(:),-1) == z ))}); end;
for z=-10:1:10; disp({z, numel(find( roundn(dat(:),0) == z ))}); end;
fmg; contour(lon,lat,dat,[0 0]); colorbar;
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
interp_dem_bathy
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
quit
%-- 5/6/2016 7:50 PM --%
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
interp_dem_bathy
newdat = roundn(dat,-1);
fmg; contour(lon,lat,newdat,[0 0]); colorbar;
ix=find(lon==15.75)
[ig,ix]=min(lat-15.75)
[ig,ix]=min(lat-15.5)
[ig,ix]=min(lat-14.5)
[ig,ix]=min(abs(lat-14.5))
[ig,ix]=min(abs(lat-15.75))
nansummary(lat)
[ig,ix]=min(abs(lat-15.25))
nansummary(lon)
[ig,jx]=min(abs(lon-145.73))
fmg; contour(lon(jx-500:jx+500),lat(ix-500:ix+500),dat(ix-500:ix+500,jx-500:jx+500),-1:1);
fmg; contour(lon(jx-500:jx+500),lat(ix-500:ix+500),dat(ix-500:ix+500,jx-500:jx+500),-1:1); colorbar;
fmg; contour(lon(jx-500:jx+500),lat(ix-500:ix+500),dat(ix-500:ix+500,jx-500:jx+500),[0 0]); colorbar;
[ig,jx]=min(abs(lon-145.75))
fmg; contour(lon(jx-500:jx+500),lat(ix-500:ix+500),dat(ix-500:ix+500,jx-500:jx+500),[0 0]); colorbar;
fmg; contour(lon(jx-500:jx+500),lat(ix-500:ix+500),dat(ix-500:ix+500,jx-500:jx+500),[-1 0]); colorbar;
fmg; contour(lon(jx-500:jx+500),lat(ix-500:ix+500),dat(ix-500:ix+500,jx-500:jx+500),[-1 0 1]); colorbar;
[LON,LAT]=meshgrid(lon,lat);
newdat = interp2(dem.lon,dem.lat,dem.dat,LON,LAT,'spline',nan);
dem = load(dempath);
dem.dat(dem.dat<eps) = nan;
bat = load(batpath);
bat.dat(bat.dat>eps) = nan;
newdat = interp2(dem.lon,dem.lat,dem.dat,LON,LAT,'spline',nan);
help interp2
newdat = interp2(dem.lon,dem.lat,dem.dat,LON,LAT,'cubic',nan);
fmg; contour(lon(jx-500:jx+500),lat(ix-500:ix+500),newdat(ix-500:ix+500,jx-500:jx+500),[-1 0 1]); colorbar;
fmg; contour(lon(jx-500:jx+500),lat(ix-500:ix+500),dat(ix-500:ix+500,jx-500:jx+500),[-1 0 1]); colorbar;
close all
newdat = interp2(lon,lat,dat,LON,LAT,'spline',nan);
newdat = interp2(lon,lat,dat,LON,LAT,'cubic',nan);
fmg; contour(lon(jx-500:jx+500),lat(ix-500:ix+500),dat(ix-500:ix+500,jx-500:jx+500),[-1 0 1]); colorbar;
fmg; contour(lon(jx-500:jx+500),lat(ix-500:ix+500),dat(ix-500:ix+500,jx-500:jx+500),[-1 0 1]); set(gca,'clim',[-2,+2]); colorbar;
reviewanim([],0,0,0)
help interp2
help griddata
all(newdat(:)==dat(:)|isnan(newdat(:)))
newdat=[]; clear newdat
LON=[]; LAT=[]; clear LON LAT
clear ans ig ix jx
ix = find(isnan(dat));
whos
ix=[]; clear ix
[iix,ix,jx] = find(isnan(dat));
help find
iix=[]; ix=[]; jx=[]; clear ix jx iix
help contourc
help griddata
help scatteredInterpolant
F = scatteredInterpolant(LON,LAT,dat);
[LON,LAT]=meshgrid(lon,lat);
F = scatteredInterpolant(LON,LAT,dat);
F = scatteredInterpolant(lon,lat,dat);
F = scatteredInterpolant([LON(:),LAT(:)],dat);
size([LON(:),LAT(:)])
F = scatteredInterpolant([LON(:),LAT(:)]',dat);
F = scatteredInterpolant([LON(:),LAT(:)],dat);
F = scatteredInterpolant([LON(:),LAT(:)],dat(:));
help scatteredInterpolant
newdat=dat;
newdat(isnan(dat)) = F(LON(isnan(dat)),LAT(isnan(dat)));
F
help filled_ts
help fill
help fillbad
nan+(i*nan)
nan
nan+(i*nan)
clear ans
help fillib
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
quit
%-- 5/8/2016 3:09 PM --%
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
interp_dem_bathy
fmg; contour(lon(jx-500:jx+500),lat(ix-500:ix+500),dat(ix-500:ix+500,jx-500:jx+500),[-1 0 1]); set(gca,'clim',[-2,+2]); colorbar;
[ig,ix]=min(abs(lat-15.25))
[ig,jx]=min(abs(lon-145.75))
fmg; contour(lon(jx-500:jx+500),lat(ix-500:ix+500),dat(ix-500:ix+500,jx-500:jx+500),[-1 0 1]); set(gca,'clim',[-2,+2]); colorbar;
[LON,LAT]=meshgrid(lon,lat);
lh = plot(LON(isnan(dat)),LAT(isnan(dat)),'.')
delete(lh)
ax=axis
axis(ax)
lh = plot(LON(isnan(dat)),LAT(isnan(dat)),'.')
delete(lh)
fmg; contour(lon(jx-500:jx+500),lat(ix-500:ix+500),dat(ix-500:ix+500,jx-500:jx+500),[-1 0 1]); set(gca,'clim',[-2,+2]); colorbar;
axis(ax)
plot(x(1,dx>=1),x(2,dx>=1),'r.');
lh = plot(LON(isnan(dat)),LAT(isnan(dat)),'.')
help colormap
help winter
colormap(winter)
colormap(parula)
colormap(hot)
colormap(bone)
colormap(flag)
colormap(pink)
colormap(parula)
cm = colormap;
colormap(cm)
set(gca,'clim',[-4,4])
fmg; contour(lon(jx-500:jx+500),lat(ix-500:ix+500),dat(ix-500:ix+500,jx-500:jx+500),[-2:2]); set(gca,'clim',[-5,+5]); colorbar;
lh = plot(LON(isnan(dat)),LAT(isnan(dat)),'.')
fmg; contour(lon(jx-500:jx+500),lat(ix-500:ix+500),dat(ix-500:ix+500,jx-500:jx+500),[-2:2]); set(gca,'clim',[-5,+5]); colorbar;
ax=axis
axis(ax)
lh = plot(LON(isnan(dat)),LAT(isnan(dat)),'.')
fmg; [c,h]=contour(lon(jx-500:jx+500),lat(ix-500:ix+500),dat(ix-500:ix+500,jx-500:jx+500),[-2:2]); set(gca,'clim',[-5,+5]); colorbar;
delete(c)
delete(h)
fmg; dx=250; [c,h]=contour(lon(jx-dx:jx+dx),lat(ix-dx:ix+dx),dat(ix-dx:ix+dx,jx-dx:jx+dx),[-2:2]); set(gca,'clim',[-5,+5]); colorbar;
axis(axis)
lh = plot(LON(isnan(dat)),LAT(isnan(dat)),'.')
help inter_field
help interp_field
dx=min(diff(unique(lon)))
dy=min(diff(unique(lat)))
dy=min(distance_wgs84(mean(lon(:)),unique(lat)))
size(mean(lon(:)))
help distance_wgs84
dy=min(distance_wgs84(lat(1:end-1),mean(lon(:)),lat(2:end),mean(lon(:))))
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
interp_dem_bathy
dx
dy
dx*1e3
dy*1e3
[ig,ix]=min(abs(lat-15.25))
[ig,jx]=min(abs(lon-145.75))
fmg; dx=250; [c,h]=contour(lon(jx-dx:jx+dx),lat(ix-dx:ix+dx),dat(ix-dx:ix+dx,jx-dx:jx+dx),[-2:2]); set(gca,'clim',[-5,+5]); colorbar;
fmg; contourf(lon,lat,dat,-1000:100:1000); colorbar;
[LON,LAT]=meshgrid(lon,lat);
lh = plot(LON(isnan(dat)),LAT(isnan(dat)),'k.')
delete(lh)
dbquit
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
help interp_field
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
interp_dem_bathy
dx
dy
npts
fix(3.1)
dbquit
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
interp_dem_bathy
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
interp_dem_bathy
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
interp_dem_bathy
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
help downsample
help interp2
help nanmean
x(1,:,:) = [1:3;1:3]
x(2,:,:) = [1:3;1:3].*2
x
x(1,:,:)
squeeze(x(1,:,:))
mean(x)
squeeze(mean(x))
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
dat=randn(5,5)
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
dat=randn(6,6)
reshape(dat,[2,3,3])
size(dat)
2*3*#
2*3*3
reshape(dat,[3,3,3])
3*3*3
reshape(dat,[4,3,3])
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
help linterp
help liftwave
help lifting
help linspace
for ix=1:36; dat(ind2sub([6,6],ix)) = ix; end; dat
help ind2sub
for x=1:36; [ix,jx]=ind2sub([6,6],x); dat(ix,jx) = x; end; dat
6*^
6*6
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
for x=1:36; [ix,jx]=ind2sub([6,6],x); dat(ix,jx) = x; end; dat
for x=1:36; [ix,jx]=ind2sub([6,6],x); dat(jx,ix) = x; end; dat
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
for x=1:36; [ix,jx]=ind2sub([6,6],x); dat(jx,ix) = x; end; dat
reshape(dat,[4,3,3])
newdat=reshape(dat,[4,3,3]); squeeze(newdat(1,:,:))
newdat=reshape(dat',[4,3,3]); squeeze(newdat(1,:,:))
newdat=reshape(dat',[4,3,3])'; squeeze(newdat(1,:,:))
newdat=reshape(dat,[4,3,3]); squeeze(newdat(1,:,:))
newdat=reshape(dat,[4,3,3]); squeeze(newdat(:,1,:))
newdat=reshape(dat,[4,3,3]); squeeze(newdat(:,:,1))
newdat=reshape(dat,[4,3,3]); squeeze(newdat(:,1,:))
newdat=reshape(dat,[3,4,3]); squeeze(newdat(:,1,:))
newdat=reshape(dat,[3,3,4]); squeeze(newdat(:,:,1))
newdat=reshape(dat,[4,3,3]); squeeze(newdat(:,:,1))
newdat=reshape(dat,[4,3,3]); squeeze(newdat(1,:,:))
newdat=reshape(dat,[4,3,3]); squeeze(newdat(2,:,:))
newdat=reshape(dat',[4,3,3]); squeeze(newdat(1,:,:))
dat
newdat=reshape(dat',[4,3,3]'); squeeze(newdat(1,:,:))
newdat=reshape(dat',[4,3,3]); squeeze(newdat(1,:,:))
help reshape
newdat=reshape(dat',[[],3,3]); squeeze(newdat(1,:,:))
newdat=reshape(dat',[],3,3); squeeze(newdat(1,:,:))
newdat=reshape(dat,[],3,3); squeeze(newdat(1,:,:))
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
interp_dem_bathy
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
interp_dem_bathy
yix
ny
yixix
yixen
nx
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
interp_dem_bathy
yixix
xixix
ny
yix-ny
xix-nx
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
interp_dem_bathy
size(newdat)
dbcont
size(newdat)
dbquit
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
interp_dem_bathy
size(newdat)
fmg; contour(newlon,newlat,newdat,[0 0]);
fmg; contour(newlon,newlat,newdat,-100:50:100);
min(diff(newlon))
max(diff(newlon))
yix-ny:yix+ny
fmg; contourf(newlon,newlat,newdat,-100:50:100);
fmg; contourf(newlon,newlat,newdat,-100:50:100); colorbar;
fmg; contourf(newlon,newlat,newdat,-1000:50:100); colorbar; contour(newlon,newlat,newdat,[0,0],'Color','k','LineWidth',2);
dbquit
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
interp_dem_bathy
clear all
pack
interp_dem_bathy
clear all; pack
interp_dem_bathy
axis([145.74,145.8,15.12,15.17])
size(newdat)
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
interp_dem_bathy
size(newdat)
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
interp_dem_bathy
size(newdat)
clear all; pack
interp_dem_bathy
axis([145.74,145.8,15.12,15.17])
fmg; contourf(newlon,newlat,newdat,-30:2:0); colorbar; contour(newlon,newlat,newdat,[0,0],'Color','k','LineWidth',2);
titlename(['Bathy/DEM at ',num2str(target_res),' m']);
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
interp_dem_bathy
clear all; pack
interp_dem_bathy
axis([145.74,145.8,15.12,15.17])
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
interp_dem_bathy
90/5
6699/18
6699/36
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
interp_dem_bathy
6699/18
nx
target_res/dx/2/1e3
dx
dx=1e3*median(distance_wgs84(mean(lat(:)),lon(1:end-1),mean(lat(:)),lon(2:end)));
dy=1e3*median(distance_wgs84(mean(lat(:)),lon(1:end-1),mean(lat(:)),lon(2:end)));
nx = ceil(target_res/dx/2);
ny = ceil(target_res/dy/2);
nx
target_res/dx/2
dx
dy
dx=1e3*median(distance_wgs84(mean(lat(:)),lon(1:end-1),mean(lat(:)),lon(2:end)));
dy=1e3*median(distance_wgs84(mean(lat(:)),lon(1:end-1),mean(lat(:)),lon(2:end)));
nx = floor(target_res/dx/2);
ny = floor(target_res/dy/2);
xixen = nx+1:(nx*2)+1:length(lon);
yixen = ny+1:(ny*2)+1:length(lat);
newlon = lon(xixen);
newlat = lat(yixen);
6699/nx
90/dx
6699/18.4
6699/18.4058
nx*2
xixen
10+9
29-9
size(datdat)
xix-nx:xix+nx
6699/21
21*dx
dx=1e3*median(distance_wgs84(mean(lat(:)),lon(1:end-1),mean(lat(:)),lon(2:end)));
dy=1e3*median(distance_wgs84(mean(lat(:)),lon(1:end-1),mean(lat(:)),lon(2:end)));
nx = floor(target_res/dx/2);
ny = floor(target_res/dy/2);
xixen = nx+1:(nx*2)+1:length(lon);
yixen = ny+1:(ny*2)+1:length(lat);
newlon = lon(xixen);
newlat = lat(yixen);
size(datdat)  dx=1e3*median(distance_wgs84(mean(lat(:)),lon(1:end-1),mean(lat(:)),lon(2:end)));
dy=1e3*median(distance_wgs84(mean(lat(:)),lon(1:end-1),mean(lat(:)),lon(2:end)));
nx = floor(target_res/dx/2);
ny = floor(target_res/dy/2);
xixen = nx+1:(nx*2)+1:length(lon);
yixen = ny+1:(ny*2)+1:length(lat);
newlon = lon(xixen);
newlat = lat(yixen);
dx=1e3*median(distance_wgs84(mean(lat(:)),lon(1:end-1),mean(lat(:)),lon(2:end)));
dy=1e3*median(distance_wgs84(mean(lat(:)),lon(1:end-1),mean(lat(:)),lon(2:end)));
nx = floor(target_res/dx/2);
ny = floor(target_res/dy/2);
xixen = nx+1:(nx*2)+1:length(lon);
yixen = ny+1:(ny*2)+1:length(lat);
newlon = lon(xixen);
newlat = lat(yixen);
for xixix=1:numel(xixen);
for yixix=1:numel(yixen);
xix = xixen(xixix);
yix = yixen(yixix);
datdat = dat(yix-ny:yix+ny,xix-nx:xix+nx);
newdat(yixix,xixix) = nanmean(datdat(:));
end;
end;
newdat=[]; clear newlon newlat newdat
dx=1e3*median(distance_wgs84(mean(lat(:)),lon(1:end-1),mean(lat(:)),lon(2:end)));
dy=1e3*median(distance_wgs84(mean(lat(:)),lon(1:end-1),mean(lat(:)),lon(2:end)));
nx = floor(target_res/dx/2);
ny = floor(target_res/dy/2);
xixen = nx+1:(nx*2)+1:length(lon);
yixen = ny+1:(ny*2)+1:length(lat);
newlon = lon(xixen);
newlat = lat(yixen);
for xixix=1:numel(xixen);
for yixix=1:numel(yixen);
xix = xixen(xixix);
yix = yixen(yixix);
datdat = dat(yix-ny:yix+ny,xix-nx:xix+nx);
newdat(yixix,xixix) = nanmean(datdat(:));
end;
end;
ny
yix
yix+ny
newdat=[]; clear newlon newlat newdat
dx=1e3*median(distance_wgs84(mean(lat(:)),lon(1:end-1),mean(lat(:)),lon(2:end)));
dy=1e3*median(distance_wgs84(mean(lat(:)),lon(1:end-1),mean(lat(:)),lon(2:end)));
nx = floor(target_res/dx/2);
ny = floor(target_res/dy/2);
xixen = nx+1:(nx*2)+1:length(lon)-nx;
yixen = ny+1:(ny*2)+1:length(lat)-ny;
newlon = lon(xixen);
newlat = lat(yixen);
for xixix=1:numel(xixen);
for yixix=1:numel(yixen);
xix = xixen(xixix);
yix = yixen(yixix);
datdat = dat(yix-ny:yix+ny,xix-nx:xix+nx);
newdat(yixix,xixix) = nanmean(datdat(:));
end;
end;
6699/18
6699/19
90/dx
dx
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
interp_dem_bathy
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
interp_dem_bathy
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
3600/3
3/3600
ans*111e3
111e3/3600
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
interp_dem_bathy
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
interp_dem_bathy
rawdat=[]; clear all; pack
interp_dem_bathy
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
interp_dem_bathy
rawdat=[]; clear all; pack
interp_dem_bathy
axis([145.74,145.8,15.12,15.17])
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
interp_dem_bathy
axis([145.74,145.8,15.12,15.17])
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
interp_dem_bathy
axis([145.74,145.8,15.12,15.17])
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
interp_dem_bathy
axis([145.74,145.8,15.12,15.17])
set_pcolor_cursor
dir ../*cursor*m
dir ../../\*cursor*m
set_surf_cursor
interp2(lon(652:699),lat(
interp2(lon(652:699),lat)
interp2(lon(652:699),lat(260:294),dat(260:294,652:699))
n=lon(652:699); t=lat(260:294); d=dat(260:294,652:699);
interp2(n,t,d,lon(278),lat(679))
interp2(n,t,d,lon(278),lat(679),'spline')
interp2(n,t,d,lon(679),lat(278),'spline')
[N,T]=meshgrid(n,t);
interp2(n,t,d,N(isnan(d)),T(isnan(d)),'spline')
interp2(n,t,d,N(isnan(d)),T(isnan(d)),'spline',nan)
nd=interp2(n,t,d,N(isnan(d)),T(isnan(d)),'spline',nan)
[m,d]=lastwarn
help warning
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
interp_dem_bathy
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
interp_dem_bathy
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
interp_dem_bathy
dir *saip*mat
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
load('pibhmc_bathy_5m_saipan.mat');
load('pibhmc_bathy_5m_saipan.dods.mat');
fmg; contour(lon,lat,dat,-100);
fmg; contour(lon,lat,dat,[-100 -100]);
[LON,LAT]=meshgrid(lon,lat);
lh = plot(LON(isnan(dat)),LAT(isnan(dat)),'k.')
min(dat(:))
nansummary(dat)
fmg; hist(dat)
fmg; hist(dat(:),100);
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
interp_dem_bathy
[ig,ix]=min(abs(lat-15.25))
[ig,jx]=min(abs(lon-145.75))
axis([145.74,145.8,15.12,15.17])
clear ig ix jx
n=lon(652:699); t=lat(260:294); d=dat(260:294,652:699);
nd=interp2(n,t,d,N(isnan(d)),T(isnan(d)),'spline',nan);
[N,T]=meshgrid(n,t);
nd=interp2(n,t,d,N(isnan(d)),T(isnan(d)),'spline',nan);
warning
warning('on','matlab:interp2:nanstrip')
warning
nd=interp2(n,t,d,N(isnan(d)),T(isnan(d)),'spline',nan);
set_surf_cursor
clear n t d N T nd
nix=652:699; tix=260:294; n=lon(nix); t=lat(tix); d=dat(tix,nix);
nd=interp2(n,t,d,N(isnan(d)),T(isnan(d)),'spline',nan);
nix=652:699; tix=260:294; n=lon(nix); t=lat(tix); d=dat(tix,nix); [N,T]=meshgrid(n,t);
nd=interp2(n,t,d,N(isnan(d)),T(isnan(d)),'spline',nan);
nix=652:699; tix=278:280; n=lon(nix); t=lat(tix); d=dat(tix,nix); [N,T]=meshgrid(n,t);
nd=interp2(n,t,d,N(isnan(d)),T(isnan(d)),'spline',nan);
help interp1
help interp2
nix=652:699; tix=277:280; n=lon(nix); t=lat(tix); d=dat(tix,nix); [N,T]=meshgrid(n,t);
nd=interp2(n,t,d,N(isnan(d)),T(isnan(d)),'spline',nan);
nix=652:699; tix=277:281; n=lon(nix); t=lat(tix); d=dat(tix,nix); [N,T]=meshgrid(n,t);
nd=interp2(n,t,d,N(isnan(d)),T(isnan(d)),'spline',nan);
nix=652:699; tix=275:283; n=lon(nix); t=lat(tix); d=dat(tix,nix); [N,T]=meshgrid(n,t);
nd=interp2(n,t,d,N(isnan(d)),T(isnan(d)),'spline',nan);
nix=652:699; tix=276:283; n=lon(nix); t=lat(tix); d=dat(tix,nix); [N,T]=meshgrid(n,t);
nd=interp2(n,t,d,N(isnan(d)),T(isnan(d)),'spline',nan);
nix=652:699; tix=277:283; n=lon(nix); t=lat(tix); d=dat(tix,nix); [N,T]=meshgrid(n,t);
nd=interp2(n,t,d,N(isnan(d)),T(isnan(d)),'spline',nan);
nix=652:699; tix=277:282; n=lon(nix); t=lat(tix); d=dat(tix,nix); [N,T]=meshgrid(n,t);
nd=interp2(n,t,d,N(isnan(d)),T(isnan(d)),'spline',nan);
nix=652:699; tix=277:281; n=lon(nix); t=lat(tix); d=dat(tix,nix); [N,T]=meshgrid(n,t);
nd=interp2(n,t,d,N(isnan(d)),T(isnan(d)),'spline',nan);
nix=652:699; tix=260:294; n=lon(nix); t=lat(tix); d=dat(tix,nix); [N,T]=meshgrid(n,t);
nd=interp2(n,t,d,N(isnan(d)),T(isnan(d)),'spline',nan);
nix=652:699; tix=260:294; n=lon(nix); t=lat(tix); d=dat(tix,nix); [N,T]=meshgrid(n,t); % THIS IS THE FULL BOX
nix=652:692; tix=260:294; n=lon(nix); t=lat(tix); d=dat(tix,nix); [N,T]=meshgrid(n,t);
nd=interp2(n,t,d,N(isnan(d)),T(isnan(d)),'spline',nan);
nix=652:691; tix=260:294; n=lon(nix); t=lat(tix); d=dat(tix,nix); [N,T]=meshgrid(n,t);
nd=interp2(n,t,d,N(isnan(d)),T(isnan(d)),'spline',nan);
nix=665:691; tix=260:294; n=lon(nix); t=lat(tix); d=dat(tix,nix); [N,T]=meshgrid(n,t);
nd=interp2(n,t,d,N(isnan(d)),T(isnan(d)),'spline',nan);
nix=666:691; tix=260:294; n=lon(nix); t=lat(tix); d=dat(tix,nix); [N,T]=meshgrid(n,t);
nd=interp2(n,t,d,N(isnan(d)),T(isnan(d)),'spline',nan);
nix=667:691; tix=260:294; n=lon(nix); t=lat(tix); d=dat(tix,nix); [N,T]=meshgrid(n,t);
nd=interp2(n,t,d,N(isnan(d)),T(isnan(d)),'spline',nan);
nix=669:691; tix=260:294; n=lon(nix); t=lat(tix); d=dat(tix,nix); [N,T]=meshgrid(n,t);
nd=interp2(n,t,d,N(isnan(d)),T(isnan(d)),'spline',nan);
nix=675:691; tix=260:294; n=lon(nix); t=lat(tix); d=dat(tix,nix); [N,T]=meshgrid(n,t);
nd=interp2(n,t,d,N(isnan(d)),T(isnan(d)),'spline',nan);
nix=675:689; tix=260:294; n=lon(nix); t=lat(tix); d=dat(tix,nix); [N,T]=meshgrid(n,t);
nd=interp2(n,t,d,N(isnan(d)),T(isnan(d)),'spline',nan);
nix=675:685; tix=260:294; n=lon(nix); t=lat(tix); d=dat(tix,nix); [N,T]=meshgrid(n,t);
nd=interp2(n,t,d,N(isnan(d)),T(isnan(d)),'spline',nan);
nix=675:685; tix=1:777; n=lon(nix); t=lat(tix); d=dat(tix,nix); [N,T]=meshgrid(n,t);
nd=interp2(n,t,d,N(isnan(d)),T(isnan(d)),'spline',nan);
nix=675:685; tix=260:294; n=lon(nix); t=lat(tix); d=dat(tix,nix); [N,T]=meshgrid(n,t);
nd=interp2(n,t,d,N(isnan(d)),T(isnan(d)),'spline',nan);
nix=652:699; tix=260:294; n=lon(nix); t=lat(tix); d=dat(tix,nix); [N,T]=meshgrid(n,t); % THIS IS THE FULL BOX
nix=675:685; tix=260:294; n=lon(nix); t=lat(tix); d=dat(tix,nix); [N,T]=meshgrid(n,t);
nd=interp2(n,t,d,N(isnan(d)),T(isnan(d)),'spline',nan);
nix=675:685; tix=255:294; n=lon(nix); t=lat(tix); d=dat(tix,nix); [N,T]=meshgrid(n,t);
nd=interp2(n,t,d,N(isnan(d)),T(isnan(d)),'spline',nan);
nix=675:685; tix=254:294; n=lon(nix); t=lat(tix); d=dat(tix,nix); [N,T]=meshgrid(n,t);
nd=interp2(n,t,d,N(isnan(d)),T(isnan(d)),'spline',nan);
nix=675:685; tix=253:294; n=lon(nix); t=lat(tix); d=dat(tix,nix); [N,T]=meshgrid(n,t);
nd=interp2(n,t,d,N(isnan(d)),T(isnan(d)),'spline',nan);
nix=675:685; tix=1:294; n=lon(nix); t=lat(tix); d=dat(tix,nix); [N,T]=meshgrid(n,t);
nd=interp2(n,t,d,N(isnan(d)),T(isnan(d)),'spline',nan);
nix=675:685; tix=227:294; n=lon(nix); t=lat(tix); d=dat(tix,nix); [N,T]=meshgrid(n,t);
nd=interp2(n,t,d,N(isnan(d)),T(isnan(d)),'spline',nan);
nix=675:685; tix=226:294; n=lon(nix); t=lat(tix); d=dat(tix,nix); [N,T]=meshgrid(n,t);
nd=interp2(n,t,d,N(isnan(d)),T(isnan(d)),'spline',nan);
nix=675:685; tix=228:294; n=lon(nix); t=lat(tix); d=dat(tix,nix); [N,T]=meshgrid(n,t);
nd=interp2(n,t,d,N(isnan(d)),T(isnan(d)),'spline',nan);
nix=675:685; tix=242:294; n=lon(nix); t=lat(tix); d=dat(tix,nix); [N,T]=meshgrid(n,t);
nd=interp2(n,t,d,N(isnan(d)),T(isnan(d)),'spline',nan);
nix=675:685; tix=247:294; n=lon(nix); t=lat(tix); d=dat(tix,nix); [N,T]=meshgrid(n,t);
nd=interp2(n,t,d,N(isnan(d)),T(isnan(d)),'spline',nan);
nix=675:685; tix=251:294; n=lon(nix); t=lat(tix); d=dat(tix,nix); [N,T]=meshgrid(n,t);
nd=interp2(n,t,d,N(isnan(d)),T(isnan(d)),'spline',nan);
axis tight
help area
lookfor area
help polyarea
lookfor contig
help grayconnected
help bwselect
which bwselect
help grayconnected
bw = grayconnected(dat,277,682,0)
bw = grayconnected(dat,277,682,0);
numel(find(bw))
nix=652:699; tix=260:294; n=lon(nix); t=lat(tix); d=dat(tix,nix); [N,T]=meshgrid(n,t); % THIS IS THE FULL BOX
bw = grayconnected(d,277,682,0);
bw = grayconnected(d,17,30,0);
bw = grayconnected(d,17,30,0)
clear bw
bw = grayconnected(d,17,30,0)
numel(find(bw))
nix=652:699; tix=260:294; n=lon(nix); t=lat(tix); d=dat(tix,nix); [N,T]=meshgrid(n,t); % THIS IS THE FULL BOX
bw = grayconnected(d,19,30,0)
numel(find(bw))
bw = grayconnected(d,19,30,0); numel(find(bw))
bw = grayconnected(d,18,30,0); numel(find(bw))
bw = grayconnected(d,17,30,0); numel(find(bw))
bw = grayconnected(d,17,29,0); numel(find(bw))
d(isnan(d)) = -99999;
bw = grayconnected(d,17,29,0); numel(find(bw))
bw = grayconnected(d,17,29,1); numel(find(bw))
help grayconnected
d
bw = grayconnected(d,17,29,1); numel(find(bw))
help grayconnected
d(isnan(d)) = -99999;
bw = grayconnected(d,17,29,1); numel(find(bw))
fmg; contourf(n,t,d);
d(d==-99999)=nan;
fmg; contourf(n,t,d);
d(isnan(d)) = -99999;
fmg; contourf(n,t,d);
bw = grayconnected(d,17,29,1); numel(find(bw))
fmg; contourf(n,t,d);
set_surf_cursor
bw = grayconnected(d,17,29,100); numel(find(bw))
help grayconnected
bw = grayconnected(d,19,29,100); numel(find(bw))
bw = grayconnected(d,19,29,1); numel(find(bw))
find(isnan(dat),1)
dat(1,1)
help grayconnected
help bwselect
help roipoly
bw = grayconnected(d,19,29,1); numel(find(bw))
find(isnan(d) & ~bw(isnan(d)))
size(d),size(bw)
find(isnan(d) & bw))
find(isnan(d) & bw)
find(isnan(d) & bw(:))
find(isnan(d(~bw)))
fmg; contourf(n,t,d);
numel(find(isnan(d(~bw))))
numel(find(d(~bw)))
numel(find(isnan(d(~bw))))
numel(find(d(~bw)==-99999))
bw = grayconnected(d,19,29,1); numel(find(bw))
d(d==-99999)=nan;
bw = grayconnected(d,19,29,1); numel(find(bw))
fmg; contourf(n,t,d);
bw = grayconnected(d,19,29,0); numel(find(bw))
bw = grayconnected(d,19,29,nan); numel(find(bw))
d(d==-99999)=nan;
d(isnan(d)) = -99999;
bw = grayconnected(d,19,29,0); numel(find(bw))
lookfor edge
help pdeexpd
which pdeexpd
help edge
numel(find(d(~bw)==-99999))
numel(find(d(~bw)<359))
numel(find(d(~bw)>0))
numel(find(d(~bw)<-359))
set_surf_cursor
0.1571-0.1568
(0.1571-0.1568)*111e3
help edge
edge(d)
bw=edge(d)
bw=edge(d); fmg; contour(n,t,bw);
d(d==-99999)=nan;
bw=edge(d); fmg; contour(n,t,bw);
numel(find(bw))
d(isnan(d)) = -99999;
bw = grayconnected(d,19,29,0); numel(find(bw))
d(bw)
help grayconnected
help encode
help rltool
help run_length_encode
type run_length_encode.m
help spline2d
help spline
which spline
help spline
help spline2d
help curve
help rcurve
lookfor fitting
which fittype
help curvefit
help splines
help spapi
help splines
lookfor edge
fmg; contourf(n,t,d);
d(d==-99999)=nan;
fmg; contourf(n,t,d);
firstix = find(~isnan(d),1);
size(~isnan(d))
firstix = find(~isnan(dat),1);
clear firstix
firstix = find(~isnan(dat));
size(~isnan(dat))
help find
[ix,jx] = find(isnan(d));
ix
[ix,jx] = find(isnan(d));
nbix = find(~isnan(d(ix-1,jx-1)) | ~isnan(d(ix+1,jx-1)) | ~isnan(d(ix-1,jx+1)) | ~isnan(d(ix+1,jx+1)));
nbix = find(~isnan(d(ix-1,jx-1)) | ~isnan(d(ix+1,jx-1)) | ~isnan(d(ix-1,jx+1)) | ~isnan(d(ix+1,jx+1))); size(nbix)
size(LON)
size(dat)
fmg; contourf(n,t,d); colorbar;
plot(N(nbix),T(nbix),'k.')
clear nbix
[nbix,nbjx] = find(~isnan(d(ix-1,jx-1)) | ~isnan(d(ix+1,jx-1)) | ~isnan(d(ix-1,jx+1)) | ~isnan(d(ix+1,jx+1))); size(nbix)
nbix
nbjx
[nbix,nbjx] = find(~isnan(d(ix-1,jx-1)) | ~isnan(d(ix+1,jx-1)) | ~isnan(d(ix-1,jx+1)) | ~isnan(d(ix+1,jx+1))); size(nbix), nbix=ix(nbix); nbjx=jx(nbjx);
plot(N(nbix),T(nbix),'k.')
N(nbix)
nbix
njix
nbjx
nbix = find(~isnan(d(ix-1,jx-1)) | ~isnan(d(ix+1,jx-1)) | ~isnan(d(ix-1,jx+1)) | ~isnan(d(ix+1,jx+1))); size(nbix), nbix=ix(nbix);
nbix
nbix = find(~isnan(d(ix-1,jx-1)) | ~isnan(d(ix+1,jx-1)) | ~isnan(d(ix-1,jx+1)) | ~isnan(d(ix+1,jx+1))); size(nbix),
nbix
[nbix,nbjx] = find(~isnan(d(ix-1,jx-1)) | ~isnan(d(ix+1,jx-1)) | ~isnan(d(ix-1,jx+1)) | ~isnan(d(ix+1,jx+1))); size(nbix), nbix=ix(nbix); nbjx=jx(nbjx);
nbix
nbjx
nansummary(nbjx)
nansummary(nbix)
lh = plot(N(nbix,nbjx),T(nbix,nbjx),'k.')
lh = plot(N(nbix,nbjx),T(nbix,nbjx),'k.');
delete(lh)
fmg; contourf(n,t,d); colorbar;
lh = plot(N(nbix,nbjx),T(nbix,nbjx),'k.');
delete(lh)
[ix,jx] = find(isnan(d));
lh = plot(N(ix,jx),T(ix,jx),'k.');
delete(lh)
[ig,ix,jx] = find(isnan(d));
lh = plot(N(ix,jx),T(ix,jx),'k.');
help find
ix = find(isnan(d));
lh = plot(N(ix),T(ix),'k.');
delete(lh)
[yix,xix] = ind2sub(size(N),ix);
size(yix)
size(ix)
yix
nbix = find(~isnan(d(yix-1,xix-1)) | ~isnan(d(yix+1,xix-1)) | ~isnan(d(yix-1,xix+1)) | ~isnan(d(yix+1,xix+1))); size(nbix), nbix=ix(nbix);
nansummary(nbix)
nbix = find(~isnan(d(yix-1,xix-1)) | ~isnan(d(yix+1,xix-1)) | ~isnan(d(yix-1,xix+1)) | ~isnan(d(yix+1,xix+1))); size(nbix),
nbix
size(yix)
nbix
nbix = find(~isnan(d(yix-1,xix-1)) | ~isnan(d(yix+1,xix-1)) | ~isnan(d(yix-1,xix+1)) | ~isnan(d(yix+1,xix+1))); size(nbix),
nbix = find(~isnan(d(yix-1,xix-1)) | ~isnan(d(yix+1,xix-1)) | ~isnan(d(yix-1,xix+1)) | ~isnan(d(yix+1,xix+1))); size(nbix), nbix=ix(nbix);
ix = find(isnan(d));
[yix,xix] = ind2sub(size(N),ix);
nbix = find(~isnan(d(yix-1,xix-1)) | ~isnan(d(yix+1,xix-1)) | ~isnan(d(yix-1,xix+1)) | ~isnan(d(yix+1,xix+1))); size(nbix), nbix=ix(nbix);
yix
xix
nbix = find(~isnan(d(yix-1,xix-1)) | ~isnan(d(yix+1,xix-1)) | ~isnan(d(yix-1,xix+1)) | ~isnan(d(yix+1,xix+1))); size(nbix), nbix=ix(nbix);
ix = find(isnan(d));
[yix,xix] = ind2sub(size(N),ix);
nbix = find(~isnan(d(yix-1,xix-1)) | ~isnan(d(yix+1,xix-1)) | ~isnan(d(yix-1,xix+1)) | ~isnan(d(yix+1,xix+1))); size(nbix), nbix=ix(nbix);
nbix = find(~isnan(d(yix-1,xix-1)) | ~isnan(d(yix+1,xix-1)) | ~isnan(d(yix-1,xix+1)) | ~isnan(d(yix+1,xix+1))); size(nbix),
nbix = find(~isnan(d(yix-1,xix-1)) | ~isnan(d(yix+1,xix-1)) | ~isnan(d(yix-1,xix+1)) | ~isnan(d(yix+1,xix+1))); size(nbix), nbix=ix(nbix);
size(xix)
nansummary(xix)
size(d)
35*48
size(~isnan(d(yix-1,xix-1)))
size(~isnan(d(yix-1,xix'-1)))
timenow
size(~isnan(d(yix,xix)))
help sub2ind
sub2ind([3,3],0,1)
edgeix = find_field_holes(d);
lh = plot(N(edgeix),T(edgeix),'k.');
edgeix = find_field_holes(d);
args={pi,[3:5],[]}
args(1) = {}
args{1} = []
args={pi,[3:5],[]}
args{1} = {}
args={pi,[3:5],[]}
args(1) = {}
args={pi,[3:5],[]}
args{1}
args(1) = []
now-(82*365.25)
datestr(now-(82*365.25))
datestr(now-(83*365.25))
edgeix = find_field_holes(d);
delete(lh)
lh = plot(N(edgeix),T(edgeix),'k.');
delete(lh)
edgeix = find_field_holes(d,true);
lh = plot(N(edgeix),T(edgeix),'k.');
delete(lh)
edgeix = find_field_holes(d,true);
args
dbquit
edgeix = find_field_holes(d,true);
delete(lh)
lh = plot(N(edgeix),T(edgeix),'k.');
delete(lh)
edgeix = find_field_holes(d,2,true);
lh = plot(N(edgeix),T(edgeix),'k.');
delete(lh)
edgeix = find_field_holes(d,2,false);
lh = plot(N(edgeix),T(edgeix),'k.');
nix=652:699; tix=260:294; n=lon(nix); t=lat(tix); d=dat(tix,nix); [N,T]=meshgrid(n,t); % THIS IS THE FULL BOX
nd=interp2(n,t,d,N(isnan(d)),T(isnan(d)),'spline',nan);
nix=652:699; tix=240:294; n=lon(nix); t=lat(tix); d=dat(tix,nix); [N,T]=meshgrid(n,t); % *DEEP* BOX
fmg; contourf(n,t,d); colorbar;
edgeix = find_field_holes(d,2,true);
dbstop if error
edgeix = find_field_holes(d,2,true);
dbup
dx
xix
max(xix)
size(d)
dbquit
edgeix = find_field_holes(d,2,true);
dbup
dx
max(xix)
size(d,1)
dbquit
edgeix = find_field_holes(d,2,true);
dbup
max(xix)
dx
dbquit
edgeix = find_field_holes(d,2,true);
fmg; contourf(n,t,d); colorbar;
lh = plot(N(edgeix),T(edgeix),'k.');
delete(lh)
edgeix = find_field_holes(d,2,true,true);
lh = plot(N(edgeix),T(edgeix),'k.');
delete(lh)
edgeix = find_field_holes(d,1,true,true);
lh = plot(N(edgeix),T(edgeix),'k.');
nd=interp2(n,t,d,N(isnan(d)),T(isnan(d)),'spline',nan);
nix=652:699; tix=200:294; n=lon(nix); t=lat(tix); d=dat(tix,nix); [N,T]=meshgrid(n,t); % *DEEPER* BOX
nd=interp2(n,t,d,N(isnan(d)),T(isnan(d)),'spline',nan);
nix=652:699; tix=100:294; n=lon(nix); t=lat(tix); d=dat(tix,nix); [N,T]=meshgrid(n,t); % *DEEPER* BOX
nd=interp2(n,t,d,N(isnan(d)),T(isnan(d)),'spline',nan);
s = warning('off','matlab:interp2:nanstrip')
nix=652:699; tix=10:294; n=lon(nix); t=lat(tix); d=dat(tix,nix); [N,T]=meshgrid(n,t); % *DEEPER* BOX
nd=interp2(n,t,d,N(isnan(d)),T(isnan(d)),'spline',nan);
dbup
dbquit
dbclear all
nix=652:699; tix=100:294; n=lon(nix); t=lat(tix); d=dat(tix,nix); [N,T]=meshgrid(n,t); % *DEEPER* BOX
nd=interp2(n,t,d,N(isnan(d)),T(isnan(d)),'spline',nan);
nix=652:699; tix=50:294; n=lon(nix); t=lat(tix); d=dat(tix,nix); [N,T]=meshgrid(n,t); % BAD
nd=interp2(n,t,d,N(isnan(d)),T(isnan(d)),'spline',nan);
nix=652:699; tix=30:294; n=lon(nix); t=lat(tix); d=dat(tix,nix); [N,T]=meshgrid(n,t); % BAD
nd=interp2(n,t,d,N(isnan(d)),T(isnan(d)),'spline',nan);
nix=652:699; tix=20:294; n=lon(nix); t=lat(tix); d=dat(tix,nix); [N,T]=meshgrid(n,t); % BAD
nd=interp2(n,t,d,N(isnan(d)),T(isnan(d)),'spline',nan);
nix=652:699; tix=15:294; n=lon(nix); t=lat(tix); d=dat(tix,nix); [N,T]=meshgrid(n,t); % BAD
nd=interp2(n,t,d,N(isnan(d)),T(isnan(d)),'spline',nan);
nix=652:699; tix=19:294; n=lon(nix); t=lat(tix); d=dat(tix,nix); [N,T]=meshgrid(n,t); % JUST OK
nd=interp2(n,t,d,N(isnan(d)),T(isnan(d)),'spline',nan);
nix=652:699; tix=20:294; n=lon(nix); t=lat(tix); d=dat(tix,nix); [N,T]=meshgrid(n,t); % JUST OK
nd=interp2(n,t,d,N(isnan(d)),T(isnan(d)),'spline',nan);
nix=652:699; tix=19:294; n=lon(nix); t=lat(tix); d=dat(tix,nix); [N,T]=meshgrid(n,t); % BAD
nd=interp2(n,t,d,N(isnan(d)),T(isnan(d)),'spline',nan);
fmg; contourf(n,t,d); colorbar;
fmg; contourf(n,t,d); daspect([1,cosd(t(1)),1]); colorbar;
nix=1:699; tix=20:294; n=lon(nix); t=lat(tix); d=dat(tix,nix); [N,T]=meshgrid(n,t); % JUST OK
nd=interp2(n,t,d,N(isnan(d)),T(isnan(d)),'spline',nan);
fmg; contourf(n,t,d); daspect([1,cosd(t(1)),1]); colorbar;
nix=100:699; tix=20:294; n=lon(nix); t=lat(tix); d=dat(tix,nix); [N,T]=meshgrid(n,t); % JUST OK
nd=interp2(n,t,d,N(isnan(d)),T(isnan(d)),'spline',nan);
nix=200:699; tix=20:294; n=lon(nix); t=lat(tix); d=dat(tix,nix); [N,T]=meshgrid(n,t); % JUST OK
nd=interp2(n,t,d,N(isnan(d)),T(isnan(d)),'spline',nan);
nix=600:699; tix=20:294; n=lon(nix); t=lat(tix); d=dat(tix,nix); [N,T]=meshgrid(n,t); % JUST OK
nd=interp2(n,t,d,N(isnan(d)),T(isnan(d)),'spline',nan);
nix=500:699; tix=20:294; n=lon(nix); t=lat(tix); d=dat(tix,nix); [N,T]=meshgrid(n,t); % JUST OK
nd=interp2(n,t,d,N(isnan(d)),T(isnan(d)),'spline',nan);
nix=300:699; tix=20:294; n=lon(nix); t=lat(tix); d=dat(tix,nix); [N,T]=meshgrid(n,t); % JUST OK
nd=interp2(n,t,d,N(isnan(d)),T(isnan(d)),'spline',nan);
nix=200:699; tix=20:294; n=lon(nix); t=lat(tix); d=dat(tix,nix); [N,T]=meshgrid(n,t); % JUST OK
nd=interp2(n,t,d,N(isnan(d)),T(isnan(d)),'spline',nan);
nix=250:699; tix=20:294; n=lon(nix); t=lat(tix); d=dat(tix,nix); [N,T]=meshgrid(n,t); % JUST OK
nd=interp2(n,t,d,N(isnan(d)),T(isnan(d)),'spline',nan);
nix=275:699; tix=20:294; n=lon(nix); t=lat(tix); d=dat(tix,nix); [N,T]=meshgrid(n,t); % JUST OK
nd=interp2(n,t,d,N(isnan(d)),T(isnan(d)),'spline',nan);
nix=287:699; tix=20:294; n=lon(nix); t=lat(tix); d=dat(tix,nix); [N,T]=meshgrid(n,t); % JUST OK
nd=interp2(n,t,d,N(isnan(d)),T(isnan(d)),'spline',nan);
nix=294:699; tix=20:294; n=lon(nix); t=lat(tix); d=dat(tix,nix); [N,T]=meshgrid(n,t); % JUST OK
nd=interp2(n,t,d,N(isnan(d)),T(isnan(d)),'spline',nan);
nix=297:699; tix=20:294; n=lon(nix); t=lat(tix); d=dat(tix,nix); [N,T]=meshgrid(n,t); % JUST OK
nd=interp2(n,t,d,N(isnan(d)),T(isnan(d)),'spline',nan);
nix=300:699; tix=20:294; n=lon(nix); t=lat(tix); d=dat(tix,nix); [N,T]=meshgrid(n,t); % JUST OK
nd=interp2(n,t,d,N(isnan(d)),T(isnan(d)),'spline',nan);
nix=299:699; tix=20:294; n=lon(nix); t=lat(tix); d=dat(tix,nix); [N,T]=meshgrid(n,t); % JUST OK
nd=interp2(n,t,d,N(isnan(d)),T(isnan(d)),'spline',nan);
fmg; contourf(n,t,d); daspect([1,cosd(t(1)),1]); colorbar;
nix=300:699; tix=20:294; n=lon(nix); t=lat(tix); d=dat(tix,nix); [N,T]=meshgrid(n,t); % JUST OK
nd=interp2(n,t,d,N(isnan(d)),T(isnan(d)),'spline',nan);
fmg; contourf(n,t,d); daspect([1,cosd(t(1)),1]); colorbar;
nix=299:699; tix=20:294; n=lon(nix); t=lat(tix); d=dat(tix,nix); [N,T]=meshgrid(n,t); % JUST OK
nd=interp2(n,t,d,N(isnan(d)),T(isnan(d)),'spline',nan);
fmg; contourf(n,t,d); daspect([1,cosd(t(1)),1]); colorbar;
set_surf_cursor
titlename('BAD')
titlename('OK')
nix=300:699; tix=20:294; n=lon(nix); t=lat(tix); d=dat(tix,nix); [N,T]=meshgrid(n,t); % JUST OK
nd=interp2(n,t,d,N(isnan(d)),T(isnan(d)),'spline',nan);
oldd=d;
d(232,400)
d(233,400)
d(233,400) = nan;
nd=interp2(n,t,d,N(isnan(d)),T(isnan(d)),'spline',nan);
d(234,400) = nan;
nd=interp2(n,t,d,N(isnan(d)),T(isnan(d)),'spline',nan);
d(235,400) = nan;
nd=interp2(n,t,d,N(isnan(d)),T(isnan(d)),'spline',nan);
timenow
d(236,400) = nan;
nd=interp2(n,t,d,N(isnan(d)),T(isnan(d)),'spline',nan);
fmg; contourf(n,t,d); daspect([1,cosd(t(1)),1]); colorbar;
fmg; contourf(n,t,d); daspect([1,cosd(t(1)),1]); colorbar; set_surf_cursor
d(237,400) = nan;
nd=interp2(n,t,d,N(isnan(d)),T(isnan(d)),'spline',nan);
fmg; contourf(n,t,d); daspect([1,cosd(t(1)),1]); colorbar; set_surf_cursor
nd=interp2(n,t,d,N(isnan(d)),T(isnan(d)),'spline',nan);
find(isnan(d(1,:))&isnan(d(:end)
find(isnan(d(:,1))&isnan(d(:,end)))
help isnan
find(isnan(d(:,1))&isnan(d(:,end)))
find(isnan(d(:,1))&isnan(d(:,end)))'
fmg; contourf(n,t,d); daspect([1,cosd(t(1)),1]); colorbar; set_surf_cursor
nd=interp2(n,t,d,N(isnan(d)),T(isnan(d)),'spline',nan);
titlename('Made just bad')
dbstop if error
nd=interp2(n,t,d,N(isnan(d)),T(isnan(d)),'spline',nan);
ncolnan,nrownan
size(X),size(V)
size(Y)
dbquit
nd=interp2(n,t,d,N(isnan(d)),T(isnan(d)),'spline',nan);
size(X),size(Y),size(V)
Xgrid
iscell(Xgrid)
fullgrid
fmg; contour(V);
fmg; contourf(V);
size(nanv)
size(sum(isnan(V)))
size(sum(isnan(V),1))
size( find( sum(isnan(V)) ) )
size( find( sum(isnan(V),1) ) )
newV=V;
clear newV
oldV=V;
size(V), V(:,nanv) = []; size(V)
fullgrid
dbquit
d=oldd;
nd=interp2(n,t,d,N(isnan(d)),T(isnan(d)),'spline',nan);
d(237,400) = nan;
nd=interp2(n,t,d,N(isnan(d)),T(isnan(d)),'spline',nan);
d=oldd;
d(233:236,400) = nan;
nd=interp2(n,t,d,N(isnan(d)),T(isnan(d)),'spline',nan);
d=oldd;
d(234:237,400) = nan;
nd=interp2(n,t,d,N(isnan(d)),T(isnan(d)),'spline',nan);
dbquit
d=oldd;
nd=interp2(n,t,d,N(isnan(d)),T(isnan(d)),'spline',nan);
d(235:237,400) = nan;
nd=interp2(n,t,d,N(isnan(d)),T(isnan(d)),'spline',nan);
dbquit
d=oldd;
d(236:237,400) = nan;
nd=interp2(n,t,d,N(isnan(d)),T(isnan(d)),'spline',nan);
dbquit
d=oldd;
d(237:237,400) = nan;
nd=interp2(n,t,d,N(isnan(d)),T(isnan(d)),'spline',nan);
d([235,237],400) = nan;
nd=interp2(n,t,d,N(isnan(d)),T(isnan(d)),'spline',nan);
d([234:235,237],400) = nan;
nd=interp2(n,t,d,N(isnan(d)),T(isnan(d)),'spline',nan);
d([233:235,237],400) = nan;
nd=interp2(n,t,d,N(isnan(d)),T(isnan(d)),'spline',nan);
dbstop interp2>stripnanwrapper
nd=interp2(n,t,d,N(isnan(d)),T(isnan(d)),'spline',nan);
size(V)
dbquit
dbstatus
dbclear interp2>stripnanwrapper
dbstatus
d=oldd;
d([233:235,237:238],400) = nan;
nd=interp2(n,t,d,N(isnan(d)),T(isnan(d)),'spline',nan);
dbquit
d([233:235,237:237],400) = nan;
nd=interp2(n,t,d,N(isnan(d)),T(isnan(d)),'spline',nan);
dbquit
d([233:235,237],400) = nan;
nd=interp2(n,t,d,N(isnan(d)),T(isnan(d)),'spline',nan);
dbquit
d=oldd;
d([233:235,237],400) = nan;
nd=interp2(n,t,d,N(isnan(d)),T(isnan(d)),'spline',nan);
oldV=V;
V(:,nanv) = []; size(V)
V=oldV;
fmg; contourf(V);
fmg; contourf(V); set(gca,'YDir','reverse');
fmg; contourf(V); set(gca,'XDir','reverse','YDir','reverse');
fmg; contourf(V); set(gca,'XDir','reverse');
fmg; contourf(V); set(gca,'YDir','reverse');
size(V)
size(nanv)
fmg; contourf(V); set(gca,'YDir','reverse');
fmg; contourf(V); set(gca,'YDir','reverse'); daspect([1,cosd(14),1]);
fmg; contourf(V'); set(gca,'YDir','reverse'); daspect([1,cosd(14),1]);
fmg; contourf(V'); set(gca,'XDir','reverse'); daspect([1,cosd(14),1]);
fmg; contourf(V'); daspect([1,cosd(14),1]);
dbquit
none
minnan
minnans
size(d)
size(dat)
minnans
nix(1)
minnans
nix(1)
minnans
dbquit
dbclear all
isodd(3)
isodd(2)
minnans
17*777/60
17*777/60/60
which interp2
checknans
checknan
tic; toc;
minnans
777/60
help tic
help toc
minnans
0.023*957
0.022*957
minnans
0.004*957
minnans
q
checknan
checknans
minnans
777/60
timenow
minnans
0.027*777
0.027*500
sum(0.027*[1:777])
sum(0.027*[1:777])/60
0.027*777
minnans
timenow
minnans
tic,
minnans
toc,
more off
tic,
minnans
toc,
more on
672/60
size(nix),size(tix)
fmg; contourf(lon(nix),lat(tix),dat(tix,nix)); daspect([1,cosd(t(1)),1]); colorbar; set_surf_cursor
fmg; contourf(lon(nix),lat(tix),dat(tix,nix)); daspect([1,cosd(t(1)),1]); colorbar; set_surf_cursor; contour(lon(nix),lat(tix),dat(tix,nix),[0,0],'Color','k','LineWidth',2);
daspect([1,cosd(lat(1)),1]);
edgeix = find_field_holes(dat,1,true,true);
lh = plot(LON(edgeix),LAT(edgeix),'k.');
delete(lh)
lh = plot(LON(edgeix),LAT(edgeix),'R.');
delete(lh)
edgeix = find_field_holes(dat,1,true,false);
lh = plot(LON(edgeix),LAT(edgeix),'R.');
delete(lh)
edgeix = find_field_holes(dat);
fmg; contourf(lon,lat,dat,-1000:50:100); colorbar;
contour(lon,lat,dat,[0,0],'Color','k','LineWidth',2);
%axis([145.74,145.8,15.12,15.17]);
titlename(['Bathy/DEM at ',num2str(target_res),' m (Laolao Bay)']);
% One arcsecond of latitude
target_res = 30;
fmg; contourf(lon,lat,dat,-1000:50:100); colorbar;
contour(lon,lat,dat,[0,0],'Color','k','LineWidth',2);
%axis([145.74,145.8,15.12,15.17]);
titlename(['Bathy/DEM at ',num2str(target_res),' m (Laolao Bay)']);
fmg; contourf(lon,lat,dat,-1000:50:100); colorbar;
contour(lon,lat,dat,[0,0],'Color','k','LineWidth',2);
daspect([1,cosd(lat(1)),1]);
titlename(['Bathy/DEM at ',num2str(target_res),' m wider bathy']);
edgeix = find_field_holes(d);
fmg; contourf(n,t,d); daspect([1,cosd(t(1)),1]); colorbar; set_surf_cursor
size(n),size(t),size(d)
size(oldd)
oldd=[]; clear oldd
nix=300:699; tix=20:294; n=lon(nix); t=lat(tix); d=dat(tix,nix); [N,T]=meshgrid(n,t); % JUST OK
fmg; contourf(n,t,d); daspect([1,cosd(t(1)),1]); colorbar; set_surf_cursor
nix=1:699; tix=20:294; n=lon(nix); t=lat(tix); d=dat(tix,nix); [N,T]=meshgrid(n,t); % JUST OK
nd=interp2(n,t,d,N(isnan(d)),T(isnan(d)),'spline',nan);
fmg; contourf(n,t,d); daspect([1,cosd(t(1)),1]); colorbar; set_surf_cursor
nix=652:699; tix=20:294; n=lon(nix); t=lat(tix); d=dat(tix,nix); [N,T]=meshgrid(n,t); % JUST OK
fmg; contourf(n,t,d); daspect([1,cosd(t(1)),1]); colorbar; set_surf_cursor
nix=652:699; tix=100:294; n=lon(nix); t=lat(tix); d=dat(tix,nix); [N,T]=meshgrid(n,t); % *DEEPER* BOX
nd=interp2(n,t,d,N(isnan(d)),T(isnan(d)),'spline',nan);
fmg; contourf(n,t,d); daspect([1,cosd(t(1)),1]); colorbar; set_surf_cursor
nix=633:699; tix=100:294; n=lon(nix); t=lat(tix); d=dat(tix,nix); [N,T]=meshgrid(n,t); % *DEEPER* BOX
fmg; contourf(n,t,d); daspect([1,cosd(t(1)),1]); colorbar; set_surf_cursor
nd=interp2(n,t,d,N(isnan(d)),T(isnan(d)),'spline',nan);
nix=500:699; tix=100:294; n=lon(nix); t=lat(tix); d=dat(tix,nix); [N,T]=meshgrid(n,t); % *DEEPER* BOX
nd=interp2(n,t,d,N(isnan(d)),T(isnan(d)),'spline',nan);
fmg; contourf(n,t,d); daspect([1,cosd(t(1)),1]); colorbar; set_surf_cursor
nix=633:799; tix=100:294; n=lon(nix); t=lat(tix); d=dat(tix,nix); [N,T]=meshgrid(n,t); % *DEEPER* BOX
nd=interp2(n,t,d,N(isnan(d)),T(isnan(d)),'spline',nan);
fmg; contourf(n,t,d); daspect([1,cosd(t(1)),1]); colorbar; set_surf_cursor
nd=interp2(n,t,d,N(isnan(d)),T(isnan(d)),'spline',nan);
fmg; contourf(n,t,nd); daspect([1,cosd(t(1)),1]); colorbar; set_surf_cursor
size(n),size(t),size(d),size(nd)
nd=d; nd(isnan(d))=interp2(n,t,d,N(isnan(d)),T(isnan(d)),'spline',nan);
size(n),size(t),size(d),size(nd)
fmg; contourf(n,t,nd); daspect([1,cosd(t(1)),1]); colorbar; set_surf_cursor
nd=d; nd(isnan(d))=interp2(n,t,d,N(isnan(d)),T(isnan(d)),'spline',nan);
fmg; contourf(n,t,nd); daspect([1,cosd(t(1)),1]); colorbar; set_surf_cursor
fmg; contourf(n,t,d); daspect([1,cosd(t(1)),1]); colorbar; set_surf_cursor
nix=633:799; tix=200:294; n=lon(nix); t=lat(tix); d=dat(tix,nix); [N,T]=meshgrid(n,t); % *DEEPER* BOX
fmg; contourf(n,t,d); daspect([1,cosd(t(1)),1]); colorbar; set_surf_cursor
nd=d; nd(isnan(d))=interp2(n,t,d,N(isnan(d)),T(isnan(d)),'spline',nan);
fmg; contourf(n,t,nd); daspect([1,cosd(t(1)),1]); colorbar; set_surf_cursor
daspect([1,cosd(lat(1)),1]);
nix=633:799; tix=200:294; n=lon(nix); t=lat(tix); d=dat(tix,nix); [N,T]=meshgrid(n,t); % *DEEPER* BOX
nd=d; nd(isnan(d))=interp2(n,t,d,N(isnan(d)),T(isnan(d)),'spline',nan);
fmg; contourf(n,t,nd); daspect([1,cosd(t(1)),1]); colorbar; set_surf_cursor
nix=633:799; tix=100:294; n=lon(nix); t=lat(tix); d=dat(tix,nix); [N,T]=meshgrid(n,t); % *DEEPER* BOX
nd=d; nd(isnan(d))=interp2(n,t,d,N(isnan(d)),T(isnan(d)),'spline',nan);
fmg; contourf(n,t,nd); daspect([1,cosd(t(1)),1]); colorbar; set_surf_cursor
nix=633:799; tix=100:294; n=lon(nix); t=lat(tix); d=dat(tix,nix); [N,T]=meshgrid(n,t); % *DEEPER* BOX
nd=d; nd(isnan(d))=interp2(n,t,d,N(isnan(d)),T(isnan(d)),'spline',nan);
fmg; contourf(n,t,d); daspect([1,cosd(t(1)),1]); colorbar; set_surf_cursor
fmg; contourf(n,t,nd); daspect([1,cosd(t(1)),1]); colorbar; set_surf_cursor
which anefkeys
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
anefkeys
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
dir('C:\Users\gramer\Documents\RSMAS\Coastal\thesis\src\data\hycom\eFKEYS_01-1-7-2012\eFKEYS_archv.2012_003_18_3zuvwts.nc
dir('C:\Users\gramer\Documents\RSMAS\Coastal\thesis\src\data\hycom\eFKEYS_01-1-7-2012\eFKEYS_archv.2012_003_18_3zuvwts.nc')
dir('C:\Users\gramer\Documents\RSMAS\Coastal\thesis\src\data\hycom\eFKEYS_01-1-7-2012\eFKEYS_archv.2012_003*')
dir('C:\Users\gramer\Documents\RSMAS\Coastal\thesis\src\data\hycom\eFKEYS_01-1-7-2012\eFKEYS*')
dir('C:\Users\gramer\Documents\RSMAS\Coastal\thesis\src\data\hycom\eFKEYS_01-1-7-2012\*')
dir('C:\Users\gramer\Documents\RSMAS\Coastal\thesis\src\data\hycom\*')
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
anefkeys
close(nc); clear nc
nc = mDataset(ncpath);
nj_info(nc),
close(nc); clear nc
anefkeys
close(nc); clear nc
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
anefkeys
nj_info(nc),
size(u)
close(nc); clear nc
dbquit
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
anefkeys
jd
hr
close(nc); clear nc
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
anefkeys
%-- 5/14/2016 9:01 PM --%
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
interp_dem_bathy
close all
nix=633:799; tix=200:294; n=lon(nix); t=lat(tix); d=dat(tix,nix); [N,T]=meshgrid(n,t); % *DEEPER* BOX
fmg; contourf(n,t,d); daspect([1,cosd(t(1)),1]); colorbar; set_surf_cursor
cm = colormap(r2b);
nd=d; nd(isnan(d))=interp2(n,t,d,N(isnan(d)),T(isnan(d)),'spline',nan);
fmg; contourf(n,t,nd); daspect([1,cosd(t(1)),1]); colorbar; set_surf_cursor
fmg; contourf(lon,lat,dat,-1000:50:100); colorbar;
contour(lon,lat,dat,[0,0],'Color','k','LineWidth',2);
daspect([1,cosd(lat(1)),1]);
titlename(['Bathy/DEM at ',num2str(target_res),' m wider bathy']);
nd=d; nd(isnan(d))=interp2(n,t,d,N(isnan(d)),T(isnan(d)),'linear',nan);
fmg; contourf(n,t,nd); daspect([1,cosd(t(1)),1]); colorbar; set_surf_cursor; titlename('Linear');
nix=633:799; tix=1:200; n=lon(nix); t=lat(tix); d=dat(tix,nix); [N,T]=meshgrid(n,t); % North reefs
fmg; contourf(n,t,nd); daspect([1,cosd(t(1)),1]); colorbar; set_surf_cursor
fmg; contourf(n,t,d); daspect([1,cosd(t(1)),1]); colorbar; set_surf_cursor
nix=633:799; tix=400:494; n=lon(nix); t=lat(tix); d=dat(tix,nix); [N,T]=meshgrid(n,t); % *DEEPER* BOX
fmg; contourf(n,t,d); daspect([1,cosd(t(1)),1]); colorbar; set_surf_cursor
nd=d; nd(isnan(d))=interp2(n,t,d,N(isnan(d)),T(isnan(d)),'spline',nan);
fmg; contourf(n,t,nd); daspect([1,cosd(t(1)),1]); colorbar; set_surf_cursor
fmg; contourf(n,t,nd); daspect([1,cosd(t(1)),1]); colorbar; set_surf_cursor; titlename('Spline');
nix=633:899; tix=300:594; n=lon(nix); t=lat(tix); d=dat(tix,nix); [N,T]=meshgrid(n,t); % *DEEPER* BOX
fmg; contourf(n,t,d); daspect([1,cosd(t(1)),1]); colorbar; set_surf_cursor
nd=d; nd(isnan(d))=interp2(n,t,d,N(isnan(d)),T(isnan(d)),'spline',nan);
fmg; contourf(n,t,nd); daspect([1,cosd(t(1)),1]); colorbar; set_surf_cursor; titlename('Spline');
nix=433:899; tix=200:694; n=lon(nix); t=lat(tix); d=dat(tix,nix); [N,T]=meshgrid(n,t); % *DEEPER* BOX
nd=d; nd(isnan(d))=interp2(n,t,d,N(isnan(d)),T(isnan(d)),'spline',nan);
nix=533:899; tix=200:694; n=lon(nix); t=lat(tix); d=dat(tix,nix); [N,T]=meshgrid(n,t); % *DEEPER* BOX
nd=d; nd(isnan(d))=interp2(n,t,d,N(isnan(d)),T(isnan(d)),'spline',nan);
nix=433:899; tix=200:694; n=lon(nix); t=lat(tix); d=dat(tix,nix); [N,T]=meshgrid(n,t); % *DEEPER* BOX
nd=d; nd(isnan(d))=interp2(n,t,d,N(isnan(d)),T(isnan(d)),'spline',nan);
nix=433:899; tix=300:694; n=lon(nix); t=lat(tix); d=dat(tix,nix); [N,T]=meshgrid(n,t); % *DEEPER* BOX
nd=d; nd(isnan(d))=interp2(n,t,d,N(isnan(d)),T(isnan(d)),'spline',nan);
fmg; contourf(n,t,d); daspect([1,cosd(t(1)),1]); colorbar; set_surf_cursor
fmg; contourf(n,t,nd); daspect([1,cosd(t(1)),1]); colorbar; set_surf_cursor; titlename('Spline');
nix=433:899; tix=300:694; n=lon(nix); t=lat(tix); d=dat(tix,nix); [N,T]=meshgrid(n,t); % *DEEPER* BOX
nd=d; nd(isnan(d))=interp2(n,t,d,N(isnan(d)),T(isnan(d)),'spline',nan);
fmg; contourf(n,t,nd); daspect([1,cosd(t(1)),1]); colorbar; set_surf_cursor; titlename('Spline');
71000/12
ans*1.26*1.51
ans*3
ans*2
fmg; contourf(lon,lat,dat,-1000:50:100); colorbar;
contour(lon,lat,dat,[0,0],'Color','k','LineWidth',2);
daspect([1,cosd(lat(1)),1]);
titlename(['Bathy/DEM at ',num2str(target_res),' m wider bathy']);
set_surf_cursor
flat
flon
lon
flon
s
bat=[]; clear bat batpath
dem=[]; clear dem dempath
datdat=[]; clear datdat xix yix
clear rawdat rawlat rawlon xixen xixix yixen yixix
clear s
LON=[]; LAT=[]; clear LON LAT
clear dx dy rawdat rawlat rawlon xixen xixix yixen yixix
d=[]; nd=[]; N=[]; T=[]; clear d nanix nix nx ny nd n t N T
d=[]; nd=[]; N=[]; T=[]; clear ans d nanix nix nx ny nd n t N T
d=[]; nd=[]; N=[]; T=[]; clear ans d nanix nix tix nx ny nd n t N T
pack
fill_dem_bathy
fdat=[];
clear Top Tops Right Rights Bottom Bottoms Left Lefts fdat flat flon ix nixen tixen
clear s
fill_dem_bathy
rawdat
fill_dem_bathy
size(nanix)
ix
nixen
tixen
fill_dem_bathy
fmg; contourf(lon,lat,newdat,-1000:50:100); colorbar;
contour(lon,lat,newdat,[0,0],'Color','k','LineWidth',2);
daspect([1,cosd(lat(1)),1]); set_surf_cursor;
titlename(['Bathy/DEM at ',num2str(target_res),' m AFTER SPLINE']);
edgeix = find_field_holes(rawdat);
[LON,LAT] = meshgrid(lon,lat);
lh = plot(LON(edgeix),LAT(edgeix),'r.','MarkerSize',0.5);
LON=[]; LAT=[]; clear LON LAT
clear ans lh
newdat=[]; clear newdat
fill_dem_bathy
newdat=[]; clear newdat
fill_dem_bathy
set_surf_cursor
fmg; contourf(lon,lat,newdat,-1000:50:100); colorbar;
contour(lon,lat,newdat,[0,0],'Color','k','LineWidth',2);
daspect([1,cosd(lat(1)),1]); set_surf_cursor;
titlename(['Bathy/DEM at ',num2str(target_res),' m AFTER SPLINE']);
edgeix = find_field_holes(rawdat);
[LON,LAT] = meshgrid(lon,lat);
lh=plot(LON(edgeix),LAT(edgeix),'r.','MarkerSize',0.5);
LON=[]; LAT=[]; clear LON LAT
plot([Lefts,Lefts,Rights,Rights],[Tops,Bottoms,Tops,Bottoms],'kp','MarkerFaceColor','w');
Tops = [653,673,761,653,394,];
Rights = [779,926,888,913,116,];
Bottoms = [255,626,691,611,375,];
Lefts = [480,872,816,869, 99,];
plot([Lefts,Lefts,Rights,Rights],[Tops,Bottoms,Tops,Bottoms],'kp','MarkerFaceColor','w');
fmg; contourf(lon,lat,newdat,-1000:50:100); colorbar;
contour(lon,lat,newdat,[0,0],'Color','k','LineWidth',2);
daspect([1,cosd(lat(1)),1]); set_surf_cursor;
titlename(['Bathy/DEM at ',num2str(target_res),' m AFTER SPLINE']);
edgeix = find_field_holes(rawdat);
[LON,LAT] = meshgrid(lon,lat);
lh=plot(LON(edgeix),LAT(edgeix),'r.','MarkerSize',0.5);
LON=[]; LAT=[]; clear LON LAT
plot([Lefts,Lefts,Rights,Rights],[Tops,Bottoms,Tops,Bottoms],'kp','MarkerFaceColor','w');
Lefts
fmg; contourf(lon,lat,newdat,-1000:50:100); colorbar;
contour(lon,lat,newdat,[0,0],'Color','k','LineWidth',2);
daspect([1,cosd(lat(1)),1]); set_surf_cursor;
titlename(['Bathy/DEM at ',num2str(target_res),' m AFTER SPLINE']);
edgeix = find_field_holes(rawdat);
[LON,LAT] = meshgrid(lon,lat);
lh=plot(LON(edgeix),LAT(edgeix),'r.','MarkerSize',0.5);
LON=[]; LAT=[]; clear LON LAT
plot(lon([Lefts,Lefts,Rights,Rights]),lat([Tops,Bottoms,Tops,Bottoms]),'kp','MarkerFaceColor','w');
fmg; contourf(lon,lat,newdat,-1000:50:100); colorbar;
contour(lon,lat,newdat,[0,0],'Color','k','LineWidth',2);
daspect([1,cosd(lat(1)),1]); set_surf_cursor;
titlename(['Bathy/DEM at ',num2str(target_res),' m AFTER SPLINE']);
edgeix = find_field_holes(rawdat);
[LON,LAT] = meshgrid(lon,lat);
lh=plot(LON(edgeix),LAT(edgeix),'r.','MarkerSize',0.5);
LON=[]; LAT=[]; clear LON LAT
plot(lon([Lefts,Lefts,Rights,Rights]),lat([Tops,Bottoms,Tops,Bottoms]),'k-p','MarkerFaceColor','w');
fmg; contourf(lon,lat,newdat,-1000:50:100); colorbar;
contour(lon,lat,newdat,[0,0],'Color','k','LineWidth',2);
daspect([1,cosd(lat(1)),1]); set_surf_cursor;
titlename(['Bathy/DEM at ',num2str(target_res),' m AFTER SPLINE']);
edgeix = find_field_holes(rawdat);
[LON,LAT] = meshgrid(lon,lat);
lh=plot(LON(edgeix),LAT(edgeix),'r.','MarkerSize',0.5);
LON=[]; LAT=[]; clear LON LAT
%plot(lon([Lefts,Lefts,Rights,Rights]),lat([Tops,Bottoms,Tops,Bottoms]),'k-p','MarkerFaceColor','w');
for ix=1:numel(Tops)
Top = Tops(ix);
Right = Rights(ix);
Bottom = Bottoms(ix);
Left = Lefts(ix);
plot(lon([Left,Left,Right,Right]),lat([Top,Bottom,Top,Bottom]),'k-p','MarkerFaceColor','w');
end;
fmg; contourf(lon,lat,newdat,-1000:50:100); colorbar;
contour(lon,lat,newdat,[0,0],'Color','k','LineWidth',2);
daspect([1,cosd(lat(1)),1]); set_surf_cursor;
titlename(['Bathy/DEM at ',num2str(target_res),' m AFTER SPLINE']);
edgeix = find_field_holes(rawdat);
[LON,LAT] = meshgrid(lon,lat);
lh=plot(LON(edgeix),LAT(edgeix),'r.','MarkerSize',0.5);
LON=[]; LAT=[]; clear LON LAT
%plot(lon([Lefts,Lefts,Rights,Rights]),lat([Tops,Bottoms,Tops,Bottoms]),'k-p','MarkerFaceColor','w');
for ix=1:numel(Tops)
Top = Tops(ix);
Right = Rights(ix);
Bottom = Bottoms(ix);
Left = Lefts(ix);
plot(lon([Left,Right,Right,Left,Left]),lat([Top,Top,Bottom,Bottom,Top]),'k-p','MarkerFaceColor','w');
end;
fmg; contourf(lon,lat,newdat,-40:2:0); colorbar;
%fmg; contourf(lon,lat,newdat,-1000:50:100); colorbar;
contour(lon,lat,newdat,[0,0],'Color','k','LineWidth',2);
daspect([1,cosd(lat(1)),1]); set_surf_cursor;
titlename(['Bathy/DEM at ',num2str(target_res),' m AFTER SPLINE']);
edgeix = find_field_holes(rawdat);
[LON,LAT] = meshgrid(lon,lat);
lh=plot(LON(edgeix),LAT(edgeix),'r.','MarkerSize',0.5);
LON=[]; LAT=[]; clear LON LAT
%plot(lon([Lefts,Lefts,Rights,Rights]),lat([Tops,Bottoms,Tops,Bottoms]),'k-p','MarkerFaceColor','w');
for ix=1:numel(Tops)
Top = Tops(ix);
Right = Rights(ix);
Bottom = Bottoms(ix);
Left = Lefts(ix);
plot(lon([Left,Right,Right,Left,Left]),lat([Top,Top,Bottom,Bottom,Top]),'k-p','MarkerFaceColor','w');
end;
fmg; contourf(lon,lat,newdat,-40:2:4); colorbar;
%fmg; contourf(lon,lat,newdat,-1000:50:100); colorbar;
contour(lon,lat,newdat,[0,0],'Color','k','LineWidth',2);
daspect([1,cosd(lat(1)),1]); set_surf_cursor;
titlename(['Bathy/DEM at ',num2str(target_res),' m AFTER SPLINE']);
edgeix = find_field_holes(rawdat);
[LON,LAT] = meshgrid(lon,lat);
lh=plot(LON(edgeix),LAT(edgeix),'r.','MarkerSize',0.5);
LON=[]; LAT=[]; clear LON LAT
%plot(lon([Lefts,Lefts,Rights,Rights]),lat([Tops,Bottoms,Tops,Bottoms]),'k-p','MarkerFaceColor','w');
for ix=1:numel(Tops)
Top = Tops(ix);
Right = Rights(ix);
Bottom = Bottoms(ix);
Left = Lefts(ix);
plot(lon([Left,Right,Right,Left,Left]),lat([Top,Top,Bottom,Bottom,Top]),'k-p','MarkerFaceColor','w');
end;
fmg; contourf(lon,lat,newdat,-40:2:6); colorbar;
%fmg; contourf(lon,lat,newdat,-1000:50:100); colorbar;
contour(lon,lat,newdat,[0,0],'Color','k','LineWidth',2);
daspect([1,cosd(lat(1)),1]); set_surf_cursor;
titlename(['Bathy/DEM at ',num2str(target_res),' m AFTER SPLINE']);
edgeix = find_field_holes(rawdat);
[LON,LAT] = meshgrid(lon,lat);
lh=plot(LON(edgeix),LAT(edgeix),'r.','MarkerSize',0.5);
LON=[]; LAT=[]; clear LON LAT
%plot(lon([Lefts,Lefts,Rights,Rights]),lat([Tops,Bottoms,Tops,Bottoms]),'k-p','MarkerFaceColor','w');
for ix=1:numel(Tops)
Top = Tops(ix);
Right = Rights(ix);
Bottom = Bottoms(ix);
Left = Lefts(ix);
plot(lon([Left,Right,Right,Left,Left]),lat([Top,Top,Bottom,Bottom,Top]),'k-p','MarkerFaceColor','w');
end;
doSave
clear Top Tops Right Rights Bottom Bottoms Left Lefts ans ix nanix nixen tixen
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
interp_dem_bathy
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
interp_dem_bathy
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
interp_dem_bathy
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
interp_dem_bathy
fill_dem_bathy
fmg; contourf(lon,lat,dat,-1000:50:100); colorbar;
contour(lon,lat,dat,[0,0],'Color','k','LineWidth',2);
daspect([1,cosd(lat(1)),1]); set_surf_cursor;
titlename(['Bathy/DEM at ',num2str(target_res),' m wider bathy']);
fmg; contourf(lon,lat,newdat,-1000:50:100); colorbar;
contour(lon,lat,newdat,[0,0],'Color','k','LineWidth',2);
daspect([1,cosd(lat(1)),1]); set_surf_cursor;
titlename(['Bathy/DEM at ',num2str(target_res),' m AFTER SPLINE']);
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
interp_dem_bathy
fill_dem_bathy
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
interp_dem_bathy
help interp2
fill_dem_bathy
ix
nansummary(fdat(nanix))
fmg; contourf(flon,flat,fdat);
fmg; contourf(lon(nixen),lat(tixen),fdat(tixen,nixen));
fmg; contourf(lon(nixen),lat(tixen),dat(tixen,nixen));
set_surf_cursor
dbcont
timenow
fill_dem_bathy
ix
size( find(FLON(yix,xix-1)<0 & FLON(yix,xix+1)<0) )
min(xix),max(xix)
size(nanlandix)
nansummary(nanlandix)
nansummary(yix)
nansummary(xix)
help ind2sub
dbquit
FLON=[]; FLAT=[]; clear FLON FLAT
fdat=[]; clear fdat flat flon
clear Top Tops Right Rights Bottom Bottoms Left Lefts ans ix nanix nixen tixen
clear xix yix
fill_dem_bathy
nansummary(yix)
nansummary(xix)
size( find(FLON(yix,xix-1)<0 & FLON(yix,xix+1)<0) )
badlandix = find(dat(yix,xix-1)<0 & fdat(yix,xix+1)<0);
dbquit
clear Top Tops Right Rights Bottom Bottoms Left Lefts ans ix nanix nixen tixen
clear xix yix
fdat=[]; clear fdat flat flon
FLON=[]; FLAT=[]; clear FLON FLAT
nansummary(badlandix)
badlandix'
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
interp_dem_bathy
fill_dem_bathy
fmg; contourf(flon,flat,fdat,-1000:50:100); colorbar;
nansummary(badlandix)
nansummary(previx)
nansummary(nextix)
badlandix = previx(badlandix);
nansummary(badlandix)
size(badlandix),size(fdat(previx))
size(badlandix),size(fdat(previx(badlandix)))
dbquit
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
interp_dem_bathy
fill_dem_bathy
fmg; contourf(flon,flat,fdat,-1000:50:100); colorbar;
size(fdat(previx(badlandix))),size(mean([fdat(previx(badlandix));fdat(nextix(badlandix))]))
size(fdat(previx(badlandix))),size(mean([fdat(previx(badlandix)),fdat(nextix(badlandix))]))
size(fdat(previx(badlandix))),size(mean([fdat(previx(badlandix)),fdat(nextix(badlandix))]'))
size(fdat(previx(badlandix))),size(mean([fdat(previx(badlandix)),fdat(nextix(badlandix))],2))
fdat(previx(badlandix)) = mean([fdat(previx(badlandix)),fdat(nextix(badlandix))],2);
fmg; contourf(flon,flat,fdat,-1000:50:100); colorbar;
previx(badlandix)
mean([fdat(previx(badlandix)),fdat(nextix(badlandix))],2);
mean([fdat(previx(badlandix)),fdat(nextix(badlandix))],2)
plot(flon(previx(badlandix)),flat(previx(badlandix)),'rp')
plot(FLON(previx(badlandix)),FLAT(previx(badlandix)),'rp')
dbquit
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
interp_dem_bathy
fill_dem_bathy
nanlandix(badlandix)
fdat(nanlandix(badlandix))
dbquit
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
interp_dem_bathy
fill_dem_bathy
timenow
fill_dem_bathy
newdat=[]; clear newdat
fill_dem_bathy
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
interp_dem_bathy
fill_dem_bathy
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
interp_dem_bathy
fill_dem_bathy
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
interp_dem_bathy
delete pibhc_bathy_5m_tiniain.dods.mat
dir('Saving C:\Users\gramer\Documents\MATLAB\ecoforecasts\coast\pibhmc_bathy_5m_tinian.dods.mat
dir('C:\Users\gramer\Documents\MATLAB\ecoforecasts\coast\pibhmc_bathy_5m_tinian.dods.mat');
delete('C:\Users\gramer\Documents\MATLAB\ecoforecasts\coast\pibhmc_bathy_5m_tinian.dods.mat');
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
quit
%-- 5/18/2016 9:51 AM --%
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
interp_dem_bathy
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
quit
%-- 5/18/2016 10:04 AM --%
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
demnm = 'usgs_dem_10m_tinian';
batnm = 'pibhmc_bathy_5m_tinian';
batpath = fullfile(get_ecoforecasts_path('coast'),[batnm,'.dods.mat']);
if ( ~exist(batpath,'file') )
disp(['Extracting ',batpath]);
url = ['http://oos.soest.hawaii.edu/thredds/dodsC/',batnm];
nc = mDataset(url);
end;
delete('C:\Users\gramer\Documents\MATLAB\ecoforecasts\coast\pibhmc_bathy_5m_tinian.dods.mat');
if ( ~exist(batpath,'file') )
disp(['Extracting ',batpath]);
url = ['http://oos.soest.hawaii.edu/thredds/dodsC/',batnm];
nc = mDataset(url);
end;
nj_info(nc),
batpath
nj_info(nc),
dat = cast(nc{'elev'}(:,:),'double');
close(nc); clear nc
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
quit
%-- 5/18/2016 10:14 AM --%
batnm = 'pibhmc_bathy_5m_tinian';
batpath = fullfile(get_ecoforecasts_path('coast'),[batnm,'.dods.mat']);
url = ['http://oos.soest.hawaii.edu/thredds/dodsC/',batnm];
nc = mDataset(url);
dat = nc{'elev'}(:,:);
close(nc); clear nc
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
quit
%-- 5/18/2016 10:57 AM --%
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
batnm = 'pibhmc_bathy_5m_tinian';
batpath = fullfile(get_ecoforecasts_path('coast'),[batnm,'.dods.mat']);
url = ['http://oos.soest.hawaii.edu/thredds/dodsC/',batnm];
nc = mDataset(url);
dat = nc{'elev'}(:,:);
close(nc); clear nc
pack
nc = mDataset(url);
dat = nc{'elev'}(1000:end-1000,1000:end-1000);
lon = cast(nc{'lon'}(1000:end-1000),'double');
lat = cast(nc{'lat'}(1000:end-1000),'double');
close(nc); clear nc
fmg; contourf(lon,lat,dat,-1000:50:100); colorbar;
save('testing.mat')
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
quit
%-- 5/18/2016 11:08 AM --%
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
load('testing.mat')
load('testing.mat');
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
batnm = 'pibhmc_bathy_5m_tinian';
batpath = fullfile(get_ecoforecasts_path('coast'),[batnm,'.dods.mat']);
load('testing.mat')
url = ['http://oos.soest.hawaii.edu/thredds/dodsC/',batnm];
nc = mDataset(url);
lon = cast(nc{'lon'}(1000:end-1000),'double');
lat = cast(nc{'lat'}(1000:end-1000),'double');
dat = cast(nc{'elev'}(1000:end-1000,1000:end-1000),'double');
close(nc); clear nc
save('tinian-test.mat')
fmg; contourf(lon,lat,dat,[0,0]); colorbar;
fmg; contourf(lon,lat,dat,[-5,-5]); colorbar;
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
batnm = 'pibhmc_bathy_5m_tinian';
batpath = fullfile(get_ecoforecasts_path('coast'),[batnm,'.dods.mat']);
batnm = 'pibhmc_bathy_5m_tinian'; batpath = fullfile(get_ecoforecasts_path('coast'),[batnm,'.dods.mat']); url = ['http://oos.soest.hawaii.edu/thredds/dodsC/',batnm];
nc = mDataset(url);
lon = cast(nc{'lon'}(1000:end-1000),'double');
lat = cast(nc{'lat'}(1:end-1000),'double');
nj_info(nc),
dat = cast(nc{'elev'}(1:end-1000,1000:end-1000),'double');
close(nc); clear nc
clear ans
save('tinian-test.mat')
pack
disp('Saved');
fmg; contourf(lon,lat,dat,[-5,-5]); colorbar;
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
batnm = 'pibhmc_bathy_5m_tinian'; batpath = fullfile(get_ecoforecasts_path('coast'),[batnm,'.dods.mat']); url = ['http://oos.soest.hawaii.edu/thredds/dodsC/',batnm];
nc = mDataset(url);
lon = cast(nc{'lon'}(1000:end-1000),'double');
lat = cast(nc{'lat'}(1000:end),'double');
dat = cast(nc{'elev'}(1000:end,1000:end-1000),'double');
close(nc); clear nc
clear ans
save('tinian-test.mat')
pack
disp('Saved');
fmg; contourf(lon,lat,dat,[-5,-5]); colorbar;
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
batnm = 'pibhmc_bathy_5m_tinian'; batpath = fullfile(get_ecoforecasts_path('coast'),[batnm,'.dods.mat']); url = ['http://oos.soest.hawaii.edu/thredds/dodsC/',batnm];
nc = mDataset(url);
lon = cast(nc{'lon'}(1000:end-1000),'double');
lat = cast(nc{'lat'}(1:end),'double');
dat = cast(nc{'elev'}(1:end,1000:end-1000),'double');
close(nc); clear nc
pack
save('tinian-big-test.mat')
disp('Saved');
fmg; contourf(lon,lat,dat,[-5,-5]); colorbar;
disp(['Saving ',batpath]);
batnm = 'pibhmc_bathy_5m_tinian'; batpath = fullfile(get_ecoforecasts_path('coast'),[batnm,'.PARTIAL.dods.mat']); url = ['http://oos.soest.hawaii.edu/thredds/dodsC/',batnm];
disp(['Saving ',batpath]);
save(batpath,'url','lon','lat','dat','-v7.3');
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
quit
%-- 5/18/2016 11:46 AM --%
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
interp_dem_bathy
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
batnm = 'pibhmc_bathy_5m_tinian'; batpath = fullfile(get_ecoforecasts_path('coast'),[batnm,'.PARTIAL.dods.mat']); url = ['http://oos.soest.hawaii.edu/thredds/dodsC/',batnm];
nc = mDataset(url);
dat = cast(nc{'elev'}(1:end,500:end-500),'double');
lon = cast(nc{'lon'}(500:end-500),'double');
lat = cast(nc{'lat'}(1:end),'double');
close(nc); clear nc
pack
disp(['Saving ',batpath]);
save(batpath,'url','lon','lat','dat','-v7.3');
fmg; contourf(lon,lat,dat,[-5,-5]); colorbar;
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
quit
%-- 5/18/2016 11:55 AM --%
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
interp_dem_bathy
corners = [ ...
480,255 , 779,255 , 779,653 , 480,653 ; ...
872,626 , 926,626 , 926,673 , 872,673 ; ...
816,691 , 888,691 , 888,761 , 816,761 ; ...
869,611 , 913,611 , 913,653 , 869,653 ; ...
99 ,375 , 116,375 , 116,394 , 99 ,394 ; ...
357,559 , 436,559 , 436,604 , 357,604 ; ...
291,145 , 362,145 , 362,223 , 291,223 ; ...
];
size(corners)
099
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
interp_dem_bathy
fill_dem_bathy
Tops
Rights
Bottoms
Lefts
ix
newdat=[]; clear newdat
FLON=[]; FLAT=[]; clear FLON FLAT
fdat=[]; clear fdat flat flon
LON=[]; LAT=[]; clear LON LAT
clear Top Tops Right Rights Bottom Bottoms Left Lefts ans ix nanix nixen tixen
clear s
clear corners rix
clear nanix nanlandix yix xix previx nextix badlandix
rawdat=[]; clear rawdat
pack
fill_dem_bathy
rawdat=[]; clear rawdat
clear edgeix lh nanlandix
newdat=[]; clear newdat
pack
fill_dem_bathy
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
interp_dem_bathy
fill_dem_bathy
Tops
Rights
Bottoms
Lefts
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
interp_dem_bathy
fill_dem_bathy
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
interp_dem_bathy
fill_dem_bathy
clear edgeix lh nanlandix
newdat=[]; clear newdat
fdat=[]; clear fdat flat flon
rawdat=[]; clear rawdat
clear nanix nanlandix yix xix previx nextix badlandix
clear corners rix
clear s
clear Top Tops Right Rights Bottom Bottoms Left Lefts ans ix nanix nixen tixen
fil
fill_dem_bathy
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
interp_dem_bathy
fill_dem_bathy
rawdat=[]; clear rawdat
newdat=[]; clear newdat
fdat=[]; clear fdat flat flon
clear Top Tops Right Rights Bottom Bottoms Left Lefts ans ix nanix nixen tixen
clear nanix nanlandix yix xix previx nextix badlandix
clear edgeix lh nanlandix
fill_dem_bathy
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
interp_dem_bathy
fill_dem_bathy
fmg; text(1,1,' \leftarrow foo');
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
interp_dem_bathy
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
quit
interp_dem_bathy
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
batnm = 'pibhmc_bathy_60m_rota'; batpath = fullfile(get_ecoforecasts_path('coast'),[batnm,'.dods.mat']); url = ['http://oos.soest.hawaii.edu/thredds/dodsC/',batnm];
nc = mDataset(url);
nj_info(nc),
lat = cast(nc{'lat'}(:),'double');
lon = cast(nc{'lon'}(:),'double');
dat = cast(nc{'elev'}(:,:),'double');
close(nc); clear nc
fmg; contourf(lon,lat,dat,[-5,-5]); colorbar;
fmg; contourf(lon,lat,dat,[0,0]); colorbar;
fmg; contourf(lon,lat,dat); colorbar;
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
interp_dem_bathy
help rectangle
help ginput
help gbm
timenow
help gruler
help gline
lookfor rectangle
help dragrect
help impositionrect
help imrect
help getrect
r = rect
r = getrect
help rectangle
rectangle(r);
help rectangle
rectangle('Position',r);
r1 = getrect; rectangle('Position',r1);
r2 = getrect; rectangle('Position',r2);
r3 = getrect; rectangle('Position',r3);
r4 = getrect; rectangle('Position',r4);
dat=[]; lat=[]; lon=[]; clear dat lat lon ans
clear region target_res url demurl demnm cmbpath batnm
clear region target_res url demurl demnm cmbpath batnm doSave
save('rota-rectangles.mat');
dir('rota-rectangles.mat');
x=dir('rota-rectangles.mat')
clear x
interp_dem_bathy
r
for cix=1:4; [xerr,xix]=min(abs(lon-r(cix)
rect2bbox(r)
b=rect2bbox(r)
clear b
help bbox2rect
help err
help erf
ix(3)=pi
clear ix
help plaid
help pid
help meshgrid
help slice
help surf
isvector([pi])
isvector([pi,3])
isvector([pi,3;1,3])
clear ans
b
help ind2sub
help meshgrid
ndims(dat)
help squeeze
help max
b=rect2bbox(r); ix=bbox2ind(b,lon,lat)
timenow
whos('rota-rectangles.mat')
x = load('rota-rectangels.mat')
x = load('rota-rectangles.mat')
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
quit
%-- 5/18/2016 4:37 PM --%
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
load('rota-rectangles.mat')
b=rect2bbox(r); ix=bbox2ind(b,lon,lat)
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
interp_dem_bathy
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
interp_dem_bathy
load('rota-rectangles.mat')
b=rect2bbox(r); ix=bbox2ind(b,lon,lat)
rehash
b=rect2bbox(r); ix=bbox2ind(b,lon,lat)
[LON,LAT] = meshgrid(lon,lat);
b=rect2bbox(r); ix=bbox2ind(b,LON,LAT)
fclose('all'); close all; clear all; clear classes; clear java; dbclear all; more on; pack
interp_dem_bathy
type cvtrects.m
cvtrects
Lefts
Rights
ig
cvtrects
Tops
fmg; contourf(lon,lat,newdat,-1000:50:100); colorbar;
contour(lon,lat,newdat,[0,0],'Color','k','LineWidth',2);
daspect([1,cosd(lat(1)),1]); set_surf_cursor;
titlename(['Bathy/DEM at ',num2str(target_res),' m AFTER SPLINE']);
fmg; contourf(lon,lat,dat,-1000:50:100); colorbar;
for ix=1:numel(Tops)
Top = Tops(ix);
Right = Rights(ix);
Bottom = Bottoms(ix);
Left = Lefts(ix);
plot(lon([Left,Right,Right,Left,Left]),lat([Top,Top,Bottom,Bottom,Top]),'k-p','MarkerFaceColor','w');
text(lon(Right),lat(Top),[' \larrow ',num2str(ix)]);
end;
ix
Tops
Rights
Bottoms
Lefts
ixes
cvtrects
fill_dem_bathy
ixen
tixen
Bottom,Top
size(lat)
clear s
clear nanix nanlandix yix xix previx nextix badlandix nbadland
FLON=[]; FLAT=[]; clear FLON FLAT
fdat=[]; clear fdat flat flon
clear Top Tops Right Rights Bottom Bottoms Left Lefts ans edgeix ix nanix nixen tixen
newdat=[]; clear newdat
clear r r1 r2 r3 r4 r5
clear rs
size(lon)
size(ixes)
Left,Right
fill_dem_bathy
Left,Right
Rights
Lefts
ixes
Tops
Bottoms
clear ans s nanix nanlandix yix xix previx nextix badlandix nbadland r r1 r2 r3 r4 r5 rs Top Tops Right Rights Bottom Bottoms Left Lefts ans edgeix ix nanix nixen tixen
clear ans s nanix nanlandix yix xix previx nextix badlandix nbadland b r r1 r2 r3 r4 r5 rs Top Tops Right Rights Bottom Bottoms Left Lefts ans edgeix ix nanix nixen tixen
LAT=[]; LON=[]; newdat=[]; rawdat=[]; clear LAT LON newdat rawdat flon flat
clear ans ig s nanix nanlandix yix xix previx nextix badlandix nbadland b r r1 r2 r3 r4 r5 rs Top Tops Right Rights Bottom Bottoms Left Lefts ans edgeix ix nanix nixen tixen
clear ans ig s nanix nanlandix yix xix previx nextix badlandix nbadland b r r1 r2 r3 r4 r5 rs ixes Top Tops Right Rights Bottom Bottoms Left Lefts ans edgeix ix nanix nixen tixen
clear ans ig s nanix nanlandix yix xix previx nextix badlandix nbadland b r r1 r2 r3 r4 r5 rs ixes Top Tops Right Rights Bottom Bottoms Left Lefts ans edgeix ix nanix nixen tixen; LAT=[]; LON=[]; newdat=[]; rawdat=[]; clear LAT LON newdat rawdat flon flat
cvtrects
lon(173),lat(147)
lon(173),lat(549)
clear ans ig s nanix nanlandix yix xix previx nextix badlandix nbadland b r r1 r2 r3 r4 r5 rs ixes Top Tops Right Rights Bottom Bottoms Left Lefts ans edgeix ix nanix nixen tixen; LAT=[]; LON=[]; newdat=[]; rawdat=[]; clear LAT LON newdat rawdat flon flat
cvtrects
help error
size(LON)
size(lon)
clear ans ig s nanix nanlandix yix xix previx nextix badlandix nbadland b r r1 r2 r3 r4 r5 rs ixes Top Tops Right Rights Bottom Bottoms Left Lefts ans edgeix ix nanix nixen tixen; LAT=[]; LON=[]; newdat=[]; rawdat=[]; clear LAT LON newdat rawdat flon flat
cvtrects
clear ans ig s nanix nanlandix yix xix previx nextix badlandix nbadland b r r1 r2 r3 r4 r5 rs ixes Top Tops Right Rights Bottom Bottoms Left Lefts ans edgeix ix nanix nixen tixen; LAT=[]; LON=[]; newdat=[]; rawdat=[]; clear LAT LON newdat rawdat flon flat
fill_dem_bathy
min(yix)
yix
min(yix)
clear ans ig s nanix nanlandix yix xix previx nextix badlandix nbadland b r r1 r2 r3 r4 r5 rs ixes Top Tops Right Rights Bottom Bottoms Left Lefts ans edgeix ix nanix nixen tixen; LAT=[]; LON=[]; newdat=[]; rawdat=[]; clear LAT LON newdat rawdat flon flat
clear ans ig s nanix nanlandix yix xix previx nextix badlandix nbadland b r r1 r2 r3 r4 r5 rs ixes Top Tops Right Rights Bottom Bottoms Left Lefts ans edgeix ix nanix nixen tixen; LAT=[]; LON=[]; FLAT=[]; FLON=[]; fdat=[]; newdat=[]; rawdat=[]; clear LAT LON FLAT FLON fdat newdat rawdat flon flat
pack
fill_dem_bathy
size(yix)
nansummary(yix)
size(lat)
size(nanlandix)
size(fdat)
nansummary(xix)
fmg; contourf(FLON,FLAT,fdat); colorbar;
fmg; contourf(FLON,FLAT,fdat); colorbar; daspect([1,cosd(FLAT(1)),1]);
fmg; contourf(lon,lat,newdat,-40:2:6); colorbar;
%fmg; contourf(lon,lat,newdat,-1000:50:100); colorbar;
contour(lon,lat,newdat,[0,0],'Color','k','LineWidth',2);
daspect([1,cosd(lat(1)),1]); set_surf_cursor;
titlename(['Bathy/DEM at ',num2str(target_res),' m AFTER SPLINE']);
edgeix = find_field_holes(rawdat);
[LON,LAT] = meshgrid(lon,lat);
plot(LON(edgeix),LAT(edgeix),'r.','MarkerSize',0.5);
LON=[]; LAT=[]; clear LON LAT
%plot(lon([Lefts,Lefts,Rights,Rights]),lat([Tops,Bottoms,Tops,Bottoms]),'k-p','MarkerFaceColor','w');
for ix=1:numel(Tops)
Top = Tops(ix);
Right = Rights(ix);
Bottom = Bottoms(ix);
Left = Lefts(ix);
plot(lon([Left,Right,Right,Left,Left]),lat([Top,Top,Bottom,Bottom,Top]),'k-p','MarkerFaceColor','w');
text(lon(Right),lat(Top),[' \larrow ',num2str(ix)]);
end;
fmg; contourf(lon,lat,newdat,-40:2:6); colorbar;
%fmg; contourf(lon,lat,newdat,-1000:50:100); colorbar;
contour(lon,lat,newdat,[0,0],'Color','k','LineWidth',2);
daspect([1,cosd(lat(1)),1]); set_surf_cursor;
titlename(['Bathy/DEM at ',num2str(target_res),' m AFTER SPLINE']);
edgeix = find_field_holes(rawdat);
[LON,LAT] = meshgrid(lon,lat);
plot(LON(edgeix),LAT(edgeix),'r.','MarkerSize',0.5);
LON=[]; LAT=[]; clear LON LAT
%plot(lon([Lefts,Lefts,Rights,Rights]),lat([Tops,Bottoms,Tops,Bottoms]),'k-p','MarkerFaceColor','w');
for ix=1:numel(Tops)
Top = Tops(ix);
Right = Rights(ix);
Bottom = Bottoms(ix);
Left = Lefts(ix);
plot(lon([Left,Right,Right,Left,Left]),lat([Top,Top,Bottom,Bottom,Top]),'k-p','MarkerFaceColor','w');
text(lon(Right),lat(Top),[' \leftarrow ',num2str(ix)]);
end;
ix
clear ans ig s nanix nanlandix yix xix previx nextix badlandix nbadland b r r1 r2 r3 r4 r5 rs ixes Top Tops Right Rights Bottom Bottoms Left Lefts ans edgeix ix nanix nixen tixen; LAT=[]; LON=[]; FLAT=[]; FLON=[]; fdat=[]; newdat=[]; rawdat=[]; clear LAT LON FLAT FLON fdat newdat rawdat flon flat
fill_dem_bathy
ix
fmg; contourf(FLON,FLAT,fdat); colorbar; daspect([1,cosd(FLAT(1)),1]);
stn
stn = get_station_from_station_name('mlrf1'); stn = load_all_ndbc_data(stn);
stn = plot_hires_bathymetry(stn,-[0:5:80]);
stn=[]; clear stn
clear ans ig s nanix nanlandix yix xix previx nextix badlandix nbadland b r r1 r2 r3 r4 r5 rs ixes Top Tops Right Rights Bottom Bottoms Left Lefts ans edgeix ix nanix nixen tixen; LAT=[]; LON=[]; FLAT=[]; FLON=[]; fdat=[]; newdat=[]; rawdat=[]; clear LAT LON FLAT FLON fdat newdat rawdat flon flat
7*24/6
7*24
dts = datenum(2012,1,1,[6:6:(7*24)],0,0);
datestr(dts)
clear dts
dts = datenum(2012,1,1,[6:6:(7*24)-eps],0,0);
datestr(dts)
dts = datenum(2012,1,1,[6:6:(7*24)-1],0,0);
datestr(dts)
clear dts
