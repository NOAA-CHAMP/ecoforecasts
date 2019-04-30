1;
%%%% SCRIPT to analyze American Samoa data


set_more off;

figspath = get_ecoforecasts_path('figs');

% ix = 1; %Airport lagoon
% ix = 2; %Alofau deep site
% ix = 3; %Alofau shallow site

for ix=1:2

if ( ~exist('stns','var') || isempty(stns) )
  stns = load_misst_region_stations('asam');
end;
if ( ~isfield(stns,'sea_t') )
  stns = load_clean_as_thermistor(stns);
end;

if ( ~isfield(stns,'oisst2_sst') )
  stns = read_oisst2(stns);
end;
if ( ~isfield(stns,'pfv5c_pentad_sst') )
  stns = read_pathfinder_pentad(stns);
end;
if ( ~isfield(stns,'pfv5c_cortad_sst') )
  stns = read_pathfinder_cortad(stns);
end;

if ( ~isfield(stns,'bleaching') )
  stns = load_bleach_data(stns);
end;

disp([stns(ix).station_name ' ' stns(ix).name]);

% Try to compare apples and clementines: 7dLP in situ, and 1d (spline-fit) SST
stns = verify_variable_multi(stns,'sea_t_7_day_lowpass');
stns = verify_variable_multi(stns,'pfv5c_cortad_sst_1_day_spline');
stns = verify_variable_multi(stns,'pfv5c_pentad_sst_1_day_spline');

if ( ~isempty(stns(ix).sea_t_7_day_lowpass.data) )
  scatter_fit_ts(stns(ix).sea_t_7_day_lowpass,stns(ix).bleaching,[],[],'In situ 7dLP','Bleaching%'),
  titlename(stns(ix).name); axis([25,32,0,100]);
  print('-dtiff',fullfile(figspath,[stns(ix).station_name '_sea_t_7dlp.tiff']));
end;

scatter_fit_ts(stns(ix).oisst2_sst,stns(ix).bleaching,[],[],'OISST v2 1d','Bleaching%'),
titlename(stns(ix).name); axis([25,32,0,100]);
print('-dtiff',fullfile(figspath,[stns(ix).station_name '_oisst2_1d.tiff']));

scatter_fit_ts(stns(ix).misst_sst,stns(ix).bleaching,[],[],'MISST v3 1d','Bleaching%'),
titlename(stns(ix).name); axis([25,32,0,100]);
print('-dtiff',fullfile(figspath,[stns(ix).station_name '_misst_1d.tiff']));

scatter_fit_ts(stns(ix).pfv5c_cortad_sst_1_day_spline,stns(ix).bleaching,[],[],'CoRTAD spline_1_d','Bleaching%'),
titlename(stns(ix).name); axis([25,32,0,100]);
print('-dtiff',fullfile(figspath,[stns(ix).station_name '_cortad_1d.tiff']));

scatter_fit_ts(stns(ix).pfv5c_pentad_sst_1_day_spline,stns(ix).bleaching,[],[],'Pentad spline_1_d','Bleaching%'),
titlename(stns(ix).name); axis([25,32,0,100]);
print('-dtiff',fullfile(figspath,[stns(ix).station_name '_pentad_1d.tiff']));


% Try to compare apples and apples: 7d max in situ and SST
stns = verify_variable_multi(stns,'sea_t_7_day_maximum');
stns = verify_variable_multi(stns,'oisst2_sst_7_day_maximum');
stns = verify_variable_multi(stns,'misst_sst_7_day_maximum');
stns = verify_variable_multi(stns,'pfv5c_cortad_sst_7_day_maximum');
stns = verify_variable_multi(stns,'pfv5c_pentad_sst_7_day_maximum');

if ( ~isempty(stns(ix).sea_t_7_day_maximum.data) )
  scatter_fit_ts(stns(ix).sea_t_7_day_maximum,stns(ix).bleaching,[],[],'In situ max_7_d','Bleaching%'),
  titlename(stns(ix).name); axis([25,32,0,100]);
  print('-dtiff',fullfile(figspath,[stns(ix).station_name '_sea_t_7d.tiff']));
end;

scatter_fit_ts(stns(ix).oisst2_sst_7_day_maximum,stns(ix).bleaching,[],[],'OISST v2 max_7_d','Bleaching%'),
titlename(stns(ix).name); axis([25,32,0,100]);
print('-dtiff',fullfile(figspath,[stns(ix).station_name '_oisst2_7d.tiff']));

scatter_fit_ts(stns(ix).misst_sst_7_day_maximum,stns(ix).bleaching,[],[],'MISST v3 max_7_d','Bleaching%'),
titlename(stns(ix).name); axis([25,32,0,100]);
print('-dtiff',fullfile(figspath,[stns(ix).station_name '_misst_7d.tiff']));

scatter_fit_ts(stns(ix).pfv5c_cortad_sst_7_day_maximum,stns(ix).bleaching,[],[],'CoRTAD max_7_d','Bleaching%'),
titlename(stns(ix).name); axis([25,32,0,100]);
print('-dtiff',fullfile(figspath,[stns(ix).station_name '_cortad_7d.tiff']));

scatter_fit_ts(stns(ix).pfv5c_pentad_sst_7_day_maximum,stns(ix).bleaching,[],[],'Pentad max_7_d','Bleaching%'),
titlename(stns(ix).name); axis([25,32,0,100]);
print('-dtiff',fullfile(figspath,[stns(ix).station_name '_pentad_7d.tiff']));

end;

set_more;
