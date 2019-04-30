1;
%% Script LOAD_ALL_SOUTH_FLORIDA_TEMPERATURE_SITES.m
%
% Create a "database" of in situ sea temperature measurements in south
% Florida waters, drawn from many different historical and ongoing projects
%
% Last Saved Time-stamp: <Thu 2019-03-21 17:25:54 EDT lew.gramer>

all_sites_matfname = fullfile(get_ecoforecasts_path('data'),'south_florida_temperature_sites.mat');

if ( exist(all_sites_matfname) )
  disp(['Load ',all_sites_matfname]);
  load(all_sites_matfname,'field','sites');

else
  more off;

  disp('****************************************');
  disp('Extracting sea temperatures and metadata');
  disp('****************************************');

  disp('****************************************');
  disp('Gathering FKNMS and SEAKEYS sites...');
  doPrint = false;  doSave = false;  doFigs = false;  doDiags = false;
  doCSV = false;  doGrid = false;  doNuts = false;  doRi = false;
  subrgn='FRT';
  use_habitat_map=false;
  allow_hard_bottom=true;
  calc_spatial_dt_hc;
  scaling='SS'; 
  comparison_plot_spatial_dt_hc;

  clear f*
  field.LAT=LAT; 
  field.LON=LON; 
  field.h=h; 
  field.bet=bet; 
  field.ang=ang; 
  field.hch=hch; 
  field.hcrng=hcrng;

  clear -regexp ^[a-eg-su-zA-Z]*

  sites=[]; clear sites;
  sites = ts;
  for tsix=1:numel(ts)
    sites(tsix).seatemp.date = sites(tsix).date;
    sites(tsix).seatemp.data = sites(tsix).data;
  end;
  sites = rmfield(sites,'date');
  sites = rmfield(sites,'data');

  ts=[];
  clear t*

  disp('****************************************');
  disp('Gathering FACE, SEFCRI, and SFOMC sites...');
  doPrint = false;  doSave = false;  doFigs = false;  doDiags = false;
  doCSV = false;  doGrid = false;  doNuts = false;  doRi = false;
  calc_upwelling_ekman_flux;
  sites = copy_T_sites(face,sites,'FACE_HWD_');
  sites = copy_T_sites(sefcri,sites,'SEFCRI_');
  sites = copy_T_sites(sfomc,sites,'SFOMC_');
  face=[]; sefcri=[]; sfomc=[]; clear face sefcri sfomc

  sites = copy_T_sites({canf1,cnnf1,fdepk,ftpf1,fwyf1,lkwf1,pvgf1},sites,'CMAN_');

  disp('****************************************');
  disp('Gathering NCORE sites...');
  doPrint = false;  doSave = false;  doFigs = false;  doDiags = false;
  doCSV = false;  doGrid = false;  doNuts = false;  doRi = false;
  ncore = get_ncore_2000;
  [ncore.klgf1,ncore.mrtf1] = an_ncore;
  sites = copy_T_sites(ncore,sites,'NCORE_');
  ncore=[]; clear ncore

  disp('****************************************');
  disp('Gathering FACE NF09 sites...');
  doPrint = false;  doSave = false;  doFigs = false;  doDiags = false;
  doCSV = false;  doGrid = false;  doNuts = false;  doRi = false;
  face = read_NF09_moorings;
  sites = copy_T_sites(face,sites,'FACE_EDDY_');

  % disp('****************************************');
  % disp('Gathering meteorology');
  % disp('... ??? ...');

  disp('****************************************');
  %%%% REDEFINE HERE??? just in case some script somewhere wiped this out...
  all_sites_matfname = fullfile(get_ecoforecasts_path('data'),'south_florida_temperature_sites.mat');
  %disp(['DO NOT Save ',all_sites_matfname]);
  disp(['Save ',all_sites_matfname]);
  save(all_sites_matfname,'field','sites');

  more on;
end;
