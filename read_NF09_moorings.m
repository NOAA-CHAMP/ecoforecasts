function stns = read_NF09_moorings

%1;
%% SCRIPT read_NF09_moorings.m:
%
% Call load_NF09_mooring_hobos (v.) to load (or reload MAT file from) raw CSV
% data for mooring deployments associated with Eddy/Internal Wave experiment
% on the NOAA R/V Nancy Foster cruise NF0914 in 2009: then fill in STNS
% struct with data for each site/location (instrument).
%
% This script is called, e.g., from within plot_NF09_moorings_upwelling.m (v.)
%
% Last Saved Time-stamp: <Thu 2018-03-01 20:15:29 Eastern Standard Time gramer>

%load('d:/coral/FACE/Eddy09/moorings/hobos.mat');
%csvpath = 'd:/coral/FACE/Eddy09/moorings/Data/Temp & Pressure';

% facepath = get_coral_path('FACE');
% moorpath = fullfile(facepath,'Eddy09','moorings');
% moordatapath = fullfile(moorpath,'Data','Temp & Pressure');
% 
% wasdir = cd;
% cd(moorpath)
% anhobos;
% cd(wasdir);

load_NF09_mooring_hobos;

if 0;
  fmg; plot_ts(locs.seatemp);
  xlim([datenum(2009,10,14,1,0,0),datenum(2009,10,17,23,30,0)]); datetick3;
  ylim([17.03,31.23]);
  jumpax;
end;

if 0;
  bath=[]; clear bath
  bath = read_hires_bathymetry(brwd40([2,1]),[5e3,12e3],[],true);
  plot_hires_bathymetry(bath,-[0:2:50,100]);
  plot_station_marker(brwd20([2,1])); plot_station_marker(brwd40([2,1])); plot_station_marker(brwd100([2,1]));
  plot_station_marker(boca20([2,1])); plot_station_marker(boca40([2,1]));
  
  %bath=[]; clear bath
  %bath = read_hires_bathymetry(brwd40([2,1]),[5e3,5e3],[],true);
  %plot_hires_bathymetry(bath,-[0:2:50,100]);
  %plot_station_marker(brwd20([2,1])); plot_station_marker(brwd40([2,1])); plot_station_marker(brwd100([2,1]));
  
  %bath = read_hires_bathymetry(brwd40([2,1]),[6e3,12e3],[],false);
  %plot_hires_bathymetry(bath,-[0:5:80]);
  %plot_station_marker(brwd20([2,1])); plot_station_marker(brwd40([2,1])); plot_station_marker(brwd100([2,1]));
  %plot_station_marker(boca20([2,1])); plot_station_marker(boca40([2,1]));
  
  crds = {...
      [brwd20(2),brwdbenth(2),brwdbenth(4),brwdbenth(6),brwdbenth(8),brwdbenth(10),brwd40(2),brwd100(2),...
       boca20(2),bocabenth(2),bocabenth(4),bocabenth(6),bocabenth(8),bocabenth(10),boca40(2),],...
      [brwd20(1),brwdbenth(1),brwdbenth(3),brwdbenth(5),brwdbenth(7),brwdbenth(9), brwd40(1), brwd100(1),...
       boca20(1),bocabenth(1),bocabenth(3),bocabenth(5),bocabenth(7),bocabenth(9), boca40(1),],...
      {'brwd20',  'brwdbenth1', 'brwdbenth2','brwdbenth3','brwdbenth4','brwdbenth5','brwd40',  'brwd100',...
       'boca20',  'bocabenth1', 'bocabenth2','bocabenth3','bocabenth4','bocabenth5','boca40', },...
         };
  
  %% RESULTS OF THIS STORED IN temp.m:
  %[sites,bath] = find_ngdc_slope_sites(crds,bath.ngdc_hires_bathy,7);
  [sites,bath] = find_ngdc_slope_sites(crds,bath,9,{@nanmedian,9,9,10});
  for ix=1:numel(sites.stnms); disp({sites.stnms{ix},sites.depths(ix),sites.betas(ix)}); end;
  % 'brwd20'    [-23.1300]    [0.0950]
  % 'brwdbenth1'    [-24.5550]    [0.0978]
  % 'brwdbenth2'    [-27.3143]    [0.0694]
  % 'brwdbenth3'    [-30.1652]    [0.0264]
  % 'brwdbenth4'    [-31.8453]    [0.0249]
  % 'brwdbenth5'    [-32.6659]    [0.0311]
  % 'brwd40'    [-38.0934]    [0.0521]
  % 'brwd100'    [-94.6936]    [0.0791]
  % 'boca20'    [-28.6443]    [0.0245]
  % 'bocabenth1'    [-39.6522]    [0.0524]
  % 'bocabenth2'    [-35.3209]    [0.0514]
  % 'bocabenth3'    [-33.9038]    [0.0409]
  % 'bocabenth4'    [-31.1383]    [0.0221]
  % 'bocabenth5'    [-29.4974]    [0.0182]
  % 'boca40'    [-41.4811]    [0.0464]
  mean(sites.depths([2:6])), mean(sites.betas([2:6])), mean(sites.depths([10:14])), mean(sites.betas([10:14]))
  % -29.3091
  %   0.0499
  % -33.9025
  %   0.0370
end;


% Range of "good" dates for ALL moorings
dtlms = [datenum(2009,10,15),datenum(2009,12,09)];

% DEPTHS/SLOPES ESTIMATED FROM 10 m "*_mhw" BATHYMETRY (code above) & DEPLOY/
% RETRIEVAL NOTES, e.g., Thermistor-Diagrams.pdf, Thermistor-placements.pdf,
% Eddy09-mooring-plan-deploy-retrieve-LGramer-logs.pdf
deps.brwd20  = 23.13;	slps.brwd20  = 0.095;
deps.brwdbe  = 29.31;	slps.brwdbe  = 0.050;
deps.brwd40  = 38.09;	slps.brwd40  = 0.052;
deps.brwd100 = 94.69;	slps.brwd100 = 0.079;
deps.boca20  = 28.64;	slps.boca20  = 0.025;
deps.bocabe  = 33.90;	slps.bocabe  = 0.037;
deps.boca40  = 41.48;	slps.boca40  = 0.046;

brwd20_top = deps.brwd20 - 15;
brwd40_top = deps.brwd40 - 30;
brwd100_top = deps.brwd100 - 85;
boca20_top = deps.boca20 - 15;
boca40_top = deps.boca40 - 30;


% Successful?,TNO (T_# thermistor serial #),SiteID,ThermID,Coords,Depth,Slope,Description ;
NF09_lines = ...
    {
% Broward_line
        1,31,'brwd20', 'brwd20_lgt', brwd20([2,1]),    brwd20_top+1,0.157,'Broward 20 m mooring near-sfc. light' ; ...
        1,01,'brwd20', 'brwd20_top', brwd20([2,1]),    brwd20_top+2,0.157,'Broward 20 m mooring near surface' ; ...
        1,02,'brwd20', 'brwd20_btm', brwd20([2,1]),    brwd20_top+15-3,0.157,'Broward 20 m mooring near bottom' ; ...
        0,03,'brwdbe', 'brwdbenth1', brwdbenth([2,1]), 24.56,0.098,'Broward benthic line 1 nearest 20 m' ; ...
        0,04,'brwdbe', 'brwdbenth2', brwdbenth([4,3]), 27.31,0.069,'Broward benthic line 2 near 20 m' ; ...
        0,05,'brwdbe', 'brwdbenth3', brwdbenth([6,5]), 30.17,0.026,'Broward benthic line 3 midline' ; ...
        1,06,'brwdbe', 'brwdbenth4', brwdbenth([8,7]), 31.85,0.025,'Broward benthic line 4 nearer 40 m (logged z=30.2)' ; ...
        1,07,'brwdbe', 'brwdbenth5', brwdbenth([10,9]),32.67,0.031,'Broward benthic line 5 nearest 40 m (logged z=30.8)' ; ...
        1,41,'brwd40', 'brwd40_prs', brwd40([2,1]),    brwd40_top+1,0.094,'Broward 40 m mooring near-sfc. pressure' ; ...
        1,08,'brwd40', 'brwd40_top', brwd40([2,1]),    brwd40_top+2,0.094,'Broward 40 m mooring near surface' ; ...
        1,09,'brwd40', 'brwd40_mid', brwd40([2,1]),    brwd40_top+15,0.094,'Broward 40 m mooring mid-line' ; ...
        1,10,'brwd40', 'brwd40_btm', brwd40([2,1]),    brwd40_top+30-3,0.094,'Broward 40 m mooring near bottom' ; ...
        1,42,'brwd100','brwd100_prs',brwd100([2,1]),   brwd100_top+1,0.079,'Broward 100 m mooring near-sfc. pressure' ; ...
        1,11,'brwd100','brwd100_sfc',brwd100([2,1]),   brwd100_top+2,0.079,'Broward 100 m mooring near surface' ; ...
        1,12,'brwd100','brwd100_upp',brwd100([2,1]),   brwd100_top+28,0.079,'Broward 100 m mooring upper' ; ...
        1,13,'brwd100','brwd100_low',brwd100([2,1]),   brwd100_top+56,0.079,'Broward 100 m mooring lower' ; ...
        1,14,'brwd100','brwd100_btm',brwd100([2,1]),   brwd100_top+85-3,0.079,'Broward 100 m mooring near bottom' ; ...
% Boca line
        1,15,'boca20','boca20_top',  boca20([2,1]),    boca20_top+2,0.017,'Boca 20 m mooring near surface' ; ...
        1,16,'boca20','boca20_btm',  boca20([2,1]),    boca20_top+15-3,0.017,'Boca 20 m mooring near bottom' ; ...
        0,17,'bocabe','bocabenth1',  bocabenth([2,1]), 39.65,0.052,'Boca benthic line 1 nearest 20 m' ; ...
        0,18,'bocabe','bocabenth2',  bocabenth([4,3]), 35.32,0.051,'Boca benthic line 2 near 20 m' ; ...
        0,19,'bocabe','bocabenth3',  bocabenth([6,5]), 33.90,0.041,'Boca benthic line 3 midline (logged z=32.5)' ; ...
        1,20,'bocabe','bocabenth4',  bocabenth([8,7]), 31.14,0.022,'Boca benthic line 4 nearer 40 m' ; ...
        0,21,'bocabe','bocabenth5',  bocabenth([10,9]),29.50,0.018,'Boca benthic line 5 nearest 40 m' ; ...
        1,43,'boca40','boca40_prs',  boca40([2,1]),    boca40_top+1,0.038,'Boca 40 m mooring near-sfc. pressure' ; ...
        1,22,'boca40','boca40_top',  boca40([2,1]),    boca40_top+2,0.038,'Boca 40 m mooring near surface' ; ...
        1,23,'boca40','boca40_mid',  boca40([2,1]),    boca40_top+15,0.038,'Boca 40 m mooring mid-line' ; ...
        1,24,'boca40','boca40_btm',  boca40([2,1]),    boca40_top+30-3,0.038,'Boca 40 m mooring near bottom' ; ...
    };


stns=[];
for ix = 1:size(NF09_lines,1)
  if ( NF09_lines{ix,1} )
    tno = NF09_lines{ix,2};
    stnm = NF09_lines{ix,3};
    lvnm = NF09_lines{ix,4};
    locix = find(tno==[locs.tno]);
    if ( ~isfield(stns,stnm) )
      stns.(stnm).station_name   = upper(stnm);
      stns.(stnm).lon            = NF09_lines{ix,5}(1);
      stns.(stnm).lat            = NF09_lines{ix,5}(2);
      stns.(stnm).depth          = deps.(stnm);
      stns.(stnm).slope          = slps.(stnm);
      stns.(stnm).idepths        = [];
      stns.(stnm).seatemps       = repmat(struct('date',[],'data',[]),[0,0]);
    end;
    stns.(stnm).(lvnm).lon       = NF09_lines{ix,5}(1);
    stns.(stnm).(lvnm).lat       = NF09_lines{ix,5}(2);
    stns.(stnm).(lvnm).depth     = NF09_lines{ix,6};
    stns.(stnm).(lvnm).slope     = NF09_lines{ix,7};

    stns.(stnm).(lvnm).tno       = locs(locix).tno;
    stns.(stnm).(lvnm).seatemp   = locs(locix).seatemp;
    stns.(stnm).(lvnm).seapres   = locs(locix).seapres;
    stns.(stnm).(lvnm).par       = locs(locix).par;
    stns.(stnm).(lvnm).light     = locs(locix).light;

    stns.(stnm).idepths(end+1)   = stns.(stnm).(lvnm).depth;
    stns.(stnm).seatemps(end+1)  = date_range_ts(stns.(stnm).(lvnm).seatemp,dtlms);
  end;
end;

stnms = fieldnames(stns);
for stix=1:numel(stnms)
  stnm = stnms{stix};
  cts = intersect_tses(stns.(stnm).seatemps);
  ts = [cts{:}]; cts=[]; clear cts
  stns.(stnm).seatemp.date = ts(1).date;
  stns.(stnm).seatemp.data = ts(1).data;
  stns.(stnm).seatemp.prof = [ts.data];
  stns.(stnm).seatemp.depth = stns.(stnm).idepths;
  ts=[]; clear ts
  stns.(stnm).seatemp.data = nanmean(stns.(stnm).seatemp.prof,2);
end;


if 0;
  fmg; contourf(stns.brwdbe.seatemp.date,-stns.brwdbe.seatemp.depth,stns.brwdbe.seatemp.prof'); colorbar; datetick3;
  fmg; contourf(stns.brwd40.seatemp.date,-stns.brwd40.seatemp.depth,stns.brwd40.seatemp.prof'); colorbar; datetick3;
  fmg; contourf(stns.brwd100.seatemp.date,-stns.brwd100.seatemp.depth,stns.brwd100.seatemp.prof'); colorbar; datetick3;

  fmg; contourf(stns.boca40.seatemp.date,-stns.boca40.seatemp.depth,stns.boca40.seatemp.prof'); colorbar; datetick3;
  fmg; surf(stns.boca40.seatemp.date,-stns.boca40.seatemp.depth,stns.boca40.seatemp.prof'); shading interp; colorbar; datetick3;

  %fmg; plot_ts(stns(ismember([stns.tno],[11:12,41:43])).seatemp); legend(num2str([stns(ismember([stns.tno],[11:12,41:43])).tno]'));
end;

if 0;
  ds = dir(fullfile(facepath,'*.csv'));
  for dix=1:numel(ds)
    fname = fullfile(facepath,ds(dix).name);
    stnm = regexprep(ds(dix).name,'_.emp_.ogger.*[.]csv$','');
    stnm = strrep(stnm,'-','_');
    disp([stnm,': ',fname]);
    tic,
      [n,t,r] = xlsread(fname);
    toc,
    hdrs = t(2,:);
    timeix = find(~cellfun(@isempty,strfind(hdrs,'Time,')));
    tempix = find(~cellfun(@isempty,strfind(hdrs,'Temp,')));
    tic,
      dts = datenum(r(3:end,timeix));
    toc,
    stns.(stnm).seatemp.date = dts;
    stns.(stnm).seatemp.data = n(:,tempix);
    keyboard;
  end;
  
  % goodix = find(cellfun(@length,r(:,2))==21);
  % stnm = num2str([r{goodix,2}]');
  % lats = [r(goodix,3)]';
  % lons = [r(goodix,4)]';
  % tdts = datenum(r(goodix,5),'yyyy-mm-ddTHH:MM:SS');
  % deps = [r(goodix,6)]';
  % dirs = [r{goodix,7}]';
  % spds = [r{goodix,8}]';
  % vspd = [r{goodix,9}]';
end;


return;
