
  if ~exist('fld','var'); 
    fld.date=[]; 
    fld.nwps_sigwavehgt.field = [];
    fld.nwps_peakwavedir.field = [];
    fld.nwps_peakwaveper.field = [];
  end; 
  

  fld.nwps_sigwavehgt.field = [];
  fld.nwps_peakwavedir.field = [];
  fld.nwps_peakwaveper.field = [];
  

    fld.nwps_sigwavehgt.field(dtixen,latixen,lonixen) = x.nwps_sigwavehgt.field;
    fld.nwps_peakwavedir.field(dtixen,latixen,lonixen) = x.nwps_peakwavedir.field;
    fld.nwps_peakwaveper.field(dtixen,latixen,lonixen) = x.nwps_peakwaveper.field;






yrmos = 
yr=2017;
for mo=1:12; 
  disp(mo);
  x=load(['d:/data/NWPS/NWPS_mfl_nwps_CG1_',sprintf('%04d_%02d',yr,mo),'.mat']);



  %sites = copy_T_sites({canf1,cnnf1,fdepk,ftpf1,fwyf1,lkwf1,pvgf1},sites,'CMAN_');
  sites = copy_T_sites(canf1,sites,'CMAN_');
  sites = copy_T_sites(cnnf1,sites,'CMAN_');
  sites = copy_T_sites(fdepk,sites,'CMAN_');
  sites = copy_T_sites(ftpf1,sites,'CMAN_');
  sites = copy_T_sites(fwyf1,sites,'CMAN_');
  sites = copy_T_sites(lkwf1,sites,'CMAN_');
  sites = copy_T_sites(pvgf1,sites,'CMAN_');






    if ( isfield(stns.(stnm),'depth') )
      sites(stix).depth = stns.(stnm).depth;
    end;



  doPrint = false;
  doSave = false;
  doFigs = false;
  doDiags = false;




      %%%%system(['c:/cygwin/bin/tar xvfz ../../archive/NWPS/NWPS_archive_',yrmo,'.tgz --exclude=''key_*''']);
      %%%% NEED TO CHANGE: PATH=c:/cygwin/home/gramer/perl5/bin:c:/cygwin/chroot/gensym/usr/gensym/g2/g2:.:c:/cygwin/home/gramer/bin:c:/cygwin/sbin:c:/cygwin/usr/sbin:c:/cygwin/bin:c:/cygwin/usr/bin:c:/cygwin/usr/local/bin:c:/cygwin/usr/local/bin:c:/cygwin/usr/bin:c:/cygwin/bin
      %system(['c:/cygwin/bin/tcsh -f -c ''(set PATH="c:/cygwin/home/gramer/perl5/bin:c:/cygwin/chroot/gensym/usr/gensym/g2/g2:.:c:/cygwin/home/gramer/bin:c:/cygwin/sbin:c:/cygwin/usr/sbin:c:/cygwin/bin:c:/cygwin/usr/bin:c:/cygwin/usr/local/bin:c:/cygwin/usr/local/bin:c:/cygwin/usr/bin:c:/cygwin/bin:$PATH" ; c:/cygwin/bin/tar --help)''']);
      system(['c:/cygwin/bin/tcsh -f -c ''(set PATH="c:/cygwin/home/gramer/perl5/bin:c:/cygwin/chroot/gensym/usr/gensym/g2/g2:.:c:/cygwin/home/gramer/bin:c:/cygwin/sbin:c:/cygwin/usr/sbin:c:/cygwin/bin:c:/cygwin/usr/bin:c:/cygwin/usr/local/bin:c:/cygwin/usr/local/bin:c:/cygwin/usr/bin:c:/cygwin/bin:$PATH" ; c:/cygwin/bin/tar xvfz ../../archive/NWPS/NWPS_archive_',yrmo,'.tgz --exclude="key_*")']);



    %if ( `ls mfl_nwps_CG1_${yr}${mo}*.grib2 >& /dev/null` ) then
    ls mfl_nwps_CG1_${yr}${mo}*.grib2 >& /dev/null
    if ( ${status} ) then


matize_NWPS_mfl_${yr}${mo}.log





disp([numel(unique(fdtix)),numel(unique(dtix))]);
disp(datestr(field.date(fdtix([1,2,end]))));
disp(datestr(dts(dtix([1,2,end]))));
disp(median(diff(field.date(fdtix([1,2])))));
disp(median(diff(dts(dtix([1,2])))));





From SCATTER_FIT.m:
    legs1 = [xname ' vs. ' yname];
    %                   B(1), B(2), (Stats.coeffcorr(2,1)^2), ...
    legs2 = sprintf('%.5g + %.5g*X, R^2~%0.2g (%0.2g;%0.2g;%0.2g), \n RMSE=%g, p=%0.5g, N=%d', ...
                    B(1), B(2), (Stats.regress_stats(1)), Stats.R2_1, Stats.R2_2, Stats.R2_2, ...
                    Stats.s, roundn(Stats.p(2),-5), Stats.N);


      %legend( legs1, legs2, 'Location','NorthWest');
      %legend( legs1, legs2, 'Location','Best');
      %legend( legs1, legs2, 'Location','SouthEast');

      %legend( legs1, legs2, 'Perfect fit', 'Location','NorthWest');
      %legend( legs1, legs2, 'Perfect fit', 'Location','Best');
      %legend( legs1, legs2, 'Perfect fit', 'Location','SouthEast');



legend(sh,'Location','Best');


  set(lh,'LineWidth',2); set(th,'FontSize',16);



stns = get_nwps_mat_stations({ncrmp.(dataset).ulon,ncrmp.(dataset).ulat,ncrmp.(dataset).ustn});



    
    for ix=1:numel(args{1})
    if ( isfield(str,'lon') )
      LONS = str.lon;
      LATS = str.lat;
    elseif ( isfield(str,'lons') )
      LONS = str.lons;
      LATS = str.lats;
    else
    end;
    
    if ( numel(args) > 1 )
      if ( isfield(str,args{2}) )
        varargout{1} = str.(args{2});
      end;
    end;
    






      tses(ix).data(dtix) = interp3(fld.lon,fld.lat,fld.field,tses(ix).lon,tses(ix).lat,'linear');



function tses = get_nwps_mat_tses(stns_or_coords,dts)
%function tses = get_nwps_mat_tses(stns_or_coords,dts)

  datapath = get_ecoforecasts_path('data');

  if ( ~exist('mdl','var') || isempty(mdl) )
    mdl = 'key_nwps_CG2';
  end;

  [lons,lats] = get_coords_from_args(stns_or_coords);

  % Convert DTS to a 3-hourly DATENUM vector with whole-day boundaries
  begix = find(get_hour(dts)==0,1);
  endix = find(get_hour(dts)>=21,1,'last');
  dts = dts(begix:endix);
  dts = dts(ismember(get_hour(dts),[0:3:24]));
  
  for ix=1:numel(lons)
    tses(ix).date = dts(:);
    tses(ix).data = repmat(nan,size(tses(ix).date));
    tses(ix).lon = lons(ix);
    tses(ix).lat = lats(ix);
  end;
  
  yrmos = unique(get_yearmonth(dts));

  for yrmo = yrmos(:)';
    yr = get_year(yrmo);
    mo = get_month(yrmo);

    matfname = fullfile(datapath,sprintf('NWPS_%s_%04d_%02d.mat',mdl,yr,mo));
    if ( ~exist(matfname,'file') )
      disp(['MISSING ',matfname]);
      %%%% EARLY CONTINUE
      continue;
    end;
    %disp(sprintf('Extracting NWPS for %s_%04d_%02d',mdl,yr,mo));
    fld = load(matfname);
    for ix=1:numel(tses)
      [lonerr,lonix] = min(abs(tses(ix).lon-fld.lon));
      [laterr,latix] = min(abs(tses(ix).lat-fld.lat));

      dtix = interp1(dts,1:numel(dts),fld.date,'nearest');
      tses(ix).data(dtix) = interp2(fld.lon,fld.lat,fld.field,tses(ix).lon,tses(ix).lat,'linear');
    end;
  end;

return;




fmg; contourf(key201610.lon,key201610.lat,squeeze(nanmean(key201610.nwps_sigwavehgt.field)),[0:0.1:1.8]); colorbar; set(gca,'CLim',[0,1.8]); titlename('key201610');
fmg; contourf(key201810.lon,key201810.lat,squeeze(nanmean(key201810.nwps_sigwavehgt.field)),[0:0.1:1.8]); colorbar; set(gca,'CLim',[0,1.8]); titlename('key201810');

% Compare at same TIME resolution
fmg; contourf(key201610.lon,key201610.lat,squeeze(nanmean(key201610.nwps_sigwavehgt.field)),[0:0.1:1.8]); colorbar; set(gca,'CLim',[0,1.8]); titlename('key201610');
dtix=find(ismember(get_yearday(key201810.date),get_yearday(key201610.date)));
fmg; contourf(key201810.lon,key201810.lat,squeeze(nanmean(key201810.nwps_sigwavehgt.field(dtix,:,:))),[0:0.1:1.8]); colorbar; set(gca,'CLim',[0,1.8]); titlename('key201810 3-hr');


  %% Scleractinian COVER_CAT_CD list
  % Good codes, not sacred cod
  goodcods = string({'ACR CERV','ACR PALM','AGA AGAR','COL NATA','DIC SPE.','DIP LABY','EUS FAST','HEL CUCU','MAD DECA','MAD MIRA','MEA MEAN','MON CAVE','MUS ANGU','MYC DANA','MYC SPE.','ORB ANNU','ORB FAVE','ORB FRAN','POR ASTR','POR DIVA','POR FURC','POR PORI','PSE STRI','SID RADI','SID SIDE','SOL BOUR','SOL HYAD','STE INTE',});
  goodix = find(ismember(ncrmp.(rgnyr).cod,goodcods));
  
  % Massive coral spp.
  masscods = string({'COL NATA','DIC SPE.','DIP LABY','MON CAVE','ORB ANNU','ORB FAVE','ORB FRAN','POR ASTR','SID RADI','SID SIDE','SOL BOUR','STE INTE',});
  massix = find(ismember(ncrmp.(rgnyr).cod,masscods));
  
  for ix=1:numel(ncrmp.(rgnyr).ulon)
    ncrmp.(rgnyr).uix{ix} = find(ncrmp.(rgnyr).lon==ncrmp.(rgnyr).ulon(ix) ...
                                 & ncrmp.(rgnyr).lat==ncrmp.(rgnyr).ulat(ix) ...
                                 & ismember(ncrmp.(rgnyr).cod,goodcods));
    
    ncrmp.(rgnyr).uix{ix} = intersect(ncrmp.(rgnyr).uix{ix},goodix);
    ncrmp.(rgnyr).mix{ix} = intersect(ncrmp.(rgnyr).uix{ix},massix);






function [dx,az,nearestCoords] = get_contour_distance(contourCoords_or_bath,fromCoords)
%function [dx,az,nearestCoords] = get_contour_distance(contourCoords_or_bath,fromCoords)
%
% Get the distance in KM between each point in FROMCOORDS (a 2xN matrix
% LON;LAT, OR a matrix of STRUCT with .lon and .lat files) and its nearest
% point in contour or locus of points CONTOURCOORDS; first arg may also be a
% bathymetry STRUCT with fields .lon,.lat,.field, or a CELL array containing
% {BATHY_STRUCT,NUMERIC_ISOBATH}. In any case, the first arg is passed
% through untouched to GET_NEAREST_POINTS (v.)
%
% Optionally also returns AZ azimuths [deg True] to shore, and NEARESTCOORDS
% the nearest points on the contour from each point. DX and AZ should be of
% size 1xN, FROMCOORDS of size 2xN.
%
% Last Saved Time-stamp: <Fri 2019-02-15 14:14:35 Eastern Standard Time gramer>

  if ( isstruct(fromCoords) && isfield(fromCoords,'lon') && isfield(fromCoords,'lat') )
    lon = [fromCoords(:).lon];
    lat = [fromCoords(:).lat];
    clear fromCoords
    fromCoords(1,1:numel(lon)) = lon;
    fromCoords(2,1:numel(lat)) = lat;
  else
    lon = fromCoords(1,:);
    lat = fromCoords(2,:);
  end;

  nearestCoords = get_nearest_points(contourCoords_or_bath,fromCoords);
  [dx,az] = distance_wgs84(lat,lon,nearestCoords(2,:),nearestCoords(1,:));

return;









sites=[]; clear sites
%[sites,bath] = find_ngdc_slope_sites({ncrmp.(dataset).ulon,ncrmp.(dataset).ulat},bath,9);





% Cut off at slopes and DEPTHS where horizontal convection and wave breaking are unlikely
goodix = find(~isnan(sites.betas) & sites.betas>=0.0075 & 1<=ncrmp.(dataset).unh & ncrmp.(dataset).unh<=20.0);
disp(['Good=',num2str(100*numel(goodix)/numel(sites.betas)),'%']);




% Cut off at slopes where horizontal convection and wave breaking are unlikely
goodix = find(~isnan(sites.betas) & sites.betas>=0.0075);
disp(['Good=',num2str(100*numel(goodix)/numel(sites.betas)),'%']);




  [bath,rad] = read_hires_bathymetry_for_field({A.ulon,A.ulat},false,[10e3,20e3]);

%goodix = find(~isnan(sites.betas) & sites.betas>=0.0025);



% xts = double(string(get(gca,'XTickLabel'))).^(1/expon);
% set(gca,'XTickLabel',num2str(xts,'%0.3f'));





expon = 2/3;
scatter_fit(sites.betas(goodix)'.^expon,A.utc(goodix)',['Seafloor Slope (\beta)^{',num2str(expon),'}'],[textize(A.basename),' Cvr%']); set(gca,'XScale','log'); ylim([0,100]); legend('Location','Best');
for bet=bets(:)'; annotline(bet.^expon,[],num2str(bet)); end; clear bet
% xts = double(string(get(gca,'XTickLabel'))).^(1/expon);
% set(gca,'XTickLabel',num2str(xts,'%0.3f'));




%goodix = find(~isnan(sites.betas));

% Cut off at slopes where horizontal convection and wave breaking are unlikely
goodix = find(~isnan(sites.betas) & sites.betas>=0.0100); %0.0075); %0.0050);

disp(['Good=',num2str(100*numel(goodix)/numel(sites.betas)),'%']);



% scatter_fit(sites.betas(goodix)',A.utc(goodix)','Seafloor Slope (\beta)',[textize(A.basename),' Cvr%']); set(gca,'XScale','log'); ylim([0,100]); legend('Location','Best');

% scatter_curve_fit(sites.betas(goodix)',A.utc(goodix)','a*(x-b)^n','Seafloor Slope (\beta)',[textize(A.basename),' Cvr%'],[],'ConfInt',0.67,'StartPoint',[1,0,+2/3],'Lower',[1/3,0,2/3],'Robust','LAR','MaxFunEvals',1e6); set(gca,'XScale','log'); ylim([0,100]); legend('Location','Best');

% scatter_curve_fit(sites.betas(goodix)',A.utc(goodix)','a*(x-b)^n','Seafloor Slope (\beta)',[textize(A.basename),' Cvr%'],[],'ConfInt',0.67,'StartPoint',[1,0,+2/3],'Lower',[1/3,0,1/3],'Robust','LAR','MaxFunEvals',1e6); set(gca,'XScale','log'); ylim([0,100]); legend('Location','Best');




scatter_fit(sites.betas(goodix)'.^(1/3),A.utc(goodix)','Seafloor Slope (\beta)^2^/^3',[textize(A.basename),' Cvr%']); set(gca,'XScale','log'); ylim([0,100]); legend('Location','Best');

xts = double(string(get(gca,'XTickLabel'))).^(3/2);
set(gca,'XTickLabel',num2str(xts,'%0.3f'));





%region = 'FGBNMS'; year = 2013;
%region = 'FGBNMS'; year = 2015;
%region = 'FLK'; year = 2014;
%region = 'FLK'; year = 2016;
%region = 'PRICO'; year = 2014;
%region = 'PRICO'; year = 2016;
%region = 'SEFCRI'; year = 2014;
region = 'SEFCRI'; year = 2016;
%region = 'TortugasMarq'; year = 2014;
%region = 'TortugasMarq'; year = 2016;
%region = 'USVI; year = 2013;
%region = 'USVI; year = 2015;
%region = 'USVI; year = 2017;





lonA = A.ulon;
latA = A.ulat;
lonB = B.ulon;
latB = B.ulat;




%{
for ix=1:numel(ulon)
  disp({ustn(ix),numel(uix{ix})});
  
  % if(year == 2014 && region == "SEFCRI" ||
  %    year == 2014 && region == "FLK") {
  %
  %   dat2 <- dat2 %>%
  %     dplyr::group_by(YEAR, REGION, SUB_REGION_NAME, PRIMARY_SAMPLE_UNIT, STATION_NR, LAT_DEGREES, LON_DEGREES,
  %                     MIN_DEPTH, MAX_DEPTH, ANALYSIS_STRATUM, STRAT, HABITAT_CD, PROT, cover_group) %>%
  %     dplyr::summarise(Percent_Cvr = sum(Percent_Cvr)) %>%
  %     dplyr::ungroup() %>%
  %     dplyr::group_by(YEAR, REGION, SUB_REGION_NAME, PRIMARY_SAMPLE_UNIT, LAT_DEGREES, LON_DEGREES,
  %                     MIN_DEPTH, MAX_DEPTH, ANALYSIS_STRATUM, STRAT, HABITAT_CD, PROT, cover_group) %>%
  %     dplyr::summarise(Percent_Cvr = mean(Percent_Cvr)) %>%
  %     dplyr::ungroup()


  % dplyr::group_by(YEAR, REGION, SUB_REGION_NAME, PRIMARY_SAMPLE_UNIT, STATION_NR, LAT_DEGREES, LON_DEGREES,
  %                 MIN_DEPTH, MAX_DEPTH, ANALYSIS_STRATUM, STRAT, HABITAT_CD, PROT, cover_group) %>%
  % dplyr::summarise(Percent_Cvr = sum(Percent_Cvr)) %>%
  %for cix=[1,14,2,3,23,24,8,7,21];
  for cix=[2,3];
    val = unique(dat(uix{ix},cix));
    disp(string({'    ',cix,val{:}}));
  end;
end;
%}





two_stage = false;

  two_stage = true;

if ( two_stage )
end;



[region_year,'_2stage_benthic_cover.csv'];


matfname = fullfile(ncrmppath,'FLK_2014_2stage_benthic_cover.csv');


for ix=1:numel(ulon)
  if ( numel(uix{ix}) == 1 && numel(unique(dat(uix{ix},3))) ~= 1 )
    keyboard;
  end;
  if ( numel(uix{ix}) > 1 && numel(unique(dat(uix{ix},3))) ~= 2 )
    keyboard;
  end;
end;




%{
for ix=1:numel(ulon)
  disp({ustn(ix),numel(uix{ix})});
  
  % if(year == 2014 && region == "SEFCRI" ||
  %    year == 2014 && region == "FLK") {
  %
  %   dat2 <- dat2 %>%
  %     dplyr::group_by(YEAR, REGION, SUB_REGION_NAME, PRIMARY_SAMPLE_UNIT, STATION_NR, LAT_DEGREES, LON_DEGREES,
  %                     MIN_DEPTH, MAX_DEPTH, ANALYSIS_STRATUM, STRAT, HABITAT_CD, PROT, cover_group) %>%
  %     dplyr::summarise(Percent_Cvr = sum(Percent_Cvr)) %>%
  %     dplyr::ungroup() %>%
  %     dplyr::group_by(YEAR, REGION, SUB_REGION_NAME, PRIMARY_SAMPLE_UNIT, LAT_DEGREES, LON_DEGREES,
  %                     MIN_DEPTH, MAX_DEPTH, ANALYSIS_STRATUM, STRAT, HABITAT_CD, PROT, cover_group) %>%
  %     dplyr::summarise(Percent_Cvr = mean(Percent_Cvr)) %>%
  %     dplyr::ungroup()


  % dplyr::group_by(YEAR, REGION, SUB_REGION_NAME, PRIMARY_SAMPLE_UNIT, STATION_NR, LAT_DEGREES, LON_DEGREES,
  %                 MIN_DEPTH, MAX_DEPTH, ANALYSIS_STRATUM, STRAT, HABITAT_CD, PROT, cover_group) %>%
  % dplyr::summarise(Percent_Cvr = sum(Percent_Cvr)) %>%
  %for cix=[1,14,2,3,23,24,8,7,21];
  for cix=[2,3];
    val = unique(dat(uix{ix},cix));
    disp(string({'    ',cix,val{:}}));
  end;
end;
%}








% Unique stations

[ulon,uix] = unique(lon);
ulat = lat(uix);
ustn = stn(uix);




cov = mean(x.data,2);



goodspp = string({'ACR CERV','ACR PALM','AGA AGAR','COL NATA','DIC SPE.','DIP LABY','EUS FAST','HEL CUCU','MAD DECA','MAD MIRA','MEA MEAN','MON CAVE','MUS ANGU','MYC DANA','MYC SPE.','ORB ANNU','ORB FAVE','ORB FRAN','POR ASTR','POR DIVA','POR FURC','POR PORI','PSE STRI','SID RADI','SID SIDE','SOL BOUR','SOL HYAD','STE INTE',});

for uix=1:numel(ulon)
  ix = find(lon==ulon(uix) & lat==ulat(uix));
  spix = find(ismember(cod(ix),goodspp));
  tco(uix) = sum(cov(ix(spix)));
end;





for uix=1:numel(ulon)
  ix = find(lon==ulon(uix) & lat==ulat(uix)) & ismember(cod,goodspp));
  tco(uix) = cov(ix);
end;


for uix=1:numel(ulon)
  ix = find(lon==ulon(uix) & lat==ulat(uix));
  spix = find(ismember(cod(ix),goodspp))
  tco = cov(ix(spix));
end;



if ( ~exist('nwpspath','var') )
  nwpspath = get_ecoforecasts_path('data/NWPS');
end;


        try, F=websave(fpath,url); catch ME, try, delete(fpath); catch ME, ...
               end; try, delete([fpath,'.html']); catch ME, end; F=[]; end;



for cnm={'plsf1','kywf1','sanf1','smkf1','vcaf1','lonf1','mlrf1','fwyf1','vakf1','pvgf1','lkwf1','41114','41009','41113'}; nm=cnm{:}; disp(nm); download_ndbc_files(nm,2018,true,false); end; clear cnm nm

for cnm={'plsf1','kywf1','sanf1','smkf1','vcaf1','lonf1','mlrf1','fwyf1', 'vakf1','pvgf1','lkwf1','41114','41009','41113'}; nm=cnm{:}; disp(nm); stn = get_station_from_station_name(nm); stn = load_all_ndbc_data(stn); stn=[]; clear stn; end; clear cnm nm


        try, F=websave(fpath,url); catch ME, F=[]; delete(fpath); end;

        try, F=websave(o_fpath,url); catch ME, F=[]; delete(o_fpath); end;


        try, F = websave(fpath,url); catch ME, F = []; end;





FROM compare_spatial_dt_hc_thermal_stress.m:
switch (scenario)
 case 1,	T = TC;	dts = datenum(2010,[1,2],1);	cmp=cldcmp;	scenario_desc = 'Mass Cold Snap';
 case 2,	T = TW;	dts = datenum(2009,[1,2],1);	cmp=cldcmp;	scenario_desc = 'Normal Winter';
 case 3,	T = TS;	dts = datenum(2009,[8,9],1);	cmp=wrmcmp;	scenario_desc = 'Normal Summer';
 case 4,	T = TB;	dts = datenum(1998,[8,9],1);	cmp=wrmcmp;	scenario_desc = 'Mass Bleaching';
end;



%alldts = datenum(2016,8,28,6,0,0):0.5:datenum(2016,8,28,18,0,0);





ndts = (numel(alldts)*4) + 35 - 1;



% begdt = datenum(2016,8,1);
% enddt = datenum(2016,8,31);
% hrs = [6,18];


  if ( hr ~= 6 && hr ~= 18 )
    % EARLY LOOP CONTINUATION
    continue;
  end;



1;

nwpspath = get_ecoforecasts_path('data/NWPS');
mdl = 'key_nwps_CG2';

begdt = datenum(2016,8,28);
enddt = datenum(2016,8,28);
%hrs = [6,18];
hrs = [6];

flds = {
    'nwps_wind_dir',    	'Wind_direction_from_which_blowing_surface' ; ...
    'nwps_wind_speed',  	'Wind_speed_surface' ; ...
    'nwps_sigwavehgt',     	'Significant_height_of_combined_wind_waves_and_swell_surface' ; ...
    'nwps_swellhgt',    	'Significant_height_of_swell_waves_surface' ; ...
    'nwps_primwavedir', 	'Primary_wave_direction_surface' ; ...
    'nwps_primwaveper', 	'Primary_wave_mean_period_surface' ; ...
    'nwps_VAR10_0_193', 	'VAR10-0-193_FROM_7-0-0_surface' ; ...
    'nwps_curr_dir',    	'Current_direction_surface' ; ...
    'nwps_curr_speed',  	'Current_speed_surface' ; ...
    'nwps_sealeveldev', 	'Deviation_of_sea_level_from_mean_surface' ; ...
    'nwps_water_depth', 	'Water_depth_surface' ; ...
       };

field=[]; clear field
field.date = [];

for dt=begdt:enddt
  yr = get_year(dt);
  mo = get_month(dt);
  dy = get_dom(dt);

  for hr=hrs(:)';
    %nc = mDataset('d:/ecoforecasts/data/NWPS/key_nwps_CG2_20160828_0600.grib2');
    ncfname = fullfile(nwpspath,sprintf('%s_%04d%02d%02d_%02d00.grib2',mdl,yr,mo,dy,hr));
    %DEBUG:
    disp(ncfname);
    nc = mDataset(ncfname);
    try,
      if ( ~isfield(field,'lon') )
        field.lon = cast(nc{'lon'}(:,:),'double');
        field.lat = cast(nc{'lat'}(:,:),'double');
        nlon = numel(field.lon);
        nlat = numel(field.lat);
      end;
      dts = dt + ( (hr + cast(nc{'time'}(:,:),'double'))/24 );
      ndts = numel(dts);
      field.date(end+1:end+ndts) = dts;

      for fldix=1:size(flds,1)
        fld = flds{fldix,1};
        %DEBUG:
        disp(fld);
        if ( ~isfield(field,fld) )
          field.(fld).field = repmat(nan,[0,nlat,nlon]);
        end;
        field.(fld).field(end+1:end+ndts,:,:) = cast(nc{flds{fldix,2}}(:,:,:),'double');
      end;
    catch,
      catchwarn;
    end;
    close(nc); clear nc

  end;

end;

fmg; contourf(field.lon,field.lat,squeeze(nanmean(field.nwps_wavehgt.field))); colorbar;








% Sea-surface heat flux (Q0 < 0 cooling the ocean)
% MLRF1: January minimum and mean, June mean, June maximum
q0s = [  -800,              -115,            +95,            +255];
% Lowest air temperature: rough floor for air-sea cooling
% 7th or 93rd percentile Air Temperature for the month of each event
Tas = [  7.5,               14.0,            30.0,           30.5];
q0str = {'2010 cold snap', 'Normal winter', 'Normal summer','1998 bleaching'};

nq0s = numel(q0s);






for cf = {'30m','30m_NANMEAN_9_9_10','92m','92m_NANMEAN_3_3_4',};
  load(['FRT_depth_and_beta_',cf{:},'_hc_range_depth.mat']);

  %for cs={'CCAF1','CNDF1','CGRF1','CPIF1'};
  for cs={'CGRF1'};
    stnm=cs{:};




clear DBG
% DBG.point = [1295,3586]; % BNPON
% DBG.point = [1297,3574]; % DBG.coords = [-80.1764,25.3641]; %BNPMI
if ( exist('stnm','var') )
  DBG.stnm = stnm;
else
  DBG.stnm = 'MLRF1';
  % DBG.stnm = 'CCAF1';
  % DBG.stnm = 'CNDF1';
  % DBG.stnm = 'CGRF1';
  % DBG.stnm = 'CPIF1';
end;
[DBG.lon,DBG.lat] = get_station_coords(DBG.stnm);







[DBG.dzdx,DBG.dzdy] = gradient(DBG.h,dx);
DBG.dzdx = -DBG.dzdx; DBG.dzdy = -DBG.dzdy; 




% DBG.latrng = 20;
% DBG.latixen = DBG.latix-DBG.latrng:DBG.latix+DBG.latrng;
% DBG.lonrng = 20;
% DBG.lonixen = DBG.lonix-DBG.lonrng:DBG.lonix+DBG.lonrng;





  if ( exist(matfname,'file') )
    disp(['Loading ',matfname]);
    load(matfname);

    % Somehow, the unaveraged bathymetry got processed differently before saving??
    if ( size(h,1) ~= numel(lat) )
      if ( size(h,2) ~= numel(lat) )
        error('Misshapen field! %s',matfname);
      end;
      disp('Transposing fields...');
      ang = ang';
      bet = bet';
      h = h';
      hch = hch';
      hcrng = hcrng';
    end;

    % No HC over land, and depths are positive
    if ( nanmin(h(:)) < -eps )
      disp('Inverting H sign...');
      h(h>0) = 0;
      h = abs(h);
    end;

  else





    % Somehow, the unaveraged bathymetry got processed differently before saving??
    if ( downstride == 1 )
      h = h';
      bet = bet';
      % No HC over land, and depths are positive
      h(h>0) = 0;
      h = abs(h);
    end;






  for ix=1:numel(ctix)
    %lon = nanmin(geodat(ctix(ix)).Lon)+0.05;
    %lat = mean([nanmin(geodat(ctix(ix)).Lat),nanmax(geodat(ctix(ix)).Lat)]);
    %th = text(lon,lat,upper(geodat(ctix(ix)).NAME_2),'HorizontalAlignment','center','FontSize',fontsz);
    lon = mean([nanmin(geodat(ctix(ix)).Lon),nanmax(geodat(ctix(ix)).Lon)]);
    lat = mean([nanmin(geodat(ctix(ix)).Lat),nanmax(geodat(ctix(ix)).Lat)]);
    %th = text(lon,lat,upper(geodat(ctix(ix)).NAME_2),'HorizontalAlignment','left','FontSize',fontsz);
    th = text(lon,lat,upper(geodat(ctix(ix)).NAME_2),'HorizontalAlignment','center','FontSize',fontsz);
  end;






function [fh,ah,lh,th,bath,isocs,isoch] = county_map(bbox,states,counties,isobaths,bath,fh,varargin)
%function [fh,ah,lh,th,bath,isocs,isoch] = county_map([bbox[,states[,counties[,isobaths[,bath[,fh[,NM1,VAL1,...]]]]]]])
%
% Load GADM36 shape file(s) and draw a high-resolution map showing U.S. State
% and County boundaries. Optional args: 'BBOX' (e.g., [MINX,MAXY,MINY,MAXY]),
% STATES (DFLT: string({"Florida"})), COUNTIES (DFLT: string({'Palm Beach',
% 'Broward', 'Miami-Dade', 'Monroe'})), ISOBATHS, BATH (STRUCT containing a
% .ngdc_hires_bathy field), FH Figure handle, and NM,VALUE pairs (e.g., for
% 'FontSize' (D: 18), 'LineWidth' (D: 1.5) for isobaths, county 'FaceColor'
% (D: cmu.colors('light gray')), 
%
% CALLS: @Map/SHAPEREAD, GEOSHOW; @Vision/BBOXOVERLAPRATIO TEXT. If ISOBATHS
% is non-empty, calls PLOT_HIRES_BATHYMETRY. If ISOBATHS is non-empty, but
% BATH is missing or empty, calls READ_HIRES_BATHYMETRY first.
%
% Last Saved Time-stamp: <Thu 2018-11-15 12:34:39 Eastern Standard Time gramer>

  if ( ~exist('bbox','var') || isempty(bbox) )
    %bbox = [-88,-78,+24,+32];
    bbox = [-88,-78,+24.5,+27.5];
  end;
  if ( ~exist('states','var') || isempty(states) )
    states = string({'Florida'});
  end;
  if ( ~exist('counties','var') || isempty(counties) )
    counties = string({'Palm Beach','Broward','Miami-Dade','Monroe'});
  end;
  if ( ~exist('isobaths','var') )
    isobaths = [];
  end;
  if ( ~exist('bath','var') )
    bath = [];
  end;
  if ( ~exist('fh','var') )
    fh = [];
  end;











function [fh,ah,lh,th,bath,isocs,isoch] = county_map(bbox,states,counties,isobaths,bath,fh,varargin)
%function [fh,ah,lh,th,bath,isocs,isoch] = county_map([bbox[,states[,counties[,isobaths[,bath[,fh[,NM1,VAL1,...]]]]]]])
%
% Load GADM36 shape file(s) and draw a high-resolution map showing U.S. State
% and County boundaries. Optional args: BBOX (e.g., [MINX,MAXY,MINY,MAXY]),
% STATES (DFLT: string({"Florida"})), COUNTIES (DFLT: string({'Palm Beach',
% 'Broward', 'Miami-Dade', 'Monroe'})), ISOBATHS, BATH (STRUCT containing a
% .ngdc_hires_bathy field), FH Figure handle, and NM,VALUE pairs (e.g., for
% 'FontSize' (D: 18), 'LineWidth' (D: 1.5) for isobaths, county 'FaceColor'
% (D: cmu.colors('light gray')), 
%
% CALLS: SHAPEREAD, BBOXOVERLAPRATIO, GEOSHOW, TEXT. If BATH is empty, but
% ISOBATHS is non-empty, calls READ_HIRES_BATHYMETRY. If ISOBATHS non-empty,
% calls PLOT_HIRES_BATHYMETRY.
%
% Last Saved Time-stamp: <Mon 2018-10-15 14:50:23 Eastern Daylight Time gramer>

  if ( ~exist('bbox','var') || isempty(bbox) )
    %bbox = [-88,-78,+24,+32];
    bbox = [-88,-78,+24.5,+27.5];
  end;
  if ( ~exist('states','var') || isempty(states) )
    states = string({'Florida'});
  end;
  if ( ~exist('counties','var') || isempty(counties) )
    counties = string({'Palm Beach','Broward','Miami-Dade','Monroe'});
  end;
  if ( ~exist('isobaths','var') )
    isobaths = [];
  end;
  if ( ~exist('bath','var') )
    bath = [];
  end;
  if ( ~exist('fh','var') )
    fh = [];
  end;


  % figure; usamap('FL')
  % floridahi = shaperead('usastatehi', 'UseGeoCoords', true,...
  %                       'Selector',{@(name) strcmpi(name,'Florida'), 'Name'});
  % geoshow(floridahi, 'FaceColor', [0.3 1.0, 0.675])
  % textm(floridahi.LabelLat, floridahi.LabelLon, floridahi.Name,...
  %       'HorizontalAlignment', 'center')
  % tightmap;
  
  % cties = shaperead('d:/ecoforecasts/coast/tl_2012_us_county');
  % flix = find(string({cties.STATEFP})=='12');
  
  % geodat = shaperead('d:/ecoforecasts/coast/gadm36_USA_2','UseGeoCoords',true,'Selector',{@(name)strcmpi('Florida',name),'NAME_1'});
  % latlims = bboxes(2,2:2:end);
  % sflix = find(latlims<28.0);
  % ctix = find(ismember(string({geodat.NAME_2}),string({'Miami-Dade','Broward'})));
  % fmg; 
  % geoshow(geodat(sflix));

  geodat = shaperead(get_ecoforecasts_path('coast/gadm36_USA_2'),'UseGeoCoords',true,'Selector',{@(name)ismember(name,states),'NAME_1'});
  
  bboxes = [geodat.BoundingBox];
  inix = find(bboxOverlapRatio(bbox2rect(bbox),bbox2rect(bboxes))>0);
  ctix = intersect(inix,find(ismember(string({geodat.NAME_2}),counties)));
  
  if ( isempty(fh) )
    fh = fmg;
  end;
  figure(fh);
  ah = gca;
  set(ah,'FontSize',18);
  %usamap(states);
  lh = geoshow(geodat(inix),'FaceColor',cmu.colors('light gray'));
  daspect([1,cosd(26),1]);
  for ix=1:numel(ctix)
    %lon = nanmin(geodat(ctix(ix)).Lon)+0.05;
    %lat = mean([nanmin(geodat(ctix(ix)).Lat),nanmax(geodat(ctix(ix)).Lat)]);
    %th = text(lon,lat,upper(geodat(ctix(ix)).NAME_2),'HorizontalAlignment','center','FontSize',18);
    lon = mean([nanmin(geodat(ctix(ix)).Lon),nanmax(geodat(ctix(ix)).Lon)]);
    lat = mean([nanmin(geodat(ctix(ix)).Lat),nanmax(geodat(ctix(ix)).Lat)]);
    %th = text(lon,lat,upper(geodat(ctix(ix)).NAME_2),'HorizontalAlignment','left','FontSize',18);
    th = text(lon,lat,upper(geodat(ctix(ix)).NAME_2),'HorizontalAlignment','center','FontSize',18);
  end;

  if ( ~isempty(isobaths) )
    if ( isempty(bath) )
      bath.station_name = 'County_Map';
      bath.lon = nanmean([geodat(inix).Lon]);
      bath.lat = nanmean([geodat(inix).Lat]);
      [radx,rady] = geographic_radius_m([geodat(inix).Lon],[geodat(inix).Lat]);
      bath = read_hires_bathymetry(bath,[radx,rady],[],false);
    end;
    [bath,isocs,isoch] = plot_hires_bathymetry(bath,isobaths,[],true,@contour,false,[],[],fh);
  end;

return;








%%%% ??? HACK HACK HACK
  ylim([24.5,27.5]);
  ah = gca;





From BBOX2RECT.m:
  for ixix=1:numel(boxix)
    ix = boxix(ixix);
    bbox = bboxes(:,[ix,ix+1]);
    if ( bbox(3)<=bbox(1) || bbox(4)<=bbox(2) )
      error('BBOX(%d) should have form [MINX,MAXX ; MINY,MAXY]',ixix);
    end;
    rect(ixix,1:4) = [bbox(1),bbox(3),bbox(2)-bbox(1),bbox(4)-bbox(3)];
  end;





  fldat = shaperead('d:/ecoforecasts/coast/gadm36_USA_2','UseGeoCoords',true,'Selector',{@(name)ismember(name,states),'NAME_1'});
  
  bboxes = [fldat.BoundingBox];
  inix = find(bboxOverlapRatio(bbox2rect(bbox),bbox2rect(bboxes))>0);

  lonlims = bboxes(:,1:2:end);
  latlims = bboxes(:,2:2:end);
  
  inix = bboxinside(lonlims,latlims,bbox);
  ctix = find(ismember(string({fldat(inix).NAME_2}),counties));
  








  if ( ~isempty(needix) )

    disp(['Loading original data from ',baseurl]);

    % Grid-point radii for time-series fields (e.g., seatemp_field)
    % Choose one axis slightly larger to guard against transposition errors
    yrad = 4;
    xrad = 5;






monn = datenum(0,1:12,1);
mons = datestr(monn,'mmm');
season_names = [mons(1:3,1)','    ';mons(4:6,1)','    ';mons(7:9,1)','    ';mons(10:12,1)','    '];
long_season_names = [mons(1,:),'-',mons(3,:);mons(4,:),'-',mons(6,:);mons(7,:),'-',mons(9,:);mons(10,:),'-',mons(12,:)];

  season_names      = ['JFM    ';'AMJ    ';'JAS    ';'OND    '];
  long_season_names = ['Jan-Mar';'Apr-Jun';'Jul-Sep';'Oct-Dec'];


  season_names      = ['JFM    ';'AMJ    ';'JAS    ';'OND    '];
  long_season_names = ['Jan-Mar';'Apr-Jun';'Jul-Sep';'Oct-Dec'];





  % Angle to rotate XLabels
  
  xrotn = nan;
xrotn = 90;


  if ( ~isnan(xrotn) )
    set(gca,'XTickLabelRotation',90);
  end;






dN = dt*seatemp_to_nutes(ts.prof(:,end));


if ( args.factfile.strip() == "" ):
    pass;
elif ( args.factfile == '-' ):
    pass;
else:
    factf.close();
    if ( len(stations_processed) > 0 ):
        if ( not os.path.isfile(full_factfile) ):
            if ( not args.quiet ):
                print(Fore.LIGHTRED_EX + "FACTFILE NOT A FILE?? " + full_factfile);
            raise IOError(-3,"Not a file",full_factfile);
        elif ( os.path.getsize(full_factfile) == 0 ):
            if ( not args.quiet ):
                print(Fore.LIGHTRED_EX + "Processed %d station(s) but EMPTY factfile %s" %(len(stations_processed),full_factfile));
        elif ( args.kedir == args.arcdir ):
            if ( not args.suppress ):
                print("Leaving facts for %d station(s) in %s" %(len(stations_processed),full_factfile));
        else:
            if ( not args.suppress ):
                print("Copying facts for %d station(s) to %s" %(len(stations_processed),args.kedir));
            shutil.copy2(full_factfile,args.kedir);
            

if ( args.write_mat_file ):
    if ( not args.suppress ):
        print("Saving time series to " + args.write_mat_file);
    scipy.io.savemat(args.write_mat_file, mdict={'stns': stns});








880x1320:

The Florida Keys region is designed to show the lower West part of Florida including Florida Bay and is bounded by these coordinates: 26°N 24°N 80°W and 83°W.

For day passes, there are nine image products produced in two different processing streams. In the 'SeaDAS' stream, products include: a chlor_a (Chlorophyll a) image, an ergb (Enhanced RGB) image, a flh (Fluorescence Line Height) image, and a sst (Sea Surface Temperature) image. In a second unique RRC processing stream, products include: a ci (Color Index) image, an efai (Enhanced Floating Algae Index) image, a fai (Floating Algae Index) image, a flh (Fluorescence Line Height) image, and a normal rgb image. During night passes, two products are producted: a sst (Sea Surface Temperature) image and sst4 (Sea Surface Temperature) image.

Meris data is available for this area, however Envisat's end of mission was declared on May 8th, 2012 after communications with the satellite were lost on April 8th, 2012. Meris products include: a chlor_a (Chlorophyll a) image, an ergb (Enhanced RGB) image, a flh (Fluorescence Line Height) image, an experimental MCI image, and a normal rgb image.

If you are viewing 'color' images, for any week but the current week, and click on a Google Earth link, Florida's FWC Karenia brevis data will be displayed as a layer. Since FWC makes the data available on Friday, current week images do not have an association to K. brevis data until Saturday. The K. brevis data displayed corresponds to the date of the images you are viewing. All images during any one week are linked to that week's FWC K. brevis data.

You will notice that the most current imagery date is displayed. If there are several passes, they will be seen inside each of the tabs.

All images are mapped to a cylindrical equidistant projection. Images are at 250 meter resolution.

All images are mapped to a cylindrical equidistant projection.

More information / references for any of the products can be found in the 'information' link located beneath every image.





if ( gfshome == '/cygdrive/d/ecoforecasts' ):
    default_cfgdir = os.path.join(os.path.abspath(gfshome));
    default_arcdir = os.path.join(os.path.abspath(gfshome),"bb");
    default_kedir = os.path.join(os.path.abspath(gfshome),"bb");
else:
    default_cfgdir = os.path.join(os.path.abspath(gfshome),"cfg","gfs");
    default_arcdir = os.path.join(os.path.abspath(gfshome),"archive","gfs");
    default_kedir = os.path.join(os.path.abspath(gfshome),"ke","gfs");





                    warnings.filterwarnings("ignore",category=ResourceWarning);
                    warnings.resetwarnings();

                    warnings.filterwarnings("ignore",category=ResourceWarning);
                    warnings.resetwarnings();


# Caller chooses whether to (also?) save data as time series in a MAT file
if ( args.write_mat_file.strip().lower() in ['','0','default','false','none'] ):
    args.write_mat_file = False;
elif ( not re.search('[.][mM][aA][tT]$', args.write_mat_file) ):
    args.write_mat_file = re.sub('[.][bB][bB]$','.mat',full_factfile);





else:
    write_mat_file = True;
    full_cfgfile = os.path.join(args.cfgdir,args.cfgfile);




                        try:
                            print("(%s %s %s %s %s00 %s u %.2f)" %(station,dtstr,yrstr,jdstr,hrstr,args.inst,mu), file=factf);
                            print("(%s %s %s %s %s00 %s v %.2f)" %(station,dtstr,yrstr,jdstr,hrstr,args.inst,mv), file=factf);
                            print("(%s %s %s %s %s00 %s speed %.2f)" %(station,dtstr,yrstr,jdstr,hrstr,args.inst,ms), file=factf);
                            print("(%s %s %s %s %s00 %s dir %.1f)" %(station,dtstr,yrstr,jdstr,hrstr,args.inst,md), file=factf);
                            print("(%s %s %s %s %s00 %s airtemp %.2f)" %(station,dtstr,yrstr,jdstr,hrstr,args.inst,mt), file=factf);
                            stations_processed[station] = stations_processed.get(station,0) + 1;
                        except:
                            #DEBUG:
                            code.interact(local=locals());





                except exception as Ex:
                    # If we're forced to give up, clean up after ourselves...
                    try: del(dap_data);
                    except: pass;
                    try: pyses.close(); del(pyses);
                    except: pass;
                    try:
                        if ( args.factfile != '-' ): factf.close();
                    except: pass;
                    print('All cleaned up!');
                    raise Ex;








            #dap_s = dap_data["Wind_speed_height_above_ground"];

                    #dap_s_fill = dap_s.attributes["_FillValue"];

                    #dap_s_fill = np.nan;

                    #gs = dap_s.array[timeidx,layers,latidx,lonidx].data;

                    #gs = dap_s.array[timeidx,layers,latidx,lonidx].data;

                    #ms = gs[gs != dap_s_fill].mean();



                    #DEBUG:                    print("MT:";) print(tix); print(dap_tm[tix]); print(dap_tm[tix][0]);
                    #DEBUG:                    print("Lat:"); print(latidx); print(all_lats[latidx]);
                    #DEBUG:                    print("Lon:"); print(lonidx); print(all_lons[lonidx]);
                    #DEBUG:                    print("X:"); print(dap_data["X"][lonidx]);
                    #DEBUG:                    print("Y:"); print(dap_data["Y"][latidx]);

                    #DEBUG:                    print("Depth:"); print(dap_data["Depth"][layers]);
                    #DEBUG:                    print("Dt:"); print(dtstr);
                    #DEBUG:                    print("U:"); print(gu);
                    #DEBUG:                    print("T:"); print(gt);






                except (AttributeError) as AE:
                    # AttributeError is usually what PyDAP throws from my coding errors
                    # Enter a debugging session a la MATLAB 'keyboard'
                    code.interact(local=locals());





            try:
            except (AttributeError) as AE:
            # AttributeError is the most likely coding failure from PyDAP





all_dates = [startdate + datetime.timedelta(days=x) for x in range(0,args.daystoget+1)]




parser.add_argument('startdate', metavar='YYYYMMDD', type=int, nargs='?',
                    help='Date from which to factize GFS data (DEFAULT: DAYSTOGET before the latest time index).');


# Default STARTDATE is normally DAYSTOGET before last full day
if ( args.startdate is None ):
    startdate = datetime.datetime(starting_time.year,starting_time.month,starting_time.day);
    startdate = startdate - datetime.timedelta(days=args.daystoget+1);
    args.startdate = startdate.strftime('%Y%m%d');
else:
    startdate = datetime.datetime.strptime(str(args.startdate),'%Y%m%d');

all_dates = [startdate + datetime.timedelta(days=x) for x in range(0,args.daystoget+1)]








all_dates = [startdate - datetime.timedelta(days=x-1) for x in range(args.daystoget+1, 0, -1)]






parser = configargparse.ArgumentParser(description='Generate Knowledge Engine facts (e.g., for PyKnow or G2) from Global Forecast System (GFS) or other weather model data extracted via OPeNDAP (PyDAP).',
                                       epilog='Read individual station MODEL_GRID lat/lon *indices* from CFGFILE, unless USE_STATIONS_FILE is set to "true" or a filename: then a stations.txt file is read with actual latitudes and latitudes (if "true", the default file to use is ' + default_stations_file + '). If a MODEL other than the default (' + default_model + ') or a GRID other than one of the defaults (e.g., ' + default_grid + ') is specified, then specify either another options file, or else at minimum another CFGFILE (or USE_STATIONS_FILE). One fact is generated per variable for each model-period in DAYSTOGET; alternatively, caller can specify OPeNDAP-style TIMEIDX (START:STRIDE:STOP) for the model and grid directly. Facts are generated in the "classic" CLIPS/G2 format for all variables specified (default - the basic five u,v,speed,dir,airtemp from GFS), and upon success, fact file is copied from ARCDIR into KEDIR. SUGGESTED USAGE for multiple models: Create an options file with defaults FOR EACH MODEL, then use -o (see above) to point to that options file.',
                                       default_config_files=[etc_options_file,default_options_file],
                                       fromfile_prefix_chars='@');




# Note: We expect a string here
parser.add_argument('-t','--timeidx', metavar='TIMEIDX#', type=str,
                    help='Alternative to DAYSTOGET: GFS time indices (STRING: start:stride:stop) to retrieve (DEFAULT: "%(default)s")');




# Parse FORTRAN-like (MATLAB-like, but 0-based) slice strings into Python slice() instances

layers = efsu.matlab_slice(args.layers);

#DEBUG:print(args.timeidx);
if ( args.timeidx is None ):
    timeidx = slice(-args.daystoget-1,None,None);
    
else:
    if ( args.daystoget != 0 ):
        if ( not args.quiet ):
            print(Fore.LIGHTRED_EX + "TIMEIDX specified: Ignoring DAYSTOGET ",args.daystoget);
    timeidx = efsu.matlab_slice(args.timeidx);





# If user specified STARTDATE (but not TIMEIDX)
if ( args.startdate is None ):
    #args.startdate = starting_time.strftime('%Y%m%d');
    args.startdate = '20150115';

if ( args.startdate is not None ):
    if ( args.timeidx is not None ):
        print(Fore.LIGHTRED_EX + "TIMEIDX specified: Ignoring STARTDATE ",args.startdate);
    else:
        startdate = datetime.datetime.strptime(str(args.startdate),'%Y%m%d');
        # #efsu.range_check(startdate,dap_tm[0],dap_tm[-1]);
        # # Give the caller a friendlier error...
        # if ( dap_tm[0] > startdate or startdate > dap_tm[-1] ):
        #     #DEBUG:            print(startdate,dap_tm[0][0]);
        #     raise ValueError("STARTDATE outside valid dates for %s_%s: %s before %s"%(args.model,args.expt,startdate,dap_tm[0][0]));
        # startidx = np.abs(dap_tm[:]-startdate).argmin();
        # if ( args.daystoget > dap_tm.shape[0]-startidx ):
        #     args.daystoget = dap_tm.shape[0]-startidx;
        #     print(Fore.LIGHTRED_EX + "Truncating DAYSTOGET to %s" %repr(args.daystoget));
        # timeidx = slice(startidx,startidx+args.daystoget+1,None);

all_dates = [startdate - datetime.timedelta(days=x-1) for x in range(args.daystoget+1, 0, -1)]





#for tix in range(*timeidx.indices(dap_tm.shape[0])):









    elif ( args.baseurl.strip() == "" ):


    if ( args.baseurl.strip() == "" ):





# Sample URLs:
#  dataurl = 'http://hycom.coaps.fsu.edu:8080/thredds/dodsC/glb_nrt_analysis.ascii';
#  dataurl = 'http://tds.hycom.org/thredds/dodsC/GLBa0.08/expt_91.2';
if ( args.expt.strip() == "" ):
    dataurl = args.baseurl + '/' + args.model;
else:
    dataurl = args.baseurl + '/' + args.model + '/' + args.expt;

if ( args.baseurl is None ):
    if ( args.model == 'gfs_0p25' ):
        args.baseurl = 'https://rda.ucar.edu/thredds/dodsC/files/g/ds084.1/';
    elif ( args.model == 'gfs_0p5' ):
        args.baseurl = 'https://www.ncei.noaa.gov/thredds/dodsC/gfs-g4-anl-files/';
    else:
        raise ValueError('Unknown URL for model code %s' %args.model);


if ( args.expt is None ):
    if ( args.model == 'gfs_0p25' or args.model == 'gfs_0p5' ):
        args.expt = '';
    else:
        raise ValueError('Which EXPT should I use for model %s?' %(args.model));

# Sample URLs:
#  dataurl = 'http://hycom.coaps.fsu.edu:8080/thredds/dodsC/glb_nrt_analysis.ascii';
#  dataurl = 'http://tds.hycom.org/thredds/dodsC/GLBa0.08/expt_91.2';
if ( args.expt.strip() == "" ):
    dataurl = args.baseurl + '/' + args.model;
else:
    dataurl = args.baseurl + '/' + args.model + '/' + args.expt;










# # Enter a debugging session a la MATLAB 'keyboard'
# code.interact(local=locals());





ds=face.shallow.depths
preix=find(ds.date<datenum(2015,01,30,15,30,0));
midix=find(ds.date>datenum(2015,01,30,15,30,0) & ds.date<datenum(2015,06,23,14,30,0));
pstix=find(ds.date>=datenum(2015,06,23,14,30,0));
ds.data(preix) = ds.data(preix) - nanmean(ds.data(preix));
ds.data(midix) = ds.data(midix) - nanmean(ds.data(midix));
ds.data(pstix) = ds.data(pstix) - nanmean(ds.data(pstix));



for ix=1:numel(stns);
  station=stns(ix);
  disp(station.station_name);
  save(fullfile(get_ecoforecasts_path('data'),['factize_gfs_',char(lower(station.station_name))]),'station','-v7.3');
  station=[]; clear station;
end;






for ix=1:numel(stns);
  stns(ix).u.date(end)=[]; stns(ix).u.data(end)=[]; 
  stns(ix).v.date(end)=[]; stns(ix).v.data(end)=[]; 
  stns(ix).speed.date(end)=[]; stns(ix).speed.data(end)=[]; 
  stns(ix).dir.date(end)=[]; stns(ix).dir.data(end)=[]; 
  stns(ix).airtemp.date(end)=[]; stns(ix).airtemp.data(end)=[]; 
  stns(ix)=verify_variable(stns(ix),'speed_3_d_avg');
  fh=fmg; boxplot_ts(stns(ix).speed_3_d_avg); waitfor(fh);
end;




for stix=1:numel(STATIONS.codes)
  stnm = lower(string(STATIONS.codes{stix}));
  stns(stix).station_name = stnm;
  stns(stix).lon = STATIONS.lons(stix);
  stns(stix).depth = STATIONS.depths(stix);
  if ( isempty(regexp(stnm(1),'[0-9]')) )
    stns(stix).u.date = repmat(nan,[numel(alldts)*numel(allhrs)*numel(allfhrs),1]);
    stns(stix).u.data = stns(stix).u.date;
    stns(stix).v = stns(stix).u;
    stns(stix).speed = stns(stix).u;
    stns(stix).dir = stns(stix).u;
    stns(stix).airtemp = stns(stix).u;
  end;
end;



STATIONS.u.date = repmat(nan,[numel(alldts)*numel(allhrs)*numel(allfhrs),1]);
STATIONS.u.data = STATIONS.u.date;
STATIONS.v = STATIONS.u;
STATIONS.speed = STATIONS.u;
STATIONS.dir = STATIONS.u;
STATIONS.airtemp = STATIONS.u;



%alldts = datenum(2004,3,2):floor(now-1);




alldts = datenum(2018,7,20,12,0,0);

% Forecast hours
%allfhrs = 0:6:18;
allfhrs = 6;

STATIONS = get_all_station_metadata;

% >> url = 'https://www.ncei.noaa.gov/thredds/dodsC/gfs-g4-anl-files/201707/20170720/gfsanl_4_20170720_1200_006.grb2';
% >>       url = sprintf('%s%04d%02d/%04d%02d%02d/gfsanl_4_%04d%02d%02d_%02d00_%03d.grb2',BASEURL,yr,mo,yr,mo,dy,yr,mo,dy,hr,fhr);
% >> url
% url =
% http://www.ncei.noaa.gov/thredds/dodsC/gfs-g4-anl-files/201707/20170720/gfsanl_4_20170720_1200_006.grb2


for dtix = 1:numel(alldts)
  dt = alldts(dtix);
  yr = get_year(dt);
  mo = get_month(dt);
  dy = get_dom(dt);
  hr = get_hour(dt);

  fname = sprintf('gfs_0p5_%04d%02d%02d.bb',yr,mo,dy);

  for fhrix = 1:numel(allfhrs)
    fhr = allfhrs(fhrix);





# >>> from pydap.client import open_url
# >>> from pydap.cas.urs import setup_session
# >>> dataset_url = 'http://server.example.com/path/to/dataset'
# >>> session = setup_session(username, password, check_url=dataset_url)
# >>> dataset = open_url(dataset_url, session=session)


# import matplotlib;
# # ImportError: Matplotlib qt-based backends require an external PyQt4, PyQt5,
# #  or PySide package to be installed, but it was not found.
# matplotlib.use('TkAgg');
# import matplotlib.pyplot as plt;
# #matplotlib.use('WxAgg');
# #matplotlib.use('Qt4Agg');
# #matplotlib.use('Qt5Agg');





    if ( ~is_valid_ts(ts{tix}) )
      error('TS{%d} is not a valid time series!',tix);
    end;




compare_slopes(lkwf1,'lkwf1'); compare_slopes(fwyf1,'fwyf1'); compare_slopes(rsm,'rsm'); compare_slopes(pvgf1,'pvgf1');  compare_slopes(ftpf1,'ftpf1'); compare_slopes(fdepk,'fdepk');  compare_slopes(canf1,'canf1');  compare_slopes(cnnf1,'cnnf1'); compare_slopes(face.shallow,'face.shallow');  compare_slopes(face.tcm1,'face.tcm1'); compare_slopes(face.tcm2,'face.tcm2');  compare_slopes(face.deep,'face.deep'); compare_slopes(face.HanS_HW,'face.HanS_HW');  compare_slopes(face.HanS_BR,'face.HanS_BR'); compare_slopes(face.AOML_HW,'face.AOML_HW');  compare_slopes(face.AOML_BR,'face.AOML_BR'); compare_slopes(face.brwd20,'face.brwd20');  compare_slopes(face.brwd40,'face.brwd40'); compare_slopes(face.brwd100,'face.brwd100');  compare_slopes(face.boca20,'face.boca20'); compare_slopes(face.boca40,'face.boca40'); compare_slopes(sefcri.dc1,'sefcri.dc1'); compare_slopes(sefcri.dc3,'sefcri.dc3'); compare_slopes(sefcri.bc1,'sefcri.bc1'); compare_slopes(sefcri.bc3,'sefcri.bc3');  compare_slopes(sefcri.pb1,'sefcri.pb1'); compare_slopes(sefcri.pb2,'sefcri.pb2');  compare_slopes(sefcri.mc2,'sefcri.mc2'); compare_slopes(sefcri.pela,'sefcri.pela');  compare_slopes(sefcri.evca,'sefcri.evca'); compare_slopes(sefcri.slib,'sefcri.slib');  compare_slopes(sefcri.updb,'sefcri.updb'); compare_slopes(sfomc.nw_w_btm,'sfomc.nw_w_btm');  compare_slopes(sfomc.ne_buoy,'sfomc.ne_buoy'); compare_slopes(sfomc.e_buoy,'sfomc.e_buoy');  compare_slopes(sfomc.c_buoy,'sfomc.c_buoy'); compare_slopes(sfomc.sw_buoy,'sfomc.sw_buoy');  compare_slopes(sfomc.se_buoy,'sfomc.se_buoy');





function [bet,ang,iso,dep,stn] = find_ngdc_slope_station(stn,rad_m,avg_m,fdif_m,method,extrapval)
%function [bet,ang,iso,dep,stn] = find_ngdc_slope_station(stn|{stn,bath}[,rad_m[,avg_m[,fdif_m[,method[,extrapval]]]]])
%
% Find seafloor slope BET from high-resolution bathymetry for station STN. If
% not present, calls READ_HIRES_BATHYMETRY to create STN.ngdc_hires_bathy,
% and calls FIND_NGDC_SLOPE to add fields .beta_deg and .beta to it.
%
% If AVG_M is given, calls INTERP_FIELD (v.) to downscale bathymetry to that
% resolution.
%
% ANG, angle of seafloor slope in degrees clockwise from True North, and
% local isobath angle ISO can also be returned (this is MOD(ANG-90,360)).
%
% Last Saved Time-stamp: <Wed 2018-07-25 13:22:52 Eastern Daylight Time gramer>

  if ( ~exist('stn','var') || (~isstruct(stn) && ~iscell(stn)) )
    error('First arg must be given, either STRUCT STN or cell {STN,BATH}');
  end;
  bath = [];
  if ( iscell(stn) )
    bath = stn{2};
    stn = stn{1};
  end;
  if ( ~exist('rad_m','var') || isempty(rad_m) || all(isnan(rad_m)) )
    rad_m = 1e3;
  end;
  if ( ~exist('avg_m','var') )
    avg_m = [];
  end;
  if ( ~exist('fdif_m','var') || isempty(fdif_m) )
    fdif_m = nan;
  end;
  if ( ~exist('method','var') )
    method = [];
  end;
  if ( ~exist('extrapval','var') )
    extrapval = [];
  end;

  minrad = nanmin([rad_m(:)',avg_m,fdif_m]);

  if ( numel(rad_m) == 1 )
    radx = nanmax([rad_m,fdif_m]);
    rady = nanmax([rad_m,fdif_m]);
  elseif ( numel(rad_m) == 2 )
    radx = nanmax([rad_m(1),fdif_m]);
    rady = nanmax([rad_m(2),fdif_m]);
  else
    error('Radius RAD_M must be a numeric 1- or 2-vector in [m]');
  end;

  if ( isempty(bath) )
    if ( ~isfield(stn,'ngdc_hires_bathy') )
      stn = read_hires_bathymetry(stn,[radx,rady]);
    end;
    bath = stn.ngdc_hires_bathy;
  end;
  res_m = min([bath.xres(:)',bath.yres(:)']);

  avgpts = ceil(avg_m / res_m);
  if ( avgpts > 1 )
    % We now know AVG_M is not NaN
    interpMethod = {@nanmean,avgpts};
    lons = bath.lon(1:avgpts:end);
    lats = bath.lat(1:avgpts:end);
    [LONS,LATS] = meshgrid(lons,lats);
    fld = interp_field(bath.lat,bath.lon,bath.field,LATS,LONS,interpMethod);
    bath.lon = lons;
    bath.lat = lats;
    lat_res_deg = min(diff(unique(LATS(:))));
    lon_res_deg = min(diff(unique(LONS(:))));
    bath.xres = distance_wgs84(stn.lat,stn.lon,stn.lat,stn.lon+lon_res_deg)*1e3;
    bath.yres = distance_wgs84(stn.lat,stn.lon,stn.lat+lat_res_deg,stn.lon)*1e3;
    bath.files = {'@nanmean',num2str(avgpts)};
    bath.field = reshape(fld,[numel(lats),numel(lons)]);
    res_m = min([bath.xres,bath.yres]);
    LONS=[]; LATS=[]; clear LONS LATS
    fld=[]; clear fld
  end;
  if ( isnan(fdif_m) )
    npts = 11;
  else
    npts = ceil(fdif_m / res_m);
    npts = interp1([2,3:2:11],[2,3:2:11],npts,'nearest','extrap');
  end;

  [bet,ang,iso,bath] = find_ngdc_slope(bath,stn.lon,stn.lat,npts,method,extrapval);
  if ( nargout > 3 )
    dep = interp2(bath.lon,bath.lat,bath.field,stn.lon,stn.lat,'nearest');
  end;
  if ( nargout > 4 )
    stn.ngdc_hires_bathy = bath;
  end;
  bath=[]; clear bath

return;











if ( ~isfield(lkwf1,'slope') )
  [lkwf1.slope,lkwf1.slope_orientation,lkwf1.isobath_orientation,lkwf1.ngdc_hires_bathy] = ...
      find_ngdc_slope(lkwf1.ngdc_hires_bathy,lkwf1.lon,lkwf1.lat,3);
end;
if ( ~isfield(fwyf1,'slope') )
  [fwyf1.slope,fwyf1.slope_orientation,fwyf1.isobath_orientation,lkwf1.ngdc_hires_bathy] = ...
      find_ngdc_slope(lkwf1.ngdc_hires_bathy,fwyf1.lon,fwyf1.lat,3);
end;
if ( ~isfield(rsm,'slope') )
  [rsm.slope,rsm.slope_orientation,rsm.isobath_orientation,lkwf1.ngdc_hires_bathy] = ...
      find_ngdc_slope(lkwf1.ngdc_hires_bathy,rsm.lon,rsm.lat,3);
end;
if ( ~isfield(pvgf1,'slope') )
  [pvgf1.slope,pvgf1.slope_orientation,pvgf1.isobath_orientation,lkwf1.ngdc_hires_bathy] = ...
      find_ngdc_slope(lkwf1.ngdc_hires_bathy,pvgf1.lon,pvgf1.lat,3);
end;
if ( ~isfield(ftpf1,'slope') )
  % [ftpf1.slope,ftpf1.slope_orientation,ftpf1.isobath_orientation,lkwf1.ngdc_hires_bathy] = ...
  %     find_ngdc_slope(lkwf1.ngdc_hires_bathy,ftpf1.lon,ftpf1.lat,3);
  % Once we fixed a bug in READ_HIRES_BATHYMETRY, suddenly station FTPF1 was
  % no longer within "120e3" meters of LKWF1!
  [ftpf1.slope,ftpf1.slope_orientation,ftpf1.isobath_orientation,ftpf1.ngdc_hires_bathy] = ...
      find_ngdc_slope(ftpf1.ngdc_hires_bathy,ftpf1.lon,ftpf1.lat,3);
end;
if ( ~isfield(fdepk,'slope') )
  [fdepk.slope,fdepk.slope_orientation,fdepk.isobath_orientation,lkwf1.ngdc_hires_bathy] = ...
      find_ngdc_slope(lkwf1.ngdc_hires_bathy,fdepk.lon,fdepk.lat,3);
end;
if ( ~isfield(canf1,'slope') )
  [canf1.slope,canf1.slope_orientation,canf1.isobath_orientation,canf1.ngdc_hires_bathy] = ...
      find_ngdc_slope(canf1.ngdc_hires_bathy,canf1.lon,canf1.lat,3);
end;
if ( ~isfield(cnnf1,'slope') )
  [cnnf1.slope,cnnf1.slope_orientation,cnnf1.isobath_orientation,canf1.ngdc_hires_bathy] = ...
      find_ngdc_slope(canf1.ngdc_hires_bathy,cnnf1.lon,cnnf1.lat,3);
end;

if ( exist('face','var') && ~isfield(face.shallow,'slope') )
  [face.shallow.slope,face.shallow.slope_orientation,face.shallow.isobath_orientation,lkwf1.ngdc_hires_bathy] = ...
      find_ngdc_slope(lkwf1.ngdc_hires_bathy,face.shallow.lon,face.shallow.lat,3);
  [face.tcm1.slope,face.tcm1.slope_orientation,face.tcm1.isobath_orientation,lkwf1.ngdc_hires_bathy] = ...
      find_ngdc_slope(lkwf1.ngdc_hires_bathy,face.tcm1.lon,face.tcm1.lat,3);
  [face.tcm2.slope,face.tcm2.slope_orientation,face.tcm2.isobath_orientation,lkwf1.ngdc_hires_bathy] = ...
      find_ngdc_slope(lkwf1.ngdc_hires_bathy,face.tcm2.lon,face.tcm2.lat,3);
  [face.deep.slope,face.deep.slope_orientation,face.deep.isobath_orientation,lkwf1.ngdc_hires_bathy] = ...
      find_ngdc_slope(lkwf1.ngdc_hires_bathy,face.deep.lon,face.deep.lat,3);

  if ( isfield(face,'HanS_HW') )
    [face.HanS_HW.slope,face.HanS_HW.slope_orientation,face.HanS_HW.isobath_orientation,lkwf1.ngdc_hires_bathy] = ...
        find_ngdc_slope(lkwf1.ngdc_hires_bathy,face.HanS_HW.lon,face.HanS_HW.lat,3);
    [face.HanS_BR.slope,face.HanS_BR.slope_orientation,face.HanS_BR.isobath_orientation,lkwf1.ngdc_hires_bathy] = ...
        find_ngdc_slope(lkwf1.ngdc_hires_bathy,face.HanS_BR.lon,face.HanS_BR.lat,3);
    [face.AOML_HW.slope,face.AOML_HW.slope_orientation,face.AOML_HW.isobath_orientation,lkwf1.ngdc_hires_bathy] = ...
        find_ngdc_slope(lkwf1.ngdc_hires_bathy,face.AOML_HW.lon,face.AOML_HW.lat,3);
    [face.AOML_BR.slope,face.AOML_BR.slope_orientation,face.AOML_BR.isobath_orientation,lkwf1.ngdc_hires_bathy] = ...
        find_ngdc_slope(lkwf1.ngdc_hires_bathy,face.AOML_BR.lon,face.AOML_BR.lat,3);
  end;

  if ( isfield(face,'brwd20') )
    [face.brwd20.slope,face.brwd20.slope_orientation,face.brwd20.isobath_orientation,lkwf1.ngdc_hires_bathy] = ...
        find_ngdc_slope(lkwf1.ngdc_hires_bathy,face.brwd20.lon,face.brwd20.lat,3);
    [face.brwd40.slope,face.brwd40.slope_orientation,face.brwd40.isobath_orientation,lkwf1.ngdc_hires_bathy] = ...
        find_ngdc_slope(lkwf1.ngdc_hires_bathy,face.brwd40.lon,face.brwd40.lat,3);
    [face.brwd100.slope,face.brwd100.slope_orientation,face.brwd100.isobath_orientation,lkwf1.ngdc_hires_bathy] = ...
        find_ngdc_slope(lkwf1.ngdc_hires_bathy,face.brwd100.lon,face.brwd100.lat,3);
    [face.boca20.slope,face.boca20.slope_orientation,face.boca20.isobath_orientation,lkwf1.ngdc_hires_bathy] = ...
        find_ngdc_slope(lkwf1.ngdc_hires_bathy,face.boca20.lon,face.boca20.lat,3);
    [face.boca40.slope,face.boca40.slope_orientation,face.boca40.isobath_orientation,lkwf1.ngdc_hires_bathy] = ...
        find_ngdc_slope(lkwf1.ngdc_hires_bathy,face.boca40.lon,face.boca40.lat,3);
  end;
end;

if ( exist('sefcri','var') && ~isfield(sefcri.updb,'slope') )
  [sefcri.dc1.slope,sefcri.dc1.slope_orientation,sefcri.dc1.isobath_orientation,lkwf1.ngdc_hires_bathy] = ...
      find_ngdc_slope(lkwf1.ngdc_hires_bathy,sefcri.dc1.lon,sefcri.dc1.lat,3);
  [sefcri.dc3.slope,sefcri.dc3.slope_orientation,sefcri.dc3.isobath_orientation,lkwf1.ngdc_hires_bathy] = ...
      find_ngdc_slope(lkwf1.ngdc_hires_bathy,sefcri.dc3.lon,sefcri.dc3.lat,3);
  [sefcri.bc1.slope,sefcri.bc1.slope_orientation,sefcri.bc1.isobath_orientation,lkwf1.ngdc_hires_bathy] = ...
      find_ngdc_slope(lkwf1.ngdc_hires_bathy,sefcri.bc1.lon,sefcri.bc1.lat,3);
  [sefcri.bc3.slope,sefcri.bc3.slope_orientation,sefcri.bc3.isobath_orientation,lkwf1.ngdc_hires_bathy] = ...
      find_ngdc_slope(lkwf1.ngdc_hires_bathy,sefcri.bc3.lon,sefcri.bc3.lat,3);
  [sefcri.pb1.slope,sefcri.pb1.slope_orientation,sefcri.pb1.isobath_orientation,lkwf1.ngdc_hires_bathy] = ...
      find_ngdc_slope(lkwf1.ngdc_hires_bathy,sefcri.pb1.lon,sefcri.pb1.lat,3);
  [sefcri.pb2.slope,sefcri.pb2.slope_orientation,sefcri.pb2.isobath_orientation,lkwf1.ngdc_hires_bathy] = ...
      find_ngdc_slope(lkwf1.ngdc_hires_bathy,sefcri.pb2.lon,sefcri.pb2.lat,3);
  [sefcri.mc2.slope,sefcri.mc2.slope_orientation,sefcri.mc2.isobath_orientation,lkwf1.ngdc_hires_bathy] = ...
      find_ngdc_slope(lkwf1.ngdc_hires_bathy,sefcri.mc2.lon,sefcri.mc2.lat,3);
  [sefcri.pela.slope,sefcri.pela.slope_orientation,sefcri.pela.isobath_orientation,lkwf1.ngdc_hires_bathy] = ...
      find_ngdc_slope(lkwf1.ngdc_hires_bathy,sefcri.pela.lon,sefcri.pela.lat,3);
  [sefcri.evca.slope,sefcri.evca.slope_orientation,sefcri.evca.isobath_orientation,lkwf1.ngdc_hires_bathy] = ...
      find_ngdc_slope(lkwf1.ngdc_hires_bathy,sefcri.evca.lon,sefcri.evca.lat,3);
  [sefcri.slib.slope,sefcri.slib.slope_orientation,sefcri.slib.isobath_orientation,lkwf1.ngdc_hires_bathy] = ...
      find_ngdc_slope(lkwf1.ngdc_hires_bathy,sefcri.slib.lon,sefcri.slib.lat,3);
  [sefcri.updb.slope,sefcri.updb.slope_orientation,sefcri.updb.isobath_orientation,lkwf1.ngdc_hires_bathy] = ...
      find_ngdc_slope(lkwf1.ngdc_hires_bathy,sefcri.updb.lon,sefcri.updb.lat,3);
end;

if ( exist('sfomc','var') && ~isfield(sfomc.nw_w_btm,'slope') )
  [sfomc.nw_w_btm.slope,sfomc.nw_w_btm.slope_orientation,sfomc.nw_w_btm.isobath_orientation,lkwf1.ngdc_hires_bathy] = ...
      find_ngdc_slope(lkwf1.ngdc_hires_bathy,sfomc.nw_w_btm.lon,sfomc.nw_w_btm.lat,3);
  [sfomc.ne_buoy.slope,sfomc.ne_buoy.slope_orientation,sfomc.ne_buoy.isobath_orientation,lkwf1.ngdc_hires_bathy] = ...
      find_ngdc_slope(lkwf1.ngdc_hires_bathy,sfomc.ne_buoy.lon,sfomc.ne_buoy.lat,3);
  [sfomc.e_buoy.slope,sfomc.e_buoy.slope_orientation,sfomc.e_buoy.isobath_orientation,lkwf1.ngdc_hires_bathy] = ...
      find_ngdc_slope(lkwf1.ngdc_hires_bathy,sfomc.e_buoy.lon,sfomc.e_buoy.lat,3);
  [sfomc.c_buoy.slope,sfomc.c_buoy.slope_orientation,sfomc.c_buoy.isobath_orientation,lkwf1.ngdc_hires_bathy] = ...
      find_ngdc_slope(lkwf1.ngdc_hires_bathy,sfomc.c_buoy.lon,sfomc.c_buoy.lat,3);
  [sfomc.sw_buoy.slope,sfomc.sw_buoy.slope_orientation,sfomc.sw_buoy.isobath_orientation,lkwf1.ngdc_hires_bathy] = ...
      find_ngdc_slope(lkwf1.ngdc_hires_bathy,sfomc.sw_buoy.lon,sfomc.sw_buoy.lat,3);
  [sfomc.se_buoy.slope,sfomc.se_buoy.slope_orientation,sfomc.se_buoy.isobath_orientation,lkwf1.ngdc_hires_bathy] = ...
      find_ngdc_slope(lkwf1.ngdc_hires_bathy,sfomc.se_buoy.lon,sfomc.se_buoy.lat,3);
  [sfomc.pier_cc.slope,sfomc.pier_cc.slope_orientation,sfomc.pier_cc.isobath_orientation,lkwf1.ngdc_hires_bathy] = ...
      find_ngdc_slope(lkwf1.ngdc_hires_bathy,sfomc.pier_cc.lon,sfomc.pier_cc.lat,3);
end;








    stn.ngdc_hires_bathy.yres = avg_m;



    annot_h(end+1) = ...
        dsuannotation('textarrow',...
                      [double(shax(2)),datenum(2015,9,4,19,20,0),double(shax(1)),datenum(2015,9,7,21,40,0)],...
                      [double(shax(2)),0.35,double(shax(1)),25.5],'String','Onset',...
                      'FontSize',fontsz,'VerticalAlignment','middle',...
                      'LineWidth',5,'LineStyle',':','HeadWidth',fontsz,'HeadLength',fontsz);
    annot_h(end+1) = ...
        dsuannotation('textarrow',...
                      [double(shax(3)),datenum(2015,9,16,0,20,0),double(shax(1)),datenum(2015,9,16,7,0,0)],...
                      [double(shax(3)),0.85,double(shax(1)),25.0],'String','Relaxation',...
                      'FontSize',fontsz,'VerticalAlignment','middle',...
                      'LineWidth',5,'LineStyle',':','HeadWidth',fontsz,'HeadLength',fontsz);






    % annot_h=[];
    % annot_h(end+1) = ...
    %     dsuannotation('textarrow',...
    %                   [double(shax(2)),datenum(2015,7,25,04,20,0),double(shax(1)),datenum(2015,7,28,0,0,0)],...
    %                   [double(shax(2)),0.35,double(shax(1)),25.5],'String','Onset',...
    %                   'FontSize',fontsz,'VerticalAlignment','middle',...
    %                   'LineWidth',5,'LineStyle',':','HeadWidth',fontsz,'HeadLength',fontsz);
    % annot_h(end+1) = ...
    %     dsuannotation('textarrow',...
    %                   [double(shax(3)),datenum(2015,8,2,18,0,0),double(shax(1)),datenum(2015,8,3,12,0,0)],...
    %                   [double(shax(3)),0.85,double(shax(1)),25.5],'String','Relaxation',...
    %                   'FontSize',fontsz,'VerticalAlignment','middle',...
    %                   'LineWidth',5,'LineStyle',':','HeadWidth',fontsz,'HeadLength',fontsz);







    legh=legend('Air (PVGF1)','Deep','TCM2','TCM1','Shallow','(Upwelling)',...
                'Location','SouthEast', 'Orientation','horiz');
    uistack([legh;annot_h(:)]);





function pers = natural_periods(maxper)
%function pers = natural_periods([maxper])
% Convenience function: returns local tidal, diurnal, ocean mesoscale,
% "weather band", biannual and annual periods in days (e.g., for DATENUM).
% If optional MAXPER is a numeric scalar, only returns PERS <= MAXPER.
% NOTE: Use INERTIAL_PERIOD (v.) to get the local inertial period.

  %pers = [0.5 1 4 7 14 42 183 365];
  
  if ( exist('tidefreq') > 1 )
    pers = sort([(1./(tidefreq)./24),1,4,7,42,183,365]);
  else
    % K2,S2,M2,N2,K1,diurnal,P1,O1,oc.mesoscale,min.weather,MF,max.weather,biannual,annual
    pers = [0.4986,0.5,0.5175,0.5274,0.9973,1.0,1.0027,1.0758,4,7,13.6604,42,183,365];
  end;
  
  if  ( exist('maxper','var') && isscalar(maxper) && isnumeric(maxper) )
    pers(pers>maxper) = [];
  end;

return;



function fh = fmg(fh)
%function fh = fmg([fh])
%
% Create new full-screen Figure with HOLD ON, GRID ON, default FONTSIZE 20.
% If optional FH is passed in, pass focus to and do the above to it. From
% 2016 Nov 28 on, also calls TIGHTEN_AXES (v.) on the resulting figure. NOTE
% that after calling COLORBAR (v.), which also inherits 20 pt. fonts, this
% can cause the figure boundary to clip the colorbar label. To avoid this,
% call FMGL (v., "L" for "Loose") in place of FMG, then call COLORBAR_TIGHT
% in place of COLORBAR.
%
% Last Saved Time-stamp: <Wed 2018-06-20 13:44:06 Eastern Daylight Time gramer>

  args = varargin;

  if ( exist('fh','var') )
    fh = figure(fh);
  else
    fh = figure;
  end;
  maximize_graph;
  hold on;
  grid on;
  grid minor;
  set(gca,'FontSize',20);  % For publication-ready fonts when printing

  % NOTE again: for right-labeled COLORBARs, figure may clip the label :(
  tighten_axes(fh);

return;






Predicted spatial variation in coral larval transport within the Florida Keys



  for ix=1:nannots
    if ( exist('xen','var') )
      [xen_n,yen_n] = ds2nfu(xax(ix,1),xen(ix,:),yen(ix,:));
      [xen2,yen2] = ds2nfu(xax(ix,2),xen(ix,:),yen(ix,:));
      annotation(antyp,xen_n,yen_n,args{:});
    elseif ( exist('pos','var') )
      pos_n = ds2nfu(ax(ix),pos(ix,:));
      annotation(antyp,pos_n,args{:});
    else
      error('Invalid second argument??');
    end;
  end;





  switch ( ncols )
   case 2, % Nx2 coordinate matrices X and Y
    xax = repmat(gca,[nannots,1]);
    yax = xax;
    xen = arg1;
    yen = arg2;
   case 4, % [AX1,X1,AX2,X2]
    xax = arg1(:,[1,3]);
    yax = arg2(:,[1,3]);
    xen = arg1(:,[2,4]);
    yen = arg2(:,[2,4]);

   case 3, % Nx3 coordinate matrices X and Y including (numerical) AXES handles
    ax = arg1(:,1);
    xen = arg1(:,2:end);
    yen = arg2(:,2:end);
   case 4, % Nx4 position matrices POS
    ax = repmat(gca,[nannots,1]);
    pos = arg1;
   case 5, % Nx3 position matrices POS including (numerical) AXES handles
    ax = arg1(:,1);
    pos = arg1(:,2:end);
   otherwise,
    error('First arg after ANTYP must be an Nx2, Nx3, Nx4, or Nx5 matrix');
  end;










    switch ( ncols )
     case 2, % Nx2 coordinate matrices X and Y
      xax = repmat(gca,[nannots,1]);
      yax = xax;
      xen = arg1;
      yen = arg2;
     case 4, % [AX1,X1,AX2,X2]
      xax = arg1(:,[1,3]);
      yax = arg2(:,[1,3]);
      xen = arg1(:,[2,4]);
      yen = arg2(:,[2,4]);
     otherwise,
      error('If X,Y are given, each must be either an Nx2 or Nx4 matrix');
    end;







  ax(1)=subplot_tight(4,1,[1:3]);
  titlename(['Upwelling Mechanisms - ',dtstr]);
  title('Temperatures by depth [\circC]');
  plot_ts(GAPFUN(sfomc.ne_buoy.at_45m.sbe_seatemp),...
          GAPFUN(sfomc.ne_buoy.at_30m.sbe_seatemp),...
          GAPFUN(sfomc.nw_w_btm.adcp_seatemp),...
          ts_nanify_gaps(fwyf1.ndbc_air_t,3/24),'k');
  xlim_datetick(xl);
  ylim([15,31]); grid on; grid minor; hold on;
  %plot(sfomc.ne_buoy.Nx.date(upwix),repmat(23,size(upwix)),'r^');
  %plot(sfomc.ne_buoy.Nx.date(dwnix),repmat(23,size(dwnix)),'bv');
  %ylabel('Temperature [\circC]');
  annotate_ts_event(upwMr,sfomc.ne_buoy.at_45m.sbe_seatemp,upwmT,-inf);
  legend('45 m','30 m','11 m','Air','(Upw.)', 'Location','SouthWest');
  panel_label_subplot(gca,'a','ur','FontSize',fontsz,'BackgroundColor','w');
  set(gca,'FontSize',fontsz); set(gca,'LineWidth',gridwid); set(get(gca,'Child'),'LineWidth',linewid);





  if ( numel(args) > 0 && (ischar(args{1}) || isnumeric(args{1}) ) )
    locs = args{1};
    args(1) = [];
  else
    locs = 'ul';
  end;




% Copied from R2016a version of this private function
if nargin < 1
    error(message('MATLAB:narginchk:notEnoughInputs'));
elseif nargin > 1
    error(message('MATLAB:narginchk:tooManyInputs'));
end



  annot_dts = [datenum([1999,1999],[7,8],[27,19]) ; datenum([1999,2000],[10,2],[25,25])];
  annot_h = textbox_range([shax,dpax,ekax,prax],annot_dts,{'Upwelling','Cold Fronts'},'uc');





    panel_label_subplot(gca,'d','ur','FontSize',fontsz,'BackgroundColor','w');




  if ( isscalar(ax) && nannots > 1 )
    ax = repmat(gca,[annots,1]);
  end;




   case 3, % Nx3 coordinate matrices X and Y including (numerical) AXES handles
    xax = arg1(:,1);
    yax = arg2(:,1);





    if ( size(arg2,2) == 3 )
      yax = arg2(:,1);
    else
      yax = xax;
    end;





      if ( blowupmap )
        annotaxis(fwc.Alpha.lat,'FWC Alpha/Rodeo \rightarrow',annotargs{:});
      else
        annotaxis(fwc.Alpha.lat+0.050,'FWC Alpha/Rodeo \rightarrow',annotargs{:});
      end;





      annotaxis(27.170278,'FDEP "K": met. \rightarrow',annotargs{:});




  if 0;
    fmg; plot_ts(sfomc.ne_buoy.intNx,sfomc.sw_buoy.intNx,sfomc.nw_w_btm.intNx); legend('Outer','Mid','Inshore','Location','East');
    axis([datenum(1999,[7,8],[15,29]),0,8]); datetick('x',2,'keeplimits'); %datetick3;
    ylabel('\int^0_-_h (u_{onshore}\bullet[NO2+NO3])dz   [kg\bulletm^-^1]'); set(gca,'FontSize',16);
    titlename('Accumulated Nitrate+Nitrite per meter of reef');
    if doPrint; print('-dpng',fullfile(figspath,['upwelling-sfomc-1999-N_cumkg.png'])); end;

    fmg; plot_ts(sfomc.ne_buoy.intPx,sfomc.sw_buoy.intPx,sfomc.nw_w_btm.intPx); legend('Outer','Mid','Inshore','Location','East');
    axis([datenum(1999,[7,8],[15,29]),0,1]); datetick('x',2,'keeplimits'); %datetick3;
    ylabel('\int^0_-_h (u_{onshore}\bullet[P])dz   [kg\bulletm^-^1]'); set(gca,'FontSize',16);
    titlename('Accumulated Phospate per meter of reef');
    if doPrint; print('-dpng',fullfile(figspath,['upwelling-sfomc-1999-P_cumkg.png'])); end;
  end;



      %if ( yr == 1999 )
      %  axes(prax(end)); axes(prax(end)); print('-dpng',fullfile(figspath,['upwelling-sfomc-',yrstr,'-partn-wide.png']));
      % end;




  % %% 1999, 2000, 2003, 2004, 2005?, 2006, 2007, 2009?, 2010?, 2011
  % modys = datenum(0,[7,9],[5,5]);
  % for yr=[1999, 2000, 2003, 2004, 2005, 2006, 2007, 2011];
  modys = datenum(0,[6,9],[5,5]);
  for yr=[1999, 2000, 2003, 2004, 2005, 2006, 2007, 2009, 2010, 2011];
    yrstr = num2str(yr);
    disp(['Focus on high-frequency variability in 5-min sea temperature in *',yrstr,'*']);

    try, delete(annot_h); catch, end;

    axes(shax(end)); xlim_datetick(shax,datenum(yr,1,0)+modys); datetick(shax(3),'x',2,'keeplimits');
    if (yr>2000); legend(sh_lhs([1,end]),{'Air (FWYF1)','Shallow (NW 11 m)'}, legargs{:}); end;
    if doPrint; axes(shax(end)); print('-dpng',fullfile(figspath,['upwelling-sfomc-',yrstr,'-shallow-wide.png'])); end;

    if ( yr == 1999 || yr == 2000 )
      axes(dpax(end)); xlim_datetick(dpax,datenum(yr,1,0)+modys); datetick(dpax(3),'x',2,'keeplimits');
      if doPrint; print('-dpng',fullfile(figspath,['upwelling-sfomc-',yrstr,'-deep-wide.png'])); end;
    end;

    axes(ekax(end)); xlim_datetick(ekax,datenum(yr,1,0)+modys); datetick(ekax(3),'x',2,'keeplimits');
    if ( yr > 1999 )
      panel_label_subplot(ekax(1),'d','ur','FontSize',fontsz,'BackgroundColor','w');
      panel_label_subplot(ekax(2),'e','ur','FontSize',fontsz,'BackgroundColor','w');
      panel_label_subplot(ekax(3),'f','ur','FontSize',fontsz,'BackgroundColor','w');
    end;
    if doPrint; axes(ekax(end)); print('-dpng',fullfile(figspath,['upwelling-sfomc-',yrstr,'-ekman-wide.png'])); end;

    % axes(wsax(end)); xlim_datetick(wsax,datenum(yr,1,0)+modys); datetick(wsax(3),'x',2,'keeplimits');
    % if doPrint; axes(wsax(end)); print('-dpng',fullfile(figspath,['upwelling-sfomc-',yrstr,'-stress-wide.png'])); end;

    % axes(vmax(end)); xlim_datetick(vmax,datenum(yr,1,0)+modys); datetick(vmax(3),'x',2,'keeplimits');
    % if doPrint; axes(vmax(end)); print('-dpng',fullfile(figspath,['upwelling-sfomc-',yrstr,'-variab-wide.png'])); end;

    axes(prax(end)); xlim_datetick(prax,datenum(yr,1,0)+modys); datetick(prax(3),'x',2,'keeplimits');
    if ( yr > 1999 )
      panel_label_subplot(prax(1),'g','ur','FontSize',fontsz,'BackgroundColor','w');
      panel_label_subplot(prax(2),'h','ur','FontSize',fontsz,'BackgroundColor','w');
      panel_label_subplot(prax(3),'i','ur','FontSize',fontsz,'BackgroundColor','w');
    end;
    if doPrint; axes(prax(end)); print('-dpng',fullfile(figspath,['upwelling-sfomc-',yrstr,'-partn-wide.png'])); end;
    % if ( yr == 1999 )
    %   axes(prax(end)); xlim_datetick(prax,datenum(yr,1,0)+modys); datetick(prax(3),'x',2,'keeplimits');
    %   if doPrint; axes(prax(end)); axes(prax(end)); print('-dpng',fullfile(figspath,['upwelling-sfomc-',yrstr,'-partn-wide.png'])); end;
    % end;

    annot_dts = sfomc_upwelling_dates(datenum(yr,1,0)+modys);
    annot_h=textbox_range([shax,dpax,ekax,prax],annot_dts);

    if doPause; pause; else; pause(0.5); end;
  end;
  %xlim([shax,dpax,ekax,wsax],datenum(1999,1,0)+modys); jumpax([shax,dpax,ekax,wsax],nan,365);
  %keyboard;










  if ( exist('dtrng','var') )
    dts = dts(dtrng(1) <= dts & dts <= dtrng(2));
  end;






  shan=textbox_range(shax(1),andts,{'Upwelling','Cold Fronts'},'uc');
  textbox_range(shax,andts);
  dpan=textbox_range(dpax(1),andts,{'Upwelling','Cold Fronts'},'uc');
  textbox_range(dpax,andts);
  ekan=textbox_range(ekax(1),andts,{'Upwelling','Cold Fronts'},'uc');
  textbox_range(ekax,andts);
  pran=textbox_range(prax(1),andts,{'Upwelling','Cold Fronts'},'uc');
  textbox_range(prax,andts);





  shan=textbox_range(shax(1),andts(1,:),'Upwelling','uc');
  textbox_range(shax(1),andts(2,:),'Cold Fronts','uc');
  textbox_range(shax,andts);
  dpan=textbox_range(dpax,andts(1,:),'Upwelling','uc');
  textbox_range(dpax,andts(2,:),'Cold Fronts','uc');
  textbox_range(dpax,andts(2,:),'Cold Fronts','uc');
  ekan=textbox_range(ekax,andts(1,:),'Upwelling','uc');
  textbox_range(ekax,andts(2,:),'Cold Fronts','uc');
  pran=textbox_range(prax,andts(1,:),'Upwelling','uc');
  textbox_range(prax,andts(2,:),'Cold Fronts','uc');






  shan=textbox_range(shax,datenum([1999,1999],[7,8],[27,19]),'Upwelling','uc');
  textbox_range(shax,datenum([1999,2000],[10,2],[25,25]),'Cold Fronts','uc');
  dpan=textbox_range(dpax,datenum([1999,1999],[7,8],[27,19]),'Upwelling','uc');
  textbox_range(dpax,datenum([1999,2000],[10,2],[25,25]),'Cold Fronts','uc');
  ekan=textbox_range(ekax,datenum([1999,1999],[7,8],[27,19]),'Upwelling','uc');
  textbox_range(ekax,datenum([1999,2000],[10,2],[25,25]),'Cold Fronts','uc');
  pran=textbox_range(prax,datenum([1999,1999],[7,8],[27,19]),'Upwelling','uc');
  textbox_range(prax,datenum([1999,2000],[10,2],[25,25]),'Cold Fronts','uc');








  if ( YnotX )
    yl = xlim(ax(axix));
  else
    yl = ylim(ax(axix));
  end;
  
  for ix = 1:size(xcoords,1)
    if ( YnotX )
      coords = [yl(1),xcoords(ix,1),(yl(2)-yl(1)),xcoords(ix,2)-xcoords(ix,1)];
    else
      coords = [xcoords(ix,1),yl(1),xcoords(ix,2)-xcoords(ix,1),(yl(2)-yl(1))];
    end;
    
    rectpos = ds2nfu(coords);
    anh(ix)=annotation('textbox',rectpos,'Color','k','BackgroundColor',...
                       [.8,.8,.8],'FaceAlpha',0.3,'EdgeColor','none', ...
                       'String',lbls{ix},'FontSize',24,'FontName','Times', ...
                       alignargs{:},args{:});
  end;




  % annotation('rectangle',[0.1809 0.0387 0.0721 0.9168],'Color','r','LineWidth',2);
  % annotation('rectangle',[0.2552 0.0387 0.0626 0.9168],'Color','b','LineWidth',2);
  % annotation('rectangle',[0.3193 0.0387 0.0831 0.9168],'Color','m','LineWidth',2);





  grid on; center_axes;



ds = [...
    0,...
    nw_dx, repmat(nw_dx,[1,numel(sfomc.nw_w_btm.(fld).depths)]), nw_dx,...
    c_dx, repmat(c_dx,[1,numel(sfomc.sw_buoy.(fld).depths)]), c_dx,...
    ne_dx, repmat(ne_dx,[1,numel(sfomc.ne_buoy.(fld).depths(nezix))]), ne_dx,...
     ];
zs = [ ...
    0;...
    0; sfomc.nw_w_btm.(fld).depths'; -20;...
    0; sfomc.sw_buoy.(fld).depths'; -20;...
    0; sfomc.ne_buoy.(fld).depths(nezix)'; -50; ...
     ];
Ts = [ ...
    endpt,...
    endpt, sfomc.nw_w_btm.(fld).prof(nwdix,:), sfomc.sw_buoy.(fld).prof(swdix,end),...
    endpt, sfomc.sw_buoy.(fld).prof(swdix,:), sfomc.sw_buoy.(fld).prof(swdix,end),...
    endpt, sfomc.ne_buoy.(fld).prof(nedix,nezix), sfomc.ne_buoy.(fld).prof(nedix,end),...
     ];






  axis([datenum(1999,[7,8],[15,29]),0,0.8]); datetick('x',2,'keeplimits'); %datetick3;





1;
%% SCRIPT plot_upwelling_section_surf.m:
%

if ( ~exist('doPrint','var') )
  doPrint = false;
end;

%dt=datenum(1999,7,28,10,22,36);

if ( ~exist('fld','var') )
  %fld = 'seatemp';	fldstr = 'Sea temperature [^oC]';
  fld = 'N';		fldstr = '[NO3+NO2] [\mumol]';
end;

if ( strcmp(fld,'seatemp') )
  dt=datenum(1999,7,28,14,46,55);
  swdt = datenum(1999,7,27,16,37,04);
  [ig,fwydix]=min(abs(fwyf1.ndbc_air_t.date-dt));
  endpt = fwyf1.ndbc_air_t.data(fwydix);
elseif ( strcmp(fld,'N') )
  dt = datenum(1999,7,28,10,40,00);
  swdt = datenum(1999,7,27,17,50,00);
  endpt = 0;
end;
[ig,nedix]=min(abs(sfomc.ne_buoy.(fld).date-dt));
[ig,swdix]=min(abs(sfomc.sw_buoy.(fld).date-swdt));
[ig,nwdix]=min(abs(sfomc.nw_w_btm.(fld).date-dt));

%[nw_dx,nw_az] = station_dist(sfomc.nw_w_btm,sfomc.c_buoy); nw_dx = -nw_dx;
[nw_dx,nw_az,nw_nearestCoords] = get_contour_distance(sfomc.c_buoy.ngdc_hires_bathy,sfomc.nw_w_btm);

[c_dx,c_az,c_nearestCoords] = get_contour_distance(sfomc.c_buoy.ngdc_hires_bathy,sfomc.c_buoy);

%[ne_dx,ne_az] = station_dist(sfomc.c_buoy,sfomc.ne_buoy);
%ne_dx = 1;
[ne_dx,ne_az,ne_nearestCoords] = get_contour_distance(sfomc.c_buoy.ngdc_hires_bathy,sfomc.ne_buoy);
ne_dx = 3.20;

% ds = [...
%     -c_dx,...
%     repmat(nw_dx,[1,2]),...
%     repmat(0,[1,numel(sfomc.sw_buoy.(fld).depths)]),...
%     repmat(ne_dx,[1,numel(sfomc.ne_buoy.(fld).depths)+1]),...
%      ]; % + c_dx;
ds = [...
%    0,...
    repmat(nw_dx,[1,numel(sfomc.nw_w_btm.(fld).depths)]),...
    repmat(c_dx,[1,numel(sfomc.sw_buoy.(fld).depths)]),...
    repmat(ne_dx,[1,numel(sfomc.ne_buoy.(fld).depths)+1]),...
     ];
zs = [ ...
%    0;...
    sfomc.nw_w_btm.(fld).depths';...
    sfomc.sw_buoy.(fld).depths';...
    sfomc.ne_buoy.(fld).depths';...
    -50
     ];
Ts = [ ...
%    endpt,...
    sfomc.nw_w_btm.(fld).prof(nwdix,:),...
    sfomc.sw_buoy.(fld).prof(swdix,:),...
    sfomc.ne_buoy.(fld).prof(nedix,:),...
    sfomc.ne_buoy.(fld).prof(nedix,end),...
     ];

% xq = -c_dx:0.01:ne_dx;
% yq = -50:0.05:0;
xq = 0:0.01:ne_dx;
yq = -50.15:0.05:0;
[Xq,Yq] = meshgrid(xq,yq);

Vq = griddata(ds,zs,Ts,Xq,Yq);

brown = [0.4,0.0,0.0];

fmg;

%[res,ig,lh]=plot_bathy_transect(sfomc.sw_buoy,ne_dx,ne_az,sfomc.c_buoy.ngdc_hires_bathy);
[res,sfomc.c_buoy,lh]=plot_bathy_transect({sfomc.c_buoy,c_nearestCoords},ne_dx+0.100,c_az+180);
badix = find(isnan(res.depths));
res.lon(badix) = []; res.lat(badix) = []; res.field(badix) = []; res.range(badix) = []; res.depths(badix) = []; 
set(lh,'LineWidth',9,'Color',brown);
legend off
%xlim([-0.100,ne_dx+0.100]); ylim([-53,3]);
xlim([nw_dx-0.100,ne_dx+0.100]); ylim([-53,3]);
axis(axis);

BWx = [res.range,res.range(end),res.range(1),res.range(1)];
BWy = [res.depths,min(res.depths),min(res.depths),res.depths(1)];
%flh = fill(BWx,BWy,brown);
% BWx = size(Vq,1)*( (BWx-min(BWx))/(max(BWx)-min(BWx)) );
% BWy = size(Vq,2)*( (BWy-min(BWy))/(max(BWy)-min(BWy)) );
% BW = poly2mask(BWx,BWy,size(Vq,1),size(Vq,2));
% Vq(BW) = nan;

%BWix=[];
if ( ~exist('BWix','var') )
  %tic, BWix = find(inside(Xq,Yq,BWx,BWy)); toc,
  tic, BWix = find(inpolygon(Xq,Yq,BWx,BWy)); toc,
end;
Vq(BWix) = nan;
sh = surf(Xq,Yq,Vq);
ws = warning('OFF','MATLAB:handle_graphics:Patch:NumColorsNotEqualNumVertsException');
shading interp
warning(ws); clear ws
set(sh,'FaceAlpha',0.8);
flh = fill(BWx,BWy,brown);
%%xlim([-0.8,1.2]);
%xlim([-1.0,1.2]);
%xlim([-1.2,1.2]);
%ylim([-53,3]);

cbh=colorbar('Location','EastOutside');
%xlabel(cbh,textize(fld));
xlabel(cbh,fldstr);
plot3(nw_dx,-sfomc.nw_w_btm.depth,1,'rp','MarkerSize',24,'LineWidth',2);
%plot3(0,-sfomc.c_buoy.depth,1,'rp','MarkerSize',24,'LineWidth',2);
plot3(c_dx,-sfomc.c_buoy.depth,1,'rp','MarkerSize',24,'LineWidth',2);
plot3(ne_dx,-sfomc.ne_buoy.depth,1,'rp','MarkerSize',24,'LineWidth',2);
xlabel('Distance from shore [km]');

if ( doPrint )
  print('-dpng',fullfile(get_coral_path,'CRCP','Upwelling','CoRIS',[mfilename,'_',fld,'.png']));
end;











%[res,ig,lh]=plot_bathy_transect(sfomc.sw_buoy,ne_dx,ne_az,sfomc.c_buoy.ngdc_hires_bathy);
[res,sfomc.c_buoy,lh]=plot_bathy_transect({sfomc.c_buoy,c_nearestCoords},ne_dx+0.100,c_az+180);


%xlim([-0.100,ne_dx+0.100]); ylim([-53,3]);
xlim([nw_dx-0.100,ne_dx+0.100]); ylim([-53,3]);


BWix=[];
if ( ~exist('BWix','var') )
  %tic, BWix = find(inside(Xq,Yq,BWx,BWy)); toc,
  tic, BWix = find(inpolygon(Xq,Yq,BWx,BWy)); toc,
end;

%flh = fill(BWx,BWy,brown);
% BWx = size(Vq,1)*( (BWx-min(BWx))/(max(BWx)-min(BWx)) );
% BWy = size(Vq,2)*( (BWy-min(BWy))/(max(BWy)-min(BWy)) );
% BW = poly2mask(BWx,BWy,size(Vq,1),size(Vq,2));
% Vq(BW) = nan;



%set(get(cbh,'Label'),'String','[NO3+NO2]');




    0,...
    0;...
    endpt,...


%BW = poly2mask(BWx,BWy,size(Vq,1),size(Vq,2));



% res=plot_station_bathy_transect(sfomc.c_buoy,dx,az,'ngdc_hires_bathy')
% fmg; plot_bathy_transect({sfomc.ne_buoy,sfomc.c_buoy},[],sfomc.ne_buoy.isobath_orientation-90); ylim([-55,5]);
% help plot_bathy_transect
axis(axis);
%plot_bathy_transect({sfomc.ne_buoy,sfomc.c_buoy},[],sfomc.ne_buoy.isobath_orientation-90);
%% function [results,stn,lh] = plot_bathy_transect(stn,rad,azs,bathynm)




[ig,nedix]=min(abs(sfomc.ne_buoy.seatemp.date-dt));
[ig,swdix]=min(abs(sfomc.sw_buoy.seatemp.date-dt));
[ig,nwdix]=min(abs(sfomc.nw_w_btm.seatemp.date-dt));




dt=datenum(1999,7,28,10,22,36);
[ig,fwydix]=min(abs(fwyf1.ndbc_air_t.date-dt));
[ig,nedix]=min(abs(sfomc.ne_buoy.seatemp.date-dt));
[ig,nedix]=min(abs(sfomc.ne_buoy.seatemp.date-dt));
[ig,swdix]=min(abs(sfomc.sw_buoy.seatemp.date-dt));
[ig,swdix]=min(abs(sfomc.sw_buoy.seatemp.date-dt));
[ig,nwdix]=min(abs(sfomc.nw_w_btm.seatemp.date-dt));
[ig,nwdix]=min(abs(sfomc.nw_w_btm.seatemp.date-dt));




  for ix = 1:size(xcoords,1)
    coords = [xcoords(ix,1),yl(1),xcoords(ix,2)-xcoords(ix,1),yl(2)];
    rectpos = ds2nfu(coords);
    if ( ~isempty(lbls{ix}) )
      anh(ix)=annotation('textbox',rectpos,'Color','k','BackgroundColor',[.8, ...
                          .8,.8],'FaceAlpha',0.3,'EdgeColor','none','String', ...
                         '','FontSize',24,'FontName','Times'); % lbls{ix}
    else
      anh(ix)=annotation('rectangle',rectpos,'Color','none','FaceColor',[.8,.8,.8],'FaceAlpha',0.3,'Color','none');
    end;
  end;





keyboard;
      anh(ix)=annotation('textbox',rectpos,,'Color','none','FaceColor',[.8,.8,.8],'FaceAlpha',0.3,'Color','none','String',lbls{ix},'FontSize',24);




  annotargs = {'Color','none','FaceColor',[.8,.8,.8],'FaceAlpha',0.3,'Color','none'};



      anh(ix)=annotation('rectangle',rectpos,'Color','none','FaceColor',[.8,.8,.8],'FaceAlpha',0.3,'Color','none');


    %anh(ix)=annotation('rectangle',rectpos,'Color','m','LineWidth',0);




25 Gmol/year of [NO2+NO3], 1.5 Gmol/year of [PO4], and 3.6 Gmol/year of [Si]




  %tbuf = 0.02;
  tbuf = 0.04; % Leave room for a LARGE title




  %xl = datenum(1999,[7,8],[15,29]); %Broader upwelling period




    annotation(fh, 'textbox', txtpos, 'String', lbl, 'FitBoxToText','on',  'verticalalignment', 'bottom');






1;

hf = figure;
for ia = 1:4
    ha(ia) = subplot(2,2,ia);
end
pos = get(ha, 'position');
dim = cellfun(@(x) x.*[1 1 0.5 0.5], pos, 'uni',0);
annotation(hf, 'textbox', dim{1}, 'String', 'test 1');
annotation(hf, 'textbox', dim{2}, 'String', 'test 2', 'FitBoxToText','on');
annotation(hf, 'textbox', dim{3}, 'String', 'test 3', 'verticalalignment', 'bottom');
annotation(hf, 'textbox', dim{4}, 'String', 'test 4', 'FitBoxToText','on',  'verticalalignment', 'bottom');






  shax(2) = subplot_tight(7,1,4:5); ylabel('Alongshore'); %title('Shallow alongshore current');







  % QC
  % Loop over each timestamp, to find the bin number just below the surface spike
  stn.adcp_sfc_height.date = dts;
  stn.adcp_sfc_height.data = repmat(max(hgts),[numel(dts),1]);

  last_good_bin = repmat(numel(hgts),[numel(dts),1]);
  for ix = 1:numel(dts)
    uprof = stn.(ufld).prof(ix,:);
    vprof = stn.(vfld).prof(ix,:);
    if ( any(~isnan(uprof)) && any(~isnan(vprof)) )
      % uspikeval = prctile(abs(uprof),99);
      % vspikeval = prctile(abs(vprof),99);
      prof = uv_to_spd(uprof,vprof);
      spikeval = prctile(prof,99);

      % NOTE: Using >=99th PRCTILE means we will ALWAYS find a spike bin.
      % (PRCTILE always returns a value within the range of its argument.)
      %spikeix = find( (abs(uprof)>=uspikeval | abs(vprof)>=vspikeval), 1 );
      spikeix = find( (prof>=spikeval), 1 );
      if ( ~isempty(spikeix) )
        stn.adcp_sfc_height.data(ix) = hgts(spikeix);

        % Beam angles 20o: side lobe contamination in last 6% + 1 bin of profile
        last_good_bin(ix) = spikeix - (floor(side_lobe_contam_frac*spikeix) + 1);
        %DEBUG:
if (last_good_bin(ix)<11 || last_good_bin(ix)==14); keyboard; end;
        stn.(ufld).prof(ix,last_good_bin(ix)+1:end) = nan;
        stn.(vfld).prof(ix,last_good_bin(ix)+1:end) = nan;
      end;
    end;
    stn.(ufld).data(ix,1) = nanmean(uprof);
    stn.(vfld).data(ix,1) = nanmean(vprof);
  end;












  [face.HanS_HW.adcp_speed,face.HanS_HW.adcp_dir] = ts_uv_to_spddir(face.HanS_HW.adcp_u,face.HanS_HW.adcp_v,false);
  [face.HanS_BR.adcp_speed,face.HanS_BR.adcp_dir] = ts_uv_to_spddir(face.HanS_BR.adcp_u,face.HanS_BR.adcp_v,false);
  [face.AOML_HW.adcp_speed,face.AOML_HW.adcp_dir] = ts_uv_to_spddir(face.AOML_HW.adcp_u,face.AOML_HW.adcp_v,false);
  [face.AOML_BR.adcp_speed,face.AOML_BR.adcp_dir] = ts_uv_to_spddir(face.AOML_BR.adcp_u,face.AOML_BR.adcp_v,false);







        % Beam angles 20o: side lobe contamination in last 6% + 1 bin of profile
        last_good_bin(ix) = spikeix - (ceil(side_lobe_contam_frac*spikeix) + 1);




  % QC
  % Loop over each timestamp, to find the bin number just below the surface spike
  stn.adcp_sfc_height.date = dts;
  stn.adcp_sfc_height.data = repmat(max(hgts),[numel(dts),1]);

  last_good_bin = repmat(numel(hgts),[numel(dts),1]);
  for ix = 1:numel(dts)
    uprof = stn.(ufld).prof(ix,:);
    vprof = stn.(vfld).prof(ix,:);
    if ( any(~isnan(uprof)) && any(~isnan(vprof)) )
      uspikeval = prctile(abs(uprof),99);
      vspikeval = prctile(abs(vprof),99);

      % NOTE: Using >=99th PRCTILE means we will ALWAYS find a spike bin.
      % (PRCTILE always returns a value within the range of its argument.)
      spikeix = find( (abs(uprof)>=uspikeval | abs(vprof)>=vspikeval), 1 );
      if ( ~isempty(spikeix) )
        stn.adcp_sfc_height.data(ix) = hgts(spikeix);

        % Beam angles 20o: side lobe contamination in last 6% + 1 bin of profile
        last_good_bin(ix) = spikeix - (ceil(side_lobe_contam_frac*spikeix) + 1);
        stn.(ufld).prof(ix,last_good_bin(ix)+1:end) = nan;
        stn.(vfld).prof(ix,last_good_bin(ix)+1:end) = nan;
      end;
    end;
    stn.(ufld).data(ix,1) = nanmean(uprof);
    stn.(vfld).data(ix,1) = nanmean(vprof);
  end;









  hgts = -deps;






    face.HanS_HW.adcp_v.depths = face.HanS_HW.depth - (face.HanS_HW.adcp_blank_dist + (face.HanS_HW.adcp_bin_size*[0:14]));


    face.AOML_BR.adcp_u.depths = face.AOML_BR.adcp_bin_depths(face.AOML_BR.valid_bins);
    face.AOML_BR.adcp_u.depths = face.AOML_BR.depth - (face.AOML_BR.adcp_blank_dist + (face.AOML_BR.adcp_bin_size*[0:14]));







    %DEBUG:
    disp('Currents');
    face.AOML_HW.adcp_bin_size = 0.5;
    face.AOML_HW.adcp_blank_dist = 0.5;

    face.HanS_HW.adcp_bin_size = 2.0;
    face.HanS_HW.adcp_blank_dist = 1.0;

    face.AOML_BR.adcp_bin_size = 0.5;
    face.AOML_BR.adcp_blank_dist = 0.5;

    face.HanS_BR.adcp_bin_size = 2.0;
    face.HanS_BR.adcp_blank_dist = 1.0;
    
    % FROM "Meta data for FACE AOML and Hazen and Sawyer current  data set  5 20 2011 .docx":
    face.AOML_BR.adcp_depths = [6.35,5.85,5.35,4.85,4.35,3.85,3.35,2.85,2.35,1.85,1.35,0.85,0.35,];
-0.15
-0.65







    plotargs = {'rp','MarkerFaceColor','w','MarkerSize',10};

    %annotargs = {'FontSize',fontsz,'BackgroundColor','w'};
    annotargs = {'FontSize',fontsz,'BackgroundColor','none'};

    % Fix the AXES bounds, so stations outside the mapped area are ignored
    axis(axis);

    if ( exist('pvgf1','var') )
      plot(pvgf1.lon,pvgf1.lat,plotargs{:});
      annotaxis(pvgf1.lat,'Port Everglades: met. \rightarrow',annotargs{:});
    end;
    if ( exist('ftpf1','var') )
      plot(ftpf1.lon,ftpf1.lat,plotargs{:});
      annotaxis(ftpf1.lat,'41114: Hs \rightarrow',annotargs{:});
    end;

    if ( exist('face','var') )
      plot(face.shallow.lon,face.shallow.lat,plotargs{:});
      plot(face.tcm1.lon,face.tcm1.lat,plotargs{:});
      plot(face.tcm2.lon,face.tcm2.lat,plotargs{:});
      plot(face.deep.lon,face.deep.lat,plotargs{:});
      annotaxis(face.deep.lat,'FACE HWD moorings: u,v \rightarrow',annotargs{:});

      if ( isfield(face,'AOML_HW') )
        % plot(face.AOML_HW.lon,face.AOML_HW.lat,plotargs{:});
        % plot(face.HanS_HW.lon,face.HanS_HW.lat,plotargs{:});
        % ah=annotaxis(face.AOML_HW.lat,'FACE/H&S HW: u,v \rightarrow',annotargs{:});
        plot(face.AOML_BR.lon,face.AOML_BR.lat,plotargs{:});
        plot(face.HanS_BR.lon,face.HanS_BR.lat,plotargs{:});
        annotaxis(face.AOML_BR.lat,'FACE/H&S BR: u,v \rightarrow',annotargs{:});
      end;
    end;

    if ( exist('sfomc','var') )
      plot(sfomc.ne_buoy.lon,sfomc.ne_buoy.lat,plotargs{:});
      plot(sfomc.c_buoy.lon,sfomc.c_buoy.lat,plotargs{:});
      plot(sfomc.sw_buoy.lon,sfomc.sw_buoy.lat,plotargs{:});
      plot(sfomc.pier_cc.lon,sfomc.pier_cc.lat,plotargs{:});
      plot(sfomc.nw_w_btm.lon,sfomc.nw_w_btm.lat,plotargs{:});
        annotaxis(sfomc.nw_w_btm.lat,'NSU NW/W,C,NE/E: Hs,u,v \rightarrow',annotargs{:});
        annotaxis(sfomc.sw_buoy.lat,'NSU SW buoy: u,v \rightarrow',annotargs{:});
    end;








      if ( isfield(face,'brwd20') )
        plot(face.brwd20.lon,face.brwd20.lat,'w.');
        plot(face.brwd40.lon,face.brwd40.lat,'w.');
        plot(face.brwd100.lon,face.brwd100.lat,'w.');
        plot(face.boca20.lon,face.boca20.lat,'w.');
        plot(face.boca40.lon,face.boca40.lat,'w.');
      end;






plot_ts(sfomc.ne_buoy.adcp_x_45m,sfomc.ne_buoy.adcp_x_30m,sfomc.nw_w_btm.adcp_x_btm,sfomc.ne_buoy.adcp_x_11m);

sfomc.nw_w_btm.adcp_x_btm,


upwmT = 25.0; % Near-bottom sea temperature below which upwelling may be occurring




if ( doFACE )
  if 0;
    fmg;
    subplot_tight(2,1,1);
     [shax,lh(1),lh(2)] = plotyy_ts(face.shallow.seatemp,face.shallow.adcp_x);
     set(shax,'FontSize',fontsz);
     legend(lh(1:2),'7m Ts','7m u','Location','SouthWest');
     grid on; grid minor;
     titlename('FACE HWD moorings - upwelling temperature and cross-shore currents');
     set(gca,'FontSize',fontsz); set(get(gca,'Child'),'LineWidth',linewid);
    subplot_tight(2,1,2);
     [dpax,lh(1),lh(2)] = plotyy_ts(face.deep.seatemp,face.deep.adcp_x);
     set(dpax,'FontSize',fontsz);
     legend(lh(1:2),'26m Ts','26m u','Location','SouthWest');
     grid on; grid minor;
     set(gca,'FontSize',fontsz); set(get(gca,'Child'),'LineWidth',linewid);
    if doPrint; print('-dpng',fullfile(figspath,'upwelling-face-MULTI.png')); end;
  end;

  if 0;
    fmg;
    [ax,lh(1),lh(2)] = plotyy_ts(face.deep.adcp_x,face.deep.seatemp);
    lh(3) = plot(ax(2),face.shallow.seatemp.date,face.shallow.seatemp.data,'-','Color',[0,0.5,0]);
    set(ax,'FontSize',fontsz);
    xlim(datenum(2015,[1,9],[1,30])); datetick3;
    ylim(ax(1),[-0.35,+0.35]);
    ylim(ax(2),[20,33]);
    legend(lh,'26m u','26m T','7m T','Location','SouthWest', 'Orientation','horiz');
    grid on; grid minor;
    set(gca,'FontSize',fontsz); set(get(gca,'Child'),'LineWidth',linewid);
    titlename('FACE HWD moorings - upwelling temperature/cross-shore currents');
    if doPrint; print('-dpng',fullfile(figspath,'upwelling-face-MULTI-2015.png')); end;
  end;


  % FACE temperature/currents (2014)

  if 1;
    %xl = datenum(2014,[7,10],[30,1]);
    xl = datenum(2014,[7,10],[18,1]);

    fmg;
    subplot_tight(7,1,1:5);
    set(gca,'FontSize',fontsz);
    % Match Deep and Shallow colors with the 2015 plot (below)
    lh=plot_ts(pvgf1.ndbc_air_t,'k-',face.deep.seatemp,'specix',1,face.shallow.seatemp,'specix',4);
    xlim_datetick(xl);
    %set(get(gca,'XAxis'),'Vis','off');
    ylim([20.0,33.0]);
    annotate_ts_event(upwMr,face.deep.seatemp,upwmT,20.0);
    %legend(lh,'PVGF1','Deep','Shallow','Location','SouthEast', 'Orientation','horiz');
    legend('Air (PVGF1)','Deep','Shallow','(Upwelling)',...
           'Location','SouthEast', 'Orientation','horiz');
    grid on; grid minor;
    set(gca,'FontSize',fontsz); set(get(gca,'Child'),'LineWidth',linewid);
    titlename('FACE HWD moorings - air/sea temperature [^oC] (2014)');

    subplot_tight(7,1,6:7);
    set(gca,'FontSize',fontsz);
    lh=plot_ts(face.deep.adcp_btm_l,face.deep.adcp_btm_x,face.shallow.adcp_btm_l,face.shallow.adcp_btm_x);
    title('Currents');
    %xlim(datenum(2014,[7,10],[30,1])); datetick3;
    xlim_datetick('Y',xl);
    ylim([-0.90,+0.90]);
    annotate_ts_event(upwMr,face.deep.seatemp,upwmT,-0.85);
    legend('Deep ALONG','Deep cross','Shallow ALONG','Shallow cross','(Upwelling)',...
           'Location','SouthWest', 'Orientation','horiz');
    grid on; grid minor;
    title('Near-bottom currents [ms^-^1]');
    set(gca,'FontSize',fontsz); set(get(gca,'Child'),'LineWidth',linewid);

    if doPrint; print('-dpng',fullfile(figspath,'upwelling-face-2014.png')); end;
  end;


  %%%% ??? COULD TCMs HAVE MIS-ORIENTED COMPONENTS??
  face.tcm1 = station_reorient_vectors(face.tcm1,-90,'x','l','alt_x','alt_l');
  face.tcm1.alt_x.data = -face.tcm1.alt_x.data;
  face.tcm2 = station_reorient_vectors(face.tcm2,-90,'x','l','alt_x','alt_l');
  face.tcm2.alt_x.data = -face.tcm2.alt_x.data;

  % FACE temperature/currents (2015)

  if 1;
    %xl = datenum(2015,[7,10],[30,1]);
    xl = datenum(2015,[7,10],[18,1]);
    fmg;
    subplot_tight(7,1,1:5);
    lh=plot_ts(pvgf1.ndbc_air_t,'k-',face.deep.seatemp,'specix',1,face.tcm2.seatemp,'specix',2,face.tcm1.seatemp,'specix',3,face.shallow.seatemp,'specix',4);
    xlim_datetick(xl); %datetick3;
    %set(get(gca,'XAxis'),'Vis','off');
    ylim([20.0,33.0]);
    annotate_ts_event(upwMr,face.deep.seatemp,upwmT,20.0);
    legend('Air (PVGF1)','Deep','TCM2','TCM1','Shallow','(Upwelling)',...
           'Location','SouthEast', 'Orientation','horiz');
    grid on; grid minor;
    set(gca,'FontSize',fontsz); set(get(gca,'Child'),'LineWidth',linewid);
    titlename('FACE HWD moorings - air/sea temperature [^oC] (2015)');

    subplot_tight(7,1,6:7);
    lh=plot_ts(face.deep.adcp_btm_l,face.deep.adcp_btm_x,face.tcm2.alt_x,face.tcm1.alt_x,face.shallow.adcp_btm_x);
    title('Currents');
    xlim_datetick('Y',xl); %datetick3;
    %ylim([-1.10,+1.10]);
    ylim([-0.90,+0.90]);
    annotate_ts_event(upwMr,face.deep.seatemp,upwmT,-0.85);
    legend('Deep ALONG','Deep cross','TCM2 cross','TCM1 cross','Shallow cross','(Upwelling)',...
           'Location','SouthEast', 'Orientation','horiz');
    grid on; grid minor;
    title('Near-bottom currents [ms^-^1]');
    set(gca,'FontSize',fontsz); set(get(gca,'Child'),'LineWidth',linewid);

    if doPrint; print('-dpng',fullfile(figspath,'upwelling-face-2015.png')); end;
  end;


  % FACE temperature/currents (Sep 2015)

  if 1;
    xl = datenum(2015,[8,10],[31,1]);

    % Ocean deep response
    fmg;
    subplot_tight(7,1,1:3);
    lh=plot_ts(pvgf1.ndbc_air_t,'k-',face.deep.seatemp,face.tcm2.seatemp,face.tcm1.seatemp,face.shallow.seatemp);
    legend(lh,'Air (PVGF1)','Deep','TCM2','TCM1','Shallow','Location','SouthEast', 'Orientation','horiz');
    %ylim([22.5,31]);
    ylim([20,33]);
    ylabel('T [^oC]');
    xlim_datetick(xl); %datetick3;
    %set(get(gca,'XAxis'),'Vis','off');
    grid on; grid minor;
    set(gca,'FontSize',fontsz); set(get(gca,'Child'),'LineWidth',linewid);
    titlename('FACE HWD moorings - temperature/lower currents (Sep 2015)');

    subplot_tight(7,1,4:5);
    lh=plot_ts(face.deep.adcp_btm_l,face.tcm2.alt_l,face.tcm1.alt_l,face.shallow.adcp_btm_l);
    legend(lh,'Deep','TCM2','TCM1','Shallow','Location','SouthEast', 'Orientation','horiz');
    %ylim([-1.00,+1.00]);
    ylim([-1.25,+1.25]);
    %ylabel('u_l [ms^-^1]');
    xlim_datetick(xl); %datetick3;
    %set(get(gca,'XAxis'),'Vis','off');
    grid on; grid minor;
    title('Near-bottom along-shore [ms^-^1]');
    set(gca,'FontSize',fontsz); set(get(gca,'Child'),'LineWidth',linewid);

    subplot_tight(7,1,6:7);
    lh=plot_ts(face.deep.adcp_btm_x,face.tcm2.alt_x,face.tcm1.alt_x,face.shallow.adcp_btm_x);
    legend(lh,'Deep','TCM2','TCM1','Shallow','Location','SouthEast', 'Orientation','horiz');
    ylim([-0.25,+0.25]);
    %ylabel('u_x [ms^-^1]');
    xlim_datetick('Y',xl); %datetick3;
    grid on; grid minor;
    title('Near-bottom cross-shore [ms^-^1]');
    set(gca,'FontSize',fontsz); set(get(gca,'Child'),'LineWidth',linewid);

    if doPrint; print('-dpng',fullfile(figspath,'upwelling-face-2015-Sep.png')); end;


    % Ocean near-surface response

    fmg;
    subplot_tight(7,1,1:3);
    lh=plot_ts(pvgf1.ndbc_air_t,'k-',face.deep.seatemp,face.tcm2.seatemp,face.tcm1.seatemp,face.shallow.seatemp);
    legend(lh,'Air (PVGF1)','Deep','TCM2','TCM1','Shallow','Location','SouthEast', 'Orientation','horiz');
    %ylim([22.5,31]);
    ylim([20,33]);
    %ylabel('T [^oC]');
    xlim_datetick(xl); %datetick3;
    %set(get(gca,'XAxis'),'Vis','off');
    grid on; grid minor;
    set(gca,'FontSize',fontsz); set(get(gca,'Child'),'LineWidth',linewid);
    %titlename('FACE HWD moorings - temperature/upper currents (Sep 2015)');
    titlename('FACE HWD moorings - air/sea temperature [^oC] (Sep 2015)');

    subplot_tight(7,1,4:5);
    lh=plot_ts(face.deep.adcp_sfc_l,face.tcm2.alt_l,face.tcm1.alt_l,face.shallow.adcp_sfc_l);
    legend(lh,'Deep','TCM2','TCM1','Shallow','Location','SouthEast', 'Orientation','horiz');
    %ylim([-1.00,+1.00]);
    ylim([-1.25,+1.25]);
    %ylabel('u_l [ms^-^1]');
    xlim_datetick(xl); %datetick3;
    %set(get(gca,'XAxis'),'Vis','off');
    grid on; grid minor;
    title('Near-surface along-shore [ms^-^1]');
    set(gca,'FontSize',fontsz); set(get(gca,'Child'),'LineWidth',linewid);

    subplot_tight(7,1,6:7);
    lh=plot_ts(face.deep.adcp_sfc_x,face.tcm2.alt_x,face.tcm1.alt_x,face.shallow.adcp_sfc_x);
    legend(lh,'Deep','TCM2','TCM1','Shallow','Location','SouthEast', 'Orientation','horiz');
    %ylim([-0.25,+0.25]);
    ylim([-0.42,+0.42]);
    %ylabel('u_x [ms^-^1]');
    xlim_datetick('Y',xl); %datetick3;
    grid on; grid minor;
    title('Near-surface cross-shore [ms^-^1]');
    set(gca,'FontSize',fontsz); set(get(gca,'Child'),'LineWidth',linewid);

    if doPrint; print('-dpng',fullfile(figspath,'upwelling-face-2015-Sep-upper.png')); end;


    if 0;
      % Forcing

      wsfh = fmg;
      wsax(1) = subplot_tight(7,1,1:3); ylabel('Alongshore');
      lh=plot_ts(ts_nanify_gaps(canf1.ndbc_bulk_windstress_ls,maxgp),ts_nanify_gaps(lkwf1.ndbc_bulk_windstress_ls,maxgp),ts_nanify_gaps(pvgf1.ndbc_bulk_windstress_ls,maxgp));
      xlim_datetick(xl); %datetick3;
                         %set(get(gca,'XAxis'),'Vis','off');
      ylim([-0.35,+0.35]);
      grid on; grid minor; center_axes;
      legh = ...
          legend(lh,'Cape Canaveral (1988-2016)','Lake Worth (2001-2016)',...
                 'Port Everglades (2009-2017)',...
                 'Orientation','horiz','Location','southeast');
      set(legh,'FontSize',fontsz);
      set(gca,'FontSize',fontsz); set(get(gca,'Child'),'LineWidth',linewid);
      titlename('Wind stress [Nm^-^2] near FACE section (Sep 2015)');

      wsax(2) = subplot_tight(7,1,4:6); ylabel('Cross-shore');
      lh=plot_ts(ts_nanify_gaps(canf1.ndbc_bulk_windstress_xs,maxgp),ts_nanify_gaps(lkwf1.ndbc_bulk_windstress_xs,maxgp),ts_nanify_gaps(pvgf1.ndbc_bulk_windstress_xs,maxgp));
      xlim_datetick(xl); %datetick3;
      %set(get(gca,'XAxis'),'Vis','off');
      ylim([-0.35,+0.35]);
      grid on; grid minor; center_axes;
      set(gca,'FontSize',fontsz); set(get(gca,'Child'),'LineWidth',linewid);

      wsax(3) = subplot_tight(7,1,7); ylabel('Barometer');
      lh=plot_ts(ts_nanify_gaps(canf1.ndbc_barom,maxgp),ts_nanify_gaps(lkwf1.ndbc_barom,maxgp),ts_nanify_gaps(pvgf1.ndbc_barom,maxgp),ts_nanify_gaps(rsm.rsmas_barom,maxgp),ts_nanify_gaps(fwyf1.ndbc_barom,maxgp));
      xlim(xl); datetick3;
      ylim([1005,1030]);
      grid on; grid minor;
      datetick3;
      set(gca,'FontSize',fontsz); set(get(gca,'Child'),'LineWidth',linewid);

      if doPrint; print('-dpng',fullfile(figspath,'upwelling-face-2015-Sep-stress.png')); end;
    end;


    % Ekman response

    ekfh = fmg;
    ekax(1) = subplot_tight(7,1,1:3); ylabel('Cross-shore');
    titlename('Ekman volumetric flux [m^2s^-^1] near FACE section (Sep 2015)');
    lh=plot_ts(ts_nanify_gaps(canf1.ndbc_ekman_flux_x_volume,maxgp),ts_nanify_gaps(lkwf1.ndbc_ekman_flux_x_volume,maxgp),ts_nanify_gaps(pvgf1.ndbc_ekman_flux_x_volume,maxgp),ts_nanify_gaps(rsm.rsmas_ekman_flux_x_volume,maxgp),ts_nanify_gaps(fwyf1.ndbc_ekman_flux_x_volume,maxgp));
    xlim_datetick(xl); %datetick3;
    ylim([-4.9,+4.9]);
    center_axes; grid on; 
    legend(lh,'Cape Can. (''88-''16)','Lake Worth (''01-''16)','Port Ev. (''09-''17)','RSMAS Roof (''03-''11)','Fowey (''91-''16)', 'Orientation','horiz','Location','SouthEast');
    set(gca,'FontSize',fontsz); set(get(gca,'Child'),'LineWidth',linewid);

    ekax(2) = subplot_tight(7,1,4:6); ylabel('Alongshore');
    lh=plot_ts(ts_nanify_gaps(canf1.ndbc_ekman_flux_l_volume,maxgp),ts_nanify_gaps(lkwf1.ndbc_ekman_flux_l_volume,maxgp),ts_nanify_gaps(pvgf1.ndbc_ekman_flux_l_volume,maxgp),ts_nanify_gaps(rsm.rsmas_ekman_flux_l_volume,maxgp),ts_nanify_gaps(fwyf1.ndbc_ekman_flux_l_volume,maxgp));
    xlim_datetick(xl); %datetick3;
    ylim([-4.9,+4.9]);
    center_axes; grid on; 
    set(gca,'FontSize',fontsz); set(get(gca,'Child'),'LineWidth',linewid);

    ekax(3) = subplot_tight(7,1,7); ylabel('Barom.');
    lh=plot_ts(ts_nanify_gaps(canf1.ndbc_barom,maxgp),ts_nanify_gaps(lkwf1.ndbc_barom,maxgp),ts_nanify_gaps(pvgf1.ndbc_barom,maxgp),ts_nanify_gaps(rsm.rsmas_barom,maxgp),ts_nanify_gaps(fwyf1.ndbc_barom,maxgp));
    %xlim(xl); %datetick3;
    %ylim([980,1040]);
    ylim([1006,1023])
    grid on;
    xlim_datetick('Y',xl); %datetick('x',2,'keeplimits'); %datetick3;
    set(gca,'FontSize',fontsz); set(get(gca,'Child'),'LineWidth',linewid);
    %set(ekax(1:2),'XTickLabel','');

    if doPrint; print('-dpng',fullfile(figspath,'upwelling-face-2015-Sep-ekman.png')); end;


    % % What was our fundamental sampling period?
    % fmg; plot(sfomc.nw_w_btm.adcp_x_btm.date(1:end-1),24*diff(sfomc.nw_w_btm.adcp_x_btm.date),'.'); datetick3; ylim([0,8]);

    fmg;
    subplot_tight(7,1,1:3);
     lh=plot_ts(face.deep.seatemp_3_h_hp_9_h_lp,face.deep.seatemp_11_h_hp_14_h_lp,face.deep.seatemp_17_h_hp_28_h_lp,face.deep.seatemp_48_h_hp_240_h_lp);
     titlename('Sea temperature variability: FACE deep mooring (Sep 2015)');
     ylabel('^oC');
     xlim_datetick(xl); %datetick3;
     center_axes; grid on; 
     legend(lh,'3-9h','11-14h','17-28h','2-10d', 'Orientation','horiz','Location','SouthEast');
     set(gca,'FontSize',fontsz); set(get(gca,'Child'),'LineWidth',linewid);
    subplot_tight(7,1,4:5);
     lh=plot_ts(face.deep.adcp_btm_x_3_h_hp_9_h_lp,face.deep.adcp_btm_x_11_h_hp_14_h_lp,face.deep.adcp_btm_x_17_h_hp_28_h_lp,face.deep.adcp_btm_x_48_h_hp_240_h_lp);
     title('Near-bottom cross-shore current var.');
     ylabel('ms^-^1');
     xlim_datetick(xl); %datetick3;
     center_axes; grid on; 
     %legend(lh,'3-9h','11-14h','17-28h','2-10d', 'Orientation','horiz');
     set(gca,'FontSize',fontsz); set(get(gca,'Child'),'LineWidth',linewid);
    subplot_tight(7,1,6:7);
     lh=plot_ts(face.deep.adcp_btm_l_3_h_hp_9_h_lp,face.deep.adcp_btm_l_11_h_hp_14_h_lp,face.deep.adcp_btm_l_17_h_hp_28_h_lp,face.deep.adcp_btm_l_48_h_hp_240_h_lp);
     title('Near-bottom along-shore current var.');
     ylabel('ms^-^1');
     center_axes; grid on; 
     %legend(lh,'3-9h','11-14h','17-28h','2-10d', 'Orientation','horiz');
     set(gca,'FontSize',fontsz); set(get(gca,'Child'),'LineWidth',linewid);
     xlim_datetick('Y',xl); %datetick('x',2,'keeplimits'); %datetick3;
    %xlim(datenum([2013,2015],[11,11],[12,12])); datetick3;

    if doPrint; print('-dpng',fullfile(figspath,'upwelling-face-2015-Sep-partn.png')); end;

  end;

end; %if ( doFACE )














    % Match Shallow color with the 2015 plot (below)
    lh=plot_ts(pvgf1.ndbc_air_t,'k-',face.deep.seatemp,'specix',1,face.shallow.seatemp,'Color',[0.5,0.7,0.2]);




wsmag = 1.50; % Magnitude of variability in wind stress
ekmag = 5.00; % Magnitude of variability in Ekman volumetric flux




  fhs=[];
  [lkwf1,ig,ig,ig,fhs(1)] = plot_hires_bathymetry(lkwf1,-[0:2:20,30:20:200]);
  set(gca,'FontSize',fontsz);
  %titlename('SE Florida shelf - USGS Bathymetry - instrumentation');

  % Draw red rectangle around close-up area (see last map below)
  if ( exist('sfomc','var') )
    rectangle('Position',bbox2rect(field_bbox(sfomc.c_buoy.ngdc_hires_bathy)),...
              'EdgeColor','r','LineWidth',2);
  end;


  plot_all_se_florida_upwelling_sites;
  set(gca,'FontSize',fontsz);
  title([]);
  print('-dpng',fullfile(figspath,['upwelling-bathymetry.png']));






      if ( blowupmap )
        annotaxis(sefcri.dc1.lat,'SEFCRI DC1 \rightarrow',annotargs{:});
        annotaxis(sefcri.dc3.lat,'SEFCRI DC3 \rightarrow',annotargs{:});
      else
        annotaxis(sefcri.dc3.lat,'SEFCRI DC1/3 \rightarrow',annotargs{:});
      end;




          'SPGF1', ...		% Settlement Point, GBI
          -78.994, ...		%SPGF1 - Settlement Point, GBI
          26.704, ...		%SPGF1 - Settlement Point, GBI
          0, ...		%SPGF1 - Settlement Point, GBI





      annotaxis(fwc.Alpha.lat,'FWC Alpha/Rodeo \rightarrow','FontSize',fontsz);









  if ( ~exist('ixes','var') || isempty(ixes) )
    ixes = [];
  elseif ( numel(ixes)== 1 )
    ixes = repmat(ixes,[1,numel(vars)]);
  elseif ( numel(ixes) ~= numel(vars) )
    error('If IXES is given, NUMEL(IXES) must equal NUMEL(VARS)');
  end;






  %lms = [datenum(2009,10,30,13,0,0),datenum(2009,10,30,23,0,0),24.7,29.0];
  lms = [datenum(2009,10,30,13,0,0),datenum(2009,10,30,23,0,0),24.0,29.0];

  % fmg; lh=plot_ts(face.brwd20.brwd20_top.seatemp,face.brwd20.brwd20_btm.seatemp,face.brwd40.brwd40_top.seatemp,face.brwd40.brwd40_btm.seatemp,face.brwd100.brwd100_upp.seatemp); legend(lh,'20@sfc','20@btm','40@sfc','40@btm','100@40m'); axis(lms); datetick3('x','mmm-dd HH:MM','keeplimits');
  % set(gca,'FontSize',fontsz); set(get(gca,'Child'),'LineWidth',linewid);
  % titlename('FACE 2009: Upwelling at Broward');
  % if doPrint; print('-dpng',fullfile(figspath,'upwelling-nf09-brwd.png')); end;





    %lh=plot_ts(face.deep.adcp_l,face.deep.adcp_x,face.tcm2.x,face.tcm1.x,face.shallow.adcp_x);


    %%%% ??? COULD TCMs HAVE MIS-ORIENTED COMPONENTS??
    lh=plot_ts(face.deep.adcp_x,face.tcm2.l,face.tcm1.l,face.shallow.adcp_x);






  legargs = { 'Orientation','horiz','Location','SouthEast' };


    % legend(lh,'PVGF1','Deep','TCM2','TCM1','Shallow',...
    %        'Location','SouthEast', 'Orientation','horiz');





    % % Hovmoeller of buoyancy frequency squared
    % fmg; surf(N2.date,N2.depths',N2.prof'); colorbar; shading interp; datetick3;
    % colormap(hot); axis([tmp.date(1),datenum(1999,11,1), -43,0, 0,1, 0,0.0035]);

    % % % Hovmoeller of buoyancy frequency squared as an inverse fraction of Coriolis frequency squared
    % % fmg; surf(t,-Pv(:,1),real(log( (sw_f(sfomc.(stnm).lat).^2)./N2 ))); colorbar; shading interp; datetick3;
    % % colormap(hot); axis([tmp.date(1),datenum(1999,11,1), -43,0, 0,1, log([1e-3,0.15])]);






      %fmg; contourf(Ri.date,Ri.depths,Ri.prof',[0:0.05:0.25,100]); datetick3; shading interp; set(gca,'CLim',[0,0.25]); colorbar;




figure(17); xlim(datenum(yr,1,0)+modys); datetick3;
figure(18); xlim(datenum(yr,1,0)+modys); datetick3;




               chl: [108×1 double]
        cruiseName: 'NF08'
               den: [108×1 double]
              dena: [108×1 double]
                dt: [108×1 double]
             fname: 'D:\coral\FACE\Master_Data_Sheet_21aug09.xlsx'
                 h: [108×1 double]
               lat: [108×1 double]
               lon: [108×1 double]
                 n: [108×1 double]
           nfrt_ix: [108×1 double]
               nh4: [108×1 double]
               no2: [108×1 double]
               no3: [108×1 double]
    non_boil_stnix: [0×1 double]
                 p: [108×1 double]
               pha: [108×1 double]
                 s: [108×1 double]
                si: [108×1 double]
              stnm: {108×1 cell}
                 t: [108×1 double]
               tss: [108×1 double]
                 z: [108×1 double]






>> face.NF08
ans = 
  struct with fields:

       below_24_ix: [61×1 double]
       below_25_ix: [93×1 double]
               chl: [108×1 double]
            chlmax: 0.6500
            chlmin: 0
        cruiseName: 'NF08'
               den: [108×1 double]
              dena: [108×1 double]
              dmax: 27
              dmin: 0
                dt: [108×1 double]
             fname: 'D:\coral\FACE\Master_Data_Sheet_21aug09.xlsx'
                 h: [108×1 double]
               lat: [108×1 double]
               lon: [108×1 double]
                 n: [108×1 double]
           nfrt_ix: [108×1 double]
            nh4max: 20
            nh4min: 0
              nmax: 40
              nmin: 0
            no2max: 1.6000
            no2min: 0
    non_boil_stnix: [0×1 double]
                 p: [108×1 double]
               pha: [108×1 double]
            phamax: 0.3500
            phamin: 0
              pmax: 2
              pmin: 0
                 s: [108×1 double]
                si: [108×1 double]
             simax: 20
             simin: 0
              smax: 36.8000
              smin: 34.9000
              stnm: {108×1 cell}
                 t: [108×1 double]
              tmax: 29.1000
              tmin: 7.9000
               tss: [108×1 double]
                 z: [108×1 double]
              zmax: 0
              zmin: -290





>> face.NF09
ans = 
  struct with fields:

       below_24_ix: [17×1 double]
       below_25_ix: [20×1 double]
               chl: [194×1 double]
            chlmax: 0.6500
            chlmin: 0
        cruiseName: 'NF09'
               den: [194×1 double]
              dena: [194×1 double]
              dmax: 27
              dmin: 0
                dt: [194×1 double]
             fname: 'D:\coral\FACE\Master_Data_Sheet_01Feb10.xls'
                 h: [194×1 double]
               lat: [194×1 double]
               lon: [194×1 double]
                 n: [194×1 double]
           nfrt_ix: [124×1 double]
               nh4: [194×1 double]
            nh4max: 20
            nh4min: 0
              nmax: 40
              nmin: 0
               no2: [194×1 double]
            no2max: 1.6000
            no2min: 0
               no3: [194×1 double]
    non_boil_stnix: [100×1 double]
                 p: [194×1 double]
               pha: [194×1 double]
            phamax: 0.3500
            phamin: 0
              pmax: 2
              pmin: 0
                 s: [194×1 double]
                si: [194×1 double]
             simax: 20
             simin: 0
              smax: 36.8000
              smin: 34.9000
              stnm: {194×1 cell}
                 t: [194×1 double]
              tmax: 29.1000
              tmin: 7.9000
                 z: [194×1 double]
              zmax: 0
              zmin: -29





    %ztop = 10;
    ztop = 0;
    %Ttop = 24;
    Ttop = 23;
    Tbtm = 12;
    %Tbtm = 0;
    Nbtm = 1;
    Pbtm = 0;
    Sibtm = 0;




  fmg; plot_ts(face.brwd20.brwd20_top.seatemp,face.brwd20.brwd20_btm.seatemp,face.brwd40.brwd40_top.seatemp,face.brwd40.brwd40_btm.seatemp,face.brwd100.brwd100_sfc.seatemp,face.brwd100.brwd100_upp.seatemp); legend('20@sfc','20@btm','40@sfc','40@btm','100@sfc','100@40m'); axis(ax); datetick3('x','mmm-dd HH:MM','keeplimits');





fmg; spt(2,1,1); boxplot_ts(face.deep.adcp_sfc_l); ylim([-1,1]); grid on; ylabel('FACE Dp. L Sfc.');
titlename('FACE Deep (26.2 m) & NSUOC NW (11m) along-shore surface');
spt(2,1,2); boxplot_ts(sfomc.nw_w_btm.adcp_l_uwc_hourly); ylim([-1,1]); grid on; ylabel('NSUOC NW L UWC');

fmg; spt(2,1,1); boxplot_ts(face.deep.adcp_sfc_x); ylim([-1,1]); grid on; ylabel('FACE Dp. X Sfc.');
titlename('FACE Deep (26.2 m) & NSUOC NW (11m) cross-shore surface');
spt(2,1,2); boxplot_ts(sfomc.nw_w_btm.adcp_x_uwc_hourly); ylim([-1,1]); grid on; ylabel('NSUOC NW X UWC');

fmg; spt(2,1,1); boxplot_ts(face.deep.adcp_btm_l); ylim([-1,1]); grid on; ylabel('FACE Dp. L Dp.');
titlename('FACE Deep (26.2 m) & NSUOC NW (11m) along-shore near-bottom');
spt(2,1,2); boxplot_ts(sfomc.nw_w_btm.adcp_l_uwc_hourly); ylim([-1,1]); grid on; ylabel('NSUOC NW L Dp.');

fmg; spt(2,1,1); boxplot_ts(face.deep.adcp_btm_x); ylim([-1,1]); grid on; ylabel('FACE Dp. X Dp.');
titlename('FACE Deep (26.2 m) & NSUOC NW (11m) cross-shore near-bottom');
spt(2,1,2); boxplot_ts(sfomc.nw_w_btm.adcp_x_btm_hourly); ylim([-1,1]); grid on; ylabel('NSUOC NW X Dp.');

fmg; spt(2,1,1); boxplot_ts(ts_op(face.deep.adcp_sfc_l,face.deep.adcp_btm_l,'-')); ylim([-1,1]); grid on; ylabel('FACE Dp. L Sfc. - Btm.');
titlename('FACE Deep site (26.2 m) surface-bottom velocity shears');
spt(2,1,2); boxplot_ts(ts_op(face.deep.adcp_sfc_x,face.deep.adcp_btm_x,'-')); ylim([-1,1]); grid on; ylabel('FACE Dp. X Sfc. - Btm.');

fmg; spt(2,1,1); boxplot_ts(ts_op(sfomc.nw_w_btm.adcp_l_uwc,sfomc.nw_w_btm.adcp_l_btm,'-')); ylim([-1,1]); grid on; ylabel('NSUOC NW L Sfc. - Btm.');
titlename('NSUOC NW (11m) site surface-bottom velocity shears');
spt(2,1,2); boxplot_ts(ts_op(sfomc.nw_w_btm.adcp_x_uwc,sfomc.nw_w_btm.adcp_x_btm,'-')); ylim([-1,1]); grid on; ylabel('NSUOC NW X Sfc. - Btm.');






PVGF1: COORD CHANGES:
-80.092  ==>  -80.109
26.093   ==>  26.092




    mbar = []; baro = []; clear mbar baro





% yImage = imread('as.jpg');
% set(handles.axes7,'Units','pixels');
% resizePos = get(handles.axes7,'Position');
% myImage= imresize(myImage, [resizePos(3) resizePos(3)]);
% axes(handles.axes7);
% imshow(myImage);
% set(handles.axes7,'Units','normalized');







fmg; plot(1:10,sin(1:10)); inset_image('CoRIS/Flag-map_of_Florida.png');




function ax = inset_image(fname,pos)
%function ax = inset_image(fname,pos)
%
% Display the image file FNAME at position POS in the current FIGURE.
% Creates a new AXES to contain the inset image.
%
% Last Saved Time-stamp: <Fri 2018-03-02 14:56:12 Eastern Standard Time gramer>

  fh = gcf;
  pax = gca;

  if ( ~exist(fname,'file') )
    error('No image file named "%s"!',fname);
  end;
  if ( ~exist('pos','var') || isempty(pos) )
    %pos = [0.5 0.5 0.1 0.1];
    pos = [0.5 0.5 0.1];
  end;

  ext = 1;
  if ( numel(pos) == 3 )
    ext = pos(3);
    pos(3) = [];
  end;
  if ( numel(pos) == 2 )
    x = pos(1);
    y = pos(2);
    imginfo = imfinfo(fname);
    if ( imginfo.Height > imginfo.Width )
      scl = ext*(1-y)/imginfo.Height;
    else
      scl = ext*(1-x)/imginfo.Width;
    end;
    disp([imginfo.Width,imginfo.Height]);
    w = imginfo.Width * scl;
    h = imginfo.Height * scl;
    disp([scl,w,h]);
    pos = [x y w h];
  end;

  disp(pos);
  ax = axes('Position',pos);
  imshow(fname,'Parent',ax);

fmg; ax=gca; nax=axes('Position',[0.5,0.5,0.5,0.5]); imshow('CoRIS/Flag-map_of_Florida.png','Parent',nax,'Border','tight'); plot(ax,1:10,sin(1:10));

  set(ax,'XTickLabel',[]);
  set(ax,'YTickLabel',[]);

% yImage = imread('as.jpg');
% set(handles.axes7,'Units','pixels');
% resizePos = get(handles.axes7,'Position');
% myImage= imresize(myImage, [resizePos(3) resizePos(3)]);
% axes(handles.axes7);
% imshow(myImage);
% set(handles.axes7,'Units','normalized');


  %axes(pax);

return;










  img = imread(fname);




    sfomc.nw_w_btm.adcp_u_hourly = ts_nanify_gaps(interp_ts(sfomc.nw_w_btm.adcp_u),10);
    sfomc.nw_w_btm.adcp_v_hourly = ts_nanify_gaps(interp_ts(sfomc.nw_w_btm.adcp_v),10);
    sfomc.nw_w_btm.adcp_u_uwc_hourly = ts_nanify_gaps(interp_ts(sfomc.nw_w_btm.adcp_u_uwc),10);
    sfomc.nw_w_btm.adcp_v_uwc_hourly = ts_nanify_gaps(interp_ts(sfomc.nw_w_btm.adcp_v_uwc),10);
    sfomc.nw_w_btm.adcp_dir_hourly.date = sfomc.nw_w_btm.adcp_u_hourly.date;
    sfomc.nw_w_btm.adcp_dir_hourly.data = uv_to_dir_curr(sfomc.nw_w_btm.adcp_u_hourly.data,sfomc.nw_w_btm.adcp_v_hourly.data);






function [newdts,newts,newprof] = window_func_prof(dts,ts,prof,varargin)
%function [newdts,newts,newprof] = window_func_prof(dts,ts,prof,func,nhrs,samplehrs,wintyp,plotResults)
%
% Call WINDOW_FUNC repeatedly on each COLUMN of time series profile PROF.
% Similarly, call WINDOW_FUNC on DTS and (vector) TS.
%
% Last Saved Time-stamp: <Tue 2018-02-27 19:26:44 Eastern Standard Time gramer>

  [newdts,newts] = window_func(dts,ts,varargin{:});
  [newdts_prof,newprof(:,1)] = window_func(dts,prof(:,1),varargin{:});
  if ( numel(newdts) ~= numel(newdts_prof) || any(newdts(:) ~= newdts_prof(:)) )
    error('Data-date mismatch with profile-date!');
  end;
  for ix = 1:size(prof,2);
    [ig,newprof(:,ix)] = window_func(dts,prof(:,ix),varargin{:});
    clear ig
  end;

return;








function stn = station_partition_periods(stn,fld,pers)
%function stn = station_partition_periods(stn,fld[,pers])
%
% Partition energy in time series STN.(FLD) into per-hour bands given by
% Nx2 vector PERS (DEFAULT: 3-11h,11-16h,17-28h,2-10d,150-210d,330-390d).
%
% Returns STN with new TS fields, e.g., STN.([FLD,'_3_h_hp_8_h_lp']).
%
% Last Saved Time-stamp: <Sun 2018-02-25 20:36:41 Eastern Standard Time gramer>

  if ( ~exist('pers','var') || isempty(pers) )
    pers = [ 3,11 ; 11,16 ; 17,28 ; 2*24,10*24 ; 150*24,210*24 ; 330*24,390*24 ];
  end;

  for perix = 1:size(pers,1)
    perfld = sprintf('%s_%g_h_hp_%g_h_lp',fld,pers(perix,1),pers(perix,2));
    %DEBUG:    disp(perfld);
    try,
      stn = verify_variable(stn,perfld);
    catch,
      %catchwarn;
    end;
  end;

return;







  % Shallow-site energy partitions
  % % What was our fundamental sampling period?
  % fmg; plot(sfomc.nw_w_btm.adcp_x_btm.date(1:end-1),24*diff(sfomc.nw_w_btm.adcp_x_btm.date),'.'); datetick3; ylim([0,8]);
  pers = [ 3,11 ; 11,16 ; 17,28 ; 2*24,10*24 ; 150*24,210*24 ; 330*24,390*24 ];
  sfomc.nw_w_btm = station_partition_periods(sfomc.nw_w_btm,'seatemp',pers);
  sfomc.nw_w_btm = station_partition_periods(sfomc.nw_w_btm,'adcp_x_btm',pers);
  sfomc.nw_w_btm = station_partition_periods(sfomc.nw_w_btm,'adcp_l_btm',pers);




    pers = [ 2,11 ; 11,16 ; 17,28 ; 2*24,10*24 ; 150*24,210*24 ; 330*24,390*24 ];
    face.deep = station_partition_periods(face.deep,'seatemp',pers);
    face.deep = station_partition_periods(face.deep,'adcp_btm_x',pers);
    face.deep = station_partition_periods(face.deep,'adcp_btm_l',pers);






  if ( ~exist('pers','var') || isempty(pers) )
    pers = [ 2,8 ; 10,14 ; 17,28 ; 2*24,10*24 ; 150*24,210*24 ; 330*24,390*24 ];
  end;



    plot_ts(face.deep.adcp_btm_x_2_h_hp_8_h_lp,face.deep.adcp_btm_x_10_h_hp_14_h_lp,face.deep.adcp_btm_x_17_h_hp_28_h_lp,face.deep.adcp_btm_x_48_h_hp_240_h_lp);






face.deep = station_partition_periods(face.deep,'seatemp');
face.deep = station_partition_periods(face.deep,'adcp_btm_x');

fmg;
subplot(2,1,1);
plot_ts(face.deep.seatemp_3_h_hp_8_h_lp,face.deep.seatemp_10_h_hp_14_h_lp,face.deep.seatemp_17_h_hp_28_h_lp);
ylabel('Sea temp.');
legend('3-8h','10-14h','17-28h', 'Orientation','horiz');
subplot(2,1,2);
plot_ts(face.deep.adcp_btm_x_3_h_hp_8_h_lp,face.deep.adcp_btm_x_10_h_hp_14_h_lp,face.deep.adcp_btm_x_17_h_hp_28_h_lp);
ylabel('Btm. u_x');
legend('3-8h','10-14h','17-28h', 'Orientation','horiz');
xlim(datenum(2015,[8,10],[31,1])); datetick3;




      nyq = 2*hrs_per_sample; % Nyqvist frequency
disp(nyq);
      [B, A] = butter(filterOrder, (nyq/nhrs), 'high');






function newts = ts_nanify_gaps(ts,maxgap)
%function newts = ts_nanify_gaps(ts,maxgap)
%
% Construct a new time series STRUCT NEWTS containing all the data in TS, but
% where ever FIND_DATE_RANGES (v.) finds a gap of longer than MAXGAP days in
% TS, insert an extra NaN at the beginning and end of that gap.
%
% This function is useful for example when calling PLOT_TS (v.) with linespec
% such as '.-': for the insertion of NaNs will cause the removal of the line
% joining points in the plot on either side of any large gap.
%
% DEFAULT MAXGAP: 3*MEDIAN(DIFF(TS.date))
%
% Last Saved Time-stamp: <Sun 2018-02-25 17:29:47 Eastern Standard Time gramer>

  if ( ~is_valid_ts(ts) )
    error('First arg must be a valid time series STRUCT');
  end;

  % "Normal" time step in time series
  dt = median(diff(ts.date));

  if ( ~exist('maxgap','var') || isempty(maxgap) )
    %maxgap = [];
    maxgap = 3*dt;
  end;
  if ( maxgap < 2*dt )
    error('MAXGAP is too short for sampling frequency of time series!');
  end;

  newts = ts;

  [rgs,ixes,allixes] = find_date_ranges(ts.date,maxgap);

  for ixix=size(ixes,2):-1:2
    % newts.date(ixes(2,ixix)+1:end+1) = newts.date(ixes(2,ixix):end);
    % newts.data(ixes(2,ixix)+1:end+1) = newts.data(ixes(2,ixix):end);
    % newts.date(ixes(2,ixix)) = newts.date(ixes(2,ixix))-dt;
    % newts.data(ixes(2,ixix)) = nan;

    newts.date(ixes(1,ixix)+1:end+1) = newts.date(ixes(1,ixix):end);
    newts.data(ixes(1,ixix)+1:end+1) = newts.data(ixes(1,ixix):end);
    newts.date(ixes(1,ixix)) = newts.date(ixes(1,ixix))-dt;
    newts.data(ixes(1,ixix)) = nan;
  end;

return;








  %DEBUG:  disp(sprintf('%s - %s ::',datestr(newts.date(1),0),datestr(newts.date(end),0)));
  %DEBUG:  for ixix=1:size(ixes,2); disp(sprintf('%s -> %s',datestr(ts.date(ixes(1,ixix)),0),datestr(newts.date(ixes(2,ixix)),0))); end;




    %DEBUG:    disp(ixes(:,ixix)');
    %DEBUG:    disp(sprintf('%s -> %s -> %s -> %s',datestr(newts.date(ixes(1,ixix)-1),0),datestr(newts.date(ixes(1,ixix))-dt,0),datestr(newts.date(ixes(1,ixix)),0),datestr(newts.date(ixes(2,ixix)),0)));
    %DEBUG:    disp(sprintf('%s -> %s -> %s',datestr(newts.date(ixes(1,ixix)-1),0),datestr(newts.date(ixes(1,ixix))-dt,0),datestr(newts.date(ixes(1,ixix)),0)));



  %DEBUG:  for ixix=1:size(ixes,2); disp(sprintf('%s -> %s',datestr(ts.date(ixes(1,ixix)),0),datestr(newts.date(ixes(2,ixix)),0))); end;











          'Orientation','horiz', 'Location','SouthWest');



          'Location','South');



    legend('Location','SouthWest');





  if ( ~exist('nw_adcp_l_uwc','var') )
tic,
    nw_adcp_l_uwc = ts_nanify_gaps(sfomc.nw_w_btm.adcp_l_uwc);
toc,
    nw_adcp_l_btm = ts_nanify_gaps(sfomc.nw_w_btm.adcp_l_btm);
    nw_adcp_x_uwc = ts_nanify_gaps(sfomc.nw_w_btm.adcp_x_uwc);
    nw_adcp_x_btm = ts_nanify_gaps(sfomc.nw_w_btm.adcp_x_btm);
  end;






  % Shallow-site currents

  if ( ~exist('nw_adcp_l_uwc','var') )
    nw_adcp_l_uwc = ts_nanify_gaps(sfomc.nw_w_btm.adcp_l_uwc);
    nw_adcp_l_btm = ts_nanify_gaps(sfomc.nw_w_btm.adcp_l_btm);
    nw_adcp_x_uwc = ts_nanify_gaps(sfomc.nw_w_btm.adcp_x_uwc);
    nw_adcp_x_btm = ts_nanify_gaps(sfomc.nw_w_btm.adcp_x_btm);
  end;

  shfh = fmg;
  shax(1) = subplot(7,1,1:3); title('SHALLOW (NW) site currents vs. sea temperatures');
  %plot_ts(fwyf1.ndbc_air_t,'k-',sfomc.ne_buoy.seatemp,sfomc.sw_buoy.seatemp,sfomc.nw_w_btm.adcp_seatemp);
  plot_ts(fwyf1.ndbc_air_t,'k-',sfomc.ne_buoy.at_30m.sbe_seatemp,sfomc.sw_buoy.at_15m.sbe_seatemp,sfomc.nw_w_btm.adcp_seatemp);
  xlim(wholedts);
  ylim([16.5,32.5]);
  grid on; 
  legend('Air (FWY)','Deep (NE 30 m)','Mid (SW 15 m)','Shallow (NW 11 m)', legargs{:});
  set(gca,'FontSize',fontsize); set(get(gca,'Child'),'LineWidth',linewid);

  shax(2) = subplot(7,1,4:5); ylabel('Alongshore'); %title('Shallow alongshore current');
  plot_ts(ts_nanify_gaps(sfomc.nw_w_btm.adcp_l_uwc),ts_nanify_gaps(sfomc.nw_w_btm.adcp_l_btm));
  ylim([-1.5,+1.5]); xlim(wholedts);
  grid on; center_axes; 
  legend('Upper Water Column (NW)','Near Seafloor (NW)', legargs{:});
  set(gca,'FontSize',fontsize); set(get(gca,'Child'),'LineWidth',linewid);

  shax(3) = subplot(7,1,6:7); ylabel('Cross-shore'); %title('Shallow cross-shore current');
  plot_ts(ts_nanify_gaps(sfomc.nw_w_btm.adcp_x_uwc),ts_nanify_gaps(sfomc.nw_w_btm.adcp_x_btm));
  ylim([-0.4,+0.4]); xlim(wholedts);
  grid on; center_axes;
  legend('Upper Water Column (NW)','Near Seafloor (NW)', legargs{:});
  datetick('x',2,'keeplimits'); %datetick3;
  set(gca,'FontSize',fontsize); set(get(gca,'Child'),'LineWidth',linewid);
  set(shax(1:2),'XTickLabel','');








   legend(['Air temp. (FDEPK)'],...
          ['Offshore (',num2str(sefcri.updb.depth),' m)'],...
          ['Near-shore (',num2str(sefcri.mc2.depth),' m)'],...
          'Location','SouthWest');

   legend(['Air temp. (LKWF1)'],...
          ['Offshore (',num2str(sefcri.pb2.depth),' m)'],...
          ['Near-shore (',num2str(sefcri.pb1.depth),' m)'],...
          'Location','SouthWest');

   legend(['Air temp. (PVGF1)'],...
          ['Offshore (',num2str(sefcri.bc3.depth),' m)'],...
          ['Near-shore (',num2str(sefcri.bc1.depth),' m)'],...
          'Location','SouthWest');








  set(gca,'FontSize',fontsize);






  legargs = { 'Orientation','horiz','Location','northeast' };





  legend('Cape Canaveral (1988-2016)','Lake Worth (2001-2016)','Port Everglades (2009-2017)','RSMAS Rooftop (2003-2011)','Fowey Rocks (1991-2016)', legargs{:});




  pause;

  disp('Broader focus on periods of high-frequency variability');
  axes(dpax(end)); xlim(datenum(1999,[7,8],[15,29])); datetick3;
  if doPrint; print('-dpng',fullfile(figspath,'upwelling-sfomc-1999-deep-Jul-Aug.png')); end;
  axes(shax(end)); xlim(datenum(1999,[7,8],[15,29])); datetick3;
  if doPrint; print('-dpng',fullfile(figspath,'upwelling-sfomc-1999-shallow-Jul-Aug.png')); end;
  axes(ekax(end)); xlim(datenum(1999,[7,8],[15,29])); datetick3;
  if doPrint; print('-dpng',fullfile(figspath,'upwelling-sfomc-1999-ekman-Jul-Aug.png')); end;
  axes(wsax(end)); xlim(datenum(1999,[7,8],[15,29])); datetick3;
  if doPrint; print('-dpng',fullfile(figspath,'upwelling-sfomc-1999-stress-Jul-Aug.png')); end;
  axes(vmax(end)); xlim(datenum(1999,[7,8],[15,29])); datetick3;
  if doPrint; print('-dpng',fullfile(figspath,'upwelling-sfomc-1999-variab-Jul-Aug.png')); end;
  pause;

  disp('Close focus on period of low-frequency (initial) variability');
  axes(dpax(end)); xlim(datenum(1999,7,[22,28])); datetick3;
  if doPrint; print('-dpng',fullfile(figspath,'upwelling-sfomc-1999-deep-Jul-lowfreq.png')); end;
  axes(shax(end)); xlim(datenum(1999,7,[22,28])); datetick3;
  if doPrint; print('-dpng',fullfile(figspath,'upwelling-sfomc-1999-shallow-Jul-lowfreq.png')); end;
  axes(ekax(end)); xlim(datenum(1999,7,[22,28])); datetick3;
  if doPrint; print('-dpng',fullfile(figspath,'upwelling-sfomc-1999-ekman-Jul-lowfreq.png')); end;
  axes(wsax(end)); xlim(datenum(1999,7,[22,28])); datetick3;
  if doPrint; print('-dpng',fullfile(figspath,'upwelling-sfomc-1999-stress-Jul-lowfreq.png')); end;
  axes(vmax(end)); xlim(datenum(1999,7,[22,28])); datetick3;
  if doPrint; print('-dpng',fullfile(figspath,'upwelling-sfomc-1999-variab-Jul-lowfreq.png')); end;
  pause;

  disp('CLOSE focus on period of high-frequency variability');
  axes(dpax(end)); xlim(datenum(1999,7,[28,32])); datetick3;
  if doPrint; print('-dpng',fullfile(figspath,'upwelling-sfomc-1999-deep-Jul-highfreq.png')); end;
  axes(shax(end)); xlim(datenum(1999,7,[28,32])); datetick3;
  if doPrint; print('-dpng',fullfile(figspath,'upwelling-sfomc-1999-shallow-Jul-highfreq.png')); end;
  axes(ekax(end)); xlim(datenum(1999,7,[28,32])); datetick3;
  if doPrint; print('-dpng',fullfile(figspath,'upwelling-sfomc-1999-ekman-Jul-highfreq.png')); end;
  axes(wsax(end)); xlim(datenum(1999,7,[28,32])); datetick3;
  if doPrint; print('-dpng',fullfile(figspath,'upwelling-sfomc-1999-stress-Jul-highfreq.png')); end;
  axes(vmax(end)); xlim(datenum(1999,7,[28,32])); datetick3;
  if doPrint; print('-dpng',fullfile(figspath,'upwelling-sfomc-1999-variab-Jul-highfreq.png')); end;
  pause;

  % disp('Whole period when all current meters were operating');
  % axes(dpax(end)); xlim(datenum(1999,[9,11],[18,12])); datetick3;
  % axes(shax(end)); xlim(datenum(1999,[9,11],[18,12])); datetick3;
  % axes(ekax(end)); xlim(datenum(1999,[9,11],[18,12])); datetick3;
  % axes(wsax(end)); xlim(datenum(1999,[9,11],[18,12])); datetick3;
  % pause;

  disp('Upwelling event when three current meters were operating');
  axes(dpax(end)); xlim(datenum(1999,7,[22,28])); datetick3;
  axes(shax(end)); xlim(datenum(1999,7,[22,28])); datetick3;
  axes(ekax(end)); xlim(datenum(1999,7,[22,28])); datetick3;
  axes(wsax(end)); xlim(datenum(1999,7,[22,28])); datetick3;
  axes(vmax(end)); xlim(datenum(1999,7,[22,28])); datetick3;
  pause;


  % %% 1999, 2000, 2003, 2004, 2005?, 2006, 2007, 2009?, 2010?, 2011
  % modys = datenum(0,[7,9],[5,5]);
  % for yr=[1999, 2000, 2003, 2004, 2005, 2006, 2007, 2011];
  modys = datenum(0,[6,9],[5,5]);
  for yr=[1999, 2000, 2003, 2004, 2005, 2006, 2007, 2009, 2010, 2011];
    yrstr = num2str(yr);
    disp(['Focus on high-frequency variability in 5-min sea temperature in *',yrstr,'*']);
    axes(dpax(end)); xlim(datenum(yr,1,0)+modys); datetick3;
    if doPrint; print('-dpng',fullfile(figspath,['upwelling-sfomc-',yrstr,'-deep-wide.png'])); end;
    axes(shax(end)); xlim(datenum(yr,1,0)+modys); datetick3;
    if doPrint; print('-dpng',fullfile(figspath,['upwelling-sfomc-',yrstr,'-shallow-wide.png'])); end;
    axes(ekax(end)); xlim(datenum(yr,1,0)+modys); datetick3;
    if doPrint; print('-dpng',fullfile(figspath,['upwelling-sfomc-',yrstr,'-ekman-wide.png'])); end;
    axes(wsax(end)); xlim(datenum(yr,1,0)+modys); datetick3;
    if doPrint; print('-dpng',fullfile(figspath,['upwelling-sfomc-',yrstr,'-stress-wide.png'])); end;
    axes(vmax(end)); xlim(datenum(yr,1,0)+modys); datetick3;
    if doPrint; print('-dpng',fullfile(figspath,['upwelling-sfomc-',yrstr,'-variab-wide.png'])); end;
    pause;
  end;
  %xlim([shax,dpax,ekax,wsax],datenum(1999,1,0)+modys); jumpax([shax,dpax,ekax,wsax],nan,365);
  %keyboard;


  disp('Done');
  pause;











%xlim(datenum(1999,[7,8],[21,1])); datetick3;
%xlim(datenum(1999,[7,8],[11,11])); datetick3;



 NO_2^-+NO_3^-



  fmgl;
  %nkax(1) = subplot(8,1,1);
  surf(nw.Nx.date,nw.Nx.depths,nw.Nx.cprof');
  view(0,270); shading interp; set(gca,'YDir','reverse');
  hold on;
  colorbar_tight;
  ylim([-sfomc.nw_w_btm.depth,0]); 
  %for dp=nw.Nx.depths(:)'; annotline([],dp); end;
  xlim(xl); datetick('x',2,'keeplimits'); %datetick3;
  ylabel('Shallow u \bullet [NO2+NO3] [kg\bulletm^-^1]');
  if doPrint; print('-dpng',fullfile(figspath,[fignm,'_shallow_N_cumprof.png'])); end;

  fmgl;
  %nkax(2) = subplot(8,1,2:3);
  %surf(sw.Nx.date,sw.Nx.depths,sw.Nx.prof');
  %surf(sw.Nx.date,sw.Nx.depths,sw.Nx.xprof');
  surf(sw.Nx.date,sw.Nx.depths,sw.Nx.cprof');
  view(0,270); shading interp; set(gca,'YDir','reverse');
  hold on;
  colorbar_tight;
  ylim([-sfomc.sw_buoy.depth,0]); 
  for dp=sw.Nx.depths(:)'; annotline([],dp); end;
  xlim(xl); datetick('x',2,'keeplimits'); %datetick3;
  ylabel('Mid (SW/C) u \bullet [NO2+NO3] [kg\bulletm^-^1]');
  if doPrint; print('-dpng',fullfile(figspath,[fignm,'_mid_N_cumprof.png'])); end;

  fmgl;
  %nkax(3) = subplot(8,1,4:8);
  % surf(ne.Nx.date,ne.Nx.depths,ne.Nx.prof');
  % surf(ne.Nx.date,ne.Nx.depths,ne.Nx.xprof');
  surf(ne.Nx.date,ne.Nx.depths,ne.Nx.cprof');
  view(0,270); shading interp; set(gca,'YDir','reverse');
  hold on;
  colorbar_tight;
  ylim([-sfomc.ne_buoy.depth,0]); 
  for dp=ne.Nx.depths(:)'; annotline([],dp); end;
  xlim(xl); datetick('x',2,'keeplimits'); %datetick3;
  ylabel('Deep (NE/E) u \bullet [NO2+NO3] [kg\bulletm^-^1]');
  if doPrint; print('-dpng',fullfile(figspath,[fignm,'_deep_N_cumprof.png'])); end;










fmgl;
%surf(sw.Nx.date,sw.Nx.depths,sw.Nx.prof');
%surf(sw.Nx.date,sw.Nx.depths,sw.Nx.xprof');
surf(sw.Nx.date,sw.Nx.depths,sw.Nx.cprof');
view(0,270); shading interp; set(gca,'YDir','reverse');
hold on;
colorbar_tight;
ylim([-sfomc.sw_buoy.depth,0]); 
for dp=sw.Nx.depths(:)'; annotline([],dp); end;
xlim(xl); datetick('x',2,'keeplimits'); %datetick3;
ylabel('Mid (SW/C) u \bullet [NO2+NO3] [kg\bulletm^-^1]');
if doPrint; print('-dpng',fullfile(figspath,[fignm,'_mid_N_prof.png'])); end;


fmgl;
% surf(ne.Nx.date,ne.Nx.depths,ne.Nx.prof');
% surf(ne.Nx.date,ne.Nx.depths,ne.Nx.xprof');
surf(ne.Nx.date,ne.Nx.depths,ne.Nx.cprof');
view(0,270); shading interp; set(gca,'YDir','reverse');
hold on;
colorbar_tight;
ylim([-sfomc.ne_buoy.depth,0]); 
for dp=ne.Nx.depths(:)'; annotline([],dp); end;
xlim(xl); datetick('x',2,'keeplimits'); %datetick3;
ylabel('Deep (NE/E) u \bullet [NO2+NO3] [kg\bulletm^-^1]');
if doPrint; print('-dpng',fullfile(figspath,[fignm,'_deep_N_prof.png'])); end;



fmg; plot_ts(ne.intNx,sw.intNx,nw.intNx); legend('NE','SW','NW','Location','East');
ylabel('\int_0^h u \bullet [NO2+NO3] [kg\bulletm^-^1]');
titlename('Accumulated NO_2^-+NO_3^- per meter of reef');
xlim(xl); datetick('x',2,'keeplimits'); %datetick3;
if doPrint; print('-dpng',fullfile(figspath,[fignm,'_cum_N.png'])); end;










  set(tpax(1:2),'XTickLabel','');
  set(xsax(1:2),'XTickLabel','');





    % On top subplots, leave room for a possible Title!
    if ( min(nix) == 1 )
      sh = sh - tbuf;
    end;



disp('Hit enter to continue...'); pause;
xlim([datenum(1999,7,21),datenum(1999,8,1)]); datetick3;
if doPrint; print('-dpng',fullfile(figspath,[fignm,'_event_2_1999.png'])); end;






%N+N
  ylim([0.0,0.5]); colorbar; %Adding useless COLORBAR just lines up the X axes nicely
%P
  ylim([0,0.1]); colorbar; %Adding useless COLORBAR just lines up the X axes nicely
%Si
  ylim([0,0.35]); colorbar; %Adding useless COLORBAR just lines up the X axes nicely



xlim(datenum(1999,[7,8],[15,29])); datetick3;




%%%
% NW/W mooring [11 m isobath) nutrient concentrations

nw.N = nwtemp;
% nw.N.data = sw.N.data.*(nwtemp.data-swtemp.prof(:,end))./(Tsfcpt-swtemp.prof(:,end));
% nw.N.data(nw.N.data<0 | nw.N.data>sw.N.data) = 0;
nw.N.data = sw.N.prof(:,end).*(1-(nwtemp.data-swtemp.prof(:,end))./(Tsfcpt-swtemp.prof(:,end)));
nw.N.data(0>nw.N.data | nw.N.data>sw.N.prof(:,end)) = 0;

nw.P = nwtemp;
nw.P.data = sw.P.prof(:,end).*(nwtemp.data-swtemp.prof(:,end))./(Tsfcpt-swtemp.prof(:,end));
nw.P.data(0>nw.P.data | nw.P.data>sw.P.prof(:,end)) = 0;

nw.Si = nwtemp;
nw.Si.data = sw.Si.prof(:,end).*(nwtemp.data-swtemp.prof(:,end))./(Tsfcpt-swtemp.prof(:,end));
nw.Si.data(0>nw.Si.data | nw.Si.data>sw.Si.prof(:,end)) = 0;














btmix = find(nw.Nx.depths<nw.N.depths(end));
nw.Nx.prof(:,btmix) = repmat(nw.Nx.prof(:,btmix(1)-1),[1,numel(btmix)]);
nw.Nx.xprof(:,btmix) = repmat(nw.Nx.xprof(:,btmix(1)-1),[1,numel(btmix)]);
nw.Nx.cprof(:,btmix) = repmat(nw.Nx.cprof(:,btmix(1)-1),[1,numel(btmix)]);






ne.adcp_onshore.prof(ne.adcp_onshore.prof>0)=0;





ix11m = 38; ix30m = 19; ix45m = 4;





sfomc.ne_buoy.adcp_x_11m.date = sfomc.ne_buoy.adcp_x.date;
sfomc.ne_buoy.adcp_x_11m.data = sfomc.ne_buoy.adcp_x.prof(:,38);
sfomc.ne_buoy.adcp_x_30m.date = sfomc.ne_buoy.adcp_x.date;
sfomc.ne_buoy.adcp_x_30m.data = sfomc.ne_buoy.adcp_x.prof(:,19);
sfomc.ne_buoy.adcp_x_45m.date = sfomc.ne_buoy.adcp_x.date;
sfomc.ne_buoy.adcp_x_45m.data = sfomc.ne_buoy.adcp_x.prof(:,4);

upwix = find(abs(sfomc.ne_buoy.at_45m.sbe_seatemp.data-23)<2e-4);

fmg;
ax(1)=spt(4,1,[1:3],'sharex');
 plot_ts(sfomc.ne_buoy.at_45m.sbe_seatemp,sfomc.ne_buoy.at_30m.sbe_seatemp,sfomc.nw_w_btm.adcp_seatemp,fwyf1.ndbc_sea_t,'k');
 ylim([15,31]); grid on;
 hold on; plot(ne.Nx.date(ne.Nx.prof(:,end)>0),repmat(23,size(ne.Nx.date(ne.Nx.prof(:,end)>0))),'r^');
ax(2)=spt(4,1,4,'sharex');
 plot_ts(sfomc.ne_buoy.adcp_x_45m,sfomc.ne_buoy.adcp_x_30m,sfomc.nw_w_btm.adcp_x_btm,sfomc.ne_buoy.adcp_x_11m);
 grid on;
 hold on; plot(ne.Nx.date(ne.Nx.prof(:,end)>0),repmat(-0.5,size(ne.Nx.date(ne.Nx.prof(:,end)>0))),'r^');
xlim(datenum(1999,[7,8],[21,1])); datetick3;





xlim(datenum(1999,[7,8],[27,1])); datetick3;



fmg; spt(4,1,[1:3],'sharexy'); plot_ts(sfomc.ne_buoy.at_45m.sbe_seatemp,sfomc.ne_buoy.at_30m.sbe_seatemp,sfomc.nw_w_btm.adcp_seatemp,sfomc.nw_w_btm.sbe_seatemp,'k'); ylim([15,31]); grid on; spt(4,1,4,'sharexy'); plot_ts(sfomc.ne_buoy.adcp_x_btm,sfomc.nw_w_btm.adcp_x_btm); grid on; xlim(datenum(1999,[7,8],[15,1])); datetick3;




  fmg;
  spt(3,1,1); plot_bathy_transect(sfomc.c_buoy,[],sfomc.c_buoy.isobath_orientation-90);
  spt(3,1,2); plot_bathy_transect({sfomc.c_buoy,sfomc.sw_buoy},[],sfomc.sw_buoy.isobath_orientation-90);
  spt(3,1,3); plot_bathy_transect({sfomc.c_buoy,face.deep},[],face.deep.isobath_orientation-90);





lh=plot([sfomc.lons,face.lons],[sfomc.lats,face.lats],'wp'); 


  % Buffers around entire figure
  lbuf = 0.08;
  rbuf = 0.02;
  tbuf = 0.02;
  bbuf = 0.05;

  % Buffer around each subplot
  buf=0.03;






    dy = ((1-(tbuf+bbuf))/n) - (2*buf);
    dx = ((1-(lbuf+rbuf))/m) - (2*buf);
    % Tom and Dick
    [mix,nix] = ind2sub([m,n],ix(:));
    if ( max(nix) > n )
      error('Ecoforecasts:subplot_tight:SubplotIndexTooLarge',...
            'Index exceeds number of subplots.');
    end;
    nm = max(mix)-min(mix)+1;
    nn = max(nix)-min(nix)+1;
    ax = subplot('position',...
                 [((min(mix)-1)*(dx+2*buf))+lbuf,((n-max(nix))*(dy+2*buf))+bbuf,(dx*nm),(dy*nn)]);







    ax = subplot('position',...
                 [((min(mix)-1)*(dx+2*buf))+lbuf,((n-max(nix)+1)*(dy+2*buf))+bbuf,(dx*nm),(dy*nn)]);





function [anomts,climts,loasdts,hiasdts,lopctts,hipctts,cumfun,per] = anomalize_ts_to_ts(ts,varargin)
%function [anomts,climts,loasdts,hiasdts,lopctts,hipctts,cumfun,per] = anomalize_ts_to_ts(ts[,smoothing][,pct][,GRP_TS args])
%
% Call ANOMALIZE_TS (v.), but return time series with time stamps identical
% to those of the ANOMTS, for climatology, climatology minus anomaly standard
% deviation, climatology plus anomaly standard deviation, PCTth (or LOPCTth)
% percentile, 100 - PCTth (or HIPCTth) percentile. Otherwise, arguments and
% return values are the same as ANOMALIZE_TS.
%
% If optional SMOOTHING is LOGICAL (v.) and true (DEFAULT), vary all returned
% time series smoothly between time stamps, e.g., for a daily climatology of
% an hourly TS, do not have a single value for all 24 values of each day.
%
% SAMPLE CALL, producing a 'year-hour' (diurnal-annual) climatology:
%  [anomts,climts,loasdts,hiasdts] = anomalize_ts_to_ts(stn.ndbc_sea_t,@get_jhour_no_leap);
%
% SAMPLE CALL, with smoothed time series and using 7th and 93rd percentiles:
%  [anomts,climts,loasdts,hiasdts] = anomalize_ts_to_ts(stn.ndbc_sea_t,true,7);
%
% NOTE that for an hourly time series TS, both of the sample calls above
% would return the same time series for climatology, "HI" and "LO" SD!
%
% Last Saved Time-stamp: <Wed 2018-02-14 13:25:18 EST lew.gramer>

  args = varargin;

  % Should climatology time series be "stepped" (have one value per CUMFUN
  % period with jumps between), or vary smoothly between timestamps? 
  if ( numel(args) > 0 && islogical(args{1}) )
    % This is safe because neither ANOMALIZE_TS nor GRP_TS can have a value
    % of type LOGICAL as their first (valid) optional arg
    smoothing = args{1};
    args(1) = [];
  else
    smoothing = true;
  end;

  [anomts,clim,tid,asd,lopct,hipct,cumfun,per] = anomalize_ts(ts,args{:});
  
  climts.date = anomts.date;
  loasdts.date = anomts.date;
  hiasdts.date = anomts.date;
  lopctts.date = anomts.date;
  hipctts.date = anomts.date;

  if ( smoothing )
    % Smoothed time series
    dt = median(diff(ts.date));
    smoothtid = [datenum(0,1,1):dt:(datenum(1,1,1)-(1.1*dt))]';
    smoothclim = interp1(tid,clim,smoothtid)';
    smoothasd = interp1(tid,asd,smoothtid)';
    smoothlopct = interp1(tid,lopct,smoothtid)';
    smoothhipct = interp1(tid,hipct,smoothtid)';
    climts.data = interp1(smoothtid,smoothclim,cumfun(anomts.date));
    loasdts.data = interp1(tid,clim-asd,cumfun(anomts.date));
    hiasdts.data = interp1(tid,clim+asd,cumfun(anomts.date));
    lopctts.data = interp1(tid,clim+lopct,cumfun(anomts.date));
    hipctts.data = interp1(tid,clim+hipct,cumfun(anomts.date));
  else
    % "Stepped" time series
    climts.data = interp1(tid,clim,cumfun(anomts.date));
    loasdts.data = interp1(tid,clim-asd,cumfun(anomts.date));
    hiasdts.data = interp1(tid,clim+asd,cumfun(anomts.date));
    lopctts.data = interp1(tid,clim+lopct,cumfun(anomts.date));
    hipctts.data = interp1(tid,clim+hipct,cumfun(anomts.date));
  end;

return;



% >> 29.25-28.27
% ans =
%     0.9800
% >> 
% >> 27.97-19.73
% ans =
%     8.2400
% >> 0.98/8.24
% ans =
%     0.1189

mixprof = repmat(0.119,size(swtemp.prof));

% DELAY: 9 hours





mixprof = (1-(swtemp.prof-repmat(Tendpt,[1,size(swtemp.prof,2)])) ...
           ./ repmat((Tsfcpt-Tendpt),[1,size(swtemp.prof,2)]));
mixprof(0>mixprof | mixprof>1) = 0;

swNendpt = repmat(Nendpt,[1,size(swtemp.prof,2)]);
sw.N = swtemp;
sw.N.prof = swNendpt .* mixprof;
sw.N.data = nansum(sw.N.prof,2);

swPendpt = repmat(Pendpt,[1,size(swtemp.prof,2)]);
sw.P = swtemp;
sw.P.prof = swPendpt .* mixprof;
sw.P.data = nansum(sw.P.prof,2);

swSendpt = repmat(Sendpt,[1,size(swtemp.prof,2)]);
sw.Si = swtemp;
sw.Si.prof = swSendpt .* mixprof;
sw.Si.data = nansum(sw.Si.prof,2);








[nenitr,nephos,nesili,netemp,swtemp,nwtemp,sfctemp] = ...
    intersect_tses(ne.N,ne.P,ne.Si,ne.seatemp,sw.seatemp,nw.adcp_seatemp,sfomc.ne_buoy.at_0_6m.sbe_seatemp);



%ne.Nx.prof(:,btmix) = nan;
%ne.Nx.xprof(:,btmix) = nan;
%ne.Nx.cprof(:,btmix) = nan;




% fmg; plot(ne.Nx.date,sum(ne.Nx.cprof,2,'omitnan')); datetick3;
% ylabel('\Sigma u \bullet [NO2+NO3] [kg\bulletm\bullets^-^1]');





ylabel('NE (Deep) u \bullet [NO2+NO3] [\mumol\bulletm\bullets^-^1]');



%fmg; plot(ne.Nx.date,sum(cumsum(ne.Nx.xprof,'omitnan'),2)); datetick3;




ne.Nx.date = ne.adcp_x.date;
ne.Nx.depths = ne.adcp_depths(end:-1:1);
btmix = find(ne.Nx.depths<ne.N.depths(end));

ne.Nx.prof = interp2(ne.N.date,ne.N.depths,ne.N.prof',ne.Nx.date,ne.Nx.depths,'linear',NaN)';
%ne.Nx.prof(:,btmix) = nan;
ne.Nx.prof(:,btmix) = repmat(ne.Nx.prof(:,btmix(1)-1),[1,numel(btmix)]);

% Nutrient mass per profile bin per meter of reef-line per sample:
%  Cross-shore-velocity x nutrient-conc. x bin-height x time-per-sample
%   [m/s] * [kg/m^3] * [m] * [s] = [kg/m]
ne.Nx.xprof = ne.adcp_x.prof .* umol_to_kgm3(ne.Nx.prof,'N') ...
    .* abs(median(diff(ne.Nx.depths))) .* median(diff(ne.Nx.date)*60*60*24);
%ne.Nx.xprof(:,btmix) = nan;
ne.Nx.xprof(:,btmix) = repmat(ne.Nx.xprof(:,btmix(1)-1),[1,numel(btmix)]);

ne.Nx.cprof = cumsum(ne.Nx.xprof,'omitnan');
%ne.Nx.cprof(:,btmix) = nan;
ne.Nx.cprof(:,btmix) = repmat(ne.Nx.cprof(:,btmix(1)-1),[1,numel(btmix)]);

% firstix = find(ne.N.prof(:,end-2)>0,1);
% [dterr,ix0] = min(abs(ne.Nx.date-ne.N.date(firstix-1)));
% [dterr,ix1] = min(abs(ne.Nx.date-ne.N.date(firstix)));
% [dterr,ix2] = min(abs(ne.Nx.date-ne.N.date(firstix+1)));

fmg;
surf(ne.Nx.date,ne.Nx.depths,ne.Nx.cprof');
view(0,270); shading interp; set(gca,'YDir','reverse');
hold on;
%caxis([0,15]);
colorbar;
tighten_axes;
ylim([-sfomc.ne_buoy.depth,0]); 
for dp=ne.Nx.depths(:)'; annotline([],dp); end;
xlim(xl); datetick3;
ylabel('NE (Deep) u \bullet [NO2+NO3] [\mumol\bulletm\bullets^-^1]');

%fmg; plot(ne.Nx.date,sum(cumsum(ne.Nx.xprof,'omitnan'),2)); datetick3;
fmg; plot(ne.Nx.date,sum(ne.Nx.cprof,2,'omitnan')); datetick3;
ylabel('\Sigma u \bullet [NO2+NO3] [\mumol\bulletm\bullets^-^1]');









function fh = fmg(fh)
%function fh = fmg([fh])
%
% Create new full-screen Figure with HOLD ON, GRID ON, default FONTSIZE 12.
% If optional FH is passed in, pass focus to and do the above to it. From
% 2016 Nov 28 to 2018 Feb 09, also called TIGHTEN_AXES (v.) on the resulting
% figure; BUT For FIGUREs with right-labeled COLORBARs and a 'FontSize' of 20
% (as set below), this causes figure boundary to clip the colorbar label! :(
%
% Last Saved Time-stamp: <Fri 2018-02-09 14:50:57 Eastern Standard Time gramer>

  if ( exist('fh','var') )
    fh = figure(fh);
  else
    fh = figure;
  end;
  maximize_graph;
  hold on;
  grid on;
  set(gca,'FontSize',20);  % For publication-ready fonts when printing

return;









function fh = fmg(fh)
%function fh = fmg([fh])
%
% Create new full-screen Figure with HOLD ON, GRID ON, default FONTSIZE 12.
% If optional FH is passed in, pass focus to and do the above to it. From
% 2016 Nov 28 to 2018 Feb 09, also called TIGHTEN_AXES (v.) on the resulting
% figure; BUT For FIGUREs with right-labeled COLORBARs and a 'FontSize' of 20
% (as set below), this causes figure boundary to clip the colorbar label! :(
%
% Last Saved Time-stamp: <Fri 2018-02-09 14:48:25 Eastern Standard Time gramer>

  if ( exist('fh','var') )
    fh = figure(fh);
  else
    fh = figure;
  end;
  maximize_graph;
  hold on;
  grid on;
  %set(gca,'FontSize',12);  % For publication-ready fonts when printing
  %set(gca,'FontSize',16);  % For publication-ready fonts when printing
  set(gca,'FontSize',20);  % For publication-ready fonts when printing
  %set(gca,'FontSize',24);  % For publication-ready fonts when printing

  % % For FIGUREs with right-labeled COLORBARs, this clips the label :(
  % tighten_axes(fh);

return;





function maximize_graph(fh)
%function maximize_graph(fh)
%
% Enlarge plot figure FH to fill available screen. If FH not specified,
% enlarge the GCF (see).
% 
% Last Saved Time-stamp: <Fri 2018-02-09 14:24:13 Eastern Standard Time gramer>

  if ( ~exist('fh','var') || isempty(fh) )
    fh = gcf;
  end;
  if ( exist('isOctave') && isOctave )
    set(fh, 'units','normalized');
    set(fh, 'position',[0 0 0.9 0.9]);
  else
    % % %%%set(fh, 'units','normalized', 'outerposition',[0 0 1 1], 'position',[0 0 1 1]);
    % % set(fh, 'units','normalized', 'outerposition',[0 0 0.99 0.94], 'position',[0 0 0.97 0.92]);
    % % %set(fh, 'units','normalized', 'Position',[0 0 1 0.92]);
    % % %%set(fh, 'units','normalized', 'outerposition',[0 0 0.95 0.90], 'position',[0 0 0.97 0.92]);

    % %set(fh, 'units','normalized', 'outerposition',[0 0 0.95 0.94], 'position',[0 0 0.93 0.92]);
    % set(fh, 'units','normalized', 'outerposition',[0 0 0.95 0.94], 'innerposition',[0 0 0.92 0.91]);

    set(fh, 'units','normalized');
    % Magic Position found by actually maximizing a Figure on laptop MANANNAN's screen
    set(fh, 'position',[0.00 0.00 1.00 0.91]);
  end;

return;








ne.Nx.xprof(:,btmix) = ne.Nx.xprof(:,btmix(1)-1);





ne.Nx.date = ne.adcp_x.date;
ne.Nx.depths = ne.adcp_depths(end:-1:1);
btmix = find(ne.Nx.depths<ne.N.depths(end));

ne.Nx.prof = interp2(ne.N.date,ne.N.depths,ne.N.prof',ne.Nx.date,ne.Nx.depths,'linear',NaN)';
%ne.Nx.prof(:,btmix) = nan;
ne.Nx.prof(:,btmix) = ne.Nx.prof(:,btmix(1)-1);

% Nutrient mass per profile bin per meter of reef-line per sample:
%  Cross-shore-velocity x nutrient-conc. x bin-height x time-per-sample
%   [m/s] * [kg/m^3] * [m] * [s] = [kg/m]
ne.Nx.xprof = ne.adcp_x.prof .* umol_to_kgm3(ne.Nx.prof,'N') ...
    .* abs(median(diff(ne.Nx.depths))) .* median(diff(ne.Nx.date)*60*60*24);
ne.Nx.xprof(:,ne.Nx.depths<ne.N.depths(end)) = nan;
%ne.Nx.cprof = cumsum(ne.Nx.xprof,2,'omitnan');
ne.Nx.cprof = cumsum(ne.Nx.xprof,'omitnan');
ne.Nx.cprof(:,ne.Nx.depths<ne.N.depths(end)) = nan;

% firstix = find(ne.N.prof(:,end-2)>0,1);
% [dterr,ix0] = min(abs(ne.Nx.date-ne.N.date(firstix-1)));
% [dterr,ix1] = min(abs(ne.Nx.date-ne.N.date(firstix)));
% [dterr,ix2] = min(abs(ne.Nx.date-ne.N.date(firstix+1)));

fmg;
surf(ne.Nx.date,ne.Nx.depths,ne.Nx.cprof');
view(0,270); shading interp; set(gca,'YDir','reverse');
hold on;
%caxis([0,15]);
colorbar;
ylim([-sfomc.ne_buoy.depth,0]); 
for dp=ne.Nx.depths(:)'; annotline([],dp); end;
xlim(xl); datetick3;
ylabel('NE (Deep) u \bullet [NO2+NO3] [\mumol\bulletm\bullets^-^1]');



%fmg; plot(ne.Nx.date,sum(cumsum(ne.Nx.xprof,'omitnan'),2)); datetick3;
fmg; plot(ne.Nx.date,sum(ne.Nx.cprof,2,'omitnan')); datetick3;
ylabel('\Sigma u \bullet [NO2+NO3] [\mumol\bulletm\bullets^-^1]');








function maximize_graph(fh)
%function maximize_graph(fh)
%
% Enlarge plot figure FH to fill available screen. If FH not specified,
% enlarge the GCF (see).
% 
% Last Saved Time-stamp: <Fri 2018-02-09 14:07:12 Eastern Standard Time gramer>

  if ( ~exist('fh','var') || isempty(fh) )
    fh = gcf;
  end;
  if ( exist('isOctave') && isOctave )
    set(fh, 'units','normalized');
    set(fh, 'position',[0 0 0.9 0.9]);
  else
    %%%set(fh, 'units','normalized', 'outerposition',[0 0 1 1], 'position',[0 0 1 1]);
    set(fh, 'units','normalized', 'outerposition',[0 0 0.99 0.94], 'position',[0 0 0.97 0.92]);
    %set(fh, 'units','normalized', 'Position',[0 0 1 0.92]);
    %%set(fh, 'units','normalized', 'outerposition',[0 0 0.95 0.90], 'position',[0 0 0.97 0.92]);
  end;

return;







stnm='SRVI2';
stn = load_portal_data(stnm,'CWD');
qcfname = [upper(stnm),'_portal_ALL_qc.mat'];
disp(['Loading ',qcfname]);
qc = load(qcfname);

% % qc.stn = qcstation(qc.stn,{'air_t_2_degc','salinity_deep_psu','salinity_shallow_psu',});
% % qc.stn = qcstation(qc.stn,{'air_t_2_degc','salinity_deep_psu','salinity_shallow_psu',});
% qc.stn = qcstation(qc.stn,'air_t_2_degc');
% fmg; plot_ts(qc.stn.air_t_2_degc);
% qc.stn = qcstation(qc.stn,'air_t_2_degc');
% fmg; plot_ts(qc.stn.air_t_2_degc,qc.stn.seatemp_shallow_degc);

qc.stn = qcstation(qc.stn,{'speed_1_kt',{'salinity_deep_psu','salinity_shallow_psu',},});
inp = input('ENTER to save and move to next field, "s"kip to next (ABANDONS edits), or "q"uit: ','s');
save('SRVI2_portal_ALL_qc.mat','-struct','qc','-v7.3');








fmg; surf(N2.date,N2.depths',N2.prof); colorbar; shading interp; datetick3;




[Ts,Zs] = meshgrid(spd.date,sfomc.(stnm).adcp_depths);
spd.prof = interp2(Ts,Zs,spd.prof,spd.date,spd.depths);
spds = interp2(Ts,Zs,spd.prof',spd.date,spd.depths','spline');



[N2,Pv] = gsw_Nsquared(S:-1:1,:)',T(end:-1:1,:)',z(end:-1:1)',sfomc.(stnm).lat);



% %t0 = datenum(1997,6,1); tN = datenum(1997,11,1);
% t0 = datenum(1997,7,1); tN = datenum(1997,10,1);
t0 = datenum(1997,7,1); tN = datenum(1997,10,1);
%%fld = 'simple_ndbc_erai_erai_30a_net_flux_term';
%fld = 'ndbc_erai_erai_30a_net_flux_term';
fld = 'ndbc_erai_erai_30a_avhrr_hc_dTdt';
lon = date_range_ts(hb.lonf1.(fld),[t0,tN]);
mlr = date_range_ts(hb.mlrf1.(fld),[t0,tN]);
[ig,lont0ix] = min(abs(hb.lonf1.ndbc_sea_t.date-t0));
lont0 = hb.lonf1.ndbc_sea_t.data(lont0ix);
[ig,mlrt0ix] = min(abs(hb.mlrf1.ndbc_sea_t.date-t0));
mlrt0 = hb.mlrf1.ndbc_sea_t.data(mlrt0ix);

fmg;
plot(lon.date,lont0+cumsum(lon.data),mlr.date,mlrt0+cumsum(mlr.data)); datetick3;
axis(axis);
plot_ts(hb.lonf1.ndbc_sea_t,'k:',hb.mlrf1.ndbc_sea_t,'k-');
legend('LONF1','MLRF1');
titlename('Temperature change');
if doPrint; print('-dpng',fullfile(figspath,[basenm,'-lonf1-mlrf1-simple_ndbc_erai_erai_30a_net_flux_1_d_sum.png'])); end;






disp('2014 is very like to 1997');
%{
fmg; plot_ts(lonf1.ndbc_sea_t,mlrf1.ndbc_sea_t,fwyf1.ndbc_sea_t);
legend('LONF1','MLRF1','FWYF1');
axis([datenum(1996,5,1),datenum(2017,11,1),27,34]); datetick3;
titlename('Hourly average sea temperature');

fmg; plot_ts(lonf1.ndbc_sea_t_24_h_avg,mlrf1.ndbc_sea_t_24_h_avg,fwyf1.ndbc_sea_t_24_h_avg);
legend('LONF1','MLRF1','FWYF1');
axis([datenum(1996,5,1),datenum(2017,11,1),27,34]); datetick3;
titlename('24-hour average sea temperature');

fmg; plot_ts(lonf1.ndbc_sea_t_7_d_avg,mlrf1.ndbc_sea_t_7_d_avg,fwyf1.ndbc_sea_t_7_d_avg);
legend('LONF1','MLRF1','FWYF1');
axis([datenum(1996,5,1),datenum(2017,11,1),27,34]); datetick3;
titlename('7-day average sea temperature');

fmg; plot_ts(lonf1.ndbc_sea_t,mlrf1.ndbc_sea_t,fwyf1.ndbc_sea_t);
legend('LONF1','MLRF1','FWYF1');
axis([datenum(2014,4,15),datenum(2014,12,15),27.5,33.5]); datetick3;
titlename('Hourly average sea temperature');

fmg; plot_ts(lonf1.ndbc_sea_t_7_d_avg,mlrf1.ndbc_sea_t_7_d_avg,fwyf1.ndbc_sea_t_7_d_avg);
legend('LONF1','MLRF1','FWYF1');
axis([datenum(2014,4,15),datenum(2014,12,15),27.5,33.5]); datetick3;
titlename('7-day average sea temperature');
%}




disp('2014 is very like to 1997');
%{
%}







fmg;
plot_ts(lonf1.ndbc_sea_t_7_d_avg,mlrf1.ndbc_sea_t_7_d_avg,fwyf1.ndbc_sea_t_7_d_avg);
legend('LONF1','MLRF1','FWYF1');
axis([datenum(1998,4,15),datenum(1998,12,15),27.5,33.5]); datetick3;
titlename('7-day average sea temperature');

fmg;
plot_ts(lonf1.ndbc_sea_t_7_d_avg,mlrf1.ndbc_sea_t_7_d_avg,fwyf1.ndbc_sea_t_7_d_avg);
legend('LONF1','MLRF1','FWYF1');
axis([datenum(2015,4,15),datenum(2015,12,15),27.5,33.5]); datetick3;
titlename('7-day average sea temperature');




doFigs = false; subrgn='FRT'; use_habitat_map=true; allow_hard_bottom=true; calc_spatial_dt_hc;



    matfname = fullfile(datapath,[lower(stns(ix).station_name) '_gom_hycom_',expt,'.mat']);


%bbox = [-80.5,-80.2,25.0,25.3];
bbox = [-80.9,-80.5,24.7,25.1];
%[lonix,latix] = bboxinside(sites.lons,sites.lats,bbox,true,0);
inix = bboxinside(sites.lons,sites.lats,bbox,false,0);

inix = 1:numel(sites.lons);

bath=[]; clear bath
%bath = read_hires_bathymetry_for_field({sites.lons,sites.lats},false,10e3);
%plot_hires_bathymetry(bath);
%axis(bbox);
bath = read_hires_bathymetry_for_field({sites.lons(inix),sites.lats(inix)},false,20e3);
%bath = read_hires_bathymetry_for_field({sites.lons(inix),sites.lats(inix)},true,10e3);
bath.ngdc_hires_bathy.field(bath.ngdc_hires_bathy.field>=-0.1) = nan;
plot_hires_bathymetry(bath);

for stix=1:1:numel(sites.stnms);
%for stix=inix(:)';
  stnm = sites.stnms(stix);
  u = stns.(stnm{:}).raw_nwps_ardhuin_surface_drift_u.data;
  v = stns.(stnm{:}).raw_nwps_ardhuin_surface_drift_v.data;
  lons = stns.(stnm{:}).lon + cumsum(cosd(25)*u*3600./111e3);
  lats = stns.(stnm{:}).lat + cumsum(v*3600./111e3);
  pvlh(stix) = plot(lons,lats,'r-','LineWidth',2);
end;








fmg;
%for stix=1:3:numel(sites.stnms);
for stix=1:1:numel(sites.stnms);
  stnm = sites.stnms(stix);
  stns.(stnm{:}).raw_nwps_ardhuin_surface_drift_speed.date = stns.(stnm{:}).raw_nwps_ardhuin_surface_drift_u.date;
  stns.(stnm{:}).raw_nwps_ardhuin_surface_drift_speed.data = ...
      uv_to_spd(stns.(stnm{:}).raw_nwps_ardhuin_surface_drift_u.data,...
                stns.(stnm{:}).raw_nwps_ardhuin_surface_drift_v.data);
  %lh = plot3(stns.(stnm{:}).raw_nwps_ardhuin_surface_drift_speed.date,repmat(stns.(stnm{:}).lon,[numel(stns.(stnm{:}).raw_nwps_ardhuin_surface_drift_speed.data),1]),stns.(stnm{:}).raw_nwps_ardhuin_surface_drift_speed.data,'.');
  %set(lh,'MarkerSize',10);

  stns.(stnm{:}) = verify_variable(stns.(stnm{:}),{'raw_nwps_ardhuin_surface_drift_u_3_d_sum','raw_nwps_ardhuin_surface_drift_v_3_d_sum'});
  stns.(stnm{:}).raw_nwps_ardhuin_surface_drift_speed_3_d_sum.date = stns.(stnm{:}).raw_nwps_ardhuin_surface_drift_u_3_d_sum.date;
  stns.(stnm{:}).raw_nwps_ardhuin_surface_drift_speed_3_d_sum.data = ...
      uv_to_spd(stns.(stnm{:}).raw_nwps_ardhuin_surface_drift_u_3_d_sum.data,...
                stns.(stnm{:}).raw_nwps_ardhuin_surface_drift_v_3_d_sum.data);

  ts = ts_nanify_gaps(stns.(stnm{:}).raw_nwps_ardhuin_surface_drift_speed_3_d_sum,3);
  lh = plot3(ts.date,repmat(stns.(stnm{:}).lon,[numel(ts.data),1]),ts.data.*3600/1e3,'-');
  set(lh,'LineWidth',2);
end;
%datetick3;
datetick3('x','mmmm');
view(75,66);
ylabel('Longitude'); zlabel('Kilometers');
if doPrint; print('-dpng','figs/fknms_stokes_vs_longitude_75_66.png'); end;

pause;
axis([datenum(2016,8,31),datenum(2016,10,1),-80.3,-80]); view(3);
if doPrint; print('-dpng','figs/fknms_waves_vs_longitude_UK_75_66_2016_Sep.png'); end;

view(0,0);
if doPrint; print('-dpng','figs/fknms_waves_vs_longitude_UK_0_0_2016_Sep.png'); end;







  %% 1999, 2000, 2003, 2004, 2005?, 2006, 2007, 2009?, 2010?, 2011
  %for yr=[2007]'




legend('20 Lgt','20 Top','20 Btm','40 Mid','40 Btm', 'Orientation','horizontal','Location','Best');



%lhs=plot_ts(axs(1),stns.brwd20.brwd20_top.seatemp,stns.brwd20.brwd20_btm.seatemp,stns.brwd40.brwd40_mid.seatemp,stns.brwd40.brwd40_btm.seatemp); set(lhs,'Marker','none');




%% BOCA (second code block)

[fh,lhs,axs] = multiplot_station(stns.brwd20.brwd20_lgt,{'seatemp','par'},...
                                 'Broward 20m & 40m Mooring Data',...
                                 [],{'T [^oC]','Light, Pressure [dbar]'},...
                                 [datenum(2009,10,27),datenum(2009,11,9)],...
                                 {[26,30],[0,40]},[],'k-');
hold(axs(2),'on');
plot_ts(axs(2),stns.brwd40.brwd40_prs.seapres,'b-',stns.boca40.boca40_prs.seapres,'m-');
hold(axs(1),'on');
lhs=plot_ts(axs(1),stns.brwd20.brwd20_btm.seatemp,stns.brwd40.brwd40_mid.seatemp,stns.brwd40.brwd40_btm.seatemp,stns.boca20.boca20_btm.seatemp,stns.boca40.boca40_mid.seatemp,stns.boca40.boca40_btm.seatemp); set(lhs,'Marker','none');
%datetick3('x',2,'keeplimits');
%print('-dtiff','broward_boca_20m_40m_snail.tiff');
disp('Figure not printed this time...');





lhs=plot_ts(axs(1),stns.brwd20.brwd20_top.seatemp,stns.brwd20.brwd20_btm.seatemp,stns.brwd40.brwd40_top.seatemp,stns.brwd40.brwd40_mid.seatemp,stns.brwd40.brwd40_btm.seatemp); set(lhs,'Marker','none');




%FROM read_NF09_moorings.m:

%% RESULTS OF THIS STORED IN temp.m:
%[sites,bath] = find_ngdc_slope_sites(crds,bath.ngdc_hires_bathy,7);
[sites,bath] = find_ngdc_slope_sites(crds,bath,9,{@nanmedian,9,9,10});
%}

% Range of "good" dates for ALL moorings
dtlms = [datenum(2009,10,15),datenum(2009,12,09)];

% DEPTHS/SLOPES ESTIMATED FROM 10 m "*_mhw" BATHYMETRY (code above) & DEPLOY/
% RETRIEVAL NOTES, e.g., Thermistor-Diagrams.pdf, Thermistor-placements.pdf,
% Eddy09-mooring-plan-deploy-retrieve-LGramer-logs.pdf
deps.brwd20 = 23.45;
deps.brwdbe = 29.23;
deps.brwd40 = 38.05;
deps.brwd100 = 94.82;
deps.boca20 = 28.70;
deps.bocabe = 33.68;
deps.boca40 = 41.27;

brwd20_top = deps.brwd20 - 15;
brwd40_top = deps.brwd40 - 30;
brwd100_top = deps.brwd100 - 85;
boca20_top = deps.boca20 - 15;
boca40_top = deps.boca40 - 30;

slps.brwd20 = 0.157;
slps.brwdbe = 0.061;
slps.brwd40 = 0.094;
slps.brwd100 = 0.079;
slps.boca20 = 0.017;
slps.bocabe = 0.048;
slps.boca40 = 0.038;


% Successful?,TNO (T_# thermistor serial #),SiteID,ThermID,Coords,Depth,Slope,Description ;
NF09_lines = ...
    {
% Broward_line
        1,31,'brwd20', 'brwd20_lig', brwd20([2,1]),    brwd20_top+1,0.157,'Broward 20 m mooring near-sfc. light' ; ...
        1,01,'brwd20', 'brwd20_top', brwd20([2,1]),    brwd20_top+2,0.157,'Broward 20 m mooring near surface' ; ...
        1,02,'brwd20', 'brwd20_btm', brwd20([2,1]),    brwd20_top+15-3,0.157,'Broward 20 m mooring near bottom' ; ...
        0,03,'brwdbe', 'brwdbe_01',  brwdbenth([2,1]), 23.59,0.134,'Broward benthic line 01 nearest 20 m' ; ...
        0,04,'brwdbe', 'brwdbe_02',  brwdbenth([4,3]), 27.84,0.114,'Broward benthic line 02 near 20 m' ; ...
        0,05,'brwdbe', 'brwdbe_03',  brwdbenth([6,5]), 30.21,0.028,'Broward benthic line 03 midline' ; ...
        1,06,'brwdbe', 'brwdbe_04',  brwdbenth([8,7]), 31.91,0.014,'Broward benthic line 04 nearer 40 m (logged z=30.2)' ; ...
        1,07,'brwdbe', 'brwdbe_05',  brwdbenth([10,9]),32.60,0.017,'Broward benthic line 05 nearest 40 m (logged z=30.8)' ; ...
        1,41,'brwd40', 'brwd40_pre', brwd40([2,1]),    brwd40_top+1,0.094,'Broward 40 m mooring near-sfc. pressure' ; ...
        1,08,'brwd40', 'brwd40_top', brwd40([2,1]),    brwd40_top+2,0.094,'Broward 40 m mooring near surface' ; ...
        1,09,'brwd40', 'brwd40_mid', brwd40([2,1]),    brwd40_top+15,0.094,'Broward 40 m mooring mid-line' ; ...
        1,10,'brwd40', 'brwd40_btm', brwd40([2,1]),    brwd40_top+30-3,0.094,'Broward 40 m mooring near bottom' ; ...
        1,42,'brwd100','brwd100_pre',brwd100([2,1]),   brwd100_top+1,0.079,'Broward 100 m mooring near-sfc. pressure' ; ...
        1,11,'brwd100','brwd100_sfc',brwd100([2,1]),   brwd100_top+2,0.079,'Broward 100 m mooring near surface' ; ...
        1,12,'brwd100','brwd100_upp',brwd100([2,1]),   brwd100_top+28,0.079,'Broward 100 m mooring upper' ; ...
        1,13,'brwd100','brwd100_low',brwd100([2,1]),   brwd100_top+56,0.079,'Broward 100 m mooring lower' ; ...
        1,14,'brwd100','brwd100_btm',brwd100([2,1]),   brwd100_top+85-3,0.079,'Broward 100 m mooring near bottom' ; ...
% Boca line
        1,15,'boca20','boca20_top',  boca20([2,1]),    boca20_top+2,0.017,'Boca 20 m mooring near surface' ; ...
        1,16,'boca20','boca20_btm',  boca20([2,1]),    boca20_top+15-3,0.017,'Boca 20 m mooring near bottom' ; ...
        0,17,'bocabe','bocabe_01',   bocabenth([2,1]), 39.80,0.034,'Boca benthic line 01 nearest 20 m' ; ...
        0,18,'bocabe','bocabe_02',   bocabenth([4,3]), 34.64,0.082,'Boca benthic line 02 near 20 m' ; ...
        0,19,'bocabe','bocabe_03',   bocabenth([6,5]), 33.85,0.087,'Boca benthic line 03 midline (logged z=32.5)' ; ...
        1,20,'bocabe','bocabe_04',   bocabenth([8,7]), 31.02,0.009,'Boca benthic line 04 nearer 40 m' ; ...
        0,21,'bocabe','bocabe_05',   bocabenth([10,9]),29.09,0.026,'Boca benthic line 05 nearest 40 m' ; ...
        1,43,'boca40','boca40_pre',  boca40([2,1]),    boca40_top+1,0.038,'Boca 40 m mooring near-sfc. pressure' ; ...
        1,22,'boca40','boca40_top',  boca40([2,1]),    boca40_top+2,0.038,'Boca 40 m mooring near surface' ; ...
        1,23,'boca40','boca40_mid',  boca40([2,1]),    boca40_top+15,0.038,'Boca 40 m mooring mid-line' ; ...
        1,24,'boca40','boca40_btm',  boca40([2,1]),    boca40_top+30-3,0.038,'Boca 40 m mooring near bottom' ; ...
    };










% Successful?,TNO (T_# thermistor serial #),SiteID,ThermID,Coords,Depth,Slope,Description ;
broward_line = ...
    { ...
        1,31,'brwd20', 'brwd20_lig', brwd20([2,1]),    brwd20_top+1,0.157,'Broward 20 m mooring near-sfc. light' ; ...
        1,01,'brwd20', 'brwd20_top', brwd20([2,1]),    brwd20_top+2,0.157,'Broward 20 m mooring near surface' ; ...
        1,02,'brwd20', 'brwd20_btm', brwd20([2,1]),    brwd20_top+15-3,0.157,'Broward 20 m mooring near bottom' ; ...
        0,03,'brwdbe', 'brwdbe_01',  brwdbenth([2,1]), 23.59,0.134,'Broward benthic line 01 nearest 20 m' ; ...
        0,04,'brwdbe', 'brwdbe_02',  brwdbenth([4,3]), 27.84,0.114,'Broward benthic line 02 near 20 m' ; ...
        0,05,'brwdbe', 'brwdbe_03',  brwdbenth([6,5]), 30.21,0.028,'Broward benthic line 03 midline' ; ...
        1,06,'brwdbe', 'brwdbe_04',  brwdbenth([8,7]), 31.91,0.014,'Broward benthic line 04 nearer 40 m (logged z=30.2)' ; ...
        1,07,'brwdbe', 'brwdbe_05',  brwdbenth([10,9]),32.60,0.017,'Broward benthic line 05 nearest 40 m (logged z=30.8)' ; ...
        1,41,'brwd40', 'brwd40_pre', brwd40([2,1]),    brwd40_top+1,0.094,'Broward 40 m mooring near-sfc. pressure' ; ...
        1,08,'brwd40', 'brwd40_top', brwd40([2,1]),    brwd40_top+2,0.094,'Broward 40 m mooring near surface' ; ...
        1,09,'brwd40', 'brwd40_mid', brwd40([2,1]),    brwd40_top+15,0.094,'Broward 40 m mooring mid-line' ; ...
        1,10,'brwd40', 'brwd40_btm', brwd40([2,1]),    brwd40_top+30-3,0.094,'Broward 40 m mooring near bottom' ; ...
        1,42,'brwd100','brwd100_pre',brwd100([2,1]),   brwd100_top+1,0.079,'Broward 100 m mooring near-sfc. pressure' ; ...
        1,11,'brwd100','brwd100_sfc',brwd100([2,1]),   brwd100_top+2,0.079,'Broward 100 m mooring near surface' ; ...
        1,12,'brwd100','brwd100_upp',brwd100([2,1]),   brwd100_top+28,0.079,'Broward 100 m mooring upper' ; ...
        1,13,'brwd100','brwd100_low',brwd100([2,1]),   brwd100_top+56,0.079,'Broward 100 m mooring lower' ; ...
        1,14,'brwd100','brwd100_btm',brwd100([2,1]),   brwd100_top+85-3,0.079,'Broward 100 m mooring near bottom' ; ...
    };

boca20_top = 28.70 - 15;
boca40_top = 41.27 - 30;

% TNO (T_# thermistor serial #),Successful?,SiteID,ThermID,Coords,Depth,Slope,Description ;
boca_line = ...
    { ...
        1,15,'boca20','boca20_top',  boca20([2,1]),    boca20_top+2,0.017,'Boca 20 m mooring near surface' ; ...
        1,16,'boca20','boca20_btm',  boca20([2,1]),    boca20_top+15-3,0.017,'Boca 20 m mooring near bottom' ; ...
        0,17,'bocabe','bocabe_01',   bocabenth([2,1]), 39.80,0.034,'Boca benthic line 01 nearest 20 m' ; ...
        0,18,'bocabe','bocabe_02',   bocabenth([4,3]), 34.64,0.082,'Boca benthic line 02 near 20 m' ; ...
        0,19,'bocabe','bocabe_03',   bocabenth([6,5]), 33.85,0.087,'Boca benthic line 03 midline (logged z=32.5)' ; ...
        1,20,'bocabe','bocabe_04',   bocabenth([8,7]), 31.02,0.009,'Boca benthic line 04 nearer 40 m' ; ...
        0,21,'bocabe','bocabe_05',   bocabenth([10,9]),29.09,0.026,'Boca benthic line 05 nearest 40 m' ; ...
        1,43,'boca40','boca40_pre',  boca40([2,1]),    boca40_top+1,0.038,'Boca 40 m mooring near-sfc. pressure' ; ...
        1,22,'boca40','boca40_top',  boca40([2,1]),    boca40_top+2,0.038,'Boca 40 m mooring near surface' ; ...
        1,23,'boca40','boca40_mid',  boca40([2,1]),    boca40_top+15,0.038,'Boca 40 m mooring mid-line' ; ...
        1,24,'boca40','boca40_btm',  boca40([2,1]),    boca40_top+30-3,0.038,'Boca 40 m mooring near bottom' ; ...
    };

both_lines = { broward_line(:) ; boca_line(:) };







% DEPTHS KENTUCKY-WINDED BASED ON DEPLOYMENT and RETRIEVAL NOTES

% TNO (T_# thermistor serial #), Successful?, SiteID, ThermID, Coords, Depth, Description ;
broward_line = ...
    { ...
        31, 1, 'brwd20',  'brwd20_lig',  brwd20,            5,    'Broward 20 m mooring near surface light' ; ...
        01, 1, 'brwd20',  'brwd20_top',  brwd20,            6,    'Broward 20 m mooring near surface' ; ...
        02, 1, 'brwd20',  'brwd20_btm',  brwd20,            18,   'Broward 20 m mooring near bottom' ; ...
        03, 0, 'brwdbe',  'brwdbe_01',   brwdbenth([1,2]),  nan,  'Broward benthic line 01 nearest 20 m' ; ...
        04, 0, 'brwdbe',  'brwdbe_02',   brwdbenth([3,4]),  nan,  'Broward benthic line 02 near 20 m' ; ...
        05, 0, 'brwdbe',  'brwdbe_03',   brwdbenth([5,6]),  nan,  'Broward benthic line 03 midline' ; ...
        06, 1, 'brwdbe',  'brwdbe_04',   brwdbenth([7,8]),  30.2, 'Broward benthic line 04 nearer 40 m' ; ...
        07, 1, 'brwdbe',  'brwdbe_05',   brwdbenth([9,10]), 30.8, 'Broward benthic line 05 nearest 40 m' ; ...
        41, 1, 'brwd40',  'brwd40_pre',  brwd40,            10,   'Broward 40 m mooring near surface pressure' ; ...
        08, 1, 'brwd40',  'brwd40_top',  brwd40,            11,   'Broward 40 m mooring near surface' ; ...
        09, 1, 'brwd40',  'brwd40_mid',  brwd40,            19,   'Broward 40 m mooring mid-line' ; ...
        10, 1, 'brwd40',  'brwd40_btm',  brwd40,            28,   'Broward 40 m mooring near bottom' ; ...
        42, 1, 'brwd100', 'brwd100_pre', brwd100,           15,   'Broward 100 m mooring near surface pressure' ; ...
        11, 1, 'brwd100', 'brwd100_sfc', brwd100,           16,   'Broward 100 m mooring near surface' ; ...
        12, 1, 'brwd100', 'brwd100_upp', brwd100,           43,   'Broward 100 m mooring upper' ; ...
        13, 1, 'brwd100', 'brwd100_low', brwd100,           71,   'Broward 100 m mooring lower' ; ...
        14, 1, 'brwd100', 'brwd100_btm', brwd100,           98,   'Broward 100 m mooring near bottom' ; ...
    };

% TNO (T_# thermistor serial #), Successful?, SiteID, ThermID, Coords, Depth, Description ;
boca_line = ...
    { ...
        15, 1, 'boca20',  'boca20_top',  boca20,            5,    'Boca 20 m mooring near surface' ; ...
        16, 1, 'boca20',  'boca20_btm',  boca20,            18,   'Boca 20 m mooring near bottom' ; ...
        17, 0, 'bocabe',  'bocabe_01',   bocabenth([1,2]),  nan,  'Boca benthic line 01 nearest 20 m' ; ...
        18, 0, 'bocabe',  'bocabe_02',   bocabenth([3,4]),  nan,  'Boca benthic line 02 near 20 m' ; ...
        19, 0, 'bocabe',  'bocabe_03',   bocabenth([5,6]),  32.5, 'Boca benthic line 03 midline' ; ...
        20, 1, 'bocabe',  'bocabe_04',   bocabenth([7,8]),  nan,  'Boca benthic line 04 nearer 40 m' ; ...
        21, 0, 'bocabe',  'bocabe_05',   bocabenth([9,10]), nan,  'Boca benthic line 05 nearest 40 m' ; ...
        43, 1, 'boca40',  'boca40_pre',  boca40,            10,   'Boca 40 m mooring near surface pressure' ; ...
        22, 1, 'boca40',  'boca40_top',  boca40,            11,   'Boca 40 m mooring near surface' ; ...
        23, 1, 'boca40',  'boca40_mid',  boca40,            19,   'Boca 40 m mooring mid-line' ; ...
        24, 1, 'boca40',  'boca40_btm',  boca40,            28,   'Boca 40 m mooring near bottom' ; ...
    };

both_lines = { broward_line ; boca_line };








%crds = {flip(brwd20),brwdbenth([2,1]),brwdbenth([4,3]),brwdbenth([6,5]),brwdbenth([8,7]),brwdbenth([10,9]),flip(brwd40),flip(brwd100),};

crds = {...
    [brwd20(2),brwdbenth(2),brwdbenth(4),brwdbenth(6),brwdbenth(8),brwdbenth(10),brwd40(2),brwd100(2),...
     boca20(2),bocabenth(2),bocabenth(4),bocabenth(6),bocabenth(8),bocabenth(10),boca40(2),],...
    [brwd20(1),brwdbenth(1),brwdbenth(3),brwdbenth(5),brwdbenth(7),brwdbenth(9), brwd40(1), brwd100(1),...
     boca20(1),bocabenth(1),bocabenth(3),bocabenth(5),bocabenth(7),bocabenth(9), boca40(1),],...
    {'brwd20',  'brwdbenth1', 'brwdbenth2','brwdbenth3','brwdbenth4','brwdbenth5','brwd40',  'brwd100',...
     'brwd20',  'brwdbenth1', 'brwdbenth2','brwdbenth3','brwdbenth4','brwdbenth5','brwd40', },...
       };

[sites,bath] = find_ngdc_slope_sites(crds,bath.ngdc_hires_bathy,7);




% Subset bathymetry (e.g., Miami, F010, or in this case, Palm Beach, Florida) to four sub-quadrants
coastpath = get_ecoforecasts_path('coast');
fname = 'pb_mhw';
[all_lons,all_lats,all_zs] = asc2mat(fname);

lat = []; lon = []; dat = []; clear lat lon dat

latquadir = 'SN';
lonquadir = 'WE';
latix = [1,ceil(numel(all_lats)/2),floor(numel(all_lats)/2),numel(all_lats)];
lonix = [1,ceil(numel(all_lons)/2),floor(numel(all_lons)/2),numel(all_lons)];
for dlat=1:2
  lat = all_lats(latix((dlat*2)-1):latix((dlat*2)));
  for dlon=1:2
    lon = all_lons(lonix((dlon*2)-1):lonix((dlon*2)));
    dat = all_zs(latix((dlat*2)-1):latix((dlat*2)), lonix((dlon*2)-1):lonix((dlon*2)));
    
    quadir = [latquadir(dlat),lonquadir(dlon)];
    matfbasename = sprintf('%s_%s.mat',fname,quadir);
    matfname = fullfile(coastpath,matfbasename);
    disp(['Saving ',matfname]);
    %save(matfname,'lat','lon','dat','fname','quadir','-v7.3');
    save(matfname,'lat','lon','dat','fname','quadir');
    disp({matfbasename,min(lon),max(lon),min(lat),max(lat)});
    lon = []; dat = []; clear lon dat
  end;
  lat = []; clear lat
end;





% DEPTHS Kentucky-winded

broward_line = ...
    { ...
        31, 1, 'brwd20',  'brwd20_lig',  brwd20,          5,    'Broward 20 m mooring near surface light' ; ...
        01, 1, 'brwd20',  'brwd20_top',  brwd20,          6,    'Broward 20 m mooring near surface' ; ...
        02, 1, 'brwd20',  'brwd20_btm',  brwd20,          18,   'Broward 20 m mooring near bottom' ; ...
        03, 0, 'brwdbe',  'brwdbe_01',   brwdbenth[1,2],  nan,  'Broward benthic line 01 nearest 20 m' ; ...
        04, 0, 'brwdbe',  'brwdbe_02',   brwdbenth[3,4],  nan,  'Broward benthic line 02 near 20 m' ; ...
        05, 0, 'brwdbe',  'brwdbe_03',   brwdbenth[5,6],  nan,  'Broward benthic line 03 midline' ; ...
        06, 1, 'brwdbe',  'brwdbe_04',   brwdbenth[7,8],  30.2, 'Broward benthic line 04 nearer 40 m' ; ...
        07, 1, 'brwdbe',  'brwdbe_05',   brwdbenth[9,10], 30.8, 'Broward benthic line 05 nearest 40 m' ; ...
        41, 1, 'brwd40',  'brwd40_pre',  brwd40,          10,   'Broward 40 m mooring near surface pressure' ; ...
        08, 1, 'brwd40',  'brwd40_top',  brwd40,          11,   'Broward 40 m mooring near surface' ; ...
        09, 1, 'brwd40',  'brwd40_mid',  brwd40,          19,   'Broward 40 m mooring mid-line' ; ...
        10, 1, 'brwd40',  'brwd40_btm',  brwd40,          28,   'Broward 40 m mooring near bottom' ; ...
        42, 1, 'brwd100', 'brwd100_pre', brwd100,         15,   'Broward 100 m mooring near surface pressure' ; ...
        11, 1, 'brwd100', 'brwd100_sfc', brwd100,         16,   'Broward 100 m mooring near surface' ; ...
        12, 1, 'brwd100', 'brwd100_upp', brwd100,         43,   'Broward 100 m mooring upper' ; ...
        13, 1, 'brwd100', 'brwd100_low', brwd100,         71,   'Broward 100 m mooring lower' ; ...
        14, 1, 'brwd100', 'brwd100_btm', brwd100,         98,   'Broward 100 m mooring near bottom' ; ...
    };

% tno (T_# thermistor serial #), Successful?, SiteID, ThermID, Depth, Description ;
boca_line = ...
    { ...
        15, 1, 'boca20',  'boca20_top',  boca20,          5,    'Boca 20 m mooring near surface' ; ...
        16, 1, 'boca20',  'boca20_btm',  boca20,          18,   'Boca 20 m mooring near bottom' ; ...
        17, 0, 'bocabe',  'bocabe_01',   bocabenth[1,2],  nan,  'Boca benthic line 01 nearest 20 m' ; ...
        18, 0, 'bocabe',  'bocabe_02',   bocabenth[3,4],  nan,  'Boca benthic line 02 near 20 m' ; ...
        19, 0, 'bocabe',  'bocabe_03',   bocabenth[5,6],  32.5, 'Boca benthic line 03 midline' ; ...
        20, 1, 'bocabe',  'bocabe_04',   bocabenth[7,8],  nan,  'Boca benthic line 04 nearer 40 m' ; ...
        21, 0, 'bocabe',  'bocabe_05',   bocabenth[9,10], nan,  'Boca benthic line 05 nearest 40 m' ; ...
        43, 1, 'boca40',  'boca40_pre',  boca40,          10,   'Boca 40 m mooring near surface pressure' ; ...
        22, 1, 'boca40',  'boca40_top',  boca40,          11,   'Boca 40 m mooring near surface' ; ...
        23, 1, 'boca40',  'boca40_mid',  boca40,          19,   'Boca 40 m mooring mid-line' ; ...
        24, 1, 'boca40',  'boca40_btm',  boca40,          28,   'Boca 40 m mooring near bottom' ; ...
    };




function sites = copy_T_sites(stns,sites,lbl)
%function sites = copy_T_sites(stns,sites,lbl)

  if ( ~exist('lbl','var') )
    lbl = '';
  end;

  %for str={'face','sefcri','sfomc'};
  cstnms = fieldnames(stns);
  for cstnmix = 1:numel(cstnms);
    stnm = cstnms{cstnmix};
    if ( ~isfield(stns.(stnm),'lon') )
      continue;
    end;

    stix = numel(sites) + 1;
    sites(stix).station_name = [lbl,upper(stnm)];
    sites(stix).lon = stns.(stnm).lon;
    sites(stix).lat = stns.(stnm).lat;
    sites(stix).depth = stns.(stnm).depth;
    if ( isfield(stns.(stnm),'ngdc_depth') )
      sites(stix).ngdc_depth = stns.(stnm).ngdc_depth;
    end;
    if ( isfield(stns.(stnm),'ngdc_beta') )
      sites(stix).ngdc_beta = stns.(stnm).ngdc_beta;
    elseif ( isfield(stns.(stnm),'slope') )
      sites(stix).ngdc_beta = stns.(stnm).slope;
    else
      keyboard;
    end;
    if ( isfield(stns.(stnm),'ngdc_beta_ang') )
      sites(stix).ngdc_beta_ang = stns.(stnm).ngdc_beta_ang;
    elseif ( isfield(stns.(stnm),'slope_orientation') )
      sites(stix).ngdc_beta_ang = stns.(stnm).slope_orientation;
    else
      keyboard;
    end;
    if ( isfield(stns.(stnm),'seatemp') )
      sites(stix).date = stns.(stnm).seatemp.date;
      sites(stix).data = stns.(stnm).seatemp.data;
      if ( isfield(stns.(stnm).seatemp,'prof') )
        sites(stix).prof = stns.(stnm).seatemp.prof;
      end;
    elseif ( isfield(stns.(stnm),'t') )
      sites(stix).date = stns.(stnm).t.date;
      sites(stix).data = stns.(stnm).t.data;
    elseif ( isfield(stns.(stnm),'adcp_seatemp') )
      sites(stix).date = stns.(stnm).adcp_seatemp.date;
      sites(stix).data = stns.(stnm).adcp_seatemp.data;
    else
      keyboard;
    end;
  end;

return;







FROM anfknms.m:
  if ( ~exist('sites','var') )
    % Get depth and slope according to bathymetric data
    % Default F010 bathymetry resolution for south Florida is 30 m: 3 pts. ~ 90 m
    sites = find_ngdc_slope_sites(stns,[],3); 
    sites7 = find_ngdc_slope_sites(stns,[],7); 
keyboard;



  tic,
  toc,

  tic,
  toc,
keyboard;
  %save(matfname);





%[X,T] = intersect_tses(sfomc.nw_w_btm.adcp_x_btm,sfomc.nw_w_btm.adcp_seatemp); Tprof=T.data;




ts = sfomc.ne_buoy.seatemp;
[X,T] = intersect_tses(sfomc.ne_buoy.adcp_x_btm,sfomc.ne_buoy.seatemp);

dt = median(diff(ts.date))*3600*24;
dN = dt*temp_to_nutes(ts.prof(:,end));

%fmg; plot(ts.date,cumsum(dN)); datetick3;





  if ( ~exist('T','var') || ~isnumeric(T) || ...
       ~exist('u','var') || ~isnumeric(u) || ...
       any(size(T) ~= size(u)) )
    error('T and u must be numerical arrays of identical size!');
  end;






  if ( any(size(T) ~= size(u) )
    error('T and u must have identical size!');
  end;





function nutes = net_nute_flux(T,u,lon,lat,sp)
%function nutes = net_nute_flux(T,u,lon,lat[,sp])

  nutes = [];

  if ( ~exist('T','var') || ~isnumeric(T) || ...
       ~exist('u','var') || ~isnumeric(u) || ...
       any(size(T) ~= size(u)) )
    error('T and u must be numerical arrays of identical size!');
  end;
  if ( ~exist('lon','var') || ~ismatrix(lon) || ~isnumeric(lon) || ...
       ~exist('lat','var') || ~ismatrix(lat) || ~isnumeric(lat) )
    error('LON and LAT must be numerical matrices!');
  end;
  if ( any(size(T) ~= size(u) )
    error('T and u must have identical size!');
  end;
  if ( ~exist('sp','var') || isempty(sp) )
    sp = 'N';
  end;

  switch (upper(sp)),
   case 'N',	a = 1.620;	b = 0.0;	eta = 2.420;
   case 'P',	a = 0.073;	b = 0.0;	eta = 0.132;
   case 'Si',	a = 0.370;	b = 0.0;	eta = 0.790;
   otherwise,	error('Unrecognized chemical species "%s"',sp);
  end;


  Tub = 23;
  Tlb = 12;

  T(T<Tlb) = Tlb;
  umol_per_T = (a.*(Tub-T) + b);
  umol_per_T(nute_per_T < 0) = 0;

  kgm3 = umol_to_kgm3(umol,sp);

  nutes = sum(u.*kgm3);

return;





function kgm3dy = umol_u_to_kgm3dy(umol,u,varargin)
%function kgm3dy = umol_u_to_kgm3dy(umol,u[,sp])

  kgm3dy = [];

  % [m/dy] = [m/s] * s/h * h/dy
  mdy = u.*3600.*24;

  % [kg/m^3/dy] = 
  kgm3dy = mdy.*umol_to_kgm3(umol,sp,varargin{:});

return;






  if ( is_ts(T) || is_ts(u) )
    if ( ~is_valid_ts(T) )
      error('T is not a valid time series!');
    end;
    if ( ~is_valid_ts(u) )
      error('u is not a valid time series!');
    end;
    [uTS,TTS] = intersect_tses(u,T);
    u = uTS.data;
    T = TTS.data;
  end;




% goodix = find( ~isnan(aExtC{1}) & aExtC{1}>0 & ~isnan(aExtC{2}) & aExtC{2}>0 & ...
%                ~isnan(aExtC{3}) & aExtC{3}>0 & ~isnan(aExtC{4}) & aExtC{4}>0 & ...
%                ~isnan(ah) & ah>0 & ah<20 & ~isnan(abet) & abet>0 );
goodix = find( ~isnan(ah) & ah>0 & ah<20 & ~isnan(abet) & abet>0 );

[sortbet,sortbetix] = sort(abet(goodix)); 
sortbetix = goodix(sortbetix);
scenix{1} = sortbetix( ~isnan(aExtC{1}(sortbetix)) & aExtC{1}(sortbetix)>0 );
scenix{2} = sortbetix( ~isnan(aExtC{2}(sortbetix)) & aExtC{2}(sortbetix)>0 );
scenix{3} = sortbetix( ~isnan(aExtC{3}(sortbetix)) & aExtC{3}(sortbetix)>0 );
scenix{4} = sortbetix( ~isnan(aExtC{4}(sortbetix)) & aExtC{4}(sortbetix)>0 );

fmg;
plot(abet(scenix{1}),aExtC{1}(scenix{1}),'mo-',abet(scenix{1}),aT{1}(scenix{1}),'ms:');
plot(abet(scenix{2}),aExtC{2}(scenix{2}),'bo-',abet(scenix{2}),aT{2}(scenix{2}),'bs:');
plot(abet(scenix{3}),aExtC{3}(scenix{3}),'ro-',abet(scenix{3}),aT{3}(scenix{3}),'rs:');
plot(abet(scenix{4}),aExtC{4}(scenix{4}),'ko-',abet(scenix{4}),aT{4}(scenix{4}),'ks:');
titlename('Seafloor slope vs. sea temperature extremes');


[sorth,sorthix] = sort(ah(goodix)); 
sorthix = goodix(sorthix);
scenix{1} = sorthix( ~isnan(aExtC{1}(sorthix)) & aExtC{1}(sorthix)>0 );
scenix{2} = sorthix( ~isnan(aExtC{2}(sorthix)) & aExtC{2}(sorthix)>0 );
scenix{3} = sorthix( ~isnan(aExtC{3}(sorthix)) & aExtC{3}(sorthix)>0 );
scenix{4} = sorthix( ~isnan(aExtC{4}(sorthix)) & aExtC{4}(sorthix)>0 );

fmg;
plot(sorth(scenix{1}),aExtC{1}(scenix{1}),'mo-',sorth(scenix{1}),aT{1}(scenix{1}),'ms:');
plot(sorth(scenix{2}),aExtC{2}(scenix{2}),'bo-',sorth(scenix{2}),aT{2}(scenix{2}),'bs:');
plot(sorth(scenix{3}),aExtC{3}(scenix{3}),'ro-',sorth(scenix{3}),aT{3}(scenix{3}),'rs:');
plot(sorth(scenix{4}),aExtC{4}(scenix{4}),'ko-',sorth(scenix{4}),aT{4}(scenix{4}),'ks:');
titlename('Seafloor depth vs. sea temperature extremes');








% goodix = find(~isnan(aExtC{1}) & ~isnan(aExtC{2}) & ~isnan(aExtC{3}) & ...
%               ~isnan(aExtC{4}) & aExtC{1}>0 & aExtC{2}>0 & aExtC{3}>0 & ...
%               aExtC{4}>0 & ~isnan(ah) & ah>0 & ah<20 & ~isnan(abet) & abet>0 );
goodix = find(~isnan(aExtC{1}) & aExtC{1}>0 & ~isnan(aExtC{2}) & aExtC{2}>0 & ...
              ~isnan(aExtC{3}) & aExtC{3}>0 & ~isnan(aExtC{4}) & aExtC{4}>0 & ...
              ~isnan(ah) & ah>0 & ah<20 & ~isnan(abet) & abet>0 );





%[ig,sortedix] = sort(aExtC(validix)); sortedix = validix(sortedix);
%[sortbet,sortix] = sort(abet);
[sortbet,sortbetix] = sort(abet(validix)); sortbetix = validix(sortbetix);

fmg;
plot(sortbet,aExtC{1}(sortbetix),'mo-',sortbet,aT{1}(sortbetix),'ms:');
plot(sortbet,aExtC{2}(sortbetix),'bo-',sortbet,aT{2}(sortbetix),'bs:');
plot(sortbet,aExtC{3}(sortbetix),'ro-',sortbet,aT{3}(sortbetix),'rs:');
plot(sortbet,aExtC{4}(sortbetix),'ko-',sortbet,aT{4}(sortbetix),'ks:');


[sorth,sorthix] = sort(ah(validix)); sorthix = validix(sorthix);

fmg;
plot(sorth,aExtC{1}(sorthix),'mo-',sorth,aT{1}(sorthix),'ms:');
plot(sorth,aExtC{2}(sorthix),'bo-',sorth,aT{2}(sorthix),'bs:');
plot(sorth,aExtC{3}(sorthix),'ro-',sorth,aT{3}(sorthix),'rs:');
plot(sorth,aExtC{4}(sorthix),'ko-',sorth,aT{4}(sorthix),'ks:');








subrgn='FRT'; use_habitat_map=false; allow_hard_bottom=true; calc_spatial_dt_hc; scaling='SS'; comparison_plot_spatial_dt_hc

subrgn='FK'; use_habitat_map=false; allow_hard_bottom=true; calc_spatial_dt_hc; doReport=1; scenario=1; scaling='SS'; compare_spatial_dt_hc


fmg; plot(sortbet,aExtC(sortix),'.-',sortbet,aT(sortix),'r.-');




      [err,ix] = min(abs(LON(:)-stns.(stnm).lon)+abs(LAT(:)-stns.(stnm).lat));
      [jix,iix] = ind2sub(size(LON),ix);
      if ( 1<=jix && jix<=size(h,2) && 1<=iix && iix<size(h,1) )
        if ( isnan(h(iix,jix)) )
          % iix = iix-1:iix+1;
          % jix = jix-1:jix+1;
        end;
        ah(stix) = nanmean(nanmean(h(iix,jix)));
        if ( ~isnan(ah(stix)) )
          validix(end+1) = stix;
        end;
        abet(stix) = nanmean(nanmean(bet(iix,jix)));
        ahcrng(stix) = nanmean(nanmean(hcrng(iix,jix)));
        ahch(stix) = nanmean(nanmean(hch(iix,jix)));

        for scenario=1:4;







          aTW(stix) = nanmean(nanmean(TW(iix,jix)));
          aQ0(stix) = nanmean(dTdtq0_simple(scenario,iix,jix)*24);
          ahc(stix) = nanmean(dTdthc_SS(scenario,iix,jix));
          adT(stix) = nanmean(dTdt_SS(scenario,iix,jix));

          aQ0(stix) = nanmean(dTdtq0_simple(scenario,iix,jix)*24);
          ahc(stix) = nanmean(dTdthc_SS(scenario,iix,jix));
          adT(stix) = nanmean(dTdt_SS(scenario,iix,jix));
          aTS(stix) = nanmean(nanmean(TS(iix,jix)));

          aTB(stix) = nanmean(nanmean(TB(iix,jix)));
          aQ0(stix) = nanmean(dTdtq0_simple(scenario,iix,jix)*24);
          ahc(stix) = nanmean(dTdthc_SS(scenario,iix,jix));
          adT(stix) = nanmean(dTdt_SS(scenario,iix,jix));









{scenario}



T{1} = TC;	






        if ( ~isnan(ah(stix)) )
          validix(end+1) = stix;
        end;




        aQ0(stix) = nanmean(dTdtq0_simple(scenario,iix,jix)*24);
        ahc(stix) = nanmean(dTdthc_SS(scenario,iix,jix));
        adT(stix) = nanmean(dTdt_SS(scenario,iix,jix));






  if ( ~isfield(ts1,'date') || ~isfield(ts2,'date') || ~isfield(ts1,'data') || ~isfield(ts2,'data') )
    error('First and second args TS1 and TS2 must be either scalar or time series structs!');
  end;




    xlim(datenum(2015,9,[1,30])); datetick3;




    set(get(gca,'XAxis'),'Vis','off');




    xlim(gca,datenum(2014,[8,10],1)); datetick('x','keeplimits'); set(gca,'xticklabel','');


    xlim(gca,datenum(2014,[8,10],1)); datetick('x','keeplimits');








    xlim(datenum(2015,[8,10],1)); datetick3;




    fmg;
    subplot(7,1,1:5);
    plot_ts(pvgf1.ndbc_air_t,'k-',face.deep.seatemp,face.tcm2.seatemp,face.tcm1.seatemp,face.shallow.seatemp);
    %xlim(datenum(2015,[8,10],1));
    ylim([20,33]);
    legend('PVGF1','Deep','TCM2','TCM1','Shallow','Location','SouthEast');
    grid on;
    titlename('FACE HWD moorings - upwelling and cross-shore currents');

    subplot(7,1,6:7);
    plot_ts(face.deep.adcp_l,face.deep.adcp_x,face.tcm2.x,face.tcm1.x,face.shallow.adcp_x);
    %plot_ts(face.deep.adcp_btm_l,face.deep.adcp_btm_x,face.tcm2.x,face.tcm1.x,face.shallow.adcp_btm_x);
    xlim(datenum(2015,[8,10],1)); datetick3;
    %ylim([-0.35,+0.35]);
    %ylim([-0.60,+0.60]);
    ylim([-1.10,+1.10]);
    ylim([-1.10,+1.10]);
    legend('Deep ALONG','Deep cross','TCM2 cross','TCM1 cross','Shallow cross',...
           'Orientation','horiz', 'Location','SouthEast');
    grid on;

    if doPrint; print('-dpng',fullfile(figspath,'upwelling-face-2015-onshore.png')); end;







if ( doPlots )

  if 1;
    fmg;
    subplot(3,1,1); title('Far North (Martin County)');
     plot_ts(sefcri.updb.hourly_t,sefcri.mc2.hourly_t); 
     ylim([12,32]);
     legend(['Offshore (',num2str(sefcri.updb.depth),' m)'],...
            ['Near-shore (',num2str(sefcri.mc2.depth),' m)'],'Location','SouthWest');
     grid on;
    subplot(3,1,2); title('Northern (Palm Beach)');
     plot_ts(sefcri.pb2.hourly_t,sefcri.pb1.hourly_t);
     ylim([12,32]);
     legend(['Offshore (',num2str(sefcri.pb2.depth),' m)'],...
            ['Near-shore (',num2str(sefcri.pb1.depth),' m)'],'Location','SouthWest');
     grid on;
    subplot(3,1,3); title('Central (Ft. Lauderdale)');
     plot_ts(sefcri.bc3.hourly_t,sefcri.bc1.hourly_t);
     ylim([12,32]);
     legend(['Offshore (',num2str(sefcri.bc3.depth),' m)'],...
            ['Near-shore (',num2str(sefcri.bc1.depth),' m)'],'Location','SouthWest');
     grid on;
    if doPrint; print('-dpng',fullfile(figspath,'upwelling-sefcri-2010-2013.png')); end;
  end;

  if 1;
    fmg;
    %spt(2,1,1);
    subplot(2,1,1);
    shax = plotyy_ts(jack.shallow.seatemp,jack.shallow.adcp_x);
    legend('11m Ts','11m u','Location','SouthWest');
    grid on;
    titlename('FACE HWD moorings - upwelling and cross-shore currents');
    %spt(2,1,2);
    subplot(2,1,2);
    dpax = plotyy_ts(jack.deep.seatemp,jack.deep.adcp_x);
    legend('27m Ts','27m u','Location','SouthWest');
    grid on;
    if doPrint; print('-dpng',fullfile(figspath,'upwelling-jack.png')); end;
  end;

  if 1;
    fmg;
    [ax,lh(1),lh(2)] = plotyy_ts(jack.deep.adcp_x,jack.deep.seatemp);
    lh(3) = plot(ax(2),jack.shallow.seatemp.date,jack.shallow.seatemp.data,'-','Color',[0,0.5,0]);
    xlim(datenum(2015,[1,10],1)); datetick3;
    ylim(ax(1),[-0.35,+0.35]);
    ylim(ax(2),[20,32]);
    legend(lh,'27m u','27m T','11m T','Location','SouthWest');
    grid on;
    titlename('FACE HWD moorings - upwelling and cross-shore currents');
    if doPrint; print('-dpng',fullfile(figspath,'upwelling-jack-2015.png')); end;
  end;

  if 1;
    fmg;
    subplot(7,1,1:5);
    plot_ts(jack.deep.seatemp,jack.tcm2.seatemp,jack.tcm1.seatemp,jack.shallow.seatemp);
    ylim([20,32]);
    xlim(datenum(2015,[8,10],1));
    legend('Deep','TCM2','TCM1','Shallow','Location','SouthEast');
    grid on;
    titlename('FACE HWD moorings - upwelling and cross-shore currents');

    subplot(7,1,6:7);
    plot_ts(jack.deep.adcp_l,jack.deep.adcp_x,jack.tcm2.x,jack.tcm1.x,jack.shallow.adcp_x);
    %ylim([-0.35,+0.35]);
    %ylim([-0.60,+0.60]);
    ylim([-1.10,+1.10]);
    xlim(datenum(2015,[8,10],1)); datetick3;
    legend('Deep ALONG','Deep cross','TCM2 cross','TCM1 cross','Shallow cross','Location','SouthEast');
    grid on;

    if doPrint; print('-dpng',fullfile(figspath,'upwelling-jack-2015-onshore.png')); end;
  end;

  if 1;
    fmg;
    subplot(7,1,1:3);
    plot_ts(jack.deep.seatemp,jack.tcm2.seatemp,jack.tcm1.seatemp,jack.shallow.seatemp);
    legend('Deep','TCM2','TCM1','Shallow','Location','SouthEast');
    ylim([27,31]);
    ylabel('T [^oC]');
    xlim(datenum(2015,9,[18,31])); datetick3;
    set(get(gca,'XAxis'),'Vis','off');
    grid on;
    titlename('FACE HWD moorings - upwelling and cross-shore currents');

    subplot(7,1,4:5);
    plot_ts(jack.deep.adcp_l,jack.tcm2.l,jack.tcm1.l,jack.shallow.adcp_l);
    legend('Deep','TCM2','TCM1','Shallow','Location','SouthEast');
    ylim([-0.80,+0.80]);
    ylabel('v [ms^-^1]');
    xlim(datenum(2015,9,[18,31])); datetick3;
    set(get(gca,'XAxis'),'Vis','off');
    grid on;

    subplot(7,1,6:7);
    plot_ts(jack.deep.adcp_x,jack.tcm2.x,jack.tcm1.x,jack.shallow.adcp_x);
    legend('Deep','TCM2','TCM1','Shallow','Location','SouthEast');
    ylim([-0.20,+0.20]);
    ylabel('u [ms^-^1]');
    xlim(datenum(2015,9,[18,31])); datetick3;
    grid on;

    if doPrint; print('-dpng',fullfile(figspath,'upwelling-jack-2015-onshore-Sep.png')); end;
  end;

end;








    if ( ~isfield(lofld,'date') )
      error('If HIXY.date exists, then LOFLD.date must also!');
    end;




    fld = reshape(fldvec,[1,numel(hixy.lat),numel(hixy.lon)]);



        %if ( isempty(RAW) || ~isempty(strfind(RAW{1},'NO DATA')) )





      res=[]; clear res
      if ( yr>=2003 && mo>=8 )
        res.seatemp.date = datenum(RAW(goodix,1));
        tcolix = 3;
        if any([RAW{goodix,tcolix}]>50); tcolix = 4; end;
        res.seatemp.data = [RAW{goodix,tcolix}]';
      elseif ( isnumeric(RAW{firstix,2}) )
        goodix = find(cellfun(@ischar,RAW(firstix:end,1))'&~isnan([RAW{firstix:end,2}]));
        goodix = goodix + firstix - 1;
        res.seatemp.date = datenum(RAW(goodix,1)) + [RAW{goodix,2}]';
        tcolix = 3;
        if any([RAW{goodix,tcolix}]>50); tcolix = 4; end;
        res.seatemp.data = [RAW{goodix,tcolix}]';
      else
        res.seatemp.date = datenum(strcat(RAW(firstix:end,1),{' '},RAW(firstix:end,2)));
        res.seatemp.data = [RAW{firstix:end,3}]';
      end;







        goodix = find(~isnan([RAW{firstix:end,2}]));




if (yr==2003&&mo==8); keyboard; end;





      if ( isempty(RAW) || strncmpi(RAW{1},'NO DATA',length('NO DATA')) )




  for dix=27:numel(ds)



        if ( ~isfield(brwd,stnm) || ~isfield(brwd.(stnm),'seatemp') )
          brwd.(stnm).seatemp = res.seatemp;
          stnms{end+1} = stnm;
          brwd.stns(numel(stnms)).station_name = stnm;
          brwd.stns(numel(stnms)).seatemp = res.seatemp;
        else
          brwd.(stnm) = merge_station_data(brwd.(stnm),res);
          stix = find(strcmp(stnms,stnm));
          brwd.stns(stix) = merge_station_data(brwd.stns(stix),res);
        end;







1;

set_more off;

if ( ~exist('facepath','var') )
  facepath = get_coral_path('FACE');
end;

browpath = fullfile(facepath,'Partners','Broward');

matfname = fullfile(facepath,'broward_thermistors.mat');

if ( exist(matfname,'file') )
  disp(['Loading ',matfname]);
  load(matfname);

else
  disp(['Extracting Broward data from ',browpath]);

  brwd = [];
  stnms = {};
  monms = {'Jan','Feb','March','April','May','June','July','Aug','Sept','Oct','Nov','Dec'};

  ds = dir(fullfile(browpath,'*.xls'));
  for dix=1:numel(ds)
    C = textscan(ds(dix).name,'%[^-]-%f.xls');
    stnm = C{1}{:};
    yr = C{2};

    %if ( ~isempty(yr) && yr >= 2003 )
    if ( ~isempty(yr) && 1980 <= yr && yr <= 2050 )

      fname = fullfile(ds(dix).folder,ds(dix).name);
      disp(fname);
      x = importdata(fname); 
      %[N,T,RAW]=xlsread(fname,'Year');

      if ( ~isfield(x,'textdata') || ~isfield(x.textdata,'Aug') ...
           || ~isfield(x,'data') || ~isfield(x.data,'Aug') )

        disp(['INVALID FORMAT: ',fname]);
        %DEBUG:
        keyboard;

      else
        %for mo=8:9
        for mo=1:12
          monm = monms{mo};
          txt = x.textdata.(monm);
          %dat = x.data.(monm);
          if ( numel(txt)>0 && strncmpi(txt{1},'NO DATA',length('NO DATA')) )
            %DEBUG:
            disp(['Skipping ',num2str(yr),' ',monm]);
            break;
          end;

          % Format of Excel files is so hashed up, this is the best we can do
          % without lots of manual/special processing for each year and site :(
          [N,T,RAW]=xlsread(fname,monm);
          firstix = find(cellfun(@ischar,RAW(:,1)),1);
          res=[]; clear res
          %res.seatemp.date = datenum(RAW(firstix:end,5));
          if ( isnumeric(RAW{firstix,2}) )
            goodix = find(~isnan([RAW{firstix:end,2}]));
            res.seatemp.date = datenum(RAW(firstix:end,1)) + [RAW{firstix:end,2}]';
          else
            res.seatemp.date = datenum(strcat(RAW(firstix:end,1),{' '},RAW(firstix:end,2)));
          end;
          res.seatemp.data = [RAW{firstix:end,3}]';

          % Do some simple fix-ups
          badix = find(diff(res.seatemp.date)<(0.9/24))+1;
          badix = union(badix,find(isnan(res.seatemp.date)));
          if ( ~isempty(badix) )
            disp(['Skipping ',num2str(numel(badix)),' bad records: ',fname,': ',monm]);
            res.seatemp.date(badix) = [];
            res.seatemp.data(badix) = [];
          end;
          badix = find(diff(res.seatemp.date)<(0.9/24))+1;
          badix = union(badix,find(isnan(res.seatemp.date)));
          if ( ~isempty(badix) )
            disp(['Skipping ',num2str(numel(badix)),' MORE bad records: ',fname,': ',monm]);
            res.seatemp.date(badix) = [];
            res.seatemp.data(badix) = [];
          end;

          if ( ~is_valid_ts(res.seatemp) )
            disp(['INVALID TIME SERIES: ',fname,': ',monm]);
            %DEBUG:
            keyboard;
          elseif ( any(5>res.seatemp.data|res.seatemp.data>35) )
            disp(['INVALID TEMP. DATA: ',fname,': ',monm]);
            %DEBUG:
            keyboard;
          else
            if ( ~isfield(brwd,stnm) || ~isfield(brwd.(stnm),'seatemp') )
              brwd.(stnm).seatemp = res.seatemp;
              stnms{end+1} = stnm;
              brwd.stns(numel(stnms)).station_name = stnm;
              brwd.stns(numel(stnms)).seatemp = res.seatemp;
            else
              brwd.(stnm) = merge_station_data(brwd.(stnm),res);
              stix = find(strcmp(stnms,stnm));
              brwd.stns(stix) = merge_station_data(brwd.stns(stix),res);
            end;
          end; %if ( ~is_valid_ts(res.seatemp) ) else
          res=[]; clear res
        end; %for mo=1:12
      end; %if ( ~isfield(x,'textdata') || ~isfield(x.textdata,'Year')... else ...
      x=[]; clear x
    end; %if ( 1980 <= yr && yr <= 2020 )
  end; %for dix=1:numel(ds)

  disp(['Saving ',matfname]);
  save(matfname,'brwd');

end;

set_more;













  monms = get_monthnames();




1;


if ( ~exist('facepath','var') )
  facepath = get_coral_path('FACE');
end;

browpath = fullfile(facepath,'Partners','Broward');

matfname = fullfile(facepath,'broward_thermistors.mat');

if ( exist(matfname,'file') )
  disp(['Loading ',matfname]);
  load(matfname);

else
  disp(['Extracting Broward data from ',browpath]);

  brwd=[];
  monms = get_monthnames();

  ds = dir(fullfile(browpath,'*.xls'));
  for dix=1:numel(ds)
    C = textscan(ds(dix).name,'%[^-]-%f.xls');
    stnm = C{1}{:};
    yr = C{2};

    if ( ~isempty(yr) && 1980 <= yr && yr <= 2050 )

      fname = fullfile(ds(dix).folder,ds(dix).name);
      disp(fname);
      x = importdata(fname); 
      %[N,T,RAW]=xlsread(fname,'Year');

      if ( ~isfield(x,'textdata') || ~isfield(x.textdata,'Year') ...
           || ~isfield(x,'data') || ~isfield(x.data,'Year') )

        disp(['INVALID FORMAT: ',fname]);
        %DEBUG:
        keyboard;

      else
        for mo=1:12
          monm = monms(mo,:);

        res.seatemp.date = datenum(yr,1,0) + x.data.Year(:,6);
        res.seatemp.data = x.data.Year(:,3);

        % Do some simple fix-ups
        badix = find(diff(res.seatemp.date)<(0.9/24))+1;
        badix = union(badix,find(isnan(res.seatemp.date)
        if ( ~isempty(badix) )
          disp(['Skipping ',num2str(numel(badix)),' bad records: ',fname]);
          res.seatemp.date(badix) = [];
          res.seatemp.data(badix) = [];
        end;

        if ( ~is_valid_ts(res.seatemp) )
          disp(['INVALID TIME SERIES: ',fname]);
          %DEBUG:
          keyboard;
        else
          if ( ~isfield(brwd,stnm) || ~isfield(brwd.(stnm),'seatemp') )
            brwd.(stnm).seatemp = res.seatemp;
          else
            brwd.(stnm) = merge_station_data(brwd.(stnm),res);
          end;
        end; %if ( ~is_valid_ts(res.seatemp) ) else
        res=[]; clear res
      end; %if ( ~isfield(x,'textdata') || ~isfield(x.textdata,'Year')... else ...
      x=[]; clear x
    end; %if ( 1980 <= yr && yr <= 2020 )
  end; %for dix=1:numel(ds)

  disp(['NOT Saving ',matfname]);
  %save(matfname,'brwd');

end;








fmg;
tpax(1) = subplot(8,1,1);
%plot_ts(fwyf1.ndbc_air_t,'k-',nw.sbe_seatemp,nw.adcp_seatemp);
plot_ts(fwyf1.ndbc_air_t,'k-',nw.adcp_seatemp);
ylim([26.5,31.5]); colorbar; %Adding useless COLORBAR just lines up the X axes nicely
xlim(xl); datetick3;
ylabel('NW');
grid on;
titlename(['Temperature Hovmuellers ',dtstr]);

tpax(2) = subplot(8,1,2:3);
%contourf(sw.seatemp.date,sw.seatemp.depths,sw.seatemp.prof');
surf(sw.seatemp.date,sw.seatemp.depths,sw.seatemp.prof');
%view(0,90); shading interp;
view(0,270); shading interp; set(gca,'YDir','reverse');
hold on;
contour(sw.seatemp.date,sw.seatemp.depths,sw.seatemp.prof',[28.8,28.8],'LineWidth',0.5,'Color','k');
%caxis([12,32]); colorbar;
caxis([25.8,31.6]); colorbar;
ylim([-sfomc.sw_buoy.depth,0]); 
for dp=sw.seatemp.depths(:)'; annotline([],dp); end;
xlim(xl); datetick3;
ylabel('SW');

tpax(3) = subplot(8,1,4:8);
%contourf(ne.seatemp.date,ne.seatemp.depths,ne.seatemp.prof');
surf(ne.seatemp.date,ne.seatemp.depths,ne.seatemp.prof'); 
%view(0,90); shading interp;
view(0,270); shading interp; set(gca,'YDir','reverse');
hold on;
[c,h]=contour(ne.seatemp.date,ne.seatemp.depths,ne.seatemp.prof',[23,23],'LineWidth',0.5,'Color','k');
caxis([12,32]); colorbar;
ylim([-sfomc.ne_buoy.depth,0]); 
for dp=ne.seatemp.depths(:)'; annotline([],dp); end;
xlim(xl); datetick3;
ylabel('NE (Deep) Temp.');

if doPrint; print('-dpng',fullfile(upwpath,[fignm,'_temps.png'])); end;


fmg;
xsax(1) = subplot(8,1,1);
%plot_ts(fwyf1.ndbc_bulk_windstress_ls,'k-',nw.adcp_x);
%ylim([-0.1,+0.1]); colorbar; %Adding useless COLORBAR just lines up the X axes nicely
%contourf(nw.adcp_x.date,nw.adcp_depths,nw.adcp_x.prof');
surf(nw.adcp_x.date,nw.adcp_depths,nw.adcp_x.prof'); 
%view(0,90); shading interp;
view(0,270); shading interp; set(gca,'YDir','reverse');
hold on;
contour(nw.adcp_x.date,nw.adcp_depths,nw.adcp_x.prof',[0,0],'LineWidth',0.5,'Color','k');
%caxis([-0.1,+0.1]);
colorbar;
ylim([-sfomc.nw_w_btm.depth,0]); 
xlim(xl); datetick3;
ylabel('NW');
titlename(['Cross-shore current Hovmuellers ',dtstr]);

xsax(2) = subplot(8,1,2:3);
%contourf(sw.adcp_x.date,sw.adcp_depths,sw.adcp_x.prof');
surf(sw.adcp_x.date,sw.adcp_depths,sw.adcp_x.prof'); 
%view(0,90); shading interp;
view(0,270); shading interp; set(gca,'YDir','reverse');
hold on;
contour(sw.adcp_x.date,sw.adcp_depths,sw.adcp_x.prof',[0,0],'LineWidth',0.5,'Color','k');
%caxis([-0.1,+0.1]);
colorbar;
ylim([-sfomc.sw_buoy.depth,0]); 
xlim(xl); datetick3;
ylabel('SW');

xsax(3) = subplot(8,1,4:8);
%contourf(ne.adcp_x.date,ne.adcp_depths,ne.adcp_x.prof');
surf(ne.adcp_x.date,ne.adcp_depths,ne.adcp_x.prof'); 
%view(0,90); shading interp;
view(0,270); shading interp; set(gca,'YDir','reverse');
hold on;
contour(ne.adcp_x.date,ne.adcp_depths,ne.adcp_x.prof',[0,0],'LineWidth',0.5,'Color','k');
%caxis([-0.1,+0.1]);
colorbar;
ylim([-sfomc.ne_buoy.depth,0]); 
xlim(xl); datetick3;
ylabel('NE (Deep) Cross');

if doPrint; print('-dpng',fullfile(upwpath,[fignm,'_across.png'])); end;


fmg;
lsax(1) = subplot(8,1,1);
%plot_ts(fwyf1.ndbc_bulk_windstress_ls,'k-',nw.adcp_l);
%ylim([-0.1,+0.1]); colorbar; %Adding useless COLORBAR just lines up the X axes nicely
%contourf(nw.adcp_l.date,nw.adcp_depths,nw.adcp_l.prof');
surf(nw.adcp_l.date,nw.adcp_depths,nw.adcp_l.prof'); 
%view(0,90); shading interp
view(0,270); shading interp; set(gca,'YDir','reverse');
hold on;
contour(nw.adcp_l.date,nw.adcp_depths,nw.adcp_l.prof',[0,0],'LineWidth',0.5,'Color','k');
%caxis([-0.1,+0.1]);
colorbar;
ylim([-sfomc.nw_w_btm.depth,0]); 
xlim(xl); datetick3;
ylabel('NW');
titlename(['Alongshore current Hovmuellers ',dtstr]);

lsax(2) = subplot(8,1,2:3);
%contourf(sw.adcp_l.date,sw.adcp_depths,sw.adcp_l.prof');
surf(sw.adcp_l.date,sw.adcp_depths,sw.adcp_l.prof'); 
%view(0,90); shading interp;
view(0,270); shading interp; set(gca,'YDir','reverse');
hold on;
contour(sw.adcp_l.date,sw.adcp_depths,sw.adcp_l.prof',[0,0],'LineWidth',0.5,'Color','k');
%caxis([-0.1,+0.1]);
colorbar;
ylim([-sfomc.sw_buoy.depth,0]); 
xlim(xl); datetick3;
ylabel('SW');

lsax(3) = subplot(8,1,4:8);
%contourf(ne.adcp_l.date,ne.adcp_depths,ne.adcp_l.prof');
surf(ne.adcp_l.date,ne.adcp_depths,ne.adcp_l.prof');
%view(0,90); shading interp;
view(0,270); shading interp; set(gca,'YDir','reverse');
hold on;
contour(ne.adcp_l.date,ne.adcp_depths,ne.adcp_l.prof',[0,0],'LineWidth',0.5,'Color','k');
%caxis([-0.1,+0.1]);
colorbar;
ylim([-sfomc.ne_buoy.depth,0]); 
xlim(xl); datetick3;
ylabel('NE (Deep) Along');

if doPrint; print('-dpng',fullfile(upwpath,[fignm,'_along.png'])); end;






ylim([ne.seatemp.depths(end),0]); 





fmg;
tpax(1) = subplot(8,1,1);
plot_ts(nw.TP);
ylim([0,0.5]); colorbar; %Adding useless COLORBAR just lines up the X axes nicely
xlim(xl); datetick3;
ylabel('NW');
titlename(['[P] Hovmoellers ',dtstr]);

tpax(2) = subplot(8,1,2:3);
surf(sw.TP.date,sw.TP.depths,sw.TP.prof');
view(0,270); shading interp; set(gca,'YDir','reverse');
hold on;
caxis([0,5]); colorbar;
ylim([sw.seatemp.depths(end),0]); 
for dp=sw.seatemp.depths(:)'; annotline([],dp); end;
xlim(xl); datetick3;
ylabel('SW');

tpax(3) = subplot(8,1,4:8);
surf(ne.TP.date,ne.TP.depths,ne.TP.prof'); 
view(0,270); shading interp; set(gca,'YDir','reverse');
hold on;
caxis([0,15]); colorbar;
ylim([ne.seatemp.depths(end),0]); 
for dp=ne.seatemp.depths(:)'; annotline([],dp); end;
xlim(xl); datetick3;
ylabel('NE (Deep) PO_4');

if doPrint; print('-dpng',fullfile(upwpath,[fignm,'_phos.png'])); end;






mixprof = (swtemp.prof-repmat(Tendpt,[1,size(swtemp.prof,2)])) ...
             ./ repmat((Tsfcpt-Tendpt),[1,size(swtemp.prof,2)]);

swNendpt = repmat(Nendpt,[1,size(swtemp.prof,2)]);
sw.TN = swtemp;
sw.TN.prof = swNendpt .* mixprof;
sw.TN.prof(0>sw.TN.prof | swtemp.prof>28.8) = 0;
sw.TN.data = nansum(sw.TN.prof,2);

swPendpt = repmat(Pendpt,[1,size(swtemp.prof,2)]);
sw.TP = swtemp;
sw.TP.prof = swPendpt .* mixprof;
sw.TP.prof(0>sw.TP.prof | swtemp.prof>28.8) = 0;
sw.TP.data = nansum(sw.TP.prof,2);

swSendpt = repmat(Sendpt,[1,size(swtemp.prof,2)]);
sw.TS = swtemp;
sw.TS.prof = swSendpt .* mixprof;
sw.TS.prof(0>sw.TS.prof | swtemp.prof>28.8) = 0;
sw.TS.data = nansum(sw.TS.prof,2);







nw.TN = nwtemp;
nw.TN.data = Nendpt.*(nwtemp.data-Tendpt)./(Tsfcpt-Tendpt);
nw.TN.data(nw.TN.data<0 | nw.TN.data>Nendpt) = 0;





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CALCULATE NUTRIENT CONCENTRATIONS and MIXING RATES!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ne.TN = ne.seatemp;
ne.TP = ne.seatemp;
ne.TS = ne.seatemp;

ne.TN.data = (23-ne.seatemp.data).*1.62;
ne.TN.data(ne.TN.data<0) = nan;
ne.TN.prof = (23-ne.seatemp.prof).*1.62;
ne.TN.prof(ne.TN.prof<0) = nan;

% Concentrations at other sites are approximated by estimated heat mixing rates
[nenitr,netemp,swtemp,nwtemp,sfctemp] = ...
    intersect_tses(ne.TN,ne.TP,ne.TS,ne.seatemp,sw.seatemp,nw.adcp_seatemp,nw.sbe_seatemp);

Tendpt = squeeze(netemp.prof(:,end));
Tsfcpt = sfctemp.data;
Nendpt = squeeze(nenitr.prof(:,end));
Pendpt = squeeze(nephos.prof(:,end));
Sendpt = squeeze(nesili.prof(:,end));


sw.TN = swtemp;
sw.TN.data = nenitr.data.*(swtemp.data-netemp.data)./(sfctemp.data-netemp.data);
sw.TN.data(sw.TN.data<0) = nan;
sw.TN.prof = nenitr.prof.*(swtemp.prof-netemp.prof)./(sfctemp.data-netemp.prof);
sw.TN.prof(sw.TN.prof<0 | swtemp.prof>28.8) = nan;

nw.TN = nwtemp;
nw.TN.data = nenitr.prof(:,end).*(nwtemp.data-netemp.prof(:,end))./(sfctemp.data-netemp.prof(:,end));
nw.TN.data(nw.TN.data<0) = nan;








        badix = find(jday < 0 | isnan(res.seatemp.data));



        res.seatemp.date = datenum(x.textdata.Year(4:end,5));



if ( ~exist('nw','var') || any(nw.xl~=xl) )


if ( ~exist('sw','var') )




titlename('Cross-shore current Hovmuellers');

titlename('Alongshore current Hovmuellers');





%%% ORIGINAL RESULT FOR PAPER:
%%% ORIGINAL RESULT FOR PAPER:
%%% ORIGINAL RESULT FOR PAPER:
dtcrit=(get_year(dts)>=1997 & (get_season(dts)==2|get_season(dts)==3)); dtcritstr = 'Spring/Summer';
ztop = 10;
Tbtm = 12;
Ttop = 24;

Nbtm = 1;
scatter_fit(T1(Tbtm<=T1&T1<=Ttop&z>ztop&N>Nbtm&dtcrit),...
            N(Tbtm<=T1&T1<=Ttop&z>ztop&N>Nbtm&dtcrit),...
            sprintf('%s T(%g<=T<=%g)',dtcritstr,Tbtm,Ttop),...
            sprintf('N(N>%g,z<-%g)',Nbtm,ztop)),
legend('Location','NorthEast');






% scatter_fit(N(Tbtm<=T1&T1<=Ttop&z>ztop&N>Nbtm&dtcrit),...
%             T1(Tbtm<=T1&T1<=Ttop&z>ztop&N>Nbtm&dtcrit),...
%             sprintf('N(N>%g,z<-%g)',Nbtm,ztop),...
%             sprintf('%s T(%g<=T<=%g)',dtcritstr,Tbtm,Ttop)),




scatter_fit(T1(Tbtm<=T1&T1<=Ttop&z>ztop&dtcrit),...
            P(Tbtm<=T1&T1<=Ttop&z>ztop&dtcrit),...
            sprintf('%s T(%g<=T<=%g)',dtcritstr,Tbtm,Ttop),...
            sprintf('P(P>%g,z<-%g)',Pbtm,ztop)),



scatter_fit(T1(Tbtm<=T1&T1<=Ttop&z>ztop&P>Pbtm&dtcrit),...
            P(Tbtm<=T1&T1<=Ttop&z>ztop&P>Pbtm&dtcrit),...
            sprintf('%s T(%g<=T<=%g)',dtcritstr,Tbtm,Ttop),...
            sprintf('P(P>%g,z<-%g)',Pbtm,ztop)),



hold on;
contour(sfomc.sw_buoy.seatemp.date,sfomc.sw_buoy.seatemp.depths,sfomc.sw_buoy.seatemp.prof',[23,23],'LineWidth',2,'Color','k');



for dj=0:7:42;
  disp(datestr(xl(1)+dj));
  xlim(ax(1),xl+dj); datetick3(ax(1)); 
  xlim(ax(4),xl+dj); datetick3(ax(4)); 
  xlim(ax(7),xl+dj); datetick3(ax(7)); 
  pause;
end;




hold on;
contour(sfomc.sw_buoy.adcp_l.date,sfomc.sw_buoy.adcp_depths,sfomc.sw_buoy.adcp_l.prof',[-0.05,+0.05],'LineWidth',2,'Color','k');
hold on;
contour(sfomc.ne_buoy.adcp_l.date,sfomc.ne_buoy.adcp_depths,sfomc.ne_buoy.adcp_l.prof',[-0.05,+0.05],'LineWidth',2,'Color','k');






%{
fmg;
subplot(7,1,1:2);
contourf(sfomc.c_buoy.seatemp.date,sfomc.c_buoy.seatemp.depths,sfomc.c_buoy.seatemp.prof');
colorbar;
xlim(datenum(2000,[7,8],[1,31])); datetick3;

subplot(7,1,3:7);
contourf(sfomc.e_buoy.seatemp.date,sfomc.e_buoy.seatemp.depths,sfomc.e_buoy.seatemp.prof');
colorbar;
xlim(datenum(2000,[7,8],[1,31])); datetick3;
%}








xlim(datenum(1999,[7,8],[15,31])); datetick3;




  tic,
  toc,   tic,
  toc,



if strcmpi(stn.station_name,'canf1'); tic, end;
if strcmpi(stn.station_name,'canf1'); toc, end;


dpfh = fmg;
dpax(1) = subplot(7,1,1:3); title('DEEP site currents vs. sea temperatures');
%plot_ts(fwyf1.ndbc_air_t,'k-',sfomc.ne_buoy.seatemp,sfomc.sw_buoy.seatemp,sfomc.nw_w_btm.adcp_seatemp);
plot_ts(fwyf1.ndbc_air_t,'k-',sfomc.ne_buoy.at_30m.sbe_seatemp,sfomc.sw_buoy._at_15m.sbe_seatemp,sfomc.nw_w_btm.adcp_seatemp);
ylim([16.5,32.5]);
grid on; set(gca,'XTickLabel','');
legend('Air (FWY)','Deep (NE 30 m)','Mid (SW 15 m)','Shallow (NW 11 m)', legargs{:});




plot_ts(fwyf1.ndbc_air_t,'k-',sfomc.se_buoy.seatemp,sfomc.ne_buoy.seatemp,sfomc.sw_buoy.seatemp,sfomc.nw_w_btm.adcp_seatemp);
xlim(sfomc.ne_buoy.seatemp.date([1,end]));
ylim([16.5,32.5]);
grid on; set(gca,'XTickLabel','');
legend('Air (FWY)','Offshore (SE)','Deep (NE)','Mid (SW)','Shallow (NW)', legargs{:});



plot_ts(fwyf1.ndbc_air_t,'k-',sfomc.se_buoy.seatemp,sfomc.ne_buoy.seatemp,sfomc.sw_buoy.seatemp,sfomc.nw_w_btm.adcp_seatemp);




ylim([-3.0,+3.0]);



ylim([19.5,31.5]);




ekfh = fmg;
ekax(1) = subplot(7,1,1:7); ylabel('Cross-shore flux [m^2s^-^1]'); %title('Shallow cross-shore current');
plot_ts(lkwf1.ndbc_ekman_flux_volume,pvgf1.ndbc_ekman_flux_volume,rsm.rsmas_ekman_flux_volume,fwyf1.ndbc_ekman_flux_volume);
ylim([-3.0,+3.0]); xlim(sfomc.ne_buoy.seatemp.date([1,end]));
titlename('Ekman volumetric flux near SFOMC section');
grid on; center_axes;
legend('Lake Worth (2001-2016)','Port Everglades (2009-2017)','RSMAS Rooftop (2003-2011)','Fowey Rocks (1991-2016)', legargs{:});
datetick3;




ylim([-0.5,+0.5]); xlim(sfomc.ne_buoy.seatemp.date([1,end]));





      if ~isempty(lnkax); lmfn(lnkax,framesize); end;
        if ~isempty(lnkax); lmfn(lnkax,[ordmin,ordmin+(framesize*(ordmax-ordmin))]); end;





  lnkax = [];
  if ( ~isscalar(ax) )
    lnkax = ax(2:end);
    ax = ax(1);
  end;



  fh=figure;
  hold on;
  grid on;
  set(gca,'FontSize',16);  % For publication-ready fonts when printing
  tighten_axes(fh);





  endix = 0;
  for drix=1:numel(drs)
    dr = drs(drix);
    n(drix) = 0;
    ix = find(B==dr);
    if ( numel(ix)>0 )
      if ( haveSpeed )
        n(drix) = round(cumfun(sts.data(ix)));
      else
        n(drix) = round(numel(ix));
      end;
    end;
    wd(endix+1:endix+n(drix)) = repmat(dr,[1,n(drix)]);
    endix = endix + n(drix);
  end;




enddt = get_yearseason(ts.date(end)+88);





            if ( ~isempty(delix) )




keyboard;
        ttl = strvcat(ttl,'And fields:',linkedfldnms{:});





  if ( nargs == 0 )
    % DEFAULT: Quality-control all time-series fields in STN
    flds = fieldnames(stn);
  else
    while ( nargs > 0 )
      arg = args{1};
      if ( ischar(arg) )
        flds{end+1} = arg;
      elseif ( iscellstr(arg) )
        %flds(end+1,end+numel(arg)) = arg{:};
        flds(end+1) = {arg};
      elseif ( isvector(arg) && isnumeric(arg) )
        if ( numel(arg) == 1 )
          frame_size = arg;
        elseif ( numel(arg) == 2 )
          jump_initial_frame = arg;
        end;
      else
        error('ecoforecasts:qcstation:UnknownArg',...
                'Optional arg must be a CHAR, CELLSTR, or numeric vector');
      end;
      args(1) = [];
      nargs = nargs - 1;
    end; %while ( nargs > 0 )
  end;







    % Cells within cells: If caller wants linked and/or display-by fields


      if ( numel(fldcel) == 2 )
           && (ischar(fld{1}) || iscellstr(fld{1})) ...
           && (ischar(fld{2}) || iscellstr(fld{2})) )
        fld2 = fld{2};
        fld = fld{1};
      else
        warning('ecoforecasts:qcstation:UnknownArg',...
                'A cell of cells must have exactly two CHAR or CELLSTR elts.');
      end;





From QCSTATION.m:
    % Special handling for linked or display-buddy field arguments...
    if ( iscellstr(fld) )

      if ( ~isempty(regexp(fld{1},'^linked')) )
        linked_tag = fld{1};
        fld(1) = [];
        if ( numel(fld) == 1 )
          fld = fld{1};
          disp([fld,' had no fields to link']);
        else
          linkedflds = fld(2:end);
          if ( ~isempty(regexp(linked_tag,'display')) )
            dispedflds = linkedflds;
          end;
          fld = fld{1};
          disp([fld,' HAS LINKED FIELDS: ']);
          disp(linkedflds);
        end;
      else
        if ( ~isempty(regexp(fld{1},'^display')) )
          fld(1) = [];
        end;
        dispedflds = fld(2:end);
        fld = fld{1};
        disp([fld,' will also display with: ']);
        disp(dispedflds);
      end; %if ( ~isempty(regexp(fld{1},'^linked')) ) else

    end;






  if ( ~exist('xyz','var') || isempty(xyz) )
    xyz = 'x';
  end;
  if ( ~exist('axfun','var') || isempty(axfun) )
    %axfun = [];
    if ( exist('datetick3') )
      axfun = @datetick3;
    elseif ( exist('datetick2_5') )
      axfun = @datetick2_5;
    elseif ( exist('datetick2') )
      axfun = @datetick2;
    elseif ( exist('datetick') )
      axfun = @datetick;
    else
      axfun = [];
    end;
  end;
  if ( ~exist('exit_at_end','var') || isempty(exit_at_end) )
    exit_at_end = true;
  end;






  % Process PERIODS string
  n = min(length(periods),2);
  periods_substr = lower(periods(1:n));
  switch ( periods_substr ),
   case {'year'},
    grpfun = @(x)(datenum(get_year(x),1,1));
    begdt = datenum(get_year(dts(1)),1,1);
    enddt = datenum(get_year(dts(end))+1,1,1);
   case {'yearseason'},
    grpfun = @get_yearseason;
    begdt = get_yearseason(dts(1));
    enddt = get_yearseason(dts(end)+88);
   case {'yearmonth'},
    grpfun = @get_yearmonth;
    begdt = get_yearmonth(dts(1));
    enddt = datenum(get_year(dts(end)),get_month(dts(end))+1,1);
   case {'yearweek'},
    grpfun = @get_yearweek;
    begdt = get_yearweek(dts(1));
    [wk,yr] = get_week(dts(end));
    % GET_WEEK always returns WK between 1 and 52
    enddt = datenum(yr,1,1) + (wk*7);
   case {'yearday'},
    grpfun = @floor;
    begdt = floor(dts(1));
    enddt = ceil(dts(end));
   case {'yearhour'},
    grpfun = @get_yearhour;
    begdt = get_yearhour(dts(1));
    enddt = floor(dts(end)) + ((get_hour(dts(end)) + 1)/24);

   otherwise,
    error('Period "%s" not yet implemented!',periods);
  end;







  d = grpfun(dts);
  ud = unique(d);

  % Maximum period resolution 
  DT = max(diff(ud));

  ud(end+1) = grpfun(ud(end) + (DT*1.001));

  alldts = ud(1):dt:ud(end);
  if ( dts(1) < alldts(1) || alldts(end) < dts(end) )
    error('GRPFUN returned values that do not encompass our dates!');
  end;
  alld = grpfun(alldts);
  ualld = unique(alld);

  % NOTE: The last period is unreliable - always include it for now
  idx = 1:numel(d);
  for uix = 1:numel(ud)-1
    dix = find(ismember(d,ud(uix)));
    allix = find(ismember(alld,ud(uix)));
    if ( (numel(dix)/numel(allix)) >= pct )
      idx(end+1:end+numel(dix)) = dix;
    end;
  end;
  dix = find(ismember(d,ud(end)));
  idx(end+1:end+numel(dix)) = dix;






  n = min(length(periods),2);
  periods_substr = lower(periods(1:n));
  switch ( periods_substr ),
   case {'yy','yr','ye'},
    begdt = datenum(get_year(dts(1)),1,1);
    enddt = datenum(get_year(dts(end))+1,1,1);
   case {'se'},
    begdt = get_yearseason(dts(1));
    enddt = get_yearmonseason(dts(end)+93);
   case {'ym','mo'},
    begdt = datenum(get_year(dts(1)),get_month(dts(1)),1,1);
    enddt = datenum(get_year(dts(end)),get_month(dts(end))+1,1);
   case {'yw','wk','we'},
   case {'yd','dy','da'},
   otherwise,
    error('Unrecognized period string "%s"',periods);
  end;

  alldts = begdt:dt:enddt;





  % CONVENIENCE if caller specified a GRPFUN that does not RETURN DATENUMs
  switch ( char(grpfun) ),
  get_jhour_num_no_leap, get_jhour_num.m
  get_year_2digit,
  get_year.m
  get_yearday.m
  get_triad.m
  get_minute.m
  get_week.m
  get_season.m
  get_pentad.m
  get_month.m
  get_jday_no_leap.m
  get_jday.m
  get_doy.m
  get_hour.m
  get_dom.m






if ( doReport )
  [ig,sortedix] = sort(aExtC(validix)); sortedix = validix(sortedix);
  disp([scenario_desc,' (',datestr(dts(1),'mmm yyyy'),')']);
  disp(sprintf('%-15s: % 6s  % 5s  % 6s  % 6s  % 7s  % 7s  % 6s % 6s ','Station',...
               'Depth','Beta','HCdep','HCrng','Q0','HC','PredT','ObsT'));

  rmse=0; gts=0; lts=0;
  for stix=sortedix(:)'
    stnm = stnms{stix};
    tdif = aT(stix)-aExtC(stix);
    if (tdif)>1; ltGT=' >'; gts=gts+1;
    elseif (tdif)<-1; ltGT='<'; lts=lts+1; else; ltGT=''; end;
    rmse = rmse + (tdif^2);
    disp(sprintf('%-15s: % 6.1f  % 5.3f % 6.1f  % 6.0f  % 7.3f  % 7.3f  % 5.1f  % 5.1f %s',...
                 stnm,ah(stix),abet(stix),ahch(stix),ahcrng(stix),aQ0(stix),ahc(stix),aT(stix),aExtC(stix),ltGT));
  end;
  rmse = sqrt(rmse/numel(sortedix));
  disp(sprintf('%g >, %g <, %g ~ : RMSE=%g',gts,lts,numel(validix)-gts-lts,rmse));
end;






tgtdt = datenum(2003,8,14,18,00,0);
tgtdt = datenum(2006,8,14,18,00,0);



% jix   = [ojix,   mjix,   ljix ;   ojix,   mjix,   ljix]; lat = lat(jix);
% iix   = [oiix,   miix,   liix ;   oiix,   miix,   liix]; lon = lon(iix)';


%{
%}
%%%%%%%%%% HACK - DEBUGGING with BNPON and BNPMI as examples
%%%%%%%%%% HACK
%        BNPON,  BNPMI,  LONF1
[LON,LAT] = meshgrid(lon,lat); 
stn=get_station_from_station_name('bnpon'); [err,ix]=min(abs(LON(:)-stn.lon)+abs(LAT(:)-stn.lat)); [ojix,oiix]=ind2sub(size(LON),ix);

jix   = [405,    406,    198 ;    405,    406,    198]; lat = lat(jix);
iix   = [754,    750,    476 ;    754,    750,    476]; lon = lon(iix)';

h     = [5.45,   6.83,   2.01 ;   5.45,   6.83,   2.01];
bet   = [0.0123, 0.0070, 0.0006 ; 0.0123, 0.0070, 0.0006];
hcrng = [736,    276,    184 ;    555,    555,    555];
hch   = [6.3,    5.6,    1.8 ;    6.3,    5.6,    1.8];

% h,bet,hcrng,hch,
% 26+(squeeze((dTdt_SS(1,:,:)*14)+(dTdt_SS(2,:,:)*16))),
%%%%%%%%%% HACK






  bigix = find(abs(u1)>=0.001 & abs(u2)>=0.001 & abs(v1)>=0.001 & abs(v2)>=0.001);



  %Vraw = cast(nc{'V'}(1:275,1:5,:,:),'double');
  %DEBUG:  V = cast(nc{'V'}(1,:,:,:),'double');



  %Uraw = cast(nc{'U'}(1:275,1:5,:,:),'double');



tix=5; ixen=295:320; jxen=245:295;
plot_lores_coastline(LON(jxen,ixen)-360,LAT(jxen,ixen));
contourfm(LAT(jxen,ixen),LON(jxen,ixen)-360,squeeze(spd(tix,1,jxen,ixen)),[0:0.05:14]);
quiverm(LAT(jxen,ixen),LON(jxen,ixen)-360,squeeze(u(tix,1,jxen,ixen)),squeeze(v(tix,1,jxen,ixen)));
quiverm(LAT(jxen,ixen),LON(jxen,ixen)-360,squeeze(u(tix,2,jxen,ixen)),squeeze(v(tix,2,jxen,ixen)));
quiverm(LAT(jxen,ixen),LON(jxen,ixen)-360,squeeze(u(tix,3,jxen,ixen)),squeeze(v(tix,3,jxen,ixen)));
legend([num2str(dep(1)),' m'],[num2str(dep(1)),' m'],[num2str(dep(2)),' m'],[num2str(dep(3)),' m'],'Location','Best');
titlename(datestr(dts(tix)));





   % double time(time=1140);
   %   :long_name = "time";
   %   :units = "days since 0000-01-01 00:00:00";
   %   :bounds = "time_bound";
   %   :calendar = "noleap";
   %
   % >> t = cast(nc{'time'}(:),'double');
   % >> dts = datenum(0,1,1,0,0,0) + t;
   % >> datestr(dts([1,end]))
   % ans =
   % 02-Oct-2004
   % 09-Aug-2099
   % datestr(dts([1,555]))
   % ans =
   % 02-Oct-2004
   % 19-Nov-2050
   % >> datestr(dts([1,275]))
   % ans =
   % 02-Oct-2004
   % 27-Jul-2027
   % >> dep(1:5)
   % ans =
   %      5
   %     15
   %     25
   %     35
   %     45




   % >> datestr(dts([1,200]))
   % ans =
   % 02-Oct-2004
   % 28-Apr-2021





%{
%%%%%%%%%% HACK - DEBUGGING with BNPON and BNPMI as examples
%%%%%%%%%% HACK: Make depth of BNPMI the same as depth of BNPON
%DEBUG:stn = get_station_from_station_name('bnpmi'); [LON,LAT] = meshgrid(lon,lat); [err,ix]=min(abs(LON(:)-stn.lon)+abs(LAT(:)-stn.lat)); [jix,iix]=ind2sub(size(LON),ix); h(iix,jix) = 5.4472;
%%%%%%%%%% HACK - DEBUG: Limit depths, ranges, and betas for testing
jix = [406,405 ; 406,405]; iix = [750,754 ; 750,754];
lon = lon(iix); lat = lat(jix);

h = [5.44,5.44 ; 6.83,6.83 ; 5.44,5.44 ; 6.83,6.83];
% bet = [0.0070,0.0123 ; 0.0123,0.0070 ; 0.0070,0.0123 ; 0.0123,0.0070];
bet = [0.0123,0.0123 ; 0.0070,0.0070 ; 0.0123,0.0123 ; 0.0070,0.0070];
hcrng = [736,736 ; 736,736 ; 276,276 ; 276,276];
%hch = [6.3,6.3 ; 6.3,6.3 ; 6.3,6.3 ; 6.3,6.3];
hch = [6.3,5.6 ; 6.3,5.6 ; 6.3,5.6 ; 6.3,5.6];
%%%%%%%%%% HACK
%%%%%%%%%% HACK
%}

%{
%%%%%%%%%% HACK
% BNPON, BNPMI, LONF1
jix = [405, 406, 198 ; 405, 406, 198]; lat = lat(jix);
iix = [754, 750, 476 ; 754, 750, 476]; lon = lon(iix)';

h     = [5.44,   6.83,   2.01 ;   5.44,   6.83,   2.01];
bet   = [0.0123, 0.0070, 0.0006 ; 0.0123, 0.0070, 0.0006];
hcrng = [736,    276,    184 ;    555,    555,    555];
hch   = [6.3,    5.6,    1.8 ;    6.3,    5.6,    1.8];

% h,bet,hcrng,hch,
% 26+(squeeze((dTdt_SS(1,:,:)*14)+(dTdt_SS(2,:,:)*16))),
%%%%%%%%%% HACK
%}






%hcrng = [736,    276,    184 ;     278,    278,    278];


% hch = [6.3,6.3 ; 5.6,5.6 ; 5.6,5.6 ; 6.3,6.3];
hch = [6.3,6.3 ; 5.6,5.6 ; 6.3,6.3 ; 5.6,5.6];







% %hcrng = [92,786 ; 92,786];
% %hcrng = [277,786 ; 277,786];
% hcrng = [277*2,786 ; 277*2,786];



      if ( err < 2*dlat )
        [jix,iix] = ind2sub(size(LON),ix);





snapix = 1;
  %DEBUG:
  if ( qix == snapix ); dx_SS_snap = dx_SS; dTdtx_SS_snap = dTdtx_SS; end;
  %DEBUG:  if ( qix == snapix ); u_SS_snap = u_SS; dTdx_SS_snap = dTdx_SS; end;




  % Temperature change due to Q0 at furthest extent of convection
keyboard;
  dTdtx_SS = dt.*q0./(rhoCp.*(h+(bet.*dx_SS)));
  dTdtx_US = dt.*q0./(rhoCp.*(h+(bet.*dx_US)));
  dTdtx_SU = dt.*q0./(rhoCp.*(h+(bet.*dx_SU)));
  dTdtx_UU = dt.*q0./(rhoCp.*(h+(bet.*dx_UU)));





TC = winter_mean + dTC;



%%%%%%%%%% HACK
  Qv_SS = uf.*h.*(bet.^(1/3));
%%%%%%%%%% HACK






    if ( ~isnan(ah(stix)) )
    end;






          ltGT='';
          if (aT(stix)-aExtC(stix))>1; ltGT=' >'; gts=gts+1;
          elseif (aT(stix)-aExtC(stix))<-1; ltGT='<'; lts=lts+1; end;
          if ( doReport )
            disp(sprintf('%-15s: % 6.1f  % 5.3f  % 6.0f  % 6.3f  % 6.3f  % 5.1f  % 5.1f %s',...
                         stnm,ah(stix),abet(stix),ahcrng(stix),aQ0(stix),ahc(stix),aT(stix),aExtC(stix),ltGT));
          end;







    if ( plot_sites )
      bnpin = get_station_from_station_name('bnpin');
      bnpon = get_station_from_station_name('bnpon');
      bnpmi = get_station_from_station_name('bnpmi');
      plot_station_marker(bnpin,~use_habitat_map,stnclr);
      plot_station_marker(bnpon,~use_habitat_map,stnclr);
      plot_station_marker(bnpmi,~use_habitat_map,stnclr);
    end;

    if ( plot_sites )
      bnpin = get_station_from_station_name('bnpin');
      bnpon = get_station_from_station_name('bnpon');
      bnpmi = get_station_from_station_name('bnpmi');
      plot_station_marker(bnpin,~use_habitat_map,stnclr);
      plot_station_marker(bnpon,~use_habitat_map,stnclr);
      plot_station_marker(bnpmi,~use_habitat_map,stnclr);
    end;







  fh=fmg;
  cbh=colorbar;
  switch (qix),
   case 1,
    %T = winter_mean + dTW*7 + dTC*3;
    T = winter_mean + dTW*23 + dTC*7;
    % T = winter_mean + (4*dTW*5) + (2*dTC*1);
    % %T = winter_mean + dTC*4;
    surf(lon,lat,zeroZ',TC'); caxis([winter_min,winter_mean]); shading interp
   case 2,
    %T = winter_mean + dTW*10;
    T = winter_mean + dTW*30;
    surf(lon,lat,zeroZ',TW'); caxis([winter_min,winter_mean]); shading interp
   case 3,
    %T = summer_mean + dTS*10;
    T = summer_mean + dTS*30;
    surf(lon,lat,zeroZ',TS'); caxis([summer_mean,summer_max]); shading interp
    cm=hot; colormap(cm(end:-1:1,:));
   case 4,
    %T = summer_mean + dTS*7 + dTB*3;
    T = summer_mean + dTS*23 + dTB*7;
    %%T = summer_mean + dTB*6;
    surf(lon,lat,zeroZ',TB'); caxis([summer_mean,summer_max]); shading interp
    cm=hot; colormap(cm(end:-1:1,:));
  end;






  winter_min = 12;
  winter_mean = 26;
  summer_mean = 26;
  summer_max = 34;









if qix>0; dTC = squeeze(dTdt_SS(1,:,:));  dTC(maxdepth<h | h<mindepth) = nan; end;
%if qix>0; dTC = squeeze(dTdt_SU(1,:,:));  dTC(maxdepth<h | h<mindepth) = nan; end;
%if qix>0; dTC = squeeze(dTdt_UU(1,:,:));  dTC(maxdepth<h | h<mindepth) = nan; end;

if qix>1; dTW = squeeze(dTdt_SS(2,:,:));  dTW(maxdepth<h | h<mindepth) = nan; end;

if qix>2; dTS = squeeze(dTdt_SS(3,:,:));  dTS(maxdepth<h | h<mindepth) = nan; end;

if qix>3; dTB = squeeze(dTdt_SS(4,:,:));  dTB(maxdepth<h | h<mindepth) = nan; end;
%if qix>3; dTB = squeeze(dTdt_SU(4,:,:));  dTB(maxdepth<h | h<mindepth) = nan; end;








)) || isnan(T(iix,jix


  %DEBUG:  if ( qix == snapix ); dTdtq0_snap = dTdtq0; end;






subrgn='UK'; use_habitat_map=true; allow_hard_bottom=true; plot_sites=true; doPrint=false; plot_spatial_dt_hc





%{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ( ~exist('use_habitat_map','var') || isempty(use_habitat_map) )
  use_habitat_map = false;
  %use_habitat_map = true;
end;
if ( ~exist('allow_hard_bottom','var') || isempty(allow_hard_bottom) )
  allow_hard_bottom = true;
  %allow_hard_bottom = false;
end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%}
%{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ( ~exist('mindepth','var') || isempty(mindepth) )
  mindepth = 1;
end;
if ( ~exist('maxdepth','var') || isempty(maxdepth) )
  maxdepth = 30;
end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%}

%{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dt = 3600;

R = (1-0.08);

%% Check seasonality of sea-surface radiative and turbulent stress forcing

%load('d:/thesis/data/mlrf1-heat_budget-erai-avhrr_weekly-ndbc-tpxo_tide-erai_DISSERT.mat');
%OR: load('d:/thesis/data/lonf1-heat_budget-erai-avhrr_weekly-ndbc-tpxo_tide-erai.mat');
%fmg; boxplot_ts(stn.simple_ndbc_erai_erai_30a_net_flux,[],'mean',true,'std',true);
%stn = verify_variable(stn,'simple_ndbc_erai_erai_30a_net_flux_1_d_sum');
%MLRF1: unique(get_month(stn.simple_ndbc_erai_erai_30a_net_flux_1_d_sum.date(stn.simple_ndbc_erai_erai_30a_net_flux_1_d_sum.data<-23000)))
%LONF1: unique(get_month(stn.simple_ndbc_erai_erai_30a_net_flux_1_d_sum.date(stn.simple_ndbc_erai_erai_30a_net_flux_1_d_sum.data<-18500)))
%MLRF1: unique(get_year(stn.simple_ndbc_erai_erai_30a_net_flux_1_d_sum.date(get_month(stn.simple_ndbc_erai_erai_30a_net_flux_1_d_sum.date)==6&stn.simple_ndbc_erai_erai_30a_net_flux_1_d_sum.data>6200))) % 1998
%LONF1: unique(get_year(stn.simple_ndbc_erai_erai_30a_net_flux_1_d_sum.date(get_month(stn.simple_ndbc_erai_erai_30a_net_flux_1_d_sum.date)>=6&stn.simple_ndbc_erai_erai_30a_net_flux_1_d_sum.data>7500))) % 1995

% % December/July "representative" means
% q0s = [-100,100];
% q0str={'Representative mean winter','Representative mean summer'};
% % MLRF: December/July means at MLRF1
% q0s = [-115,95];
% q0str={'Reef-crest mean winter','Reef-crest mean summer'};
% % LONF1: February 25th-percentile / May 75th-percentile
% q0s = [-140,285];
% q0str={'Back-reef flats cool winter','Back-reef flats warm summer'};
% % MLRF1: December 25th-percentile / May 75th-percentile
% q0s = [-240,320];
% q0str={'Reef-crest cool winter','Reef-crest warm summer'};

% %BASED ON 1_d_sum
% % MLRF1: December lower std. dev. and mean, June mean, June maximum
% q0s = [-240,-115,95,255];
% q0str={'Cold winter','Normal winter','Normal summer','Hot summer'};

% MLRF1: January minimum and mean, June mean, June maximum
q0s = [-800,-115,95,255];
%q0s = [-550,-115,95,255];
q0str={'2010 cold snap','Normal winter','Normal summer','1998 bleaching'};
snapix = 1;

nq0s = numel(q0s);

switch (subrgn),
 case '',
  error('No sub-region code (SUBRGN) was specified!');
 case 'SE',
  subrgnstr = 'SE Florida';
  subrgnbox = [-80.20,-79.95,25.45,27.00];%27.30];
 case 'UK',
  subrgnstr = 'Upper Keys';
  %subrgnbox = [-80.45,-80.05,24.95,25.55];
  subrgnbox = [-80.60,-80.05,24.90,25.55];
 case 'MK',
  subrgnstr = 'Middle Keys';
  %subrgnbox = [-81.05,-80.35,24.55,25.05];
  subrgnbox = [-81.05,-80.30,24.55,25.15];
 case 'LK',
  subrgnstr = 'Lower Keys';
  subrgnbox = [-82.05,-80.95,24.40,24.75];
 case 'DT',
  subrgnstr = 'Dry Tortugas';
  subrgnbox = [-82.95,-81.95,24.35,24.75];
 otherwise,
  error('Region %s not yet implemented!',subrgn);
end;


% Assume Horizontal Convection affects sea temperatures within about
% 270 m radius: Downsample bathymetry to the appropriate resolution
if ( ~exist('hc_radius','var') || isempty(hc_radius) )
  hc_radius = 270;
  %hc_radius = 92;
end;

% % 10-m resolution bathymetry
% bath_res = 10;
% % 30-m resolution bathymetry
% bath_res = 30;
% 92-m resolution bathymetry
bath_res = 92;

downstride = ceil(hc_radius/bath_res);
if ( downstride > 1 )
  downfunc = {@nanmean,downstride,downstride,downstride+1}; 
  downfuncstr = ['_',interpmethod_to_text(downfunc)];
else
  disp('No downsampling needed');
  downfuncstr = '';
end;

if ( plot_sites && ~exist('stns','var') )
  anfknms;
end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%}






%{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
matbasename = sprintf('FRT_depth_and_beta_%dm',bath_res);
if ( use_habitat_map )
  if ( allow_hard_bottom )
    habfname = '_coral_and_hard_bottom';
  else
    habfname = '_coral';
  end;
else
  habfname = '';
end;
matbasefname = [matbasename,'_hc_range',habfname,'.mat'];
matfname = [matbasename,downfuncstr,'_hc_range',habfname,'.mat'];

if ( ~exist('lon','var') || ~exist('lat','var')  || ...
     ~exist('h','var') || ~exist('bet','var') || ~exist('hcrng','var') )

  if ( exist(matfname,'file') )
    disp(['Loading ',matfname]);
    load(matfname);

    % Somehow, the unaveraged bathymetry got processed differently before saving??
    if ( downstride == 1 )
      h = h';
      bet = bet';
      % No HC over land, and depths are positive
      h(h>0) = 0;
      h = abs(h);
    end;

  else
    %DEBUG:
    disp(['Not found: ',matfname]);    keyboard;
    %DEBUG:    if ( use_habitat_map );      keyboard;    end;

    lon=[]; lat=[]; h=[]; bet=[]; clear lon lat h bet
    if ( exist(matbasefname,'file') )
      disp(['Loading ',matbasefname]);
      load(matbasefname);

    elseif ( use_habitat_map )
      error(['USE_HABITAT_MAP: Unable to find ',matbasefname]);

    else
      disp('Extracting bathymetry');
      fld.lon = [-82.95,-79.95];
      fld.lat = [+24.30,+27.30];
      [x,rad] = read_hires_bathymetry_for_field(fld,false);
      bath = x.ngdc_hires_bathy; x=[]; clear x
      [ig,ig,ig,bath] = find_ngdc_slope(bath,[],[],3);
      lon = bath.lon;
      lat = bath.lat;
      h = bath.field;
      bet = bath.beta;
      disp(['Saving ',matbasefname]);
      save(matbasefname,'lon','lat','h','bet');
    end; %if ( exist(matbasefname,'file') ) else

    % No HC over land, and depths are positive
    h(h>0) = 0;
    h = abs(h);

    % Downsample to appropriate resolution for HC effects
    if ( downstride > 1 )
      disp(['Downsampling H and BETA to ~',num2str(hc_radius),' m']);
      [LAT,LON] = meshgrid(lat(1:downstride:end),lon(1:downstride:end));
      H = interp_field(lat,lon,h,LAT,LON,downfunc);
      H = reshape(H,size(LAT));
      BET = interp_field(lat,lon,bet,LAT,LON,downfunc);
      BET = reshape(BET,size(LAT));

      lon=[]; lat=[]; h=[]; bet=[]; clear lon lat h bet
      lon = unique(LON(:));
      lat = unique(LAT(:));
      h = H;
      bet = BET;
      LON=[]; LAT=[]; H=[]; BET=[]; clear LON LAT H BET

      disp(['Saving ',matfname]);
      save(matfname,'lon','lat','h','bet');
    end; %if ( downstride > 1 )

  end; %if ( exist(matfname,'file') ) else

end; %if ( ~exist('lon','var') || ...

lon_rmix = find(lon<subrgnbox(1) | subrgnbox(2)<lon);
lat_rmix = find(lat<subrgnbox(3) | subrgnbox(4)<lat);
if ( ~isempty(lon_rmix) || ~isempty(lat_rmix) )
  disp(['Subsetting to sub-region ',subrgn]);
  lat(lat_rmix) = [];		lon(lon_rmix) = [];
  h(:,lat_rmix) = [];		h(lon_rmix,:) = [];
  bet(:,lat_rmix) = [];		bet(lon_rmix,:) = [];
  ang(:,lat_rmix) = [];		ang(lon_rmix,:) = [];
  hcrng(:,lat_rmix) = [];	hcrng(lon_rmix,:) = [];
end;

if ( isempty(h) || isempty(bet) )
  error('Subsetting resulted in empty set! Try reloading...');
end;


if ( ~exist('rhoCp','var') )
  disp('Calculating rhoCp');
  s = repmat(35,size(h));
  t = repmat(25,size(h));
  %h = 10.96;

  g = 9.79;						%[m/s^2]
  alph = sw_alpha(s,t,h);				%[1/K]
  rho = sw_dens(s,t,h);					%[kg/m^3]
  Cp = sw_cp(s,t,h);					%[J/kg*K]
  rhoCp = rho.*Cp;					%[J/K*m^3]
  s=[]; t=[]; rho=[]; Cp=[]; clear s t Cp rho
end;

%clear dTdthc_SS dTdthc_US dTdthc_SU dTdthc_UU

[betSz1,betSz2] = size(bet);
bix1=1:betSz1;
bix2=1:betSz2;

for qix=1:nq0s;
  % Net sea-surface heat flux, Q_0
  q0 = q0s(qix);
  %DEBUG:  disp(['Q0 = ',num2str(q0)]);

  % Net sea-surface buoyancy flux B_0
  B0 = (g.*alph.*abs(q0))./(rhoCp);
  % "Characteristic" convective velocity (Sturmann)
  uf = (B0.*h).^(1/3);
  B0=[]; clear B0

  % Volumetric flow [m^2/s]: Panel letter refers to Fig. 10, Monismith et al (2006)
  % Advective inertial (steady) momentum balance, steady thermal balance: Panel (c)
  Qv_SS = uf.*h./(bet.^(1/3));
  % Viscous (unsteady) inertial balance, balanced thermal forcing: Panel (a)
  Qv_US = sqrt((uf.^3) .* (24.*dt) .* h);
  % Advective inertial balance, unbalanced thermal forcing: Panel (f)
  Qv_SU = (bet.^(2/3)).*uf.*((uf.*(24.*dt)./h).^(3/2));
  % Viscous (unsteady) inertia, unbalanced thermal forcing: Panel (d)
  Qv_UU = bet.*(uf.^3).*((24.*dt)^2)./h;
  uf=[]; clear uf

  dTdtq0 = dt.*q0./(rhoCp.*h);
  %DEBUG:  if ( qix == snapix ); dTdtq0_snap = dTdtq0; end;

  % Convective flow rate [m/s]
  u_SS = (5.0.*Qv_SS./h) - 0.05;	u_SS(u_SS<0) = 0;
  Qv_SS=[]; clear Qv_SS

  u_US = (3.0.*Qv_US./h) - 0.026;	u_US(u_US<0) = 0;
  Qv_US=[]; clear Qv_US

  u_SU = (2.7.*Qv_SU./h) - 0.0224;	u_SU(u_SU<0) = 0;
  Qv_SU=[]; clear Qv_SU

  u_UU = (0.1.*Qv_UU./h);		u_UU(u_UU<0) = 0;
  Qv_UU=[]; clear Qv_UU


  % Furthest distance traveled by convective flow in an hour [m]
  dx_SS = u_SS .* dt;
  dx_US = u_US .* dt;
  dx_SU = u_SU .* dt;
  dx_UU = u_UU .* dt;
  % REALITY CHECK: convective flow is limited in spatial extent!
  dx_SS = min(dx_SS,hcrng,'includenan');
  dx_US = min(dx_US,hcrng,'includenan');
  dx_SU = min(dx_SU,hcrng,'includenan');
  dx_UU = min(dx_UU,hcrng,'includenan');

  % Temperature change due to Q0 at furthest extent of convection
  dTdtx_SS = dt.*q0./(rhoCp.*(h+(bet.*dx_SS)));
  dTdtx_US = dt.*q0./(rhoCp.*(h+(bet.*dx_US)));
  dTdtx_SU = dt.*q0./(rhoCp.*(h+(bet.*dx_SU)));
  dTdtx_UU = dt.*q0./(rhoCp.*(h+(bet.*dx_UU)));


  % Static temperature gradient due to Q0 over depth difference between
  % observation point, and point of further extent of convection
  dTdx_SS=(dTdtq0-dTdtx_SS)./dx_SS;	dTdx_SS(~isfinite(dTdx_SS)) = 0;
  %DEBUG:
  if ( qix == snapix ); dx_SS_snap = dx_SS; dTdtx_SS_snap = dTdtx_SS; end;
  dx_SS=[]; dTdtx_SS=[]; clear dx_SS dTdtx_SS

  dTdx_US=(dTdtq0-dTdtx_US)./dx_US;	dTdx_US(~isfinite(dTdx_US)) = 0;
  dx_US=[]; dTdtx_US=[]; clear dx_US dTdtx_US

  dTdx_SU=(dTdtq0-dTdtx_SU)./dx_SU;	dTdx_SU(~isfinite(dTdx_SU)) = 0;
  dx_SU=[]; dTdtx_SU=[]; clear dx_SU dTdtx_SU

  dTdx_UU=(dTdtq0-dTdtx_UU)./dx_UU;	dTdx_UU(~isfinite(dTdx_UU)) = 0;
  dx_UU=[]; dTdtx_UU=[]; clear dx_UU dTdtx_UU


  % Rayleigh Benard instability may dampen HC during warming (Mao Lei Patterson)
  RB = repmat(1.00,size(dTdtq0));
  RB(dTdtq0 > 0) = 0.66;

  % Temperature change due to horizontal convection at observation point
  dTdthc_SS(qix,bix1,bix2) = -RB.*R.*24.*dt.*u_SS.*dTdx_SS;
  %DEBUG:  if ( qix == snapix ); u_SS_snap = u_SS; dTdx_SS_snap = dTdx_SS; end;
  u_SS=[]; dTdx_SS=[]; clear u_SS dTdx_SS

  dTdthc_US(qix,bix1,bix2) = -RB.*R.*24.*dt.*u_US.*dTdx_US;
  u_US=[]; dTdx_US=[]; clear u_US dTdx_US

  dTdthc_SU(qix,bix1,bix2) = -RB.*R.*24.*dt.*u_SU.*dTdx_SU;
  u_SU=[]; dTdx_SU=[]; clear u_SU dTdx_SU

  dTdthc_UU(qix,bix1,bix2) = -RB.*R.*24.*dt.*u_UU.*dTdx_UU;
  u_UU=[]; dTdx_UU=[]; clear u_UU dTdx_UU
  RB=[]; clear RB

  % Net temperature change (HC+Q0) at observation point
  dTdt_SS(qix,bix1,bix2) = squeeze(dTdthc_SS(qix,:,:)) + (dTdtq0.*24);
  dTdt_US(qix,bix1,bix2) = squeeze(dTdthc_US(qix,:,:)) + (dTdtq0.*24);
  dTdt_SU(qix,bix1,bix2) = squeeze(dTdthc_SU(qix,:,:)) + (dTdtq0.*24);
  dTdt_UU(qix,bix1,bix2) = squeeze(dTdthc_UU(qix,:,:)) + (dTdtq0.*24);

  % REALITY CHECK: this convective flow model cannot reverse temperatures!
  if ( q0 > 0 )
    dTdt_SS(qix,:,:) = max(0,dTdt_SS(qix,:,:),'includenan');
    dTdt_US(qix,:,:) = max(0,dTdt_US(qix,:,:),'includenan');
    dTdt_SU(qix,:,:) = max(0,dTdt_SU(qix,:,:),'includenan');
    dTdt_UU(qix,:,:) = max(0,dTdt_UU(qix,:,:),'includenan');
  else
    dTdt_SS(qix,:,:) = min(0,dTdt_SS(qix,:,:),'includenan');
    dTdt_US(qix,:,:) = min(0,dTdt_US(qix,:,:),'includenan');
    dTdt_SU(qix,:,:) = min(0,dTdt_SU(qix,:,:),'includenan');
    dTdt_UU(qix,:,:) = min(0,dTdt_UU(qix,:,:),'includenan');
  end;


  dTdtq0_simple(qix,bix1,bix2) = dTdtq0;
  dTdtq0=[]; clear dTdtq0

end; %for qix=1:nq0s;


alph=[]; rhoCp=[]; clear alph rhoCp
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%}















TC = winter_mean + dTC;

% %TW = winter_mean + dTW*10;
% TW = winter_mean + dTW*30;
TW = winter_mean + dTW;

% %TS = summer_mean + dTS*10;
% TS = summer_mean + dTS*30;
TS = summer_mean + dTS;

% %TB = summer_mean + dTS*7 + dTB*3;
% TB = summer_mean + dTS*23 + dTB*7;
% %%TB = summer_mean + dTB*6;
TB = summer_mean + dTB;







      if ( isnan(nanmean(h(iix,jix))) )
        iix = iix(1)-1:iix(end)+1;
        jix = jix(1)-1:jix(end)+1;
      end;





 case 'FK',
  subrgnstr = 'Florida Keys';
  subrgnbox = [-82.05,-80.05,24.35,25.60];
 case 'FRT',
  subrgnstr = 'Florida Reef Tract';
  subrgnbox = [-82.05,-79.95,24.35,27.00];





axis([-80.21,-80.135,25.33,25.40]); daspect([1,cosd(25),1]); view(2); colorbar('off'); colorbar('Location','EastOutside');





  % REALITY CHECK: convective flow is limited in spatial extent!
  dx_SS(dx_SS > hc_radius) = hc_radius;
  dx_US(dx_US > hc_radius) = hc_radius;
  dx_SU(dx_SU > hc_radius) = hc_radius;
  dx_UU(dx_UU > hc_radius) = hc_radius;






  % DEBUG: REALITY CHECK: convective flow cannot reverse temperatures!
  if ( q0 > 0 )
    dTdthc_SS(qix,:,:) = min(-eps,dTdthc_SS(qix,:,:));
    dTdthc_US(qix,:,:) = min(-eps,dTdthc_US(qix,:,:));
    dTdthc_SU(qix,:,:) = min(-eps,dTdthc_SU(qix,:,:));
    dTdthc_UU(qix,:,:) = min(-eps,dTdthc_UU(qix,:,:));
  else
    dTdthc_SS(qix,:,:) = max(+eps,dTdthc_SS(qix,:,:));
    dTdthc_US(qix,:,:) = max(+eps,dTdthc_US(qix,:,:));
    dTdthc_SU(qix,:,:) = max(+eps,dTdthc_SU(qix,:,:));
    dTdthc_UU(qix,:,:) = max(+eps,dTdthc_UU(qix,:,:));
  end;




  % DEBUG: REALITY CHECK: convective flow cannot reverse temperatures!
  if q0 > 0;	dTdthc_SS(qix,:,:) = min(-eps,dTdthc_SS(qix,:,:));
  else		dTdthc_SS(qix,:,:) = max(+eps,dTdthc_SS(qix,:,:)); end;
  % DEBUG: REALITY CHECK: convective flow cannot reverse temperatures!
  if q0 > 0;	dTdthc_US(qix,:,:) = min(-eps,dTdthc_US(qix,:,:));
  else		dTdthc_US(qix,:,:) = max(+eps,dTdthc_US(qix,:,:)); end;
  % DEBUG: REALITY CHECK: convective flow cannot reverse temperatures!
  if q0 > 0;	dTdthc_SU(qix,:,:) = min(-eps,dTdthc_SU(qix,:,:));
  else		dTdthc_SU(qix,:,:) = max(+eps,dTdthc_SU(qix,:,:)); end;
  % DEBUG: REALITY CHECK: convective flow cannot reverse temperatures!
  if q0 > 0;	dTdthc_UU(qix,:,:) = min(-eps,dTdthc_UU(qix,:,:));
  else		dTdthc_UU(qix,:,:) = max(+eps,dTdthc_UU(qix,:,:)); end;




  dTdt_SS(sign(dTdt_SS(qix,:) == sign(dTdtq0))) = 



if ( ~exist('FRT_depth_and_beta_92m_NANMEAN_3_3_4_hc_range_coral.mat','file') )
disp('NANMEAN - Reef - Hit Enter to continue...'); pause;
tic,
  h=[]; bet=[]; LON=[]; LAT=[]; clear lon lat h bet LON LAT
  load('FRT_depth_and_beta_92m_NANMEAN_3_3_4_hc_range.mat');
  [LON,LAT] = meshgrid(lon,lat);
  inix=[];
  disp(numel(shps));
  dlat = 1.1*abs(lat(1)-lat(2));

  w = warning('OFF','MATLAB:triangulation:EmptyTri2DWarnId'); % We don't care about these
  for shpix=1:numel(shps)
    if (mod(shpix,checkin_count) == 0); disp(shpix); toc; tic; end;
    if ( strcmp(shps(shpix).DESCRIPT,'Coral Reef') )
      [lonix,latix] = fieldfind(lon,lat,shps(shpix).X(1:end-1),shps(shpix).Y(1:end-1));
      goodix = find(~isnan(lonix) & ~isnan(latix));
      lonix = lonix(goodix);
      latix = latix(goodix);
      if ( ~isempty(lonix) && ~isempty(latix) )
        tr = delaunayTriangulation(lonix',latix');
        ix = sub2ind(size(LON),tr.Points(:,2),tr.Points(:,1));
        ix = find(inpolygon(LON,LAT,LON(ix),LAT(ix))==1);
        inix = union(inix,ix);
      end;
    end;
  end;
  newht = repmat(nan,size(LON));
  ht=h'; newht(inix)=ht(inix); ht=[]; clear ht
  h=[]; h=newht'; newht=[]; clear newht

  newbett = repmat(nan,size(LON));
  bett=bet'; newbett(inix)=bett(inix); bett=[]; clear bett
  bet=[]; bet=newbett'; newbett=[]; clear newbett

  newangt = repmat(nan,size(LON));
  angt=ang'; newangt(inix)=angt(inix); angt=[]; clear angt
  ang=[]; ang=newangt'; newangt=[]; clear newangt

  newhcrngt = repmat(nan,size(LON));
  hcrngt=hcrng'; newhcrngt(inix)=hcrngt(inix); hcrngt=[]; clear hcrngt
  hcrng=[]; hcrng=newhcrngt'; newhcrngt=[]; clear newhcrngt

  disp(numel(find(~isnan(h))));
  save('FRT_depth_and_beta_92m_NANMEAN_3_3_4_hc_range_coral.mat',...
       'lon','lat','h','bet','ang','hcrng');
  clear inix ix shpix;
toc,
end;






%DEBUG:
fmg; surf(lon,lat,repmat(0,size(h')),-h'); colorbar; shading interp; caxis([-80,0]);
axis([-80.21,-80.135,25.33,25.40]); daspect([1,cosd(25),1]); view(2); colorbar('off'); colorbar('Location','EastOutside');
fmg; surf(lon,lat,repmat(0,size(hcrng')),hcrng'); colorbar; shading interp; caxis([0 2e3]);
axis([-80.21,-80.135,25.33,25.40]); daspect([1,cosd(25),1]); view(2); colorbar('off'); colorbar('Location','EastOutside');
[cc,ch]=contour(lon,lat,-h',-[0:2:10]);
keyboard;




        %DEBUG:        keyboard;
        %DEBUG:        [w,d]=lastwarn; if strcmp(d,'MATLAB:triangulation:EmptyTri2DWarnId'); keyboard; lastwarn(''); end;
        %DEBUG:        [w,d]=lastwarn; if strcmp(d,'MATLAB:inpolygon:ModelingWorldLower'); keyboard; lastwarn(''); end;




%DEBUG:
fmg; surf(lon,lat,repmat(0,size(h')),-h'); colorbar; shading interp; caxis([-80,0]);
axis([-80.21,-80.135,25.33,25.40]); daspect([1,cosd(25),1]); view(2); colorbar('off'); colorbar('Location','EastOutside');

fmg; surf(lon,lat,repmat(0,size(hcrng')),hcrng'); colorbar; shading interp; caxis([0 2e3]);
axis([-80.21,-80.135,25.33,25.40]); daspect([1,cosd(25),1]); view(2); colorbar('off'); colorbar('Location','EastOutside');
[cc,ch]=contour(lon,lat,-h',-[0:2:10]);





w=warning('OFF','MATLAB:inpolygon:ModelingWorldLower');
warning(w);





  for shpix=1:numel(shps)
    if (mod(shpix,checkin_count) == 0); disp(shpix); toc; tic; end;
    [lonix,latix] = fieldfind(lon,lat,shps(shpix).X(1:end-1),shps(shpix).Y(1:end-1));
    goodix = find(~isnan(lonix) & ~isnan(latix));
    lonix = lonix(goodix);    latix = latix(goodix);
    lix = unique([lonix;latix]','rows')';
    if ( ~isempty(lix) )
      lonix = lix(2,:);
      latix = lix(1,:);
      if ( numel(lonix) == 1 || numel(lonix) == 2 || shps(shpix).Shapearea < 1e5 )
        ix = sub2ind(size(LON),latix,lonix);
        %DEBUG:        keyboard;
        inix = union(inix,ix);
      elseif ( numel(lonix) > 2 )
        tr = delaunayTriangulation(lonix',latix');
        %DEBUG:
        [w,d] = lastwarn; if ( strcmp(d,'MATLAB:delaunayTriangulation:DupPtsWarnId') ); keyboard; lastwarn(''); end;
        ix = sub2ind(size(LON),tr.Points(:,2),tr.Points(:,1));
        ix = find(inpolygon(LON,LAT,LON(ix),LAT(ix))==1);
        %DEBUG:
        [w,d] = lastwarn; if ( strcmp(d,'MATLAB:inpolygon:ModelingWorldLower') ); keyboard; lastwarn(''); end;
        inix = union(inix,ix);
      end;
    end;
  end;






%DEBUG:
if 0;
  jix=3638; iix=587;
  jixen = jix-20:jix+20;
  jixen(1>jixen|jixen>sz1) = [];
  iixen = iix-20:iix+20;
  iixen(1>iixen|iixen>sz2) = [];
  fmg;
  contourf(lon(iixen),lat(jixen),h(jixen,iixen));
  colorbar('Location','East');
end;







oris = [ 315   0  45 ;...
         270 nan  90 ;...
         225 180 135 ];



maxsteps = 0;
    if ( n > maxsteps ); maxsteps = n; end;





for jix=sz1:-1:1
  if (mod(jix,checkin_count) == 0); disp([jix,sz2]); toc; tic; end;
  orig_jix = jix;
  for iix=1:sz2
    n = 0;
    jix = orig_jix;
    oldang = ang(jix,iix);
    while ( iix>0 && iix<=sz2 && jix>0 && jix<=sz1 ...
            && h(jix,iix)<0 && h(jix,iix)>-100 ...
            && bet(jix,iix) >= 0.0020 ...
            && cosd(oldang - ang(jix,iix)) >= 0 )
      niix = round(iix+u(jix,iix));
      njix = round(jix+v(jix,iix));
      hcrng(jix,iix) = hcrng(jix,iix) + dx;
      oldang = ang(jix,iix);
      iix=niix; jix=njix;
      n = n+1;
      if (mod(n,checkin_count) == 0); keyboard; end;
    end;
  end;
end;






if ( ~use_habitat_map && ~exist('rfhd','var') )
  % Paradoxically, if USE_HABITAT_MAP=True, we do not need to load these data
  read_coral_and_hard_bottom;
end;






.
.
.
if ( ~exist('FRT_depth_and_beta_92m_NANMEAN_3_3_4_coral.mat','file') )
disp('NANMEAN - Reef - Hit Enter to continue...'); pause;
tic,
  h=[]; bet=[]; LON=[]; LAT=[]; clear lon lat h bet LON LAT
  load('FRT_depth_and_beta_92m_NANMEAN_3_3_4.mat');
  [LON,LAT] = meshgrid(lon,lat);
  inix=[];
  disp(numel(shps));
  tic,
  dlat = 1.1*abs(lat(1)-lat(2));
  for shpix=1:numel(shps)
    if (mod(shpix,checkin_count) == 0); disp(shpix); toc; tic; end;
    if ( strcmp(shps(shpix).DESCRIPT,'Coral Reef') )
      ix = find(inpolygon(LON,LAT,shps(shpix).X,shps(shpix).Y)==1);
      if ( isempty(ix) )
        shpctr = mean(shps(shpix).BoundingBox);
        % Find the ONE gridpoint that is closest to this polygon
        [err,ix] = min(abs(LON(:)-shpctr(1)) + abs(LAT(:)-shpctr(2)));
        if ( err > dlat )
          ix = [];
        end;
      end;
      inix = union(inix,ix);
    end;
  end;
  toc,
  newht = repmat(nan,size(LON));
  ht=h'; newht(inix)=ht(inix); ht=[]; clear ht
  h=[]; h=newht'; newht=[]; clear newht
  newbett = repmat(nan,size(LON));
  bett=bet'; newbett(inix)=bett(inix); bett=[]; clear bett
  bet=[]; bet=newbett'; newbett=[]; clear newbett
  disp(numel(find(~isnan(h))));
  save('FRT_depth_and_beta_92m_NANMEAN_3_3_4_coral.mat','lon','lat','h','bet');
  clear inix ix shpix;
toc,
end;



if ( ~exist('FRT_depth_and_beta_92m_coral_and_hard_bottom.mat','file') )
disp('Full - RfHd - Hit Enter to continue...'); pause;
tic,
  h=[]; bet=[]; LON=[]; LAT=[]; clear lon lat h bet LON LAT
  load('FRT_depth_and_beta_92m.mat');
  [LON,LAT] = meshgrid(lon,lat);
  inix=[];
  disp(numel(shps));
  tic,
  dlat = 1.1*abs(lat(1)-lat(2));
  for shpix=1:numel(shps)
    if (mod(shpix,checkin_count) == 0); disp(shpix); toc; tic; end;
    ix = find(inpolygon(LON,LAT,shps(shpix).X,shps(shpix).Y)==1);
    if ( isempty(ix) )
      shpctr = mean(shps(shpix).BoundingBox);
      % Find the ONE gridpoint that is closest to this polygon
      [err,ix] = min(abs(LON(:)-shpctr(1)) + abs(LAT(:)-shpctr(2)));
      if ( err > dlat )
        ix = [];
      end;
    end;
    inix = union(inix,ix);
  end;
  toc,
  newht = repmat(nan,size(LON));
  ht=h'; newht(inix)=ht(inix); ht=[]; clear ht
  h=[]; h=newht'; newht=[]; clear newht
  newbett = repmat(nan,size(LON));
  bett=bet'; newbett(inix)=bett(inix); bett=[]; clear bett
  bet=[]; bet=newbett'; newbett=[]; clear newbett
  disp(numel(find(~isnan(h))));
  save('FRT_depth_and_beta_92m_coral_and_hard_bottom.mat','lon','lat','h','bet');
  clear inix ix shpix;
toc,
end;

if ( ~exist('FRT_depth_and_beta_92m_coral.mat','file') )
disp('Full - Reef - Hit Enter to continue...'); pause;
tic,
  h=[]; bet=[]; LON=[]; LAT=[]; clear lon lat h bet LON LAT
  load('FRT_depth_and_beta_92m.mat');
  [LON,LAT] = meshgrid(lon,lat);
  inix=[];
  disp(numel(shps));
  tic,
  dlat = 1.1*abs(lat(1)-lat(2));
  for shpix=1:numel(shps)
    if (mod(shpix,checkin_count) == 0); disp(shpix); toc; tic; end;
    if ( strcmp(shps(shpix).DESCRIPT,'Coral Reef') )
      ix = find(inpolygon(LON,LAT,shps(shpix).X,shps(shpix).Y)==1);
      if ( isempty(ix) )
        shpctr = mean(shps(shpix).BoundingBox);
        % Find the ONE gridpoint that is closest to this polygon
        [err,ix] = min(abs(LON(:)-shpctr(1)) + abs(LAT(:)-shpctr(2)));
        if ( err > dlat )
          ix = [];
        end;
      end;
      inix = union(inix,ix);
    end;
  end;
  toc,
  newht = repmat(nan,size(LON));
  ht=h'; newht(inix)=ht(inix); ht=[]; clear ht
  h=[]; h=newht'; newht=[]; clear newht
  newbett = repmat(nan,size(LON));
  bett=bet'; newbett(inix)=bett(inix); bett=[]; clear bett
  bet=[]; bet=newbett'; newbett=[]; clear newbett
  disp(numel(find(~isnan(h))));
  save('FRT_depth_and_beta_92m_coral.mat','lon','lat','h','bet');
  clear inix ix shpix;
toc,
end;

set_more;



%{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       TOO DOGGONE SLOW!      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

timenow;

h=[]; bet=[]; LON=[]; LAT=[]; clear lon lat h bet LON LAT
load('FRT_depth_and_beta_92m_NANMEAN_3_3_4.mat');
[LON,LAT] = meshgrid(lon,lat);
tic, inix = inpolygon(LON,LAT,rfhd.polylon,rfhd.polylat); toc,
h(~inix) = NaN;
bet(~inix) = NaN;
save('FRT_depth_and_beta_92m_NANMEAN_3_3_4_coral_and_hard_bottom.mat','lon','lat','h','bet');


h=[]; bet=[]; LON=[]; LAT=[]; clear lon lat h bet LON LAT
load('FRT_depth_and_beta_92m_NANMEAN_3_3_4.mat');
[LON,LAT] = meshgrid(lon,lat);
tic, inix = inpolygon(LON,LAT,reef.polylon,reef.polylat); toc,
h(~inix) = NaN;
bet(~inix) = NaN;
save('FRT_depth_and_beta_92m_NANMEAN_3_3_4_coral.mat','lon','lat','h','bet');


h=[]; bet=[]; LON=[]; LAT=[]; clear lon lat h bet LON LAT
load('FRT_depth_and_beta_92m.mat');
[LON,LAT] = meshgrid(lon,lat);
tic, inix = inpolygon(LON,LAT,rfhd.polylon,rfhd.polylat); toc,
h(~inix) = NaN;
bet(~inix) = NaN;
save('FRT_depth_and_beta_92m_coral_and_hard_bottom.mat','lon','lat','h','bet');


h=[]; bet=[]; LON=[]; LAT=[]; clear lon lat h bet LON LAT
load('FRT_depth_and_beta_92m.mat');
[LON,LAT] = meshgrid(lon,lat);
tic, inix = inpolygon(LON,LAT,reef.polylon,reef.polylat); toc,
h(~inix) = NaN;
bet(~inix) = NaN;
save('FRT_depth_and_beta_92m_coral.mat','lon','lat','h','bet');

timenow;
%}













    [lonix,latix] = fieldfind(lon,lat,shps(shpix).X,shps(shpix).Y);
    goodix = find(~isnan(lonix) & ~isnan(latix));
    tr = delaunayTriangulation(lonix(goodix)',latix(goodix)');
    ix = sub2ind(size(LON),tr.Points(:,2),tr.Points(:,1));




    shpixen = sub2ind(size(LON),tr.Points(:,1),tr.Points(:,2));
    if ( numel(shpixen) < 3 )
      ix = shpixen;
    else
      ix = find(inpolygon(LON,LAT,tr.Points(:,1),tr.Points(:,2))==1);
      [w,d]=lastwarn; if ( strcmp(d,'MATLAB:inpolygon:ModelingWorldLower') ); keyboard; end;
    end;
    if ( isempty(ix) )
      shpctr = mean(shps(shpix).BoundingBox);
      % Find the ONE gridpoint that is closest to this polygon
      [err,ix] = min(abs(LON(:)-shpctr(1)) + abs(LAT(:)-shpctr(2)));
keyboard;
      if ( err > dlat )
        ix = [];
      end;
    end;





  for shpix=1:numel(shps)
    if (mod(shpix,checkin_count) == 0); disp(shpix); toc; tic; end;
    [lonix,latix] = fieldfind(lon,lat,shps(shpix).X,shps(shpix).Y);
    goodix = find(~isnan(lonix) & ~isnan(latix));
    tr = delaunayTriangulation(lon(lonix(goodix)),lat(latix(goodix)));
    if ( numel(shpixen) < 3 )
      ix = shpixen;
    else
      ix = find(inpolygon(LON,LAT,tr.Points(:,1),tr.Points(:,2))==1);
      [w,d]=lastwarn; if ( strcmp(d,'MATLAB:inpolygon:ModelingWorldLower') ); keyboard; end;
    end;
    if ( isempty(ix) )
      shpctr = mean(shps(shpix).BoundingBox);
      % Find the ONE gridpoint that is closest to this polygon
      [err,ix] = min(abs(LON(:)-shpctr(1)) + abs(LAT(:)-shpctr(2)));
keyboard;
      if ( err > dlat )
        ix = [];
      end;
    end;





    w=warning('off'); tr = delaunayTriangulation(shplonix(goodix)',shplatix(goodix)'); warning(w);
keyboard;
    shpixen = sub2ind(size(LON),shplatix(goodix),shplonix(goodix));
    shpixen(diff(shpixen)==0) = [];




    w=warning('off'); tr = delaunayTriangulation(LON(shpixen)',LAT(shpixen)'); warning(w);
keyboard;




    ix = find(inpolygon(LON,LAT,LON(shpixen),LAT(shpixen))==1);
    ix = find(inpolygon(LON,LAT,LON(shplatix,shplonix),LAT(shplatix,shplonix))==1);





if ( ~exist('FRT_depth_and_beta_92m_NANMEAN_3_3_4_coral_and_hard_bottom.mat','file') )
disp('NANMEAN - RfHd - Hit Enter to continue...'); pause;
tic,
  h=[]; bet=[]; LON=[]; LAT=[]; clear lon lat h bet LON LAT
  load('FRT_depth_and_beta_92m_NANMEAN_3_3_4.mat');
  [LON,LAT] = meshgrid(lon,lat);
  inix=[];
  disp(numel(shps));
  tic,
  dlat = 1.1*abs(lat(1)-lat(2));
  for shpix=1:numel(shps)
    if (mod(shpix,checkin_count) == 0); disp(shpix); toc; tic; end;
    ix = find(inpolygon(LON,LAT,shps(shpix).X,shps(shpix).Y)==1);
    if ( isempty(ix) )
      shpctr = mean(shps(shpix).BoundingBox);
      % Find the ONE gridpoint that is closest to this polygon
      [err,ix] = min(abs(LON(:)-shpctr(1)) + abs(LAT(:)-shpctr(2)));
      if ( err > dlat )
        ix = [];
      end;
    end;
    inix = union(inix,ix);
  end;
  toc,
  newht = repmat(nan,size(LON));
  ht=h'; newht(inix)=ht(inix); ht=[]; clear ht
  h=[]; h=newht'; newht=[]; clear newht
  newbett = repmat(nan,size(LON));
  bett=bet'; newbett(inix)=bett(inix); bett=[]; clear bett
  bet=[]; bet=newbett'; newbett=[]; clear newbett
  disp(numel(find(~isnan(h))));
  save('FRT_depth_and_beta_92m_NANMEAN_3_3_4_coral_and_hard_bottom.mat','lon','lat','h','bet');
  clear inix ix shpix;
toc,
end;







  % If we got vectors, make sure their number of elements match
  if ( numel(lon)>1 || numel(lat)>1 || iscellstr(stnm) )
    if ( ~isvector(lon) || ~isvector(lat) || (~isempty(stnm) && ~iscellstr(stnm)) )
      error('If STN is a cell containing vectors, it must be only vectors');
    elseif ( numel(lon) ~= numel(lat) || (~isempty(stnm) && numel(lon) ~= numel(stnm)) )
      error('If STN is a cell containing vectors, their sizes must match');
    end;
  end;






  [cc,ch]=contour(lon,lat,h',[0 0],'-','Color',[.8,.8,.8],'LineWidth',2);
  [cc,ch]=contour(lon,lat,h',[10,30,maxdepth],'-','Color',[.8,.8,.8],'LineWidth',1.5);
  [cc2,ch2]=contour(lon,lat,h',[5,15],':','Color',[.8,.8,.8],'LineWidth',1);

  if ( ~use_habitat_map )
    % Coral reef habitat contour lines (use REEF or RFHD)
    plot(rfhd.polylon,rfhd.polylat,'-','Color',[.4,.4,.4],'LineWidth',0.5);
  elseif ( plot_coastline )
    plot_hires_coastline(bath);
  end;




  % Basic bathymetry contour lines
  if ( plot_coastline )
    [cc,ch] = contour(lon,lat,h',[0 0],'-','Color',[.8,.8,.8],'LineWidth',2);
    [cc,ch] = contour(lon,lat,h',[10 30 maxdepth],'-','Color',[.8,.8,.8],'LineWidth',1.5);
    [cc2,ch2] = contour(lon,lat,h',[5 15],':','Color',[.8,.8,.8],'LineWidth',1);
  end;

  if ( ~use_habitat_map )
    % Coral reef habitat contour lines (use REEF or RFHD)
    plot(rfhd.polylon,rfhd.polylat,'-','Color',[.4,.4,.4],'LineWidth',0.5);
  elseif ( plot_coastline )
    plot_hires_coastline(bath);
  end;








if ( ~exist('FRT_depth_and_beta_92m_NANMEAN_3_3_4_coral_and_hard_bottom.mat','file') )
disp('NANMEAN - RfHd - Hit Enter to continue...'); pause;
tic,
  h=[]; bet=[]; LON=[]; LAT=[]; clear lon lat h bet LON LAT
  load('FRT_depth_and_beta_92m_NANMEAN_3_3_4.mat');
  [LON,LAT] = meshgrid(lon,lat);
  inix=[];
  disp(numel(shps));
  tic,
  dlat = 1.1*abs(lat(1)-lat(2));
  for shpix=1:numel(shps)
    if (mod(shpix,checkin_count) == 0); disp(shpix); toc; tic; end;
    ix = find(inpolygon(LON,LAT,shps(shpix).X,shps(shpix).Y)==1);
    if ( isempty(ix) )
      if ( shps(shpix).Shapearea < 1e5 )
        shpctr = mean(shps(shpix).BoundingBox);
        % Find the ONE gridpoint that is closest to this polygon
        [err,ix] = min(abs(LON(:)-shpctr(1)) + abs(LAT(:)-shpctr(1,2)));
        if ( err > dlat )
          ix = [];
        end;
      end;
    end;
    inix = union(inix,ix);
  end;
  toc,
  newht = repmat(nan,size(LON));
  ht=h'; newht(inix)=ht(inix); ht=[]; clear ht
  h=[]; h=newht'; newht=[]; clear newht
  newbett = repmat(nan,size(LON));
  bett=bet'; newbett(inix)=bett(inix); bett=[]; clear bett
  bet=[]; bet=newbett'; newbett=[]; clear newbett
  disp(numel(find(~isnan(h))));
  save('FRT_depth_and_beta_92m_NANMEAN_3_3_4_coral_and_hard_bottom.mat','lon','lat','h','bet');
  clear inix ix shpix;
toc,
end;







keyboard;
      ix = interp2(LON,LAT,shps(shpix).X,shps(shpix).Y,'nearest');



      plot_station_marker(bnpin,~use_habitat_map);
      plot_station_marker(bnpon,~use_habitat_map);
      plot_station_marker(bnpmi,~use_habitat_map);






        hdrs = x.textdata(1,2:end);




      x = importdata(fname);
      if ( iscellstr(x) )
        % Weirdo special case - usually when there is just one data column
        cols = split(x,',');
        if ( ~strcmpi(cols{1,1},'datetime' )
        error('Invalid CHAMP Portal file download format?? %s',fname);
        dts = datenum(char(cols(2:end,1)));
        flds = cols(1,2:end);
        data = 
      end;
      if ( ~isfield(x,'textdata') || isempty(x.textdata) || ~strcmpi(x.textdata{1,1},'datetime') )
        error('Invalid CHAMP Portal file download format?? %s',fname);
      end;

      dts = datenum(x.textdata(2:end,1));
      flds = x.textdata(1,2:end);

      if ( numel(flds) ~= size(x.data,2) )
        error('Ecoforecasts:Portal:ColumnMismatch',...
              'Expected %d columns but found %d in %s',numel(flds),size(x.data,2),fname);
      end;
      if ( numel(dts) ~= size(x.data,1) )
        error('Ecoforecasts:Portal:RowMismatch',...
              'Expected %d rows but found %d in %s',numel(dts),size(x.data,1),fname);
      end;

      for fldix = 1:numel(flds)
        fld = flds{fldix};
        if ( ~strncmpi(fld,stnm_prefix,length(stnm_prefix)) )
          warning('Ecoforecasts:Portal:SuspectField',...
                  'Station name prefix %s not found in field %s',...
                  stnm_prefix,fld);
        else
          fld = lower(fld(length(stnm_prefix)+1:end));
        end;

        fres.(fld).date = dts;
        fres.(fld).data = x.data(:,fldix);
      end;






  dx_SU=[]; clear dx_SU





  switch (qix),
   case 1,
    %T = winter_mean + dTW*8 + dTC*2;
    T = winter_mean + (4*dTW*5) + (2*dTC*1);
    surf(lon,lat,zeroZ',T'); caxis([winter_min,winter_mean]); shading interp
    % %contourf(lon,lat,-T',-[winter_min:0.50:winter_mean]); caxis(-[winter_mean,winter_min]);
    % surf(lon,lat,zeroZ',-T'); caxis(-[winter_mean,winter_min]); shading interp
    % cm=parula; colormap(cm(end:-1:1,:));
    % %cm=cool; colormap(cm(end:-1:1,:))
    % %cm=cool; colormap(cm);
    % rescale_colorbar(cbh,@(x)(-x),'%g');
   case 2,
    T = winter_mean + dTW*10;
    surf(lon,lat,zeroZ',T'); caxis([winter_min,winter_mean]); shading interp
    % %contourf(lon,lat,-T',-[winter_min:0.50:winter_mean]); caxis(-[winter_mean,winter_min]);
    % surf(lon,lat,zeroZ',-T'); caxis(-[winter_mean,winter_min]); shading interp
    % cm=parula; colormap(cm(end:-1:1,:));
    % %cm=cool; colormap(cm(end:-1:1,:))
    % %cm=cool; colormap(cm);
    % rescale_colorbar(cbh,@(x)(-x),'%g');
   case 3,
    T = summer_mean + dTS*10;
    %contourf(lon,lat,T',[summer_mean:0.50:summer_max]); caxis([summer_mean,summer_max]);
    surf(lon,lat,zeroZ',T'); caxis([summer_mean,summer_max]); shading interp
    cm=hot; colormap(cm(end:-1:1,:));
    cbh=colorbar;
   case 4,
    T = summer_mean + dTS*8 + dTB*2;
    %contourf(lon,lat,T',[summer_mean:0.50:summer_max]); caxis([summer_mean,summer_max]);
    surf(lon,lat,zeroZ',T'); caxis([summer_mean,summer_max]); shading interp
    cm=hot; colormap(cm(end:-1:1,:));
    cbh=colorbar;
  end;






  % When using SURF instead of CONTOURF... The following incantations are
  % needed to make station locations and contours visible on warm maps.
  if ( q0 > 0 )
    [az,el] = view;
    view(az,270);
    set(gca,'XDir','Reverse');
    %zlim([-34,34]);
  end;






doLog = true;

for qix=1:nq0s;

  % Net sea-surface heat flux, Q_0
  q0 = q0s(qix);
  if (doLog)
    disp(['Log warming: Q0 = ',num2str(q0)]);
  else
    disp(['3 day del T: Q0 = ',num2str(q0)]);
  end;
  dT = squeeze(dTdt_SS(qix,:,:));

  fh=fmg;
  if (doLog)
    dTlog = log10(abs(dT));
    dTlog(maxdepth<h | h<mindepth) = nan;
    %contourf(lon,lat,dTlog',[-1.00:0.05:0.50]);  caxis([-1,0]);
    contourf(lon,lat,dTlog',[-1.30:0.05:0.50]);  caxis([-1.3,0.2]);
    set_pcolor_cursor(fh,@xyz_logc_id_select_cb);
  else
    if ( q0 > 0 )
      % MLRF1 July mean 28C
      T = 28 + dT*3;
      contourf(lon,lat,T',[28:0.50:33]);  caxis([28,33]);
    else
      % MLRF1 February mean 23C
      T = 23 + dT*3;
      contourf(lon,lat,T',[12:0.50:23]);  caxis([12,23]);
    end;
  end;

  if ( q0 > 0 )
    phenomstr = 'warming';
    % Reversed HOT shows warmest spots as darkest red
    cm=hot; colormap(cm(end:-1:1,:))
  else
    phenomstr = 'cooling';
    cm=cool; colormap(cm);
  end;

  cbh=colorbar;
  if (doLog)
    % Make a logarithmic color bar
    rescale_colorbar(cbh);
  end;

  daspect([1,cosd(lat(1)),1]); 

  % Adjust graph appearance based on shape of sub-region
  switch (subrgn),
   case 'SE',
    view(90,90);
    set(cbh,'Orientation','horizontal','Position',[0.30 0.30 0.50 0.02]);
   case 'UK',
    view(90,90);
    set(cbh,'Location','East');
   case 'MK',
    view(25,90);
    %set(gca,'Position',[-0.5,-0.5,2,2]);
    %set(cbh,'Orientation','horizontal','Position',[0.30 0.10 0.50 0.02]);
    set(gca,'Position',[-0.52,-0.30,2.0,1.9]);
    set(get(gca,'Title'),'Units','normalized','Position',[0.5,0.4]);
    set(cbh,'Orientation','horizontal','Position',[0.25 0.40 0.50 0.02]);
   case 'LK',
    set(cbh,'Orientation','horizontal','Position',[0.50,0.30,0.45,0.03]);
   case 'DT',
    %set(cbh,'Location','East');
  end;

  % Basic bathymetry contour lines
  [cc,ch]=contour(lon,lat,h',[0 0],'-','Color',[.8,.8,.8],'LineWidth',2);
  [cc,ch]=contour(lon,lat,h',[10,30,maxdepth],'-','Color',[.8,.8,.8],'LineWidth',1.5);
  [cc2,ch2]=contour(lon,lat,h',[5,15],':','Color',[.8,.8,.8],'LineWidth',1);

  if ( ~use_habitat_map )
    axis(axis);
    % Coral reef habitat contour lines (use REEF or RFHD)
    plot(rfhd.polylon,rfhd.polylat,'-','Color',[.4,.4,.4],'LineWidth',0.5);
  end;

  if (doLog)
    titlename(sprintf('%s: %s (%+d W/m^2) daily reef %s',subrgnstr,q0str{qix},q0,phenomstr));
    if ( doPrint )
      print('-dpng',fullfile(figspath,sprintf('%s_%s_%s_%dW.png',figbasename,subrgn,phenomstr,abs(q0s(qix)))));
    end;
  else
    titlename(sprintf('%s: temperature after 3 d %s (%+d W/m^2)',subrgnstr,q0str{qix},q0));
    if ( doPrint )
      print('-dpng',fullfile(figspath,sprintf('%s_%s_3_days_%s_%dW.png',figbasename,subrgn,phenomstr,abs(q0s(qix)))));
    end;
  end;

end; %for qix=1:nq0s;









    % Coral reef habitat contour lines (was RFHD)
    plot(reef.polylon,reef.polylat,'-','Color',[.4,.4,.4],'LineWidth',0.5);



    % Reef and hard-bottom habitat contour lines
    plot(rfhd.polylon,rfhd.polylat,'-','Color',[.8,.8,.8],'LineWidth',2);







if ( ~exist('allow_hard_bottom','var') || isempty(allow_hard_bottom) )
  if ( strcmpi(use_habitat_map,'rfhd') )
    allow_hard_bottom = true;
  else
    allow_hard_bottom = false;
  end;
end;



  if ( use_habitat_map )
    % Basic bathymetry contour lines
    [cc,ch]=contour(lon,lat,h',[0 0],'-','Color',[.8,.8,.8],'LineWidth',2);
    [cc,ch]=contour(lon,lat,h',[10,30,maxdepth],'-','Color',[.8,.8,.8],'LineWidth',1.5);
    [cc2,ch2]=contour(lon,lat,h',[5,15],':','Color',[.8,.8,.8],'LineWidth',1);
  else
    % Reef and hard-bottom habitat contour lines
    axis(axis);
    plot(rfhd.lon,rfhd.lat,'-','Color',[.8,.8,.8],'LineWidth',2);
  end;
















    if ( shps(shpix).Shapearea > 1e7 )
      ix = find(inpolygon(LON,LAT,shps(shpix).X,shps(shpix).Y)==1);
      inix = union(inix,ix);
    end;




  newh = repmat(nan,size(h));
  newh(inix) = h(inix);
  newbet = repmat(nan,size(bet));
  newbet(inix) = bet(inix);
  h=[]; h=newh; newh=[]; clear newh
  bet=[]; bet=newbet; newbet=[]; clear newbet





  newht = repmat(nan,size(LON));
  ht=h'; newht(inix)=ht(inix); ht=[]; clear ht; newh=newht'; newht=[]; clear newht;
  newbett = repmat(nan,size(LON));
  bett=bet'; newbett(inix)=bett(inix); bett=[]; clear bett; newbet=newbett'; newbett=[]; clear newbett;
  h=[]; h=newh; newh=[]; clear newh
  bet=[]; bet=newbet; newbet=[]; clear newbet




    if (mod(shpix,1000) == 0); toc; disp(shpix); tic; end;
    ix = find(inpolygon(LON,LAT,shps(shpix).X,shps(shpix).Y));



    if ( shps(shpix).Shapearea > 5e6 )
      if (numel(ix) > 500); disp(['FOUND >500: ',num2str(shpix),',',num2str(numel(ix))]); end;
    end;



FROM plot_range_vs_all_controls.m:

1;

%doPrint = false;
doPrint = true;

% Linear R2 ~ 0.58, Exp(-1/3) R2 ~ 0.92! (N=30)
doAccDist = true; HsRange=[]; fillBathGaps = false; plot_range_vs_control
 bat=[]; bath=[]; kwbat=[]; kwbath=[]; stns=[]; clear -regexp ^[a-cegi-z]*; pack

% For 7-31 to 10-1, linear R2 ~ 0.74! (N=44)
doAccDist = true; HsRange=[]; fillBathGaps = false; ctlvar = 'raw_nwps_sigwavehgt'; plot_range_vs_control
 bat=[]; bath=[]; kwbat=[]; kwbath=[]; stns=[]; clear -regexp ^[a-cegi-z]*; pack

% For 7-31 to 10-1, linear R2 ~ 0.56, exp(-1/3) R2 ~ 0.91! (N=28)
HsRange=[0,0.5]; fillBathGaps = false; plot_range_vs_control
 bat=[]; bath=[]; kwbat=[]; kwbath=[]; stns=[]; clear -regexp ^[a-cegi-z]*; pack


doAccDist = true; HsRange=[0,0.5]; fillBathGaps = false; plot_range_vs_control

persub = @ts_jas; HsRange=[]; fillBathGaps = false; ctlvar = 'raw_nwps_sigwavehgt'; plot_range_vs_control
 bat=[]; bath=[]; kwbat=[]; kwbath=[]; stns=[]; clear -regexp ^[a-cegi-oq-z]*; pack

persub = @ts_jas; HsRange=[0,0.4]; fillBathGaps = false; plot_range_vs_control
persub = @ts_jas; HsRange=[0,0.5]; fillBathGaps = false; plot_range_vs_control
 bat=[]; bath=[]; kwbat=[]; kwbath=[]; stns=[]; clear -regexp ^[a-cegi-oq-z]*; pack


%{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

doPrint = false;

% For 7-31 to 10-1, linear R2 ~ 0.49, exp(-1/3) R2 ~ 0.44! (N=32)
% (Setting fillBathGaps adds a "Channel" site! But also a "Shoals", and...?)
HsRange=[0,0.5]; fillBathGaps = true; plot_range_vs_control
 bat=[]; bath=[]; kwbat=[]; kwbath=[]; stns=[]; clear -regexp ^[a-cegi-z]*; pack

% For 7-31 to 10-31, linear R2 ~ 0.77! (N=10)
HsRange=[0,0.45]; fillBathGaps = false; plot_range_vs_control
 bat=[]; bath=[]; kwbat=[]; kwbath=[]; stns=[]; clear -regexp ^[a-cegi-z]*; pack

% For 7-31 to 10-1, linear R2 ~ 0.79! (N=36)
HsRange=[0,0.5]; fillBathGaps = true; ctlvar = 'raw_nwps_sigwavehgt'; plot_range_vs_control
 bat=[]; bath=[]; kwbat=[]; kwbath=[]; stns=[]; clear -regexp ^[a-cegi-z]*; pack

% For 7-31 to 10-1, linear R2 ~ 0.74! (N=44)
HsRange=[]; fillBathGaps = true; ctlvar = 'raw_nwps_sigwavehgt'; plot_range_vs_control
 bat=[]; bath=[]; kwbat=[]; kwbath=[]; stns=[]; clear -regexp ^[a-cegi-z]*; pack


%% Experiment: Look at HIGH wave-energy sites instead...
HsRange = [0.3,+inf];
%fillBathGaps = true;
fillBathGaps = false;
 ctlvar = 'ngdc_depth'; plot_range_vs_control
 bat=[]; bath=[]; kwbat=[]; kwbath=[]; stns=[]; clear -regexp ^[a-cegi-z]*; pack
 ctlvar = {'ngdc_depth',@nanmean,4,4,5}; plot_range_vs_control
 bat=[]; bath=[]; kwbat=[]; kwbath=[]; stns=[]; clear -regexp ^[a-cegi-z]*; pack
 ctlvar = 'ngdc_beta'; plot_range_vs_control
 bat=[]; bath=[]; kwbat=[]; kwbath=[]; stns=[]; clear -regexp ^[a-cegi-z]*; pack
 ctlvar = 'raw_nwps_sigswellhgt'; plot_range_vs_control
 bat=[]; bath=[]; kwbat=[]; kwbath=[]; stns=[]; clear -regexp ^[a-cegi-z]*; pack
 % ctlvar = {'ngdc_beta',@nanmean,3,3,5}; plot_range_vs_control
 % bat=[]; bath=[]; kwbat=[]; kwbath=[]; stns=[]; clear -regexp ^[a-cegi-z]*; pack
 ctlvar = {'ngdc_beta',@nanmean,4,4,5}; plot_range_vs_control
 bat=[]; bath=[]; kwbat=[]; kwbath=[]; stns=[]; clear -regexp ^[a-cegi-z]*; pack
 % ctlvar = {'ngdc_beta',@nanmean,5,5,5}; plot_range_vs_control
 % bat=[]; bath=[]; kwbat=[]; kwbath=[]; stns=[]; clear -regexp ^[a-cegi-z]*; pack

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%}














    set(fh, 'units','normalized', 'outerposition',[0 0 1.00 1.00], 'position',[0.0 0.0 0.92 0.92]);





    contourf(lon,lat,T',[winter_min:0.50:winter_mean]); caxis([winter_min,winter_mean]);




    % Make a logarithmic color bar
    th=get(cbh,'TickLabels');
    for thix=1:numel(th);
      th{thix}=num2str(power(10,str2num(th{thix})),'%5.2f');
    end;
    set(cbh,'TickLabels',th);




  [cc2,ch2]=contour(lon,lat,h',[5,15],':','Color',[.8,.8,.8],'LineWidth',1);




    contourf(lon,lat,T',[5:0.50:winter_mean]);  caxis([10,winter_mean]);

    contourf(lon,lat,T',[summer_mean:0.50:33]);  caxis([summer_mean,33]);





    contourf(lon,lat,T',[12:0.50:25]);  caxis([12,25]);

    contourf(lon,lat,T',[25:0.50:33]);  caxis([25,33]);






  [cc,ch]=contour(lon,lat,h',[0 0],'k-','LineWidth',2);





 case 'UK',
  subrgnbox = [-80.45,-80.00,24.95,25.55];
  subrgnstr = 'Upper Keys';

 case 'LK',
  subrgnbox = [-82.05,-80.95,24.35,24.75];
  subrgnstr = 'Lower Keys';



%BASED ON 24_h_sum
% MLRF1: December lower std. dev. and mean, June mean, June maximum
q0s = [-240,-115,95,255];
q0str={'Cold winter','Normal winter','Normal summer','Hot summer'};




% MLRF1: December 25th-percentile and mean, July mean, May 75th-percentile
q0s = [-240,-115,95,320];




  switch (subrgn),
   case 'SE',
    view(90,90);
    set(cbh,'Orientation','horizontal','Position',[0.30 0.30 0.50 0.02]);
   case 'MK',
    view(30,90);
    set(gca,'Position',[-0.5,-0.5,2,2])
    %axis([-81.0,-80.4,24.6,25.0]);
    %set(cbh,'Position',[0.30 0.10 0.02 0.50]);
    set(cbh,'Orientation','horizontal','Position',[0.30 0.10 0.50 0.02]);
   case 'LK',
    %set(cbh,'Orientation','horizontal','Location','South');
    set(cbh,'Orientation','horizontal','Position',[0.50,0.30,0.45,0.03]);
  end;




  % Basic bathymetry contour lines
  [cc,ch]=contour(lon,lat,h',[0 0],'k-','LineWidth',2);
  [cc,ch]=contour(lon,lat,h',[10:10:maxdepth],'-','Color',[.8,.8,.8],'LineWidth',1.5);
  [cc2,ch2]=contour(lon,lat,h',[5:10:maxdepth],':','Color',[.8,.8,.8],'LineWidth',1);



   case 'LK',
    %set(cbh,'Orientation','horizontal','Location','South');
    set(cbh,'Orientation','horizontal','Position',[0.0497,0.2704,0.9372,0.0344]);



switch (subrgn),
 case 'SE',
  subrgnbox = [-80.20,-79.95,26.00,27.30];
  subrgnstr = 'Southeast Florida';
 case 'UK',
  subrgnbox = [-80.40,-80.00,25.00,26.00];
  subrgnstr = 'Upper Florida Keys';
 case 'MK',
  subrgnbox = [-81.00,-80.40,24.50,25.00];
  subrgnstr = 'Middle Florida Keys';
 case 'LK',
  subrgnbox = [-82.00,-81.00,24.30,24.70];
  subrgnstr = 'Lower Florida Keys';
 case 'DT',
  subrgnbox = [-82.95,-82.00,24.30,24.70];
  subrgnstr = 'Dry Tortugas Area';
 otherwise,
  error('Region %s not yet implemented!',subrgn);
end;




switch (subrgn),
 case 'SE',
  bbox = [-80.20,-79.95,26.00,27.30];
  subrgnstr = 'Southeast Florida';
 case 'UK',
  bbox = [-80.40,-80.00,25.00,26.00];
  subrgnstr = 'Upper Florida Keys';
 case 'MK',
  bbox = [-81.00,-80.40,24.50,25.00];
  subrgnstr = 'Middle Florida Keys';
 case 'LK',
  bbox = [-82.00,-81.00,24.30,24.70];
  subrgnstr = 'Lower Florida Keys';
 case 'DT',
  bbox = [-82.95,-82.00,24.30,24.70];
  subrgnstr = 'Dry Tortugas Area';
 otherwise,
  error('Region %s not yet implemented!',subrgn);
end;




% % 30-m resolution bathymetry
% downfunc = {@nanmean,9,9,10}; downstride = 9;
% downfuncstr = 'NANMEAN_9_9_10';
% 92-m resolution bathymetry
downfunc = {@nanmean,3,3,4}; downstride = 3;
downfuncstr = 'NANMEAN_3_3_4';



    % Downsample to appropriate resolution (HC affects water over ~120 m)
    disp('Downsampling H and BETA to ~120 m');
    [LAT,LON] = meshgrid(lat(1:downstride:end),lon(1:downstride:end));
    H = interp_field(lat,lon,h,LAT,LON,downfunc);
    H = reshape(H,size(LAT));
    BET = interp_field(lat,lon,bet,LAT,LON,downfunc);
    BET = reshape(BET,size(LAT));






  %DEBUG:  disp(['Imag{B0} ',num2str(numel(find(imag(B0) ~= 0)))]);
  %DEBUG:  disp(['Imag{uf} ',num2str(numel(find(imag(uf) ~= 0)))]);
  %DEBUG:  disp(['Imag{Qv_SS} ',num2str(numel(find(imag(Qv_SS) ~= 0)))]);
  %DEBUG:  disp(['Imag{dTdtq0} ',num2str(numel(find(imag(dTdtq0) ~= 0)))]);
  %DEBUG:  disp(['Imag{u_SS} ',num2str(numel(find(imag(u_SS) ~= 0)))]);
  %DEBUG:  disp(['Imag{dTdtx_SS} ',num2str(numel(find(imag(dTdtx_SS) ~= 0)))]);
  %DEBUG:  disp(['Imag{dTdx_SS} ',num2str(numel(find(imag(dTdx_SS) ~= 0)))]);







% %nq0s = 60;
% nq0s = 4;
% %q0s = linspace(-500,0,nq0s);				% dissertation figure
% %q0s = linspace(-200,0,nq0s);
% %q0s = linspace(-250,0,nq0s);				% proposal figure
% %q0s = linspace(-250,+250,nq0s);			% smart proposal figure
% %q0s = linspace(-300,+300,nq0s);			% smarter proposal figure
% q0s = linspace(-300,+300,nq0s);			% smartest proposal figure

%q0s = [-800,400];
%q0s = [-1400,-800,-300,300,400,800];
%q0s = [-1400,800];
%q0s = [-400,400];






FROM quilt_frt_bathymetry.m:
    %pvgf1 = read_hires_bathymetry(pvgf1,[19e3,100e3],[],true);


    % % %sanf1 = plot_hires_bathymetry(sanf1,-[0:10:120],[24e3,25e3]);
    % % sanf1 = plot_hires_bathymetry(sanf1,-[0:10:120],[54e3,25e3]);
    % sanf1 = read_hires_bathymetry(sanf1,[54e3,25e3],[],true);




if 0;
  fmg; contourf(lon,lat,-h',-[0:10:120]); colorbar; daspect([1,cosd(lat(1)),1]);
  ax=axis;
  fmg; contourf(patch_lon,patch_lat,-patch_h',-[0:10:120]); colorbar; daspect([1,cosd(lat(1)),1]);
  axis(ax);
  fmg; contourf(lon,lat,-quilt_h',-[0:10:120]); colorbar; daspect([1,cosd(lat(1)),1]);
  axis(ax);
end;








for qix=1:nq0s;
  % Net sea-surface heat flux, Q_0
  q0 = q0s(qix);

  dT = squeeze(dTdt_SS(qix,:,:));
  dTlog = log10(abs(dT));
  dTlog(30<h | h<0.5) = nan;
  fmg; contourf(lon,lat,dTlog',[-1.00:0.05:0.50]);
  caxis([-1,0]);
  %camroll(-45);
  %colormap(summer);
  %colormap(parula(8));
  if ( q0 > 0 )
    phenomstr = 'warming';
    cm=hot; colormap(cm(end:-1:1,:))
  else
    phenomstr = 'cooling';
    cm=cool; colormap(cm);
  end;

  cbh=colorbar;
  th=get(cbh,'TickLabels');
  for thix=1:numel(th);
    th{thix}=num2str(power(10,str2num(th{thix})),'%5.2f');
  end;
  set(cbh,'TickLabels',th);
  daspect([1,cosd(lat(1)),1]); 
  if ( strcmp(rgn,'MK') )
    %view(45,90);
    view(30,90);
    set(gca,'Position',[-0.5,-0.5,2,2])
    axis([-81.0,-80.4,24.6,25.0]);
    set(cbh,'Position',[0.3 0.1 0.02 0.5]);
  end;

  % Basic bathymetry
  [cc,ch]=contour(lon,lat,h',[0 0],'k-','LineWidth',2);
  [cc,ch]=contour(lon,lat,h',[10:10:30],'-','Color',[.8,.8,.8],'LineWidth',1.5);
  [cc2,ch2]=contour(lon,lat,h',[5:10:30],':','Color',[.8,.8,.8],'LineWidth',1);

  titlename(sprintf('%+d W/m^2 daily reef %s',q0,phenomstr));
  if ( doPrint )
    print('-dpng',fullfile(figspath,sprintf('%s_%s_%s_%dW.png',mfilename,rgn,phenomstr,abs(q0s(2)))));
  end;
end;








  if ( strcmp(rgn,'MK') )
    %view(45,90);
    view(30,90);
    set(gca,'Position',[-0.5,-0.5,2,2])
    axis([-81.0,-80.4,24.6,25.0]);
    set(cbh,'Position',[0.3 0.1 0.02 0.5]);
  end;






  fmg; contourf(lon,lat,dTlog',[-1.00:0.05:0.50]);
  caxis([-1,0]);
  %camroll(-45);
  %colormap(summer);
  %colormap(parula(8));
  if ( q0 > 0 )
    phenomstr = 'warming';
  else
    phenomstr = 'cooling';
  end;
  cm=hot; colormap(cm(end:-1:1,:))
  cm=cool; colormap(cm);





  [cc,ch]=contour(lon,lat,h',[0:10:30],'-','Color',[.2,.2,.2],'LineWidth',2);
  [cc2,ch2]=contour(lon,lat,h',[5:10:30],':','Color',[.2,.2,.2],'LineWidth',1.5);




dT1 = squeeze(dTdt_SS(1,:,:));
dT1log = log10(abs(dT1));
dT1log(30<h | h<0.5) = nan;
fmg; contourf(lon,lat,dT1log',[-2.0:0.1:0.5]);
cbh=colorbar;
th=get(cbh,'TickLabels');
for thix=1:numel(th);
  th{thix}=num2str(power(10,str2num(th{thix})),'%5.2f');
end;
set(cbh,'TickLabels',th);
daspect([1,cosd(lat(1)),1]); 
set(gca,'CLim',[-2,1.5]);
titlename('Reef cooling per day [log(K)/d] for avg. sfc. cooling Q_0 = -100 W/m^2');
%camroll(-45);
%colormap(winter);
colormap(hot(8));
if ( strcmp(rgn,'MK') ); view(45,90); end;
if ( doPrint )
  print('-dpng',fullfile(figspath,sprintf('%s_%s_cooling_%dW.png',mfilename,rgn,abs(q0s(1)))));
end;





dT1 = squeeze(dTdt_SS(1,:,:));
dT1log = log10(abs(dT1));
dT1log(30<h | h<0.5) = nan;
fmg; contourf(lon,lat,dT1log',[-2.0:0.1:0.5]);
cbh=colorbar;
th=get(cbh,'TickLabels');
for thix=1:numel(th);
  th{thix}=num2str(power(10,str2num(th{thix})),'%5.2f');
end;
set(cbh,'TickLabels',th);
daspect([1,cosd(lat(1)),1]); 
set(gca,'CLim',[-2,1.5]);
titlename('Reef cooling per day [log(K)/d] for avg. sfc. cooling Q_0 = -100 W/m^2');
%camroll(-45);
%colormap(winter);
colormap(hot(8));
if ( strcmp(rgn,'MK') ); view(45,90); end;
if ( doPrint )
  print('-dpng',fullfile(figspath,sprintf('%s_%s_cooling_%dW.png',mfilename,rgn,abs(q0s(1)))));
end;


dT2 = squeeze(dTdt_SS(2,:,:));
dT2log = log10(abs(dT2));
dT2log(30<h | h<0.5) = nan;
fmg; contourf(lon,lat,dT2log',[-2.0:0.1:0.5]);
cbh=colorbar;
th=get(cbh,'TickLabels');
for thix=1:numel(th);
  th{thix}=num2str(power(10,str2num(th{thix})),'%5.2f');
end;
set(cbh,'TickLabels',th);
daspect([1,cosd(lat(1)),1]); 
set(gca,'CLim',[-1,1.5]);
titlename('Reef warming per day [log(K)/d] for avg. daily warming Q_0 = +100 W/m^2');
%camroll(-45);
%colormap(summer);
colormap(hot(8));
if ( strcmp(rgn,'MK') )
  view(45,90);
  set(gca,'Position',[-0.5,-0.5,2,2])
  axis([-81.0,-80.4,24.6,25.0]);
end;
if ( doPrint )
  print('-dpng',fullfile(figspath,sprintf('%s_%s_warming_%dW.png',mfilename,rgn,abs(q0s(2)))));
end;










dT1 = squeeze(dTdt_SS(1,:,:));
dT1log = log(abs(dT1));
%dT1log(30<h | h<0.5) = nan;
fmg; contourf(lon,lat,dT1log',[-6.50:0.50:5.50]);
cbh=colorbar;
th=get(cbh,'TickLabels');
for thix=1:numel(th);
  th{thix}=num2str(exp(str2num(th{thix})),'%5.2f');
end;
set(cbh,'TickLabels',th);
daspect([1,cosd(lat(1)),1]); 
set(gca,'CLim',[-2,1.5]);
titlename('Reef cooling per day [log(K)/d] for avg. sfc. cooling Q_0 = -100 W/m^2');
if ( doPrint )
  print('-dpng',fullfile(figspath,sprintf('%s_%s_cooling_%dW.png',mfilename,rgn,abs(q0s(1)))));
end;

dT2 = squeeze(dTdt_SS(2,:,:));
dT2log = log(abs(dT2));
%dT2log(30<h | h<0.5) = nan;
fmg; contourf(lon,lat,dT2log',[-6.50:0.50:5.50]);
cbh=colorbar;
th=get(cbh,'TickLabels');
for thix=1:numel(th);
  th{thix}=num2str(exp(str2num(th{thix})),'%5.2f');
end;
set(cbh,'TickLabels',th);
daspect([1,cosd(lat(1)),1]); 
set(gca,'CLim',[-2,1.5]);
titlename('Reef warming per day [log(K)/d] for avg. daily warming Q_0 = +100 W/m^2');
if ( doPrint )
  print('-dpng',fullfile(figspath,sprintf('%s_%s_warming_%dW.png',mfilename,rgn,abs(q0s(2)))));
end;





  print('-dpng',fullfile(figspath,sprintf('%s_cooling_%dW.png',mfilename,abs(q0s(1)))));




%%fmg; contourf(lon,lat,7*squeeze(dTdt_SS(1,:,:)),7*[-0.15:0.01:0.15]);
%fmg; contourf(lon,lat,7*squeeze(dTdt_SS(1,:,:)),7*[-10:0.10:0.00]);
%colorbar; daspect([1,cosd(lat(1)),1]);   
%titlename('Cooling per week [K/7d] when avg. sfc. cooling Q_0 = -100 W/m^2');






axis([-81.65,-80,24.4,26]); 





lonrmix = find(bbox(1)<=lon & lon<=bbox(2));





    sites.hs(isnan(sites.hs)) = sitehs(isnan(sites.hs)); sitehs=[]; clear sitehs
    sites.wu(isnan(sites.wu)) = sitewu(isnan(sites.wu)); sitewu=[]; clear sitewu
    sites.wv(isnan(sites.wv)) = sitewv(isnan(sites.wv)); sitewv=[]; clear sitewv
    sites.tp(isnan(sites.tp)) = sitetp(isnan(sites.tp)); sitetp=[]; clear sitetp





    hs = oversample_attenuate_field({ww3.lon,ww3.lat,ww3.hs.date,ww3.hs.field},4,[],bat,[],true);
    wu = oversample_attenuate_field({ww3.lon,ww3.lat,ww3.hs.date,ww3.u},4,[],bat,[],true);
    wv = oversample_attenuate_field({ww3.lon,ww3.lat,ww3.hs.date,ww3.v},4,[],bat,[],true);
    tp = oversample_attenuate_field({ww3.lon,ww3.lat,ww3.tp.date,ww3.tp.field},4,[],bat,[],true);






    else
      % ... OR, Interpolate using METHOD, to attenuate field near the coastline.
      newfldsz = [length(histr.date),length(deplat),length(deplon)];
      x = memory; avail_doubles = x.MaxPossibleArrayBytes/8; clear x;
      warning('OFF','MATLAB:chckxy:IgnoreNaN');
      if ( avail_doubles < (prod(newfldsz) / 4) )
        disp('Insufficient memory for true 3D interpolation: this may take a while...');
        %newfld = repmat(nan,size(histr.field));
        newfld = repmat(nan,newfldsz);
        for lonix=1:numel(deplon)
          newfld(:,:,lonix) = interp3(flat,fdts,flon,fldfld,...
                                      deplat,histr.date,deplon(lonix),method);
        end;
      else
        % We have enough memory to save the caller some time
        newfld = interp3(flat,fdts,flon,fldfld,deplat,histr.date,deplon,method);
      end;
      warning('ON','MATLAB:chckxy:IgnoreNaN');
    end;







    tp = oversample_attenuate_field({ww3.lon,ww3.lat,ww3.tp.date,ww3.tp.field},...
                                    {bat.lon,bat.lat},bat.field,[],true);





if ( ~doWaves || 1 )
  disp(['NOT using WW3 wave data']);
else
  disp(['Using WW3 wave data']);
  sites.hs = repmat(nan,[1,numel(sites.lons)]);
  sites.wu = repmat(nan,[1,numel(sites.lons)]);
  sites.wv = repmat(nan,[1,numel(sites.lons)]);
  sites.tp = repmat(nan,[1,numel(sites.lons)]);
  % Process three WW3 regions to get wave data for all sites of interest 
  %for ww3rgn = {'fks','bsc','sef'};
  for ww3rgn = {'fks'};
    % Climatologize NOAA Wave Watch III output surrounding a region RGN
    ww3 = seasonalize_ww3_region(ww3rgn{:},false,false,false);

    % Use bilinear interpolation to increase field resolution, ...
    % ... First over- or undersample the depth field to match our resolution, ...
    % ... Then mark land gridpoints in result as zero, ...
    % ... AND EITHER...
    % ... Linearly attenuate (scale) field for depths in range 0 to -20 m ...
    % ... OR, Interpolate using METHOD, to attenuate field near the coastline.
    % ... And finally, land-mask gridpoints in result as NaN
    hs = oversample_attenuate_field({ww3.lon,ww3.lat,ww3.hs.date,ww3.hs.field},...
                                    {sites.lons,sites.lats,ww3.hs.date},bat.field,[],true);
    wu = oversample_attenuate_field({ww3.lon,ww3.lat,ww3.hs.date,ww3.u},...
                                    {sites.lons,sites.lats,ww3.hs.date},bat.field,[],true);
    wv = oversample_attenuate_field({ww3.lon,ww3.lat,ww3.hs.date,ww3.v},...
                                    {sites.lons,sites.lats,ww3.hs.date},bat.field,[],true);
    tp = oversample_attenuate_field({ww3.lon,ww3.lat,ww3.tp.date,ww3.tp.field},...
                                    {sites.lons,sites.lats,ww3.hs.date},bat.field,[],true);
    interp_field(,sites.lons,sites.lats,DATE???);
    ww3=[]; clear ww3
    hs=[]; wu=[]; wv=[]; tp=[]; clear hs wu wv tp
  end;
end;






OLD SCRIPT CODE from OVERSAMPLE_ATTENUATE_FIELD.m:
  flat=[]; fdts=[]; flon=[]; clear flat fdts flon
  fldlat=[]; fldlon=[]; fldfld=[]; flddep=[]; clear fldlat fldlon fldfld flddep




  if ( dyrs(stix,find(yrs==1996)) > 200 ); flgstr = '*'; else flgstr = ' '; end;
  if ( tgap(stix,find(yrs==1996)) < 30 ); flgstr = '*'; else flgstr = ' '; end;




clear tyrs
yrs = 1991:2005;
str = [sprintf('%4d ',yrs),' STATION NAME'];
disp(str);
disp(repmat('#',[1,length(str)]));
for stix=1:numel(stnms); 
  str = '';
  for yrix=1:numel(yrs)
    yr = yrs(yrix);
    tyrs(stix,yrix) = numel(find(get_year(tses(stix).date)==yr));
    str = [str,sprintf('%02.0f%% ',(100*tyrs(stix,yrix)/(24*365)))];
  end;
  disp([sprintf('%4d ',tyrs(stix,:)),stnms{stix}]);
end;




scatter_fit(ctlvars,resvars,tctlvar,tresvar);
text(ctlvars,resvars,trim_stnms);
legend('Location','NorthEast');
%%axis([0,0.06,4,14]);
%axis([minctl,maxctl,minvar,maxvar]);
if ( doPrint )
  print('-dpng',fullfile(figspath,[printname,'.png']));
end;




try,
  %scatter_curve_fit(ctlvars',resvars','a*(x-b)^n',tctlvar,tresvar,[],'StartPoint',[1,0,-1/3]); %SS
  scatter_curve_fit(ctlvars',resvars','a*(x-b)^n',tctlvar,tresvar,[],'ConfInt',0.67,'StartPoint',[1,0,-1],'Robust','LAR','MaxFunEvals',5e4,'MaxIter',5e3,'TolFun',1e-7); %SS
  th = text(ctlvars,resvars,trim_stnms);
  axis([minctl,maxctl,minvar,maxvar]);
  if ( doPrint )
    print('-dpng',fullfile(figspath,[printname,'_hc_SS.png']));
  end;
catch ME,
  catchwarn(ME,'Exponential -1/3 fit failure');
end;








    resstr = strtrim( genresults(FM,Stats,Output,[],[],[],clev) );





% If caller does not specify USEHIGHESTRES, always CONTOUR [0,0] (or if no
% 0-contour is found, CONTOUR [-1,-1] instead); otherwise, behavior depends
% on range of BATH.lon,BATH.lat: if field is entirely within south Florida,
% load SOFLA_COAST.dat and call FILL to plot its contours.
%
% Optionally returns AX, and handle CH returned by FILL or CONTOUR.





jx=1059:1084; ix=1042:1061; 
fmg; quiver(squeeze(u(1,jx,ix)),squeeze(v(1,jx,ix))); quiver(squeeze(u(2,jx,ix)),squeeze(v(2,jx,ix)));


fmg; quiver(u(1,1,1:1000,1:1001),v(1,1,1:1000,1:1001)); quiver(u(1,2,1:1000,1:1001),v(1,2,1:1000,1:1001));

fmg; quiver(squeeze(u(1,1001:1100,1001:1101)),squeeze(v(1,1001:1100,1001:1101))); quiver(squeeze(u(2,1001:1100,1001:1101)),squeeze(v(2,1001:1100,1001:1101)));





          'ARTO1', ...
          -60.5206, ...         % ARTO1
          11.3009, ...          % ARTO1
          6.0, ...              % ARTO1





    if ( nargout > 2 )
      iso = ang - 90;
      iso(iso<0) = 360 + iso(iso<0);
    end;




      tic,
      for ix = 1:numel(uts.data)
        if ( mod(ix,1000) == 1 ); toc, disp(ix); tic, end;





  if ( ~exist('start_lon','var') || isempty(start_lon) )
    start_lon = mean(bath.lon);
  end;
  if ( ~exist('start_lat','var') || isempty(start_lat) )
    start_lat = mean(bath.lat);
  end;
  if ( ~exist('method','var') || isempty(method) )
    method = 1;
  end;
  if ( ~exist('maxiso','var') || isempty(maxiso) )
    maxiso = 30;
  end;
  if ( ~exist('miniso','var') || isempty(miniso) )
    miniso = 0;
  end;
  if ( ~exist('tgtpts','var') || isempty(tgtpts) )
    tgtpts = 1e4;
  end;
  if ( ~exist('maxres_t','var') || isempty(maxres_t) )
    maxres_t = 60;
  end;






        if ( secix == 1 )
          all_res_t = res_t;
        else
          n = numel(res_t.date);
          all_res_t.date(end+1:end+n) = res_t.date;
          all_res_t.data(end+1:end+n) = res_t.data;
        end;




        if ( halves == 1 )
          thishalf = 1:floor(numel(all_uts.date)/2);
        else
          thishalf = floor(numel(all_uts.date)/2)+1:numel(all_uts.date);
        end;





            % NOTE: If existing field has no timestamps, we *TRASH* it here!
            % (Also, this only trashes a field if its name happens to match a
            % loaded COLUMN HEADER. User-defined fields should be preserved.)
            elseif ( ~isfield(stn.(fld), 'date') || isempty(stn.(fld).date) || ...
                     ~isfield(stn.(fld), 'data') || isempty(stn.(fld).data) )
                if ( ~strcmp(fld,'station_name') && ~strcmp(fld,'lon') && ~strcmp(fld,'lat') && ~strcmp(fld,'depth') )
                  if ( isfield(stn.(fld), 'date') && isempty(stn.(fld).date) && ...
                       isfield(stn.(fld), 'data') && isempty(stn.(fld).data) )
                    warning('Ecoforecasts:mergedEmptyTS',...
                            'Field %s was empty: Replacing.', fld);
                  elseif ( isnumeric(stn.(fld)) && isnumeric(result.(fld)) && ...
                           ndims(stn.(fld)) == ndims(result.(fld)) && ...
                           all(size(stn.(fld)) == size(result.(fld))) && ...
                           all(stn.(fld)(:) == result.(fld)(:)) )
                    warning('Ecoforecasts:mergedNonTS',...
                            'Numeric fields %s are identical: replacing old value.', fld);
                  elseif ( ischar(stn.(fld)) && ischar(result.(fld)) )
                    warning('Ecoforecasts:mergedNonTS',...
                            'Merging two CHAR fields %s.', fld);
                    result.(fld) = {stn.(fld), result.(fld)};
                  elseif ( iscellstr(stn.(fld)) && ischar(result.(fld)) )
                    warning('Ecoforecasts:mergedNonTS',...
                            'Merging CELLSTR and CHAR fields %s.', fld);
                    result.(fld) = {stn.(fld){:}, result.(fld)};
                  elseif ( iscellstr(stn.(fld)) && iscellstr(result.(fld)) )
                    warning('Ecoforecasts:mergedNonTS',...
                            'Merging CELLSTR and CELLSTR fields %s.', fld);
                    result.(fld) = {stn.(fld){:}, result.(fld){:}};
                  else
                    warning('Ecoforecasts:mergedNonTS',...
                            'Field %s cannot be merged: REPLACING old value!', fld);
                  end;
                  stn = rmfield(stn, fld);
                  stn.(fld) = result.(fld);
                end;

            % Otherwise, merge all values - overwriting old data with new
            else






% If for any reason an existing field does not contain both '.data' and
% '.date' fields, a "best efforts" attempt is made to merge the two fields:
% e.g., if they're arrays of identical values, they are not merged; if they
% are CHAR (or a CHAR and a CELLSTR), they are appended together; if they are
% numeric or cell arrays of compatible dimensions, the RESULT field is
% appended to the STN field. If these best efforts fail, however, *the old
% field in STN is replaced by that in RESULT*, with a warning.


            % Handle numeric arrays specially, e.g., ADCP bin height vectors
            elseif ( isnumeric(stn.(fld)) )
              if ( isnumeric(result.(fld)) && ...
                   ndims(stn.(fld))==ndims(result.(fld)) && 
                   all(size(stn.(fld))==size(result.(fld))) )
                warning('Ecoforecasts:mergedReplaceNonTS',...
                        'Numeric field %s cannot be merged: REPLACING old value!', fld);
              else
                warning('Ecoforecasts:mergedReplaceNonTS',...
                        'Numeric field %s cannot be merged: REPLACING old value!', fld);
                stn = rmfield(stn, fld);
                stn.(fld) = result.(fld);
              end;

            % Special handling of straight-up numeric arrays and cell arrays,
            % e.g., vectors of ADCP bin heights, CELLSTRs of filenames, etc.
            elseif ( iscell(stn.(fld)) )
              if ( ischar(result.(fld)) )
              elseif ( iscell(result.(fld)) )
              else
              end;
              if ( 1 )
                warning('Ecoforecasts:mergedCell',...
                        'Appending to CELL Field %s.', fld);
                stn.(fld)(end+1:end+numel(result.(fld))) = result.(fld);
              else
                warning('Ecoforecasts:mergedReplaceCell',...
                        'Cell %s cannot be merged: REPLACING old value!', fld);
                stn = rmfield(stn, fld);
                stn.(fld) = result.(fld);
              end;







      %clh = clabel(cs,h);


      %clh = clabel(cs,h,doLabels);





% If DOLABELS (DEFAULT: empty) logical TRUE,
% pass *two* return values from @CONTOURFUN to CLABEL (v.); if DOLABELS is a
% numeric vector, then only label those contour levels; if nonempty, call
% also returns [CS,H]=@CONTOURFUN(...) and CLH=CLABEL(...). Finally, if
% DOLABELS is a cell, DOLABELS{1} must be logical for whether to do INLINE
% labels (DEFAULT: False); if DOLABELS{1}, call CLABEL(CS,DOLABELS{2:end}).
%
% If DOLABELS is True (DEFAULT: []), call CLABEL. If DOLABELS is a cell,
% first elt. must be a logical INLINELABELS (see CLABEL; DEFAULT: False);
% second arg may be logical, or a numerical array of depths to label.




%bathyrng = [120e3,95e3]; bath.lat = +25.00; bath.lon = -80.50;
%bathyrng = [95e3,95e3]; bath.lat = +25.00; bath.lon = -80.50;
%bathyrng = [80e3,80e3]; bath.lat = +25.00; bath.lon = -81.00;
%bathyrng = [50e3,50e3]; bath.lat = +25.00; bath.lon = -80.50;
%bathyrng = [20e3,20e3]; bath.lat = +25.00; bath.lon = -80.50;




From READ_HIRES_BATHYMETRY.m:

    nlats = numel(all_lats);
    nlons = numel(all_lons);
    lat_res_deg = min(diff(unique(all_lats(:))));
    lon_res_deg = min(diff(unique(all_lons(:))));
    lat_res = distance_wgs84(stn.lat,stn.lon,stn.lat+lat_res_deg,stn.lon)*1e3;
    lon_res = distance_wgs84(stn.lat,stn.lon,stn.lat,stn.lon+lon_res_deg)*1e3;
    dlat = ceil(rad(2)/lat_res);
    dlon = ceil(rad(1)/cosd(stn.lat)/lon_res);

    [lonerr,lonix]=min(abs(all_lons-stn.lon));
    [laterr,latix]=min(abs(all_lats-stn.lat));
    % NOTE: If we're stitching together multiple files, be careful around
    % file boundaries to adjust the radius of the rectangle we subset to!
    if ( lonerr > (lon_res_deg+eps) )
      ldlon = dlon; rdlon = dlon;
disp(mfilename);
keyboard;
    else
      ldlon = dlon; rdlon = dlon;
    end;
    if ( laterr > (lat_res_deg+eps) )
      tdlat = dlat; bdlat = dlat;
disp(mfilename);
keyboard;
    else
      tdlat = dlat; bdlat = dlat;
    end;
    latixen = latix-bdlat:latix+tdlat;
    lonixen = lonix-ldlon:lonix+rdlon;
    latixen(1 > latixen | latixen > nlats) = [];
    lonixen(1 > lonixen | lonixen > nlons) = [];

    if ( lonerr > (2*lon_res) || laterr > (2*lat_res) || isempty(latixen) || isempty(lonixen) )
      error('Location %f,%f lies outside domain of "%s"',stn.lon,stn.lat,bathfile);
    end;

    if ( ~isempty(latixen) && ~isempty(lonixen) )
      if ( doAsc )
        z = all_zs(latixen,lonixen)';
      else
        % Non-standard variable names from GEBCO
        if ( ~isempty(strfind(bathfile,'GEBCO')) )
          z = cast(nc{'elevation'}(latixen,lonixen),'double')';
        else
          z = cast(nc{'z'}(latixen,lonixen),'double')';
        end; %if ( ~isempty(strfind(bathfile,'GEBCO')) )
        close(nc);
      end; %if ( doAsc )
    end; %if ( ~isempty(latixen) && ~isempty(lonixen) )






FROM BATHYSTATS.m:

if 0;
  % Save this code for later, like chewing gum under the desk...
  % Find nearest shore point & distance/heading from this point to our center
  [nearestShore,coast,nearestShoreIx,dist2Shore0,az2Shore0] = ...
      get_nearest_points(bath.ngdc_hires_bathy,[bath.lon;bath.lat]);

  if ( exist('findNearestShore','var') && findNearestShore )
    onshore_hdg = az2Shore0;
  else
    onshore_hdg = azimuth_wgs84(bath.lat,bath.lon,stns(5).lat,stns(5).lon);
  end;

  disp([az2Shore0,onshore_hdg]);

  offshore_hdg = mod(onshore_hdg+180,360);
  [xlons,xlats] = transect_wgs84(bath.lon,bath.lat,-dist2Shore0:0.05:dist2Shore0,offshore_hdg);
  [xlons2,xlats2] = transect_wgs84(stns(6).lon,stns(6).lat,-dist2Shore0:0.05:dist2Shore0,offshore_hdg);
  [dist2Shore,az2Shore] = distance_wgs84(nearestShore(2),nearestShore(1),xlats,xlons);

  zs = interp2(bath.ngdc_hires_bathy.lon,bath.ngdc_hires_bathy.lat,bath.ngdc_hires_bathy.field,xlons,xlats);

  [ig,ig,ig,stnsd,stnsaz] = get_nearest_points(bath.ngdc_hires_bathy,[stns.lon;stns.lat]);

  [xlons_f,xlats_f] = transect_wgs84(bath.lon,bath.lat,-3:0.05:1,offshore_hdg);
  [dist2Shore_f,az2Shore_f] = distance_wgs84(nearestShore(2),nearestShore(1),xlats_f,xlons_f);
  zs_f = interp2(bath.ngdc_hires_bathy.lon,bath.ngdc_hires_bathy.lat,bath.ngdc_hires_bathy.field,xlons_f,xlats_f);

end;




%{
disp('Pausing for reflection...');
pause(5);
Cmaxz = repmat(nan,[1,size(Ciso,2)]);
%Cmaxzix = repmat(nan,[1,size(Ciso,2)]);
Cmaxbeta = repmat(nan,[1,size(Ciso,2)]);
[LON,LAT] = meshgrid(bath.ngdc_hires_bathy.lon,bath.ngdc_hires_bathy.lat);
disp(size(Ciso,2));
for ix = 1:size(Ciso,2)
  if ( mod(ix,1000) == 1 ); toc, disp(ix); tic, end;
  [ig,coordix] = min(abs(Ciso(1,ix)-LON(:))+abs(Ciso(2,ix)-LAT(:)));
  [yix,xix] = ind2sub(size(bath.ngdc_hires_bathy.field),coordix);
  yixen = yix-sampleix:yix+sampleix;
  xixen = xix-sampleix:xix+sampleix;
  try,
    dat = bath.ngdc_hires_bathy.field(yixen,xixen);
    [Cmaxz(ix),ixix] = nanmax(dat(:));
    %???    Cmaxzix(ix) = ixix;
    dat = bath.ngdc_hires_bathy.beta(yixen,xixen);
    Cmaxbeta(ix) = nanmax(dat(:));
  catch,
  end;
end;
LON=[]; LAT=[]; clear LON LAT
toc,
disp('interp_field');
Cmaxbeta_ALT = interp_field(bath.ngdc_hires_bathy.lat,bath.ngdc_hires_bathy.lon,bath.ngdc_hires_bathy.beta,Ciso(2,:),Ciso(1,:),{@nanmax,sampleix});
tic,
%}





  [maxlon,rix] = max(LONS,[],2);
  [minlat,bix] = min(LATS,[],1);
  [maxlat,tix] = max(LATS,[],1);




  [ig,minlatix] = min(lats(:));
  [ig,maxlatix] = max(lats(:));
  [ig,minlonix] = min(lons(:));
  [ig,maxlonix] = max(lons(:));





[anomts,clim,tid,asd,lopct,hipct] = anomalize_ts(fwyf1.ndbc_sea_t,@get_jhour_no_leap);



    % if ( stix ~= ix )
      x = res.ps(stix);
      d = station_dist(x,x0);
      if ( d < 10.0 && abs(x.depth-x0.depth)<abs(x0.depth*0.20) )
        fprintf(1,'P: % 5s,% 5.1f,% 5.1f,%.3f,%.3f  % 5s,% 5.1f,% 5.1f,%.3f,%.3f  %.0f\n',...
                x0.station_name,x0.depth,x0.ngdc_depth,x0.ngdc_beta,x0.very_hires_smooth_beta,...
                x.station_name, x.depth, x.ngdc_depth, x.ngdc_beta, x.very_hires_smooth_beta,d*1e3);
      end;






    % if ( stix ~= ix )





        fprintf(1,'P: %s,%f,%f,%f,%f \t %s,%f,%f,%f,%f\n',...
                x0.station_name,x0.depth,x0.ngdc_depth,x0.ngdc_beta,x0.very_hires_smooth_beta,...
                x.station_name, x.depth, x.ngdc_depth, x.ngdc_beta, x.very_hires_smooth_beta);









        disp({'P: ',x0.station_name,x.station_name,x0.depth,x.depth,...
              x0.ngdc_beta,x.ngdc_beta,x0.very_hires_smooth_beta,x.very_hires_smooth_beta,});






      latres = min(diff(x.ngdc_hires_bathy.lat))*111e3;
      switch ( roundn(latres,1) ),
       case 10:
        N = 7;
        smoothopt = {@nanmean,12,12,18}; % Minimum 10%
       case 20:
        N = 5;
        smoothopt = {@nanmean,6,6,5}; % Minimum 14%
       case 30:
        N = 3;
        smoothopt = {@nanmean,4,4,3}; % Minimum 19%
       case 60:
        N = 3;
        smoothopt = {@nanmean,3,3,2}; % Minimum 22%
       case 90:
        N = 2;
        smoothopt = 'linear';
       otherwise:
        error('Unknown bathymetry resolution %f [m]',latres);
      end;






    if ( ~strncmpi(recalcSlope,'w',1) )
    else
    end; %if ( ~strncmpi(recalcSlope,'w',1) ) else




  if ( ~matLoaded )
  end;





        x.lon = stn.ms(stix).lon;
        x.lat = stn.ms(stix).lat;
        x = read_hires_bathymetry(x,[1e3,1e3],[],true);
        [sLAT,sLON] = meshgrid(x.lat,x.lon);
        sbat.field = bat.field(latixen,lonixen);
        sbat.beta = bat.beta(latixen,lonixen);
        sbat.(smoothvar) = interp_field(sbat.lat,sbat.lon,sbat.(fldvar),sLAT,sLON,smoothparms);
        sbat.(smoothvar) = reshape(sbat.(smoothvar),size(sbat.beta'))';
        sbat.(smoothvar)(sbat.field >= -1.0) = nan;
        stns.(fld).(smoothvar) = interp2(sbat.lon,sbat.lat,sbat.(smoothvar),stns.(fld).lon,stns.(fld).lat,'nearest');
        sbat=[]; sLAT=[]; sLON=[]; clear sbat sLAT sLON lonix lonixen latix latixen 
        stn.ms(stix).ngdc_hires_bathy = x.ngdc_hires_bathy;
        x=[]; clear x




function flds = read_hycom_gom_reanalysis()
%function flds = read_hycom_gom_reanalysis()
% CREATE data/GoMHYCOM.nc file with latitude and longitude from DODS (expt_32.5)

  datapath = get_ecoforecasts_path('data');

  % RESULTS FROM: http://tds.hycom.org/thredds/dodsC/GOMl0.04/expt_32.5/hrly.ascii?Latitude[0:1:384],Longitude[0:1:540]
  lats = [ ...
      18.091648, 18.129667, 18.167677, 18.205679, 18.243671, 18.281658, 18.319633, 18.357603, 18.395563, 18.433516, 18.471458, 18.509394, 18.54732, 18.585238, 18.623148, 18.661049, 18.698942, 18.736826, 18.774702, 18.81257, 18.85043, 18.888279, 18.92612, 18.963955, 19.00178, 19.039595, 19.077402, 19.115202, 19.15299, 19.190773, 19.228546, 19.26631, 19.304066, 19.341812, 19.379549, 19.417278, 19.455, 19.492712, 19.530415, 19.56811, 19.605795, 19.64347, 19.681139, 19.718798, 19.756447, 19.794088, 19.83172, 19.869343, 19.906958, 19.944563, 19.98216, 20.019747, 20.057325, 20.094894, 20.132456, 20.170008, 20.20755, 20.245083, 20.282608, 20.320122, 20.357628, 20.395126, 20.432613, 20.470093, 20.507563, 20.545023, 20.582474, 20.619915, 20.657349, 20.694773, 20.732187, 20.769592, 20.806988, 20.844374, 20.881752, 20.91912, 20.956478, 20.993828, 21.031168, 21.068499, 21.10582, 21.143133, 21.180435, 21.217728, 21.255013, 21.292286, 21.329552, 21.366806, 21.404053, 21.441288, 21.478516, 21.515734, 21.55294, 21.59014, 21.627329, 21.664507, 21.701677, 21.738836, 21.775988, 21.81313, 21.85026, 21.88738, 21.924494, 21.961596, 21.998688, 22.03577, 22.072844, 22.109907, 22.146961, 22.184006, 22.221039, 22.258064, 22.295078, 22.332083, 22.369078, 22.406063, 22.443039, 22.480003, 22.51696, 22.553905, 22.590841, 22.627768, 22.664682, 22.70159, 22.738485, 22.775372, 22.812248, 22.849113, 22.88597, 22.922815, 22.959652, 22.996479, 23.033295, 23.0701, 23.106897, 23.143682, 23.180458, 23.217224, 23.25398, 23.290726, 23.327461, 23.364185, 23.4009, 23.437605, 23.4743, 23.510984, 23.547659, 23.584324, 23.620977, 23.65762, 23.694254, 23.730877, 23.767488, 23.804092, 23.840683, 23.877266, 23.913837, 23.950397, 23.986948, 24.023489, 24.060019, 24.096539, 24.133047, 24.169546, 24.206034, 24.242512, 24.27898, 24.315437, 24.351883, 24.388319, 24.424746, 24.461159, 24.497564, 24.533958, 24.570341, 24.606714, 24.643076, 24.679428, 24.715769, 24.7521, 24.78842, 24.824728, 24.861027, 24.897314, 24.933592, 24.969858, 25.006115, 25.04236, 25.078594, 25.114817, 25.151031, 25.187233, 25.223425, 25.259605, 25.295774, 25.331934, 25.368082, 25.40422, 25.440348, 25.476463, 25.512568, 25.548662, 25.584745, 25.620817, 25.65688, 25.69293, 25.72897, 25.765, 25.801016, 25.837023, 25.87302, 25.909004, 25.944979, 25.980942, 26.016893, 26.052835, 26.088764, 26.124683, 26.160591, 26.19649, 26.232374, 26.26825, 26.304113, 26.339966, 26.375807, 26.411638, 26.447458, 26.483265, 26.519062, 26.554848, 26.590624, 26.626387, 26.662138, 26.69788, 26.73361, 26.76933, 26.805037, 26.840733, 26.876417, 26.91209, 26.947754, 26.983404, 27.019045, 27.054672, 27.09029, 27.125896, 27.161491, 27.197075, 27.232647, 27.268206, 27.303755, 27.339294, 27.37482, 27.410336, 27.445839, 27.48133, 27.516811, 27.55228, 27.587738, 27.623184, 27.658619, 27.694044, 27.729456, 27.764856, 27.800245, 27.835623, 27.870987, 27.906342, 27.941685, 27.977016, 28.012335, 28.047644, 28.082941, 28.118225, 28.153498, 28.18876, 28.22401, 28.259249, 28.294476, 28.32969, 28.364893, 28.400085, 28.435266, 28.470434, 28.50559, 28.540735, 28.575869, 28.61099, 28.646101, 28.681198, 28.716284, 28.751358, 28.78642, 28.821472, 28.856512, 28.891539, 28.926554, 28.96156, 28.99655, 29.03153, 29.066498, 29.101456, 29.136398, 29.171331, 29.206253, 29.241161, 29.276058, 29.310944, 29.345816, 29.380678, 29.415527, 29.450363, 29.48519, 29.520002, 29.554804, 29.589594, 29.62437, 29.659136, 29.69389, 29.72863, 29.763361, 29.798077, 29.832783, 29.867476, 29.902157, 29.936827, 29.971483, 30.006128, 30.04076, 30.075382, 30.109991, 30.144587, 30.17917, 30.213743, 30.248304, 30.282852, 30.317387, 30.351912, 30.386423, 30.420921, 30.455408, 30.489883, 30.524345, 30.558796, 30.593235, 30.62766, 30.662075, 30.696476, 30.730865, 30.765242, 30.799606, 30.83396, 30.8683, 30.902626, 30.936943, 30.971245, 31.005537, 31.039816, 31.074081, 31.108335, 31.142576, 31.176805, 31.211023, 31.245228, 31.279419, 31.313599, 31.347765, 31.38192, 31.416063, 31.450193, 31.48431, 31.518415, 31.552507, 31.586588, 31.620657, 31.65471, 31.688755, 31.722784, 31.756802, 31.790808, 31.8248, 31.858782, 31.892748, 31.926704, 31.960648 ...
         ];
  lons = [ ...
      -98.0, -97.96002, -97.91998, -97.880005, -97.839966, -97.79999, -97.76001, -97.71997, -97.67999, -97.640015, -97.599976, -97.56, -97.52002, -97.47998, -97.44, -97.400024, -97.359985, -97.32001, -97.28003, -97.23999, -97.20001, -97.160034, -97.119995, -97.08002, -97.03998, -97.0, -96.96002, -96.91998, -96.880005, -96.839966, -96.79999, -96.76001, -96.71997, -96.67999, -96.640015, -96.599976, -96.56, -96.52002, -96.47998, -96.44, -96.400024, -96.359985, -96.32001, -96.28003, -96.23999, -96.20001, -96.160034, -96.119995, -96.08002, -96.03998, -96.0, -95.96002, -95.91998, -95.880005, -95.839966, -95.79999, -95.76001, -95.71997, -95.67999, -95.640015, -95.599976, -95.56, -95.52002, -95.47998, -95.44, -95.400024, -95.359985, -95.32001, -95.28003, -95.23999, -95.20001, -95.160034, -95.119995, -95.08002, -95.03998, -95.0, -94.96002, -94.91998, -94.880005, -94.839966, -94.79999, -94.76001, -94.71997, -94.67999, -94.640015, -94.599976, -94.56, -94.52002, -94.47998, -94.44, -94.400024, -94.359985, -94.32001, -94.28003, -94.23999, -94.20001, -94.160034, -94.119995, -94.08002, -94.03998, -94.0, -93.96002, -93.91998, -93.880005, -93.839966, -93.79999, -93.76001, -93.71997, -93.67999, -93.640015, -93.599976, -93.56, -93.52002, -93.47998, -93.44, -93.400024, -93.359985, -93.32001, -93.28003, -93.23999, -93.20001, -93.160034, -93.119995, -93.08002, -93.03998, -93.0, -92.96002, -92.91998, -92.880005, -92.839966, -92.79999, -92.76001, -92.71997, -92.67999, -92.640015, -92.599976, -92.56, -92.52002, -92.47998, -92.44, -92.400024, -92.359985, -92.32001, -92.28003, -92.23999, -92.20001, -92.160034, -92.119995, -92.08002, -92.03998, -92.0, -91.96002, -91.91998, -91.880005, -91.839966, -91.79999, -91.76001, -91.71997, -91.67999, -91.640015, -91.599976, -91.56, -91.52002, -91.47998, -91.44, -91.400024, -91.359985, -91.32001, -91.28003, -91.23999, -91.20001, -91.160034, -91.119995, -91.08002, -91.03998, -91.0, -90.96002, -90.91998, -90.880005, -90.839966, -90.79999, -90.76001, -90.71997, -90.67999, -90.640015, -90.599976, -90.56, -90.52002, -90.47998, -90.44, -90.400024, -90.359985, -90.32001, -90.28003, -90.23999, -90.20001, -90.160034, -90.119995, -90.08002, -90.03998, -90.0, -89.96002, -89.91998, -89.880005, -89.839966, -89.79999, -89.76001, -89.71997, -89.67999, -89.640015, -89.599976, -89.56, -89.52002, -89.47998, -89.44, -89.400024, -89.359985, -89.32001, -89.28003, -89.23999, -89.20001, -89.160034, -89.119995, -89.08002, -89.03998, -89.0, -88.96002, -88.91998, -88.880005, -88.839966, -88.79999, -88.76001, -88.71997, -88.67999, -88.640015, -88.599976, -88.56, -88.52002, -88.47998, -88.44, -88.400024, -88.359985, -88.32001, -88.28003, -88.23999, -88.20001, -88.160034, -88.119995, -88.08002, -88.03998, -88.0, -87.96002, -87.91998, -87.880005, -87.839966, -87.79999, -87.76001, -87.71997, -87.67999, -87.640015, -87.599976, -87.56, -87.52002, -87.47998, -87.44, -87.400024, -87.359985, -87.32001, -87.28003, -87.23999, -87.20001, -87.160034, -87.119995, -87.08002, -87.03998, -87.0, -86.96002, -86.91998, -86.880005, -86.839966, -86.79999, -86.76001, -86.71997, -86.67999, -86.640015, -86.599976, -86.56, -86.52002, -86.47998, -86.44, -86.400024, -86.359985, -86.32001, -86.28003, -86.23999, -86.20001, -86.160034, -86.119995, -86.08002, -86.03998, -86.0, -85.96002, -85.91998, -85.880005, -85.839966, -85.79999, -85.76001, -85.71997, -85.67999, -85.640015, -85.599976, -85.56, -85.52002, -85.47998, -85.44, -85.400024, -85.359985, -85.32001, -85.28003, -85.23999, -85.20001, -85.160034, -85.119995, -85.08002, -85.03998, -85.0, -84.96002, -84.91998, -84.880005, -84.839966, -84.79999, -84.76001, -84.71997, -84.67999, -84.640015, -84.599976, -84.56, -84.52002, -84.47998, -84.44, -84.400024, -84.359985, -84.32001, -84.28003, -84.23999, -84.20001, -84.160034, -84.119995, -84.08002, -84.03998, -84.0, -83.96002, -83.91998, -83.880005, -83.839966, -83.79999, -83.76001, -83.71997, -83.67999, -83.640015, -83.599976, -83.56, -83.52002, -83.47998, -83.44, -83.400024, -83.359985, -83.32001, -83.28003, -83.23999, -83.20001, -83.160034, -83.119995, -83.08002, -83.03998, -83.0, -82.96002, -82.91998, -82.880005, -82.839966, -82.79999, -82.76001, -82.71997, -82.67999, -82.640015, -82.599976, -82.56, -82.52002, -82.47998, -82.44, -82.400024, -82.359985, -82.32001, -82.28003, -82.23999, -82.20001, -82.160034, -82.119995, -82.08002, -82.03998, -82.0, -81.96002, -81.91998, -81.880005, -81.839966, -81.79999, -81.76001, -81.71997, -81.67999, -81.640015, -81.599976, -81.56, -81.52002, -81.47998, -81.44, -81.400024, -81.359985, -81.32001, -81.28003, -81.23999, -81.20001, -81.160034, -81.119995, -81.08002, -81.03998, -81.0, -80.96002, -80.91998, -80.880005, -80.839966, -80.79999, -80.76001, -80.71997, -80.67999, -80.640015, -80.599976, -80.56, -80.52002, -80.47998, -80.44, -80.400024, -80.359985, -80.32001, -80.28003, -80.23999, -80.20001, -80.160034, -80.119995, -80.08002, -80.03998, -80.0, -79.96002, -79.91998, -79.880005, -79.839966, -79.79999, -79.76001, -79.71997, -79.67999, -79.640015, -79.599976, -79.56, -79.52002, -79.47998, -79.44, -79.400024, -79.359985, -79.32001, -79.28003, -79.23999, -79.20001, -79.160034, -79.119995, -79.08002, -79.03998, -79.0, -78.96002, -78.91998, -78.880005, -78.839966, -78.79999, -78.76001, -78.71997, -78.67999, -78.640015, -78.599976, -78.56, -78.52002, -78.47998, -78.44, -78.400024, -78.359985, -78.32001, -78.28003, -78.23999, -78.20001, -78.160034, -78.119995, -78.08002, -78.03998, -78.0, -77.96002, -77.91998, -77.880005, -77.839966, -77.79999, -77.76001, -77.71997, -77.67999, -77.640015, -77.599976, -77.56, -77.52002, -77.47998, -77.44, -77.400024, -77.359985, -77.32001, -77.28003, -77.23999, -77.20001, -77.160034, -77.119995, -77.08002, -77.03998, -77.0, -76.96002, -76.91998, -76.880005, -76.839966, -76.79999, -76.76001, -76.71997, -76.67999, -76.640015, -76.599976, -76.56, -76.52002, -76.47998, -76.44, -76.400024 ...
         ];

  fname = fullfile(datapath,'GoMHYCOM.nc');

  %netcdf.open(fname,'WRITE');
  %netcdf.reDef(ncid);
  ncid = netcdf.create(fname,'CLOBBER');

  latDimID = netcdf.defDim(ncid,'latitude',numel(lats));
  lonDimID = netcdf.defDim(ncid,'longitude',numel(lons));

  latVarID = netcdf.defVar(ncid,'latitude','double',latDimID);
  lonVarID = netcdf.defVar(ncid,'longitude','double',lonDimID);

  netcdf.endDef(ncid);
  netcdf.putVar(ncid,latVarID,lats);
  netcdf.putVar(ncid,lonVarID,lons);

  netcdf.close(ncid);

return;







1;
%% SCRIPT alexanders_mit_problem.m

N = 12;
light = repmat(0,[1,N]);

facs = {};
for ix=1:N
  facs{ix} = unique([1,factor(ix),ix]);
end;

tic,
for ix=1:N
  for ixix=ix:N;
    % if ( ismember(ix,facs{ixix}) )
    %   light(ixix) = not(light(ixix));
    % end;
    if ( isint(ixix/ix) )
      light(ixix) = not(light(ixix));
    end;
  end;
  %DEBUG:
  disp(light);
end;
toc,

disp([N,numel(find(light))]);





%{
  % Field-wide spatial mean Stokes/HYCOM percentage of each time step
  mspdpct = nanmean(spdpct(:,:),2);

  % Time steps where Stokes was the greatest field-wide percentage of HYCOM
  %[mx,mxix] = max(mspdpct);
  cutoff_mspdpct = prctile(mspdpct,99);
  hiix = find(mspdpct > cutoff_mspdpct);

  % Interpolate 4 km fields to 90 (or 30, or 10??) m bathymetry
  disp('Interpolating current fields');
  % HYCOM West-to-East
  zhu = interp3(hlat,hdts(hiix),hlon,hu(hiix,:,:),bat.lat,hdts(hiix),bat.lon);
  % HYCOM South-to-North
  zhv = interp3(hlat,hdts(hiix),hlon,hv(hiix,:,:),bat.lat,hdts(hiix),bat.lon);
  zhd = uv_to_dir_curr(zhu,zhv);
  % Stokes West-to-East
  zwu = interp3(hlat,hdts(hiix),hlon,wu(hiix,:,:),bat.lat,hdts(hiix),bat.lon);
  % Stokes South-to-North
  zwv = interp3(hlat,hdts(hiix),hlon,wv(hiix,:,:),bat.lat,hdts(hiix),bat.lon);
  zwd = uv_to_dir_curr(zwu,zwv);
  [d,n,m] = size(zhu);

  % Dot-product of A (bathymetric x/y gradient) and B (current u/v field)
  A(1,1:n,1:m) = beta_x;
  A(2,1:n,1:m) = beta_y;

  B = repmat(nan,[2,n,m]);
  hxs = repmat(nan,[d,n,m]);
  wxs = repmat(nan,[d,n,m]);

  disp(['Calculating ',num2str(d),' cross-shore current fields']);
  for tix=1:d
    disp(['Bh1 ',num2str(tix)]);
    B(1,:,:) = zhu(tix,:,:);
    disp(['Bh2 ',num2str(tix)]);
    B(2,:,:) = zhv(tix,:,:);
    disp(['hxs ',num2str(tix)]);
    hxs(tix,:,:) = squeeze(dot(A,B));
    disp(['Bw1 ',num2str(tix)]);
    B(1,:,:) = zwu(tix,:,:);
    disp(['Bw2 ',num2str(tix)]);
    B(2,:,:) = zwv(tix,:,:);
    disp(['wxs ',num2str(tix)]);
    wxs(tix,:,:) = squeeze(dot(A,B));
  end;

  A = []; clear A
  B = []; clear B
%}

%{
  % Find a good (individual) testing point
  ix = find(squeeze(~isnan(wxs(1,:,:))) & squeeze(~isnan(hxs(1,:,:))) & aspect_deg == 180,1)
  [yix,xix] = ind2sub(size(bat.field),ix), yixen=yix-2:yix+2; xixen=xix-2:xix+2;

  % Look at the lay of the land (or, um, seafloor)
  bat.field(yixen,xixen),
  beta_deg(yixen,xixen),

  % Go with the flow (look at the HYCOM and Stokes currents)
  % squeeze(zhu(1,yixen,xixen)),
  % squeeze(zhv(1,yixen,xixen)),
  squeeze(zhd(1,yixen,xixen)),
  % squeeze(zwu(1,yixen,xixen)),
  % squeeze(zwv(1,yixen,xixen)),
  % squeeze(zwd(1,yixen,xixen)),

  % How did we do on the dot products?
  squeeze(hxs(1,yixen,xixen)),
  % squeeze(wxs(1,yixen,xixen)),
%}









%{
  % Dot-product of A (bathymetric x/y gradient) and B (current u/v field)
  A(1,1:n,1:m) = ubeta_x;
  A(2,1:n,1:m) = ubeta_y;

  B = repmat(nan,[2,n,m]);
  uhxs = repmat(nan,[d,n,m]);
  uwxs = repmat(nan,[d,n,m]);

  disp(['Calculating ',num2str(d),' cross-shore current fields']);
  for tix=1:d
    if ( mod(tix,floor(d/10)) == 0 ); disp(tix); end;
    %disp(['Bh1 ',num2str(tix)]);
    B(1,:,:) = zhu(tix,:,:);
    %disp(['Bh2 ',num2str(tix)]);
    B(2,:,:) = zhv(tix,:,:);
    %disp(['uhxs ',num2str(tix)]);
    %uhxs(tix,:,:) = squeeze(dot(A,B));
    uhxs(tix,reef_line_ix) = squeeze(dot(A(:,reef_line_ix),B(:,reef_line_ix)));

    %disp(['Bw1 ',num2str(tix)]);
    B(1,:,:) = zwu(tix,:,:);
    %disp(['Bw2 ',num2str(tix)]);
    B(2,:,:) = zwv(tix,:,:);
    %disp(['uwxs ',num2str(tix)]);
    %uwxs(tix,:,:) = squeeze(dot(A,B));
    uwxs(tix,reef_line_ix) = squeeze(dot(A(:,reef_line_ix),B(:,reef_line_ix)));
  end;

  A = []; clear A
  B = []; clear B


  % Result 1: Median of cross-shore Stokes magnitude as percentage of
  % cross-shore HYCOM magnitude: 36%!
  nansummary(abs(uwxs)./(abs(uhxs))),
  % Result 2: 8% of the time, Stokes is onshore while HYCOM is offshore!
  numel(find(uwxs(:)<-0.0&uhxs(:)>0.0))./numel(uwxs),
  tix = find(nanmean(uwxs(:,:),2)<-0.000065,1);
  fmg;
  contourf(wlon,wlat,squeeze(wspd(tix,:,:))); colorbar;
  h(1) = quiver(wlon,wlat,squeeze(wu(tix,:,:)),squeeze(wv(tix,:,:))); set(h(1),'Color','r');
  h(2) = quiver(wlon,wlat,squeeze(hu(tix,:,:)),squeeze(hv(tix,:,:))); set(h(2),'Color','k');
  axis([-81.3269 -80.0855 24.5793 25.5640]); set(gca,'clim',[-0.01,0.01]);
  legend(h,'Stokes','HYCOM');
  [WLON,WLAT] = meshgrid(wlon,wlat);
  WLON = [WLON;WLON];
  WLAT = [WLAT;WLAT];
  l.wu = wu; l.wv = wv; l.hu = hu; l.hv = hv;
  maxspd = 0.10;
  bigix = find(abs(l.wu)>maxspd|abs(l.wv)>maxspd|abs(l.hu)>maxspd|abs(l.hv)>maxspd);
  l.wu(bigix) = nan; l.wv(bigix) = nan; l.hu(bigix) = nan; l.hv(bigix) = nan;
  U = [squeeze(l.wu(tix,:,:));squeeze(l.hu(tix,:,:))];
  V = [squeeze(l.wv(tix,:,:));squeeze(l.hv(tix,:,:))];
  fmg;
  contourf(hlon,hlat,squeeze(hspd(tix,:,:))); colorbar;
  h = quiver(WLON,WLAT,U,V,0); set(h,'Color','r');
  axis([-81.3269 -80.0855 24.5793 25.5640]);
  titlename(datestr(wdts(tix)));
%}










tic,
  zhx = repmat(nan,size(zhu));
  zhl = repmat(nan,size(zhu));
  zwx = repmat(nan,size(zhu));
  zwl = repmat(nan,size(zhu));
toc,
tic,
  disp(['Calculating ',num2str(d),' cross-shore current fields']);
  for tix=1:d
    if ( mod(tix,floor(d/10)) == 0 ); disp(tix); end;
    [zhx(tix,:,:),zhl(tix,:,:)] = reorient_vectors(ubeta_deg,squeeze(zhu(tix,:,:)),squeeze(zhv(tix,:,:)));
    [zwx(tix,:,:),zwl(tix,:,:)] = reorient_vectors(ubeta_deg,squeeze(zwu(tix,:,:)),squeeze(zwv(tix,:,:)));
  end;
toc,




tic,
  % [zhx,zhl] = reorient_vectors(ubeta_deg,zhu,zhv);
  [zhx,zhl] = reorient_vectors(ubeta_deg,squeeze(zhu(1,:,:)),squeeze(zhv(1,:,:)));
toc,


    [hx,hl] = reorient_vectors(ubeta_ang,hu,hv);
    [wx,wl] = reorient_vectors(ubeta_ang,wu,wv);


if 1;
tic,
  % [hx,hl] = reorient_vectors(ubeta_deg,hu,hv);
  [hx,hl] = reorient_vectors(ubeta_deg,squeeze(hu(1,:,:)),squeeze(hv(1,:,:)));
toc,
tic,
  [wx,wl] = reorient_vectors(ubeta_deg,wu,wv);
toc,
end;







  U(abs(U) > 0.20) = nan;
  V(abs(V) > 0.20) = nan;



  for rowix=kr:kr:n;
    tmpl = fld(rowix-kr+1:rowix,:);

for colix=kc:kc:m; newfld(rowix/kr,colix/kc) = method(tmpl(:,colix-kc+1:colix)); end;
keyboard;
    newfld(rowix/kr,:) = method(tmpl(:));
  end;



        linkedfldnms = cellstr(strrep(dispedflds,'_','\_'));



      ttl = ['BRUSH points and ''Del'' to remove from STN.',fldnm];
      titlename(ttl);
      if ( numel(linkedflds) > 0 )
        ttl = strvcat(ttl,'And fields:',linkedflds{:});
        % ttl = strcat(ttl,', ',linkedflds);
        % for linkedix = 1:numel(linkedflds);
        %   ttl = strcat(ttl,linkedflds{linkedix});
        % end;
        title(ttl);
      end;



        %ttl = [ttl,sprintf('\n And fields: %s,',linkedflds{:})];





    if ( ischar(fld) )
      % Do nothing
    elseif ( iscellstr(fld) )
      if ( ~isempty(regexp(fld{1},'^linked')) )
        linked_tag = fld{1};
        fld(1) = [];
        fldnm(1) = [];
        if ( numel(fld) == 1 )
          disp([fld,' had no fields to link']);
          fld = fld{1};
        else
          linkedflds = fld(2:end);
          if ( ~isempty(regexp(linked_tag,'display')) )
            dispedflds = linkedflds;
          end;
          fld = fld{1};
          disp([fld,' HAS LINKED FIELDS: ']);
          disp(linkedflds);
        end;
      else
        dispedflds = fld(2:end);
        fld = fld{1};
      end;
    else
      error('Invalid field sub-argument (#%d)??',fldix);
    end;








  linkedflds = {};



  disp('You may now wish to save your work, e.g., save(''SRVI2_portal_ALL_qc.mat'',''-struct'',''qc'',''-v7.3'')');




  % Calculate bathymetric gradient (i.e., to get cross-shore direction)
  [LON,LAT] = meshgrid(bat.lon,bat.lat);
  % [aspect_deg,slope_deg,beta_y,beta_x] = gradientm(LAT,LON,bat.field,wgs84Ellipsoid);
  % aspect_deg=[]; slope_deg=[]; clear aspect_deg slope_deg
  hx = 1e3 .* distance_wgs84(bat.lat(1),bat.lon(1),bat.lat(1),bat.lon(2));
  hy = 1e3 .* distance_wgs84(bat.lat(1),bat.lon(1),bat.lat(2),bat.lon(1));
  [beta_x,beta_y] = gradientn(bat.field,3,hx,hy);
  % hx=[]; hy=[]; clear hx hy
  % We want direction of steepest DESCENT - so invert the gradient
  %beta_x = -beta_x; beta_y = -beta_y;
  beta_x = beta_x; beta_y = -beta_y;
  beta = uv_to_spd(beta_x,beta_y);
  beta_deg = uv_to_dir_curr(beta_x,beta_y);








  [aspect,slope,dFdyNorth,dFdxEast] = gradientm(bat.lat,bat.lon,bat.field,wgs84Ellipsoid);



%m = stokes_drift(wflds.windspeed.field(:),wflds.sigwavehgt.field(:),wflds.primwaveper.field(:),'monismith',20);
%m = stokes_drift(wflds.windspeed.field(:),wflds.sigwavehgt.field(:),wflds.primwaveper.field(:),'monismith',5);
m = stokes_drift(wflds.windspeed.field(:),wflds.sigwavehgt.field(:),wflds.primwaveper.field(:),'monismith',3);



%fmg; contourf(wlon,wlat,squeeze(prctile(wspd,93)),[0:0.005:0.050]); colorbar;



u(isnan(u(:,:))) = [];
wu(isnan(wu(:,:))) = [];





v = wv + hflds.v.field;





wu = interp3(wlat,wdts,wlon,mu,hflds.u.lat,hflds.u.date,hflds.u.lon);





wv = interp3(wflds.windspeed.lat,wflds.windspeed.date,wflds.windspeed.lon,wflds.monismith_surface_drift_v.field,hflds.v.lat,hflds.v.date,hflds.v.lon);




    save(matfname,'flds');



          flds.u.date(tix,1) = dt + ((hrix-1)/24);
          flds.v.date(tix,1) = dt + ((hrix-1)/24);




read_gom_hycom_reanalysis




    flds.u.date = dts';



u = interp3(wflds.windspeed.date,wflds.windspeed.lat,wflds.windspeed.lon,wflds.monismith_surface_drift_u.field,uflds.u.date,uflds.u.lat,uflds.u.lon);



          catchwarn(miniME,...
                    ['REMOVING FLDS.(',fldnm,'): not found while processing ',ncfname]);



0,0.042,0.083,0.125,0.167,0.208,0.250,0.292, ...
         0.333,0.375,0.417,0.458,0.500,0.542,0.583,0.625, ...
         0.667,0.708,0.750,0.792,0.833,0.875,0.917,0.958];




function flds = read_hycom_gom_reanalysis()
%function flds = read_hycom_gom_reanalysis()

  datapath = get_ecoforecasts_path('data');

  chrs = {'0','0.042','0.083','0.125','0.167','0.208','0.250','0.292', ...
         '0.333','0.375','0.417','0.458','0.500','0.542','0.583','0.625', ...
         '0.667','0.708','0.750','0.792','0.833','0.875','0.917','0.958'};

  dts = datenum(2016,7,31):datenum(2016,9,30);
  for dt = dts(:)'
    [y,m,d] = datevec(dt);
    for chr = chrs;
      hr = chr{:};
      % SAMPLE FILE: ${ECOFORECASTS}/data/HYCOM_reGOM_0m_1993(08-1h0).nc
      fname = sprintf('HYCOM_reGOM_0m_%04d(%02d-%02d.%3f).nc',y,m,d,hr);

      nc = mDataset(fullfile(datapath,'hycom',fname));
      close(nc); clear nc
  end;

return;







  fld = [];
  if ( nargs >= 1 )
    if ( ~isempty(args{1})
      fld = args{1};
    end;
    args(1) = [];
    nargs = nargs - 1;
  end;
  if ( isempty(fld) )
    fld = fldstr.field;
  elseif ( ischar(fld) )
    fld = fldstr.(fld);
  elseif ( isnumeric(fld) )
    % Caller convenience: Remove singleton dimension(s)
    fld = squeeze(fld);
  end;
  % At this point, FLD needs to be a numeric matrix
  if ( ~isnumeric(fld) )
    error('FLD should either be empty, a CHAR fieldname, or a numeric matrix');
  end;





    if ( ~iscellstr(fnames) )
      error('Could not determine list of filenames from arguments!');
    end;




% RESponse VARiable
if ( ~exist('resvar','var') || isempty(resvar) )
  %resvar = 'fknms_seatemp';
  resvar = 'seatemp';
end;





  if ( ~exist('msg','var') )
    msg = '';
  end;




% M'eh - load waves no matter what
if ( ~isempty(findstr(ctlvar,'nwps')) )
  w = warning('OFF','Ecoforecasts:NWPS:NoFile');
  stns = get_nwps_stations(stns,[],[],false);
  %stns = get_nwps_stations(stns,[],[],true);
  warning(w);
end;





      for fldix = 1:numel(flds)
        fld = flds{fldix};
        nval = nc{fld}(:,:,:);
        val = interp_field(lat,lon,nval,stlats,stlons,'*linear');
        badix = find(all(isnan(val)));
        val(:,badix) = interp_field(lon,lat,nval,stlons(badix),stlats(badix),@nanmean);
        for tix = 1:numel(t)
          val(tix,1:numel(stlons)) = interp2(lon,lat,squeeze(nval(tix,:,:)),stlons,stlats,'*linear')';
          badix = find(all(isnan(val(tix,:))));
          val(tix,badix) = interp2(lon,lat,squeeze(nval(tix,:,:)),stlons(badix),stlats(badix),'nearest')';
        end;
        % Careful with our memory
        nval = []; clear nval;






        if ( numel(t) == 1 )
          % This should never happen, but...
          val = interp2(lon,lat,squeeze(nval),stlons,stlats,'*linear')';
        else
          % Weirdo dimension-ordering in MATLAB...
          val = interp3(lat,t,lon,nval,stlats,t,stlons,'linear');
          % UnMESH the result
          val = squeeze(val(:,1,:));
        end;





    Cs = interp2(lon,lat,nc{'currspeed'}(:,:,:),stlons,stlats,'*linear');
    Cd = interp2(lon,lat,nc{'currdir'}(:,:,:),stlons,stlats,'*linear');
    Ws = interp2(lon,lat,nc{'windspeed'}(:,:,:),stlons,stlats,'*linear');
    Wd = interp2(lon,lat,nc{'winddir'}(:,:,:),stlons,stlats,'*linear');
    Hs = interp2(lon,lat,nc{'sigwavehgt'}(:,:,:),stlons,stlats,'*linear');
    Sw = interp2(lon,lat,nc{'sigswellhgt'}(:,:,:),stlons,stlats,'*linear');
    Pp = interp2(lon,lat,nc{'primwaveper'}(:,:,:),stlons,stlats,'*linear');
    Pd = interp2(lon,lat,nc{'primwavedir'}(:,:,:),stlons,stlats,'*linear');
    Lu = interp2(lon,lat,nc{'ardhuin_surface_drift_u'}(:,:,:),stlons,stlats,'*linear');
    Lv = interp2(lon,lat,nc{'ardhuin_surface_drift_v'}(:,:,:),stlons,stlats,'*linear');






  n.lon = [];
  n.lat = [];
  for dtix = 1:numel(dts)
    dt = dts(dtix);
    [y,m,d,H,M,S] = datevec(dt);
    basefile = sprintf('%s_%04d%02d%02d_%02d00',dataset,y,m,d,H);
    clear y m d H M S;

    ncfname = fullfile(datapath,[basefile,'.nc']);

    datafile = ncfname;
    nc = mDataset(datafile);
    if ( isempty(lat) )
      n.lat = nc{'lat'}(:); 
      n.lon = nc{'lon'}(:); 
      nlat = numel(n.lat);
      nlon = numel(n.lon);
    end;
    hrs = nc{'time'}(:);
    nhrs = numel(hrs);
    n.t(end+1:end+nhrs,1) = datenum(1,1,1) + (hrs/24); 
    n.Cs(end+1:end+nhrs,1:n.nlat,1:n.nlon) = nc{'currspeed'}(:,:,:);
    n.Ws(end+1:end+nhrs,1:n.nlat,1:n.nlon) = nc{'windspeed'}(:,:,:); 
    n.Hs(end+1:end+nhrs,1:n.nlat,1:n.nlon) = nc{'sigwavehgt'}(:,:,:); 
    n.Sw(end+1:end+nhrs,1:n.nlat,1:n.nlon) = nc{'sigswellhgt'}(:,:,:); 
    n.Lg(end+1:end+nhrs,1:n.nlat,1:n.nlon) = nc{'ardhuin_surface_drift'}(:,:,:); 
    n.Lgm(end+1:end+nhrs,1:n.nlat,1:n.nlon) = nc{'monismith_surface_drift'}(:,:,:); 
    close(nc); clear nc
  end;





rot_stnms = randn(numel(stnms),1)*30;
text(ctlvars,resvars,trim_stnms,'Rotation',rot_stnms);
th = text(ctlvars,resvars,trim_stnms); set(th,'Rotation',rot_stnms);




    ts = subset_ts(stns.(stnm).(resvar),persub);
    if ( strcmpi(ctlvar,'depth') )
      ctlvars(end+1) = -stns.(stnm).(ctlvar);
    else
      ctlvars(end+1) = stns.(stnm).(ctlvar);
%%%% ??? HACK
if ( ~isnan(ctlvars(end)) ); ctlvars(end) = stns.(stnm).ngdc_beta; end;
%%%% ??? HACK
    end;





  %ctlvar = {'ngdc_beta',@nanmax,4,4,5};
  %ctlvar = {'ngdc_beta',@(x)(prctile(x,93)),4,4,5};



0:0.0050:0.30



  ntms = numel(nc{'time'}(:));




  ntms = numel(nc{'time'}(:));


titlename(['netCDF MAX Lg: ',datestr(floor(n.t(1)))]);





    mdl = 'Mex';
    cmp = 'M2';
    % cmp = 'K1';


    mydir = pwd;
    cd(tide_mpath);
    [umj,umn,uph,uinc] = ellipse(fullfile('DATA','Model_Mex'),stn.lat,stn.lon,cmp);
    % [umj,umn,uph,uinc] = ellipse(fullfile('DATA','Model_tpxo7.2'),stn.lat,stn.lon,cmp);
    % [umj,umn,uph,uinc] = ellipse(fullfile('DATA','Model_AO'),stn.lat,stn.lon,cmp);
    cd(mydir);




% dlon = distance_wgs84(stn.lat,stn.lon,stn.lat,stn.lon+1)*1e3;
% dlat = distance_wgs84(stn.lat,stn.lon,stn.lat+1,stn.lon)*1e3;
dlon = 111e3*cosd(stn.lat);
dlat = 111e3;





    if ( ~doAccDist )
      p07d = prctile(ts.data, 7) - mean(ts.data);
      p93d = prctile(ts.data,93) - mean(ts.data);
      if ( strcmpi(ctlvar,'depth') )
        ctlvars(end+1) = -stns.(fld).(ctlvar);
      else
        ctlvars(end+1) = stns.(fld).(ctlvar);
      end;
      resvars(end+1) = p93d - p07d;
    else
      p07d = []; p93d = [];
      pers = unique(peracc(ts.date));
      for ix=1:numel(pers);
        per = pers(ix);
        perix = find(peracc(ts.date)==per);
        p07d(ix) = prctile(ts.data(perix), 7) - mean(ts.data(perix));
        p93d(ix) = prctile(ts.data(perix),93) - mean(ts.data(perix));
      end; %for ix=1:numel(pers);
      clear ix
      if ( strcmpi(ctlvar,'depth') )
        ctlvars(end+1) = -stns.(fld).(ctlvar);
      else
        ctlvars(end+1) = stns.(fld).(ctlvar);
      end;
      resvars(end+1) = nanmax(p93d - p07d);
      pers = []; clear per perix pers
    end; %if ( ~doAccDist ) else




1;
% SCRIPT PLOT_RANGE_VS_CONTROL.m
%
% Compare scalar control variable CTLVAR (e.g., seafloor slope, 'ngdc_beta')
% vs. time series response variable RESVAR (e.g., 'fknms_seatemp') for all
% sites in dataset STRUCT STNS whose site names match PATT (v. GREPSTRUCT).
% DEFAULTS: Compare stns.FKNMS_.*.beta vs. whole-distribution Jun-Sep 93rd
% minus 7th PRCTILE ranges in time series stns.FKNMS_*.fknms_seatemp. If
% STRUCT STNS does not exist, call STNS = GET_LANGDON_THERMISTORS (v.)
%
% Last Saved Time-stamp: <Fri 2016-11-11 15:04:37 Eastern Standard Time gramer>


% Load dataset to process
if ( ~exist('stns','var') )
  stns = get_langdon_thermistors;
end;





  % If we are STILL missing betas, try the old tried and true...
  nanix = find(isnan(sites.depths) | isnan(sites.betas));
  if ( ~isempty(nanix) )
    lobath.lon = bath.lon;
    lobath.lat = bath.lat;
    lobath.rx = bath.rx;
    lobath.ry = bath.ry;
    lobath = read_hires_bathymetry(lobath,[lobath.rx,lobath.ry],[],false);
    lodepths = interp2(lobath.ngdc_hires_bathy.lon,lobath.ngdc_hires_bathy.lat,...
                       lobath.ngdc_hires_bathy.field,sites.lons,sites.lats,'*linear',nan);
    % Bathymetry resolution near Key West is 10 m: 9 pts. ~ 90 m
    [lobetas,lobeta_angs,loiso_angs,lobath.ngdc_lores_bathy] = ...
        find_ngdc_slope(lobath.ngdc_hires_bathy,sites.lons,sites.lats,2);
    
    sites.depths(nanix) = lodepths(nanix);
    sites.betas(nanix) = lobetas(nanix);
    sites.beta_angs(nanix) = lobeta_angs(nanix);
    sites.iso_angs(nanix) = loiso_angs(nanix);

    % lobath=[]; clear lobath
    clear ans lobetas lobeta_angs lodepths loiso_angs nanix 
  end;







sites.stnms(isnan(sites.depths))









      if ( ~exist('ctlvars','var') )
        ctlvars(1) = stns.(fld).(ctlvar);
        resvars(1) = nanmax(p93d - p07d);
        stnms{1} = fld;
      else
        ctlvars(end+1) = stns.(fld).(ctlvar);
        resvars(end+1) = nanmax(p93d - p07d);
        stnms{end+1} = fld;
      end;




      sresvar = ['MAX(',sresvar,'_{',char(peracc),'}'];







  clr = {'k','r',[0.0,0.5,0.0],'b','m'};
  mrk = {'.','o','x','+','*','s','d','p','^','v'};
  for ix=prmix(:)'; 
    clrix = mod(ix,numel(clr));
    mrkix = 1 + floor(ix/numel(clr));
    % Left-hand (7th percentile) cluster
    lhs(end+1) = plot(p07d(:,ix),ctlprms(ix)+([1:numel(pers)]./peroff),'.','Color',clr{ix});
    % Right-hand (93rd percentile) cluster
    plot(p93d(:,ix),ctlprms(ix)+([1:numel(pers)]./peroff),'.','Color',clr{ix});
    legs{end+1} = upper(tsnms{ix});
  end;





  % Subset first TS to the desired period
  tses(1) = subset_ts(tses(1),persub);

  % Subset all TSES to one another and the desired period
  tses = intersect_tses(tses);







  peracc = @get_jday; stracc='jday';
  %peracc = @get_year; stracc='year';

  persub = @ts_isfinite; strsub='all';
  %persub = @ts_jas; strsub='summer';
  %persub = @ts_jfm; strsub='winter';
  %persub = @ts_boreal_warm; strsub='rainy';
  %persub = @ts_boreal_cool; strsub='dry';


  ts1 = subset_ts(ts1,persub);

  tses = intersect_tses(ts1,ts2,ts3,ts4,ts5);








v = ver;
if ( strcmpi(v(1).Name,'Octave') )
  isOctave = true;
else
  isOctave = false;
end;
clear v






v = ver;
if ( strcmpi(v(1).Name,'Octave') )
  isOctave = true;





  if ( islogical(args{1}) || ( iscell(args{1}) && islogical(args{1}{1}) ) )
    if ( islogical(args{1}) )
      addqc = 'qc_';
    elseif ( args{1} && ischar(args{1}{2}) )
      addqc = args{1}{2};
    else
      error('If specified, ADDQC must be a Logical or a cell array with TRUE 1st element and CHAR 2nd element');
    end;
    args(1) = [];
    nargs = nargs - 1;
  end;






            if ( isnan(ngdc_deps(stnix)) )
              stns.(fld).ngdc_depth = -hires_deps(stnix);
              stns.(fld).ngdc_beta = hibetas(stnix);
              stns.(fld).ngdc_beta_ang = hibeta_angs(stnix);
              stns.(fld).ngdc_iso_ang = hiiso_angs(stnix);
            else
              stns.(fld).ngdc_depth = -ngdc_deps(stnix);
              stns.(fld).ngdc_beta = betas(stnix);
              stns.(fld).ngdc_beta_ang = beta_angs(stnix);
              stns.(fld).ngdc_iso_ang = iso_angs(stnix);
            end;





    [betas,beta_angs,beta_isos,bath] = find_ngdc_slope(bath,lons,lats,5);





clear ans C dfiles dts dtstr fid fnm ix ixen latc latst lgnms lonc lonst metafnm stnix



%{
  stubnm = strrep(stnms{ix},'FKNMS_','');
  dirnm = fullfile(langpath,stubnm,'data','0-data');
  if ( ~exist(dirnm,'dir') )
    warning('No data directory for %s',dirnm);
  else
    datafnms = dir(fullfile(dirnm,'FKNMS_*xls*'));
    for datix = 1:numel(datafnms)
      stnmix = find(strcmp(stnms,stnm));
      if ( isempty(stnmix) )
        warning('Unknown station for data file: %s',datafnms(datix).name);
      else
        fname = fullfile(dirnm,datafnms(datix).name);
        x = importdata(fname);
keyboard;
      end;
    end;
  end;
end;
%}



    %warning('Unknown station for data file: %s',fnm);





1;

stns = [];

langpath = fullfile('data','langdon');

% Process metadata precompiled by Lew.Gramer@noaa.gov
metafnm = fullfile(langpath,'metadata.xlsx');
x = importdata(metafnm);
ixen = 2:size(x.textdata,1);
stnms = x.textdata(ixen,1);
lgnms = x.textdata(ixen,2);
latst = x.textdata(ixen,3);
lonst = x.textdata(ixen,4);
% Calculated depth in [m]
deps = x.data(:,2);

% Parse stupid 'DDD-MM.ff' coordinate-string format
latc = double(split(latst,'-'));
lats = latc(:,1)+(latc(:,2)/60);

lonc = double(split(lonst,'-'));
% And convert to longitude West
lons = -lonc(:,1)+(lonc(:,2)/60);


% Process each data directory for which we had metadata
dfiles = dir(fullfile(langpath,'*','data','0-data','FKNMS_*.csv'));

for ix = 1:numel(dfiles)
  stnm = regexprep(dfiles(ix).name,'FKNMS_(.*)_WQDATA.*','$1');
  disp(stnm);
end;
%{
  stubnm = strrep(stnms{ix},'FKNMS_','');
  dirnm = fullfile(langpath,stubnm,'data','0-data');
  if ( ~exist(dirnm,'dir') )
    warning('No data directory for %s',dirnm);
  else
    datafnms = dir(fullfile(dirnm,'FKNMS_*xls*'));
    for datix = 1:numel(datafnms)
      stnmix = find(strcmp(stnms,stnm));
      if ( isempty(stnmix) )
        warning('Unknown station for data file: %s',datafnms(datix).name);
      else
        fname = fullfile(dirnm,datafnms(datix).name);
        x = importdata(fname);
keyboard;
        if ( ~isfield(stns,stnm) )
          stns.(stnm).lon = lons(stnmix);
          stns.(stnm).lat = lats(stnmix);
          stns.(stnm).depth = deps(stnmix);
        end;
      end;
    end;
  end;
end;
%}







% Process each data directory for which we had metadata
for ix = 1:numel(stnms)
  stubnm = strrep(stnms{ix},'FKNMS_','');
  dirnm = fullfile(langpath,stubnm,'data','0-data');
  if ( ~exist(dirnm,'dir') )
    warning('No data directory for %s',dirnm);
  else
    datafnms = dir(fullfile(dirnm,'FKNMS_*xls*'));
    for datix = 1:numel(datafnms)
      stnm = regexprep(datafnms(datix).name,'FKNMS_\(.*\)[.]xls.*','\1');
      stnmix = find(strcmp(stnms,stnm));
      if ( isempty(stnmix) )
        warning('Unknown station for data file: %s',datafnms(datix).name);
      else
        fname = fullfile(dirnm,datafnms(datix).name);
        x = importdata(fname);
keyboard;
        if ( ~isfield(stns,stnm) )
          stns.(stnm).lon = lons(stnmix);
          stns.(stnm).lat = lats(stnmix);
          stns.(stnm).depth = deps(stnmix);
        end;
      end;
    end;
  end;
end;






for ix = 1:numel(stnms)
  stnm = stnms{ix};
  stndirnm = strrep(stnm,'FKNMS_','');
  dirnm = fullfile(langpath,stndirnm,'data','0-data');
  if ( exist(dirnm) )
    metafnms = dir(fullfile(dirnm,'tblStation*xls');
    if ( isempty(metafnms) )
      warning('No metadata for station %s',stnm);
      fnamestubs = {['FKNMS_',stnm]};
    else
      metadata = importdata(fullfile(dirnm,metafnms(1).name));
    end;
    for metaix = 1:numel(fnamestubs);
      datafnms = dir(fullfile(dirnm,fnamestub));
      if ( ~isempty(datafnms) )
        if ( ~isfield(stns,stnm) )
          stns.(stnm).lon_str = lons{metaix};
          stns.(stnm).lat_str = lats{metaix};
        end;
        x = importdata(fullfile(dirnm,f.name));
      end;
    end;
  end;
end;







1;

stns = [];

langpath = fullfile('data','langdon');
for d = dir(langpath)'
  dirnm = fullfile(langpath,d.name,'data','0-data');
  if ( exist(dirnm) )
    stnm = d.name;
    metafnms = dir(fullfile(dirnm,'tblStation*xls');
    if ( isempty(metafnms) )
      warning('No metadata for station %s',stnm);
      fnamestubs = {['FKNMS_',stnm]};
    else
      metadata = importdata(fullfile(dirnm,metafnms(1).name));
    end;
    for metaix = 1:numel(fnamestubs);
      datafnms = dir(fullfile(dirnm,fnamestub));
      if ( ~isempty(datafnms) )
        stns.(stnm)
    for f = dir(dirnm)'
      if ( ~isempty(regexp(f.name,'FKNMS_.*.csv')) )




  h = evalc('!hostname');





if ( ~exist('n','var') || ~isfield(n,'Hs') )


  %g.t = cast(nc{'time'}(:),'double'); 







fmg; contourf(n.lon,n.lat,squeeze(nanmax(n.Lg)),0:0.005:0.035); colorbar; 




axis(bbox(:),[0,maxHs],[0,maxHs]); daspect([1,cosd(n.lat(1)),1]);
axis(bbox(:)); daspect([1,cosd(n.lat(1)),1]);




    datrows = x(2:end);



  if ( ~isempty(strfind(axis,'x')) )
    % Leave only the bottom AXES unchanged
    botx = min(pos(:,2));
    chgx = find( (pos(:,2) - botx) > 1e-6 );
    set(axs(chgx),'XTickLabel',[]);
  end;

  if ( ~isempty(strfind(axis,'y')) )
    % Leave only the leftmost AXES unchanged
    lefy = min(pos(:,1));
    chgy = find( (pos(:,1) - lefy) > 1e-6 );
    set(axs(chgy),'YTickLabel',[]);
  end;




    saxs = axs(ix(:,1));
    for axix = 2:numel(saxs)
      set(saxs(axix),'XTickLabel',[]);
    end;

    saxs = axs(ix(:,1));
    for axix = 2:numel(saxs)
      set(saxs(axix),'YTickLabel',[]);
    end;





dpfh = fmg;
dpax(1) = subplot(7,1,1:3); title('Sea temperatures');
plot_ts(pvgf1.ndbc_air_t,'k:',jack.deep.seatemp,jack.tcm2.seatemp,jack.tcm1.seatemp,jack.shallow.seatemp);
ylim([19.5,31.5]); xlim(jack.deep.seatemp.date([1,end])); 
grid on; set(gca,'XTickLabel','');
legend('Air','Deep','TCM2','TCM1','Shallow', legargs{:});

dpax(2) = subplot(7,1,4:5); title('Deep alongshore current');
plot_ts(jack.deep.adcp_v_sfc,jack.deep.adcp_v_mid,jack.deep.adcp_v_btm);
ylim([-1.5,+1.5]); xlim(jack.deep.seatemp.date([1,end])); 
grid on; center_axes; set(gca,'XTickLabel','');
legend('Surface','Mid-column','Bottom', legargs{:});

dpax(3) = subplot(7,1,6:7); title('Deep cross-shore current');
plot_ts(jack.deep.adcp_u_sfc,jack.deep.adcp_u_mid,jack.deep.adcp_u_btm);
ylim([-0.4,+0.4]); xlim(jack.deep.seatemp.date([1,end]));
grid on; center_axes;
legend('Surface','Mid-column','Bottom', legargs{:});








ylim([-7000,+7000]);





    res = {[],[],[],[],[],[],[],[]};





    %N = (62*24)+1;
    N = (31*24)+1;
    % Ensure N is odd
    N = N + ~mod(N,2);
    %smoothargs = {N,'loess'};
    smoothargs = {N,'moving'};

    % N = 31*24;
    % idx = linspace(clm.clim.date(1),clm.clim.date(end),N);
    % smth = interp1(clm.clim.date,clm.clim.data,idx);

    dt = median(diff(ts.date));
    idx = datenum(0,1,1):dt:(datenum(1,1,1)-dt+(dt/10));
    smth = interp1(clm.clim.date,clm.clim.data,idx);












    % Fastest
    smoothargs = {'lowess'};
    smoothargs = {'loess'};
    %% Too slow
    %smoothargs = {'rloess'};
    %% Too weird
    %smoothargs = {'sgolay'};






for hrs = 12:12:36
  for cderf={'lp','sum','avg'};
lkwf1 = verify_variable(lkwf1,'ndbc_ekman_flux_volume_12_h_lp');
lkwf1 = verify_variable(lkwf1,'ndbc_ekman_flux_volume_24_h_lp');
lkwf1 = verify_variable(lkwf1,'ndbc_ekman_flux_volume_36_h_lp');
lkwf1 = verify_variable(lkwf1,'ndbc_ekman_flux_volume_12_h_sum');
lkwf1 = verify_variable(lkwf1,'ndbc_ekman_flux_volume_24_h_sum');
lkwf1 = verify_variable(lkwf1,'ndbc_ekman_flux_volume_36_h_sum');

pvgf1 = verify_variable(pvgf1,'ndbc_ekman_flux_volume_12_h_lp');
pvgf1 = verify_variable(pvgf1,'ndbc_ekman_flux_volume_24_h_lp');
pvgf1 = verify_variable(pvgf1,'ndbc_ekman_flux_volume_36_h_lp');
pvgf1 = verify_variable(pvgf1,'ndbc_ekman_flux_volume_12_h_sum');
pvgf1 = verify_variable(pvgf1,'ndbc_ekman_flux_volume_24_h_sum');
pvgf1 = verify_variable(pvgf1,'ndbc_ekman_flux_volume_36_h_sum');

jack.deep = verify_variable(jack.deep,'adcp_l_36_h_lp');
jack.deep = verify_variable(jack.deep,'adcp_x_36_h_lp');
jack.deep = verify_variable(jack.deep,'adcp_x_36_h_avg');
jack.deep = verify_variable(jack.deep,'adcp_x_36_h_sum');




    % These can be back-calculated: Raw direction is hard to work with anyway
    res.nw_w_btm = rmfield(res.nw_w_btm,{'adcp_u','adcp_v','adcp_dir'});
    res.ne_buoy = rmfield(res.ne_buoy,{'adcp_u','adcp_v','adcp_dir'});
    res.nw_w_btm = rmfield(res.nw_w_btm,grepstruct(res.nw_w_btm,'_([uv]|dir)$'));
    res.ne_buoy = rmfield(res.ne_buoy,grepstruct(res.ne_buoy,'_([uv]|dir)$'));
    res.nw_w_btm = rmfield(res.nw_w_btm,grepstruct(res.nw_w_btm,'_[uv]_'));
    res.ne_buoy = rmfield(res.ne_buoy,grepstruct(res.ne_buoy,'_[uv]_'));




    res.nw_w_btm = rmfield(res.nw_w_btm,grepstruct(res.nw_w_btm,'_var'));





% Coast-line data
if ( ~exist('fwyf1','var') )
  fwyf1 = get_station_from_station_name('fwyf1'); fwyf1 = load_all_ndbc_data(fwyf1);
  fwyf1 = station_spddir_to_uv(fwyf1,'ndbc_wind1_speed','ndbc_wind1_dir');
  fwyf1 = station_reorient_vectors(fwyf1,'isobath_orientation','ndbc_wind1_u','ndbc_wind1_v','ndbc_wind1_x','ndbc_wind1_l');
end;
if ( ~exist('lkwf1','var') )
  lkwf1 = get_station_from_station_name('lkwf1'); lkwf1 = load_all_ndbc_data(lkwf1);
  lkwf1 = station_spddir_to_uv(lkwf1,'ndbc_wind1_speed','ndbc_wind1_dir');
  lkwf1 = station_reorient_vectors(lkwf1,'isobath_orientation','ndbc_wind1_u','ndbc_wind1_v','ndbc_wind1_x','ndbc_wind1_l');
end;
if ( ~exist('pvgf1','var') )
  pvgf1 = get_station_from_station_name('pvgf1'); pvgf1 = load_all_ndbc_data(pvgf1);
  pvgf1 = station_spddir_to_uv(pvgf1,'ndbc_wind1_speed','ndbc_wind1_dir');
  pvgf1 = station_reorient_vectors(pvgf1,'isobath_orientation','ndbc_wind1_u','ndbc_wind1_v','ndbc_wind1_x','ndbc_wind1_l');
end;
if ( ~exist('ftpf1','var') )
  ftpf1 = get_station_from_station_name('41114'); ftpf1 = load_all_ndbc_data(ftpf1);
  ftpf1 = station_spddir_to_uv(ftpf1,'ndbc_sigwavehgt','ndbc_avgwavedir','ndbc_wave_u','ndbc_wave_v');
  ftpf1 = station_reorient_vectors(ftpf1,'isobath_orientation','ndbc_wave_u','ndbc_wave_v','ndbc_wave_x','ndbc_wave_l');
end;
if ( ~exist('fdepk','var') )
  fdepk = read_fdep_stevens_data('k');
  fdepk = station_spddir_to_uv(fdepk,'fdep_wind_speed','fdep_wind_dir');
  fdepk = station_reorient_vectors(fdepk,'isobath_orientation','ndbc_wind1_u','ndbc_wind1_v','ndbc_wind1_x','ndbc_wind1_l');
end;
if ( ~exist('cnnf1','var') )
  cnnf1 = get_station_from_station_name('41113'); cnnf1 = load_all_ndbc_data(cnnf1);
  cnnf1 = station_spddir_to_uv(cnnf1,'ndbc_sigwavehgt','ndbc_avgwavedir','ndbc_wave_u','ndbc_wave_v');
  cnnf1 = station_reorient_vectors(cnnf1,'isobath_orientation','ndbc_wave_u','ndbc_wave_v','ndbc_wave_x','ndbc_wave_l');
end;
if ( ~exist('canf1','var') )
  canf1 = get_station_from_station_name('41009'); canf1 = load_all_ndbc_data(canf1);
  canf1 = station_spddir_to_uv(canf1,'ndbc_wind1_speed','ndbc_wind1_dir');
  canf1 = station_reorient_vectors(canf1,'isobath_orientation','ndbc_wind1_u','ndbc_wind1_v','ndbc_wind1_x','ndbc_wind1_l');
end;








  % High detail of Jack's FACE section
  if ( exist('jack','var') )
    jack.tcm2 = plot_hires_bathymetry(jack.tcm2,-[0:1:30,40:10:120],[3e3,9e3]);
    plot_all_se_florida_upwelling_sites;
    if doPrint; print('-dpng',fullfile(figspath,['upwelling-bathymetry-jack.png'])); end;
  end;








% Code BEFORE I wrote the convenience function ANNOTAXIS...
if ( doPlots )
  fhs=[];
  [lkwf1,ig,ig,ig,fhs(1)] = plot_hires_bathymetry(lkwf1,-[0:2:20,30:20:200]);
  titlename('SE Florida shelf - USGS Bathymetry - instrumentation');

  [sefcri.dc3,ig,ig,ig,fhs(end+1)] = plot_hires_bathymetry(sefcri.dc3,-[0:2:20,30:10:50]);
  titlename('SE Florida shelf (southern) - USGS Bathymetry - instrumentation');

  [fdepk,ig,ig,ig,fhs(end+1)] = plot_hires_bathymetry(fdepk,-[0:2:20,30:10:50]);
  titlename('SE Florida shelf (northern) - USGS Bathymetry - instrumentation');

  % [canf1,ig,ig,ig,fhs(end+1)] = plot_hires_bathymetry(canf1,-[0:2:20,30:10:50]);
  % titlename('Central Florida shelf - USGS Bathymetry - instrumentation');

  for fhix = 1:numel(fhs)
    figure(fhs(fhix));
    axis(axis);
    % textlat = -81.25;
    xl = get(gca,'XLim');
    textlat = xl(1) - 0.2;

    plot(fwyf1.lon,fwyf1.lat,'rp','MarkerFaceColor','w');
    text(textlat,fwyf1.lat,'FWYF1: met. \rightarrow','Horizontal','right');
    plot(lkwf1.lon,lkwf1.lat,'rp','MarkerFaceColor','w');
    text(textlat,lkwf1.lat,'LKWF1: met. \rightarrow','Horizontal','right');
    % plot(pvgf1.lon,pvgf1.lat,'rp','MarkerFaceColor','w');
    % text(textlat,pvgf1.lat,'PVGF1 \rightarrow','Horizontal','right');
    plot(fdepk.lon,fdepk.lat,'rp','MarkerFaceColor','w');
    text(textlat,fdepk.lat,'FDEP "K": met \rightarrow','Horizontal','right');
    plot(ftpf1.lon,ftpf1.lat,'rp','MarkerFaceColor','w');
    text(textlat,ftpf1.lat,'41114: Hs \rightarrow','Horizontal','right');

    plot(sefcri.dc1.lon,sefcri.dc1.lat,'rp','MarkerFaceColor','w');
    plot(sefcri.dc3.lon,sefcri.dc3.lat,'rp','MarkerFaceColor','w');
    %text(textlat,sefcri.dc3.lat,'SEFCRI DC3 \rightarrow','Horizontal','right');
    text(textlat,sefcri.dc3.lat,'SEFCRI DC1/3 \rightarrow','Horizontal','right');
    plot(sefcri.bc1.lon,sefcri.bc1.lat,'rp','MarkerFaceColor','w');
    plot(sefcri.bc3.lon,sefcri.bc3.lat,'rp','MarkerFaceColor','w');
    %text(textlat,sefcri.bc3.lat,'SEFCRI BC3 \rightarrow','Horizontal','right');
    text(textlat,sefcri.bc3.lat,'SEFCRI BC1/3 \rightarrow','Horizontal','right');

    plot(jack.shallow.lon,jack.shallow.lat,'rp','MarkerFaceColor','w');
    if ( fhix > 1 )
      plot(jack.tcm1.lon,jack.tcm1.lat,'rp','MarkerFaceColor','w');
      plot(jack.tcm2.lon,jack.tcm2.lat,'rp','MarkerFaceColor','w');
    end;
    plot(jack.deep.lon,jack.deep.lat,'rp','MarkerFaceColor','w');
    text(textlat,jack.deep.lat,'FACE HWD moorings: u,v \rightarrow','Horizontal','right');

    plot(sefcri.pb1.lon,sefcri.pb1.lat,'rp','MarkerFaceColor','w');
    plot(sefcri.pb2.lon,sefcri.pb2.lat,'rp','MarkerFaceColor','w');
    %text(textlat,sefcri.pb2.lat,'SEFCRI PB2 \rightarrow','Horizontal','right');
    text(textlat,sefcri.pb2.lat,'SEFCRI PB1/2 \rightarrow','Horizontal','right');

    plot(sefcri.mc2.lon,sefcri.mc2.lat,'rp','MarkerFaceColor','w');
    text(textlat,sefcri.mc2.lat,'SEFCRI MC2 \rightarrow','Horizontal','right');
    if ( fhix > 1 )
      plot(sefcri.pela.lon,sefcri.pela.lat,'rp','MarkerFaceColor','w');
      plot(sefcri.evca.lon,sefcri.evca.lat,'rp','MarkerFaceColor','w');
      plot(sefcri.slib.lon,sefcri.slib.lat,'rp','MarkerFaceColor','w');
      %text(textlat,sefcri.pela.lat,'SEFCRI PeLa \rightarrow','Horizontal','right');
      %text(textlat,sefcri.pela.lat,'SEFCRI MC2/PeLa \rightarrow','Horizontal','right');
    end;
    plot(sefcri.updb.lon,sefcri.updb.lat,'rp','MarkerFaceColor','w');
    text(textlat,sefcri.updb.lat,'SEFCRI UpDB \rightarrow','Horizontal','right');

    if ( exist('sfomc','var') )
      plot(sfomc.ne_buoy.lon,sfomc.ne_buoy.lat,'rp','MarkerFaceColor','w');
      %text(textlat,sfomc.ne_buoy.lat,'NSU E buoy \rightarrow','Horizontal','right');
      if ( fhix > 1 )
        plot(sfomc.c_buoy.lon,sfomc.c_buoy.lat,'rp','MarkerFaceColor','w');
        plot(sfomc.sw_buoy.lon,sfomc.sw_buoy.lat,'rp','MarkerFaceColor','w');
        plot(sfomc.se_buoy.lon,sfomc.se_buoy.lat,'rp','MarkerFaceColor','w');
        plot(sfomc.pier_cc.lon,sfomc.pier_cc.lat,'rp','MarkerFaceColor','w');
      end;
      plot(sfomc.nw_w_btm.lon,sfomc.nw_w_btm.lat,'rp','MarkerFaceColor','w');
      %text(textlat,sfomc.nw_w_btm.lat,'NSU NW/W mooring \rightarrow','Horizontal','right');
      text(textlat,sfomc.nw_w_btm.lat,'NSU moorings: u,v \rightarrow','Horizontal','right');
    end;

    if doPrint; print('-dpng',fullfile(figspath,['upwelling-bathymetry-',num2str(fhix),'.png'])); end;
  end; %for fhix = 1:numel(fhs)

end; %if ( doPlots )











  % [az,el] = view(ax);
  % if ( el == 90 )
  %   dim3 = false;
  % else
    dim3 = true;
  % end;





  if ( dim3 )
    text(xcoord,ycoord,zcoord,txt,align{:},args{:});
  else
    text(xcoord,ycoord,txt,align{:},args{:});
  end;














  switch ( lower(x_y_or_z) )
   case 'x',
disp('x');
    pos = get(get(ax,'XLabel'),'Pos')
    lim = get(ax,'YLim')
    xcoord = coord;
    ycoord = pos(1) - ( (lim(1) - pos(1)) / 2 )
    align = {'Vertical','top'};
   case 'y',
    pos = get(get(ax,'YLabel'),'Pos');
    lim = get(ax,'XLim')';
    xcoord = pos(1) - ( (lim(1) - pos(1)) / 2 );
    ycoord = coord;
    align = {'Horizontal','right'};
   case 'z',
   otherwise,
    error('X_Y_OR_Z should specify which axis to annotate (''x'',''y'', or ''z'')');
  end;









    lbl = 'YLabel';
    lim = 'XLim';







  if ( ischar(args{1}) )
    x_y_or_z = args{1};
    args(1) = [];
    if ( lower(x_y_or_z) ~= 'x' && lower(x_y_or_z) ~= 'y' && lower(x_y_or_z) ~= 'z' )
      error('X_Y_OR_Z should specify which axis to annotate (''x'',''y'', or ''z'')');
    end;
  else
    % DEFAULT: Annotate the Y-axis
    x_y_or_z = 'y';
  end;










  if 1;
    fmg;
    %spt(6,1,-2:5);
    %subplot(6,1,1:5);
    subplot(7,1,1:5);
    plot_ts(jack.deep.seatemp,jack.tcm2.seatemp,jack.tcm1.seatemp,jack.shallow.seatemp);
    legend('Deep','TCM2','TCM1','Shallow','Location','SouthEast');
    ylim([20,32]);
    xlim(datenum(2015,[8,10],1));

    %spt(6,1,6);
    %subplot(6,1,6);
    subplot(7,1,6:7);
    plot_ts(jack.deep.adcp_x,jack.tcm2.x,jack.tcm1.x,jack.shallow.adcp_x);
    legend('Deep','TCM2','TCM1','Shallow','Location','SouthEast');
    % ylim([-0.35,+0.35]);
    ylim([-0.50,+0.50]);
    xlim(datenum(2015,[8,10],1)); datetick3;

    if doPrint; print('-dpng',fullfile(figspath,'upwelling-jack-2015.png')); end;








    % textlat = -81.25;
    xl = get(gca,'XLim');
    textlat = xl(1) - diff(xl)*1.3;







  [lkwf1,ig,ig,ig,fhs(1)] = plot_hires_bathymetry(lkwf1,-[0:2:20,30:10:150],[30e3,120e3]);
  [sefcri.dc3,ig,ig,ig,fhs(end+1)] = plot_hires_bathymetry(sefcri.dc3,-[0:2:20,30:10:50],[30e3,50e3]);
  [fdepk,ig,ig,ig,fhs(end+1)] = plot_hires_bathymetry(fdepk,-[0:2:20,30:10:50],[30e3,50e3]);
  [canf1,ig,ig,ig,fhs(end+1)] = plot_hires_bathymetry(canf1,-[0:2:20,30:10:50],[30e3,50e3]);








jack.deep = verify_variable(jack.deep,'adcp_u_36_h_lp');
jack.deep = verify_variable(jack.deep,'adcp_u_36_h_avg');
jack.deep = verify_variable(jack.deep,'adcp_u_36_h_sum');







    fh = fmg;
    subplot(1,3,1:2);
    lkwf1 = plot_hires_bathymetry(lkwf1,-[0:1:20,30:10:150],[30e3,120e3],[],[],[],[],[],fh);
    subplot(1,3,3);
    sefcri.dc3 = plot_hires_bathymetry(sefcri.dc3,-[0:1:20,30:10:150],[10e3,35e3],[],[],[],[],[],fh);









  if 1;
    fmg;
    subplot(3,1,1); title('Far North (Martin County)');
     plot_ts(sefcri.updb.hourly_t,sefcri.mc2.hourly_t); 
     ylim([12,32]); grid on; 
     legend(['Offshore (',num2str(sefcri.updb.depth),' m)'],...
            ['Near-shore (',num2str(sefcri.mc2.depth),' m)']);
    subplot(3,1,2); title('Northern (Palm Beach)');
     plot_ts(sefcri.pb2.hourly_t,sefcri.pb1.hourly_t);
     ylim([12,32]); grid on; legend('Offshore','Near-shore');
     legend(['Offshore (',num2str(sefcri.updb.depth),' m)'],...
            ['Near-shore (',num2str(sefcri.mc2.depth),' m)']);
    subplot(3,1,3); title('Central (Ft. Lauderdale)');
     plot_ts(sefcri.bc3.hourly_t,sefcri.bc1.hourly_t);
     ylim([12,32]); grid on; legend('Offshore','Near-shore');
    if doPrint; print('-dpng',fullfile(figspath,'upwelling-sefcri-2010-2013.png')); end;
  end;






  {varargin{:}},
  disp('VS.');

  disp(newargs);




(global-unset-key [C-mouse-1])
(global-unset-key [C-M-mouse-1])
(global-unset-key [C-S-mouse-1])
(global-unset-key [M-mouse-1])
(global-unset-key [M-S-mouse-1])
(global-unset-key [S-mouse-1])

(global-unset-key [C-down-mouse-1])
(global-unset-key [C-M-down-mouse-1])
(global-unset-key [C-S-down-mouse-1])
(global-unset-key [M-down-mouse-1])
(global-unset-key [M-S-down-mouse-1])
(global-unset-key [S-down-mouse-1])

(setq mode-line-coding-system-map nil mode-line-column-line-number-mode-map nil mode-line-input-method-map nil)

(fset mouse-buffer-menu nil)

plotspec_check({'Color',[0,0,0],'ro-.'});
plotspec_check('Color',[0,0,0],'ro-.');
plotspec_check('Color',[0,0,0],'ro-.r');
plotspec_check('Color',[0,0,0],'o-.r');
plotspec_check('Color',[0,0,0],'r-.r');
plotspec_check('Color',[0,0,0],'-.r');
plotspec_check('Color',[0,0,0],'r-.');
plotspec_check('Color',[0,0,0],'r-..');
plotspec_check('Color',[0,0,0],'r-.r');
plotspec_check('Color',[0,0,0],'k-.');
plotspec_check('Color',[0,0,0],'k.');





      % If caller gave a plot spec with a line-style in it
      ls = regexprep(arg,'[^-]*(-[.]|--|-|:)[^-]*','$1','once');
      if ( ~isempty(ls) )
        newargs(end+1:end+2) = {'LineStyle',ls};
        % Delete line style, so '.' in '-.' doesn't mistakenly also become a Marker 
        arg = strrep(arg,ls,'');
      end;







  if ( isempty(args) )
    error('First argument must either be a CHAR or a cell array, e.g., VARARGIN');
  end;



  if ( ~iscell(args) )
    error('Argument should be a cell array, e.g., VARARGIN');
  end;
  






      % Special handling of stupid dash-dot
      lsix = strfind(arg,'-.');
      if ( ~isempty(lsix) )
        arg(lsix:lsix+1) = [];
        ts1_args = { 'LineStyle','-.',ts1_args{:} };
      end;








%%%%%%%%%%%%%%%%%%%%
% INTERNAL FUNCTION
%%%%%%%%%%%%%%%%%%%%

function baddts = read_fdep_stevens_data_qc_fld(res,fld,minv,maxv)
  if ( isfield(res,'fdep_battv') )
    baddts = union(baddts,res.fdep_battv.date(11>res.fdep_battv.data | res.fdep_battv.data>14));
  end;
return;






plot_ts(ax(1),jack.shallow.seatemp,'-','Color',[0,0.5,0]);



 || ~isfield(jack,'shallow') 




function savedirs

  global myDirStack;

%   % If there's nothing in there, load our last save first...
%   if ( ~exist('myDirStack', 'var') || isempty(myDirStack) )
%     loaddirs;
%   end;

  % This may be set in the call to INITDIRS (v.), but if not, default to
  % saving the directory stack in a .MAT file in the same directory where we
  % keep this code, i.e., in the Ecoforecasts toolkit directory
  if ( ~exist('myDirStackFilename','var') || isempty(myDirStackFilename) )
    % myDirStackFilename = 'c:/Documents and Settings/lew.gramer/My Documents/MATLAB/savedirs.mat';
    %[pathroot, ig, ig, ig] = fileparts(mfilename('fullpath'));
    [pathroot, ig, ig] = fileparts(mfilename('fullpath'));
    myDirStackFilename = fullfile(pathroot,'savedirs.mat');
  end;

  % dirsfname = 'c:/Documents and Settings/lew.gramer/My Documents/MATLAB/savedirs.mat';
  %[pathroot, ig, ig, ig] = fileparts(mfilename('fullpath'));
  [pathroot, ig, ig] = fileparts(mfilename('fullpath'));
  dirsfname = fullfile(pathroot,'savedirs.mat');
  %save(dirsfname, 'myDirStack');
  %save(dirsfname, 'myDirStack','-v7.3');
  save(dirsfname, 'myDirStack','-v7');

return;







%%%%%%%%%%
%%% FROM my_startup.m circa 2005


% PHOME = '\\cygnus\gramer\home\matlab\';

% path(path,[PHOME]);

% % Like UNIQUE(), but more powerful and with fuzzy tolerance
% path(path,[PHOME,'consolidator/']);
% path(path,[PHOME,'ecoforecasts/']);
% path(path,[PHOME,'ecoforecasts/coast/']);
% path(path,[PHOME 'EOFs/']);


% path(path,[PHOME,'ADCP/']);
% path(path,[PHOME,'Advstats']);
% path(path,[PHOME,'AirSea/']);
% path(path,[PHOME,'Arfit/']);
% path(path,[PHOME,'Bsplines']);
% path(path,[PHOME,'Claudia/']);
% path(path,[PHOME,'climatology']);
% path(path,[PHOME,'CODAS3/']);
% path(path,[PHOME,'compstats']);
% path(path,[PHOME,'Colormap']);
% path(path,[PHOME,'coral/']);
% path(path,[PHOME,'CTDCalib/']);
% path(path,[PHOME,'DataSet/']);
% path(path,[PHOME,'DateTime']);
% path(path,[PHOME,'Dsatbx/']);
% path(path,[PHOME,'Dvt/']);
% path(path,[PHOME,'DynModes/']);
% path(path,[PHOME,'EarthSci/']);
% path(path,[PHOME,'Epic/']);
% path(path,[PHOME,'Eps/']);
% path(path,[PHOME,'Even/']);
% path(path,[PHOME,'Fmex/']);
% path(path,[PHOME,'Geodetics/']);
% path(path,[PHOME,'Geography']);
% path(path,[PHOME,'gkslib/']);
% path(path,[PHOME,'Glmlab/']);
% path(path,[PHOME,'Graphics']);
% path(path,[PHOME,'Hist']);
% path(path,[PHOME,'HtmlTool']);
% path(path,[PHOME,'ifmObana/']);
% path(path,[PHOME,'Imputation/']);
% path(path,[PHOME,'Integration']);
% path(path,[PHOME,'IO/']);
% path(path,[PHOME,'MClasses']);
% path(path,[PHOME,'MClasses/Generic']);
% path(path,[PHOME,'MClasses/Interval']);
% path(path,[PHOME,'MClasses/Matfile']);
% path(path,[PHOME,'MClasses/Mode']);
% path(path,[PHOME,'MClasses/Mooring']);
% path(path,[PHOME,'MClasses/Vos']);
% path(path,[PHOME,'MClasses/batch']);
% path(path,[PHOME,'MClasses/presto']);
% path(path,[PHOME,'MClasses/time_series']);
% path(path,[PHOME,'MClasses/Vivace/']);
% path(path,[PHOME,'MClasses/Profile']);
% path(path,[PHOME,'Makehelp/']);
% path(path,[PHOME,'mat2html/']);
% path(path,[PHOME,'matlab_netcdf_5_0/']);
% path(path,[PHOME,'mexcdf60/']);
% path(path,[PHOME,'MexEPS/']);

% %%%% ??? MAYBE NEEDED AFTER ALL
% path(path,[PHOME,'mexnc']);

% path(path,[PHOME,'Misc/']);
% path(path,[PHOME,'m_map/']);
% path(path,[PHOME,'m_map/m_namebox/']);
% path(path,[PHOME,'MySQL']);
% path(path,[PHOME,'netcdf']);
% path(path,[PHOME,'netcdf/ncfiles']);
% path(path,[PHOME,'netcdf/nctype']);
% path(path,[PHOME,'netcdf/ncutility']);

% %%%% ??? MAYBE NEEDED AFTER ALL
% path(path,[PHOME,'netcdf_toolbox']);
% path(path,[PHOME,'netcdf_toolbox/netcdf']);
% % NOTE: This path differs from that mentioned in README... (LGramer, 2005-06-29)
% path(path,[PHOME,'netcdf_toolbox/netcdf/ncsource']);
% path(path,[PHOME,'netcdf_toolbox/netcdf/nctype']);
% path(path,[PHOME,'netcdf_toolbox/netcdf/ncutility']);

% path(path,[PHOME,'NumMethods/']);
% path(path,[PHOME,'Nurbs/']);
% path(path,[PHOME,'Oceanography']);
% path(path,[PHOME,'Oceanography/tsplot']);
% path(path,[PHOME,'Oceanography/Hydrobase']);
% path(path,[PHOME,'OMP2/']);
% path(path,[PHOME,'omviz/']);
% path(path,[PHOME,'OPNML/']);
% path(path,[PHOME,'PCCruiseCD/']);
% path(path,[PHOME,'RDADCP/']);
% path(path,[PHOME,'Regularization/']);
% path(path,[PHOME,'R12/']);
% path(path,[PHOME,'Smooth/']);
% path(path,[PHOME,'snackbar']);
% % Old tools for netCDF reading and creation
% path(path,[PHOME,'snctools']);
% path([PHOME,'Statistics'],path);
% path(path,[PHOME,'Strings']);
% path(path,[PHOME,'SaGA/']);
% path(path,[PHOME,'Statbx40']);
% path(path,[PHOME,'Spatial']);
% path(path,[PHOME,'Splines/']);
% path(path,[PHOME,'Staplot/']);
% path(path,[PHOME,'Steger/']);

%%%% ??? NOTE: This 'Sw' is older, but includes additional components from
%%%% ??? NOTE:  URI for sound-velocity and deep-ocean density calculations.
% path(path,[PHOME,'Sw/']);

% path(path,[PHOME,'timeplt5/']);
% path(path,[PHOME,'TimeSeries/']);
% path(path,[PHOME,'Todd/']);
% path(path,[PHOME,'transports']);
% path(path,[PHOME,'WaveCov/']);
% path(path,[PHOME,'WaveLab802/']);
% path(path,[PHOME,'Wavelets/']);
% path(path,[PHOME,'wetcdf/']);
% path(path,[PHOME,'XmlStuff/']);

% addutils;
%%%% ??? NOTE: Commenting out ADDUTILS above actually keeps the following
%%%% ??? NOTE: *TWELVE* utility directories from being added to the path.
% datautil
% fileutil
% graphutil
% imgutil
% mathutil
% matutil
% numutil
% polyutil
% statutil
% strutil
% sysutil
% timeutil

clear PHOME


% %%%% ??? MAYBE NEEDED AFTER ALL
% %
% % netCDF initialization
% %
% global nctbx_options;
% nctbx_options.theAutoscale = 1;
% nctbx_options.theAutoNaN = 1;


% %
% % Keep tabs on the current path on this host...
% %
% fid = fopen('~/.matlab.path', 'w+');
% if ( fid > 0 )
%     fprintf(fid, '%s\n', path);
%     fclose(fid);
% end
% clear fid;



%%% FROM my_startup.m circa 2005
%%%%%%%%%%










    text(WS(end,floor(end/2)),HS(end,floor(end/2)),z(end,floor(end/2)),[num2str(per),' s per.']);





      if ( ~exist(fpath,'file') || doOverwrite )







%function stn = station_ekman_flux(stn_or_stnm,orifld,afld,qfld,wfld,dfld,pfld,taufld,ffld)







  % Add datestrings to X labels
  datetick(axen(1),'x',2,'keeplimits','keepticks');
  if ( exist('set_datetick_cursor','file') )
    set_datetick_cursor;
  end;












HACK???? ABORTED EDITING ATTEMPT!


function lh = plotyy_ts(varargin)
%function lh = plotyy_ts([ax,][ts1_1,...],[ts2_1,...],...)
%
% PLOTYY_TS is identical to PLOTYY (v.) except that values to be plotted are
% not given as separate X and Y arguments, but are instead drawn from *time
% series structs* TS1,... which are expected to have .date and .data fields.
% First arg may specify the AXES to plot in, as for PLOTYY; otherwise, any
% argument which is NOT a time series struct (or optionally a *vector* of
% time series structs) is passed on as an arg to PLOTYY.
%
% Calls DATETICK (DATETICK3 if present) to display X labels as date strings.
%
% SAMPLE CALL #1:
%  >> % Plot Molasses Reef air temperature time series on one AXES, thick blue line,
%  >> % Sombrero Key sea temperature on a separate AXES on the same figure panel.
%  >> plotyy_ts(mlrf1.ndbc_air_t,smkf1.ndbc_sea_t);
%
% SAMPLE CALL #2:
%  >> % Plot MISST sea temperature in all stations in struct vector SITES;
%  >> % Rely on PLOTYY_TS to automatically distinguish by Color and Marker
%  >> plotyy_ts(sites.misst); legend({sites.station_name});
%
% Last Saved Time-stamp: <Fri 2016-08-12 12:45:03 Eastern Daylight Time gramer>

  argix = 1;
  if ( ishandle(varargin{1}) )
    ax = varargin{1};
    argix = argix + 1;
  else
    ax = gca;
  end;

  co = get(ax,'ColorOrder');
  ms = '.ox+*sdv^<>ph';
  ncos = size(co,1);

  % Extract time series structs to plot, and plot args for each TS
  ts = {};
  pargs = {};
  nts = [0,0];
  plotfuns = [0,0];
  while ( argix <= nargin )
    arg = varargin{argix};
    if ( isstruct(arg) && isfield(arg,'date') && isfield(arg,'data') )
      if ( nts(1) == 0 )
        tsix = 1;
      elseif ( nts(2) == 0 )
        tsix = 2;
      else
        error('Only TWO args may be TSes or arrays of TSes!');
      end;
      % Gracefully deal with vectors OR MATRICES of structs passed as a single argument
      for strcix = 1:numel(arg)
        strc = arg(strcix);
        if ( ~is_valid_ts(strc) )
          error('STRUCT args must be valid time series!');
        end;
        nts(tsix) = nts(tsix) + 1;
        ts{tsix,nts(tsix)} = strc;
        pargs{tsix,nts(tsix)} = {};
      end;
    elseif ( nts(1) == 0 )
      error('First arg after optional AXES handle must be a time series struct or matrix of TS structs!');
    elseif ( isa(arg,'function_handle') )
      if ( nts(1) == 0 )
        error('Specify at least one TS before specifying a plotting function!');
      else ( ~is(plotfuns(tsix),'function_handle') )
        plotfuns(tsix) = arg;
      else
        error('More than one plotting function given for TS #%g!',tsix);
      end;
    else
      pargs{tsix,nts(tsix)} = {pargs{tsix,nts}{:} arg};
    end;
    argix = argix + 1;
  end;

  % Plot the first of each of the two arrays of time series
  if ( ~is(plotfuns(1),'function_handle') )
    plotfuns(1) = @plot;
  end;
  if ( ~is(plotfuns(2),'function_handle') )
    plotfuns(2) = @plot;
  end;
  [axen,lh(1,1),lh(2,1)] = plotyy(ax,ts{1,1}.date,ts{2,1}.data,plotfuns(1),plotfuns(2));


  % Try to ensure that every plot has a distinct plot spec
  for tsix = 1:2
    for tix = 1:nts(tsix)
      if ( isempty(pargs{tsix,tix}) )
        % Use default color order if user did not give a plot spec (see PLOT)
        % If user passed in >ncos Time Series, also distinguish them by Marker


LEFT OFF EDITING HERE!


        pargs{tsix,tix} = {pargs{tsix,tix}{:} 'Color',co(mod(tix-1,ncos)+1,:), 'Marker',ms(ceil(tix/ncos))};
    else
      % Plot spec, if given, must be the first plot arg for this TS
      arg = pargs{tix}{1};
      % The following mess is due to PLOT's ridiculously flexible calling sequence!
      if ( ~ischar(arg) || length(arg) > 4 )
        % User did not give a plot spec
        pargs{tix} = {'Color',co(mod(tix-1,ncos)+1,:) pargs{tix}{:}};
      elseif ( ~isempty(regexp(arg,'^[.ox+*sdv^<>ph:-]*[bgrcmykw][.ox+*sdv^<>ph:-]*$')) )
        % User gave a plot spec with a color in it - do nothing
      elseif ( ~isempty(regexp(arg,'^[.ox+*sdv^<>ph:-]*$')) )
        % User gave a plot spec WITHOUT a color in it
        pargs{tix} = {pargs{tix}{1} 'Color',co(mod(tix-1,ncos)+1,:) pargs{tix}{2:end}};
      else
        % User did not give a plot spec
        pargs{tix} = {'Color',co(mod(tix-1,ncos)+1,:) pargs{tix}{:}};
      end;
    end;
  end;

  % Turn on HOLD
  hold_state = ishold(ax);
  hold('on');
  % Plot each line
  %DEBUG:  disp(nts);
  for tix = 2:nts
    lh(tix) = plot(ax,ts{tix}.date,ts{tix}.data,pargs{tix}{:});
  end;
  % Restore previous HOLD state
  if ( ~hold_state )
    hold('off');
  end;

  % Add datestrings to X labels
  datetick('x',2,'keeplimits','keepticks');
  if ( exist('set_datetick_cursor','file') )
    set_datetick_cursor;
  end;

  % Work around bug in DATETICK2/DATETICK3 that shuffles current AXES
  axes(ax);

return;













    elseif ( nts(1) == 0 )
      error('First arg after optional AXES handle must be a time series struct or matrix of TS structs!');
    elseif ( isa(arg,'function_handle') )
      if ( nts(2) == 0 )
        error('Specify at least two TSes before specifying a plotting function!');
      elseif ( ~exist(plotfun1,'var') )
        plotfun1 = arg;
      elseif ( ~exist(plotfun2,'var') )
        plotfun2 = arg;
      else
        error('More than two plotting functions were given!');
      end;
    else
      pargs{tsix,nts(tsix)} = {pargs{tsix,nts}{:} arg};
    end;












        if ( tsix == 1 )
          nts1 = nts1 + 1;
        else
          nts2 = nts2 + 1;
        end;








stnm = 'lkwf1';

stn = get_station_from_station_name(stnm); stn = load_all_ndbc_data(stn);

stn = station_dewp_to_relhumid(stn,'ndbc_air_t','ndbc_dew_t','ndbc_relhumid');
stn = station_relhumid_to_spechumid(stn,'ndbc_air_t','ndbc_relhumid','ndbc_spechumid');

stn = station_bulk_windstress(stn,'ndbc_windstress','ndbc_wind1_speed',[],...
                              'ndbc_wind1_dir','ndbc_air_t','ndbc_relhumid','ndbc_barom');

stn.ekman_flux.date = stn.ndbc_windstress.date;
stn.ekman_flux.data = stn.ndbc_windstress.data ./ sw_f(stn.lat);







      res = [];

      % These query results are generally time-inverted
      [dts,sortix] = sort(dts);
      % Easier to do this check for each month individually
      dtix = find(mindt <= dts & dts <= maxdt);
      dts = dts(dtix);
      sortix = sortix(dtix);

      for fldix=2:numel(flds)
        fld = flds{fldix};
        if ( ~isempty(fld) )
          res.(fld).date = dts;
          res.(fld).data = C{fldix}(sortix);
        end;
      end; %for fldix=2:numel(flds)





function stn = station_ekman_flux(stn_or_stnm,afld,qfld,wfld,dfld,pfld,taufld,ffld,orifld)
%function stn = station_ekman_flux(stn_or_stnm,afld,qfld,wfld,dfld,pfld,taufld,ffld,orifld)
%
% Calculate cross-shore component of Ekman flux from surface wind stress.
%

  stn = get_station_from_station_name(stn_or_stnm);
  clear stn_or_stnm;

  if ( ~exist('afld','var') )	afld = 'ndbc_air_t';			end;
  if ( ~exist('qfld','var') )	qfld = {'ndbc_dew_t','ndbc_relhumid'};	end;
  if ( ~exist('wfld','var') )	wfld = 'ndbc_wind1_speed';		end;
  if ( ~exist('dfld','var') )	dfld = 'ndbc_wind1_dir';		end;
  if ( ~exist('pfld','var') )	pfld = 'ndbc_barom';			end;

  if ( ~exist('taufld','var') )
    taufld = 'ndbc_bulk_windstress';
  end;
  tauxfld = [taufld,'_x'];
  tauyfld = [taufld,'_y'];
  tauxsfld = [taufld,'_xs'];
  taulsfld = [taufld,'_ls'];

  if ( ~exist('ffld','var') )
    ffld = 'ndbc_ekman_flux';
  end;
  if ( ~exist('orifld','var') )
    orifld = 'isobath_orientation';
  end;

  dewfld = '';
  if ( iscell(qfld) )
    dewfld = qfld{1};
    qfld = qfld{2};
  end;










  if ( ~exist('dts','var') || isempty(dts) )
    [ig,dts] = get_fdep_stevens_metadata(stn);
    ig = []; clear ig
  end;






stn = station_relhumid_to_spechumid(stn,'ndbc_air_t','ndbc_relhumid','ndbc_spechumid');





  if ( ischar(stnm_or_stn_or_loc) )
  elseif ( isfield(stnm_or_stn_or_loc,'lon') )
  else
  end;



    %%%% ??? HACK:
    %dts = datenum(2011,4,1):datenum(2010,1,0,23,59,59);
    dts = datenum(2011,4,1):now;




    %s = strrep(s,'"1,0','"10');




    %%%% ??? HACK:    dts = datenum(2012,11,1):datenum(2013,1,1)-(1/24/60/60);







      % Meshuggahs because Stevens interface doesn't seem to support specifying
      % query dates in UTC: We want an exact month worth in UTC, so for EST/EDT
      % transition months, we will get slightly more or less than 24 hours...???
      begdati = datetime(begdt,'TimeZone','America/New_York','ConvertFrom','datenum');
      begdt = begdt + ( hour(tzoffset(begdati)) / 24);
      enddati = datetime(enddt,'TimeZone','America/New_York','ConvertFrom','datenum');
disp(datestr(enddt));
disp(enddati);
      enddt = enddt + ( hour(tzoffset(enddati)) / 24);
disp(datestr(enddt));






  OneSecond = (1/24/60/60);



  if ( isempty(dts) )
    dts = fdep_stns{stix,3}:fdep_stns{stix,4};
    %%%% ??? HACK
    dts = datenum(2012,11,1):datenum(2013,1,1) - OneSecond;
  end;



  for begdt = yrmos(:)';
    yr = get_year(begdt);
    mo = get_month(begdt);
    enddt = datenum(yr,mo+1,1) - OneSecond;








      hdrs = textscan(lns{1},'%q','Delimiter',',','EndOfLine','\r\n');
      C = textscan(s,'%q','Delimiter',',','EndOfLine','\r\n','Headerlines',1);




  %{
  Select a Station:
  value="8720757" A - Bing's Landing
  value="8722213" B - Binney Dock
  value="8728744" C - Dry Bar
  value="8728603" D - East Bay
  value="8725081" E - Gordon River Inlet
  value="8728703" F - Little St. Marks River
  value="8721843" G - Melbourne
  value="8725114" H - Naples Bay
  value="8728732" I - Pilot's Cove
  value="8721147" J - Ponce de Leon South
  value="8722375" K - St. Lucie Inlet
  value="8720494" L - Tolomato River
  value="8722125" M - Vero Beach
  value="8720554" N - Vilano Beach

  http://www.fldep-stevens.com/export-8722375.php?t=c&txtFDate=02/01/2013&selFTime=00&txtTDate=02/28/2013&selTTime=23&selChannel=&rdoDateTime=1&rdoWTemp=1&rdoWSpeed=2&rdoATemp=1&rdoPressure=1&rdoRainfall=1&rdoWLevel=1
  %}






str(strfind(str,char([13 10]))) = '';
% make remaining \r into \n
str(str==char(13)) = char(10);








    begdati = datetime(begdt,'TimeZone','America/New_York','ConvertFrom','datenum');
    begdt = begdt - hour(tzoffset(begdati));






    enddy = get_dom(enddt);





    begdati = datetime(begdt,'TimeZone','America/New_York','ConvertFrom','datenum');
    begdt = begdt - hour(tzoffset(begdati));

    t = hour(tzoffset(t1));






for yr = 2012
  for mo=12
    begdy = 1;
    enddy = get_dom(datenum(yr,mo+1,0));
    url = sprintf('http://www.fldep-stevens.com/export-8722375.php?t=c&txtFDate=%02d/%02d/%04d&selFTime=00&txtTDate=%02d/%02d/%04d&selTTime=23&selChannel=&rdoDateTime=1&rdoWTemp=1&rdoWSpeed=2&rdoATemp=1&rdoPressure=1&rdoRainfall=1&rdoWLevel=1',mo,begdy,yr,mo,enddy,yr);
    fname = fullfile(datpath,sprintf('%s%s-%04d%02d.csv'));
    [s,status] = urlread(url);
    if ( ~status || numel(s) < 128 )
      warning('No data from %s',url);
      continue;
    end;
    lnends = strfind(s,sprintf('\n'));
keyboard;
  end;
end;






  if ( ~isfield(bat,'beta') )
    if ( N < 3 )
      % "Aspect Angle" is actually (v. GRADIENTM) the "direction of steepest
      % descent expressed as an azimuth measured clockwise from north".
      [aspect_deg,slope_deg,beta_y,beta_x] = gradientm(LAT,LON,bat.field,wgs84Ellipsoid);
      % We are interested in direction
      beta_x = -beta_x; beta_y = -beta_y;
      bat.beta = uv_to_spd(beta_x,beta_y);
      bat.beta_deg = uv_to_dir_curr(beta_x,beta_y);
      aspect_deg=[]; slope_deg=[]; beta_y=[]; beta_x=[]; clear aspect_deg slope_deg beta_y beta_x

    else
      hx = 1e3 .* distance_wgs84(bat.lat(1),bat.lon(1),bat.lat(1),bat.lon);
      hy = 1e3 .* distance_wgs84(bat.lat(1),bat.lon(1),bat.lat,bat.lon(1));
      [beta_x,beta_y] = gradientn(bat.field,N,hx,hy);
      hx=[]; hy=[]; clear hx hy

      bat.beta = uv_to_spd(beta_x,beta_y);
      bat.beta_deg = uv_to_dir_curr(beta_x,beta_y);
      beta_y=[]; beta_x=[]; clear beta_y beta_x
    end;
  end;









    if ( N < 3 )
      [aspect_deg,bat.beta_deg,beta_y,beta_x] = gradientm(LAT,LON,bat.field,wgs84Ellipsoid);
      bat.beta = uv_to_spd(beta_x,beta_y);
      aspect_deg=[]; beta_y=[]; beta_x=[]; clear aspect_deg beta_y beta_x

    else
      hx = 1e3 .* distance_wgs84(bat.lat(1),bat.lon(1),bat.lat(1),bat.lon);
      hy = 1e3 .* distance_wgs84(bat.lat(1),bat.lon(1),bat.lat,bat.lon(1));
      if ( N == 3 )
        [beta_x,beta_y] = gradient(bat.field,hx,hy);
      else
        [beta_x,beta_y] = gradientn(bat.field,N,hx,hy);
      end;
      hx=[]; hy=[]; clear hx hy

      bat.beta = uv_to_spd(beta_x,beta_y);
      bat.beta_deg = uv_to_dir_curr(beta_x,beta_y);
      beta_y=[]; beta_x=[]; clear beta_y beta_x
    end;








          ... % NCORE experiment 2007-2008, courtesy TN Lee and S Sponaugle, by K Shulzitski in 2011.
          , ... % klgf1,25.0313833,-80.3480167,23 # Key Largo Current/Temperature 7/23 m (offshore French Reef)
          , ... % mrtf1,24.7403000,-80.7763667,23 # Marathon Current/Temperature 8/23 m (offshore Tennessee Reef)
          ... % NCORE studies 2000-2002, provided courtesy of TN Lee, S Sponaugle, and E Williams
          , ... % ncora,25.1090,-80.3803,4.1
          , ... % ncorb,25.0932,-80.3550,7.0
          , ... % ncorc,25.0673,-80.3183,21.8
          , ... % ncor1,25.0740,-80.3178,26.0
          , ... % ncor2,25.0733,-80.3183,21.8
          , ... % ncot1,25.06900,-80.31950,9.75
          , ... % ncot2,25.07383,-80.32450,7.0
          , ... % ncot3,25.07950,-80.33433,4.0
          , ... % ncot4,25.08783,-80.34717,5.6
          ...





      bat.beta = sqrt((beta_y.^2)+(beta_x.^2));




1,450,250,092


function [bet,bat] = find_ngdc_slope(bat,lon,lat,N)
%function [bet,bat] = find_ngdc_slope(bat,lon,lat,N)
%
% Find seafloor slope BET from bathymetry struct BAT at location(s) LON and
% LAT.  If necessary, calls GRADIENTM to calculate 2-point finite difference
% and add fields .aspect_deg,.beta_deg,.beta_y,.beta_x to BAT. Calls INTERP2
% to linearly interpolate slope to coordinates LON,LAT; if N>2, then finite
% difference with an N+1 point template instead of a 2-point mean.
%
% Last Saved Time-stamp: <Wed 2016-07-13 16:33:25 Eastern Daylight Time gramer>

  if ( ~exist('N','var') || isempty(N) )
    N = 0;
  end;
  if ( ~isscalar(N) || ~isnumeric(N) || (floor(N) ~= N) )
    error('If specified, N must be an integer scalar!');
  end;
  if ( ~isfield(bat,'lon') || ~isfield(bat,'lat') || ~isfield(bat,'field') )
    error('BAT must be a bathymetry struct with fields .lon,.lat,.field');
  end;

  [LON,LAT] = meshgrid(bat.lon,bat.lat);

  if ( ~isfield(bat,'beta') )
    % [bat.aspect_deg,bat.beta_deg,bat.beta_y,bat.beta_x] = gradientm(LAT,LON,bat.field,wgs84Ellipsoid);
    % bat.beta = sqrt((bat.beta_y.^2)+(bat.beta_x.^2));
    hx = distance_wgs84(bat.lat(1),bat.lon(1),bat.lat(1),bat.lon);
    hy = distance_wgs84(bat.lat(1),bat.lon(1),bat.lat,bat.lon(1));
    [bat.beta_x,bat.beta_y] = gradientn(bat.field,hx,hy);
    bat.aspect_deg = 
    bat.beta_deg = 
    bat.beta = sqrt((bat.beta_y.^2)+(bat.beta_x.^2));
  end;

  if ( N < 1 )
    bet = interp2(LAT,LON,bat.beta,lat,lon,'linear',NaN);
    if ( any(isnan(bet)) )
      error('One or more points fall outside bathymetry region!');
    end;

  else
keyboard;
    lonix = interp2(LAT,LON,repmat(1:numel(bat.lon),[numel(bat.lat),1]),lat,lon,'nearest',NaN);
    latix = interp2(LAT,LON,repmat(1:numel(bat.lat),[1,numel(bat.lon)]),lat,lon,'nearest',NaN);
    if ( any(isnan(lonix)) || any(isnan(latix)) )
      error('One or more points fall outside bathymetry region!');
    end;

    x.bath = bat;
    for ix=1:numel(lonix)
      d = distance_wgs84(bat.lat(latix(ix)),bat.lon(lonix(1)),bat.lat(latix(1)),bat.lon(lonix(1)+N+1));
      ang = bat.beta_deg(latix(ix),lonix(ix));
      [lons,lats] = transect_wgs84(x,[-d,d],ang,'bath')
      switch (N),
       case 3,
        bet(ix,1) = 
       case 5,
       case 7,
       case 9,
       otherwise,
      end; %switch (N),
      clear lons lats
    end; %for ix=1:numel(lonix)
    x=[]; clear x
  end;

  LON=[]; LAT=[]; clear LON LAT

return;









  if ( useHighestRes )
    fld = 'ngdc_hires_bathy';
  else
    fld = 'ngdc_lores_bathy';
  end;








    %disp({min(rixen),max(rixen),min(cixen),max(cixen)});
    %d = griddata(LON(rixen,cixen),LAT(rixen,cixen),DAT(rixen,cixen),MLON(rixen,cixen),MLAT(rixen,cixen));
    %d = griddata(LON(rixen,cixen),LAT(rixen,cixen),DAT(rixen,cixen),MLON(rixen,cixen),MLAT(rixen,cixen),'natural');
    %d = griddata(LON(rixen,cixen),LAT(rixen,cixen),DAT(rixen,cixen),MLON(rixen,cixen),MLAT(rixen,cixen),'cubic');






% Linearly interpolate back onto a latitude-longitude grid
lon = linspace(min(LON(:)),max(LON(:)),size(LON,2));
lat = linspace(min(LAT(:)),max(LAT(:)),size(LAT,1));

[MLON,MLAT] = meshgrid(lon,lat);

dat=[]; clear dat
%dat = repmat(nan,size(DAT));
dat = repmat(nan,[numel(lat),numel(lon)]);









[y1,x1] = ind2sub(size(LON),ix(1));
[ye,xe] = ind2sub(size(LON),ix(end));







%ix = find(-81.2<=LON&LON<=-80.8 & 24.5<=LAT&LAT<=24.8); clear ix

ul = find(LON<=-81.2&LAT<=24.5,1,'last'), [ulix,uljx] = ind2sub(size(LON),ul),
lr = find(LON>=-81.2&LAT>=24.8,1,'first'), [lrix,lrjx] = ind2sub(size(LON),lr),

LON = LON(ulix:lrix,uljx:lrjx); LAT = LAT(ulix:lrix,uljx:lrjx); nansummary(LON), nansummary(LAT),

clear ans ix lr lrix lrjx ul ulix uljx









    disp({min(rixen),max(rixen),min(cixen),max(cixen)});
    d = griddata(LON(rixen(rix),cixen(cix)),LAT(rixen(rix),cixen(cix)),DAT(rixen(rix),cixen(cix)),...
                 MLON(rixen(rix),cixen(cix)),MLAT(rixen(rix),cixen(cix)));
    dat(allrixen(roverlap+1:end-roverlap-1),allcixen(coverlap+1:end-coverlap-1)) = ...

    dat(nestrix,nestcix) = d(gridrix,gridcix)








    fullcix = begcix-coverlap:begcix+cstrd-1+coverlap;
    cix = find(1 <= fullcix & fullcix <= size(LAT,2));
    nestcix = begcix:begcix+cstrd-1;
    gridcix = begcix-coverlap:begcix+cstrd-1+coverlap;
    cix = find(1 <= cixen & cixen <= size(LAT,2));








for begrixix=1:numel(begrixen);
  fprintf(1,'%g of %g rows\n',begrixix,numel(begrixen)); 
  tic,

  begrix = begrixen(begrixix);
  rixen = begrix-roverlap:begrix+rstrd-1+roverlap;
  rix = find(1 <= rixen & rixen <= size(LAT,1));

  for begcixix=1:numel(begcixen); 
    begcix=begcixen(begcixix);
    fullcix = begcix-coverlap:begcix+cstrd-1+coverlap;
    cix = find(1 <= fullcix & fullcix <= size(LAT,2));
    nestcix = begcix:begcix+cstrd-1;
    gridcix = begcix-coverlap:begcix+cstrd-1+coverlap;
    cix = find(1 <= cixen & cixen <= size(LAT,2));

    disp({min(rixen),max(rixen),min(cixen),max(cixen)});
    d = griddata(LON(rixen(rix),cixen(cix)),LAT(rixen(rix),cixen(cix)),DAT(rixen(rix),cixen(cix)),...
                 MLON(rixen(rix),cixen(cix)),MLAT(rixen(rix),cixen(cix)));
    dat(allrixen(roverlap+1:end-roverlap-1),allcixen(coverlap+1:end-coverlap-1)) = ...

    dat(nestrix,nestcix) = d(gridrix,gridcix)
  end; 
  clear begcixix begcix cixen

  toc,
end;
clear begrixix begrix rixen









for begrixix=1:numel(begrixen);
  fprintf(1,'%g of %g rows\n',begrixix,numel(begrixen)); 
  tic,

  begrix = begrixen(begrixix);
  allrixen = begrix-roverlap:begrix+rstrd-1+roverlap;
  rixen = allrixen(1 <= allrixen & allrixen <= size(LAT,1));
  for begcixix=1:numel(begcixen); 
    begcix=begcixen(begcixix); 
    allcixen = begcix-coverlap:begcix+cstrd-1+coverlap;
    cixen = allcixen(1 <= allcixen & allcixen <= size(LAT,2));
    disp({min(rixen),max(rixen),min(cixen),max(cixen)});
    d = griddata(LON(rixen,cixen),LAT(rixen,cixen),DAT(rixen,cixen),MLON(rixen,cixen),MLAT(rixen,cixen));
    dat(allrixen(roverlap+1:end-roverlap-1),allcixen(coverlap+1:end-coverlap-1)) = ...
  end; 
  clear begcixix begcix cixen

  toc,
end;
clear begrixix begrix rixen







for begrixix=1:numel(begrixen);
  fprintf(1,'%g of %g rows\n',begrixix,numel(begrixen)); 
  tic,

  begrix = begrixen(begrixix);
  if ( begrixix == 1 )
    rixen = begrix:begrix+rstrd-1+roverlap;
  elseif ( begrixix == numel(begrixen) )
    rixen = begrix-roverlap:begrix+rstrd-1;
  else
    rixen = begrix-roverlap:begrix+rstrd-1+roverlap;
  end;
  for begcixix=1:numel(begcixen); 
    begcix=begcixen(begcixix); 
    if ( begcixix == 1 )
      cixen = begcix:begcix+cstrd-1+coverlap;
    elseif ( begcixix == numel(begcixen) )
      cixen = begcix-coverlap:begcix+cstrd-1;
    else
      cixen = begcix-coverlap:begcix+cstrd-1+coverlap;
    end;
    disp({min(rixen),max(rixen),min(cixen),max(cixen)});
    dat(rixen(roverlap+1:end-roverlap-1),allcixen(coverlap+1:end-coverlap-1)) = ...
        griddata(LON(rixen,cixen),LAT(rixen,cixen),DAT(rixen,cixen),MLON(rixen,cixen),MLAT(rixen,cixen));
  end; 
  clear begcixix begcix cixen

  toc,
end;
clear begrixix begrix rixen







for begrixix=numel(begrixen):-1:1;
  for begcixix=numel(begcixen):-1:1;






%for begrixix=1:numel(begrixen);
for begrixix=numel(begrixen):-1:1;
  fprintf(1,'%g of %g rows\n',begrixix,numel(begrixen)); 
  tic,

  begrix = begrixen(begrixix);
  allrixen = begrix-roverlap:begrix+rstrd-1+roverlap;
  rixen = allrixen(1 <= allrixen & allrixen <= size(LAT,1));
  %for begcixix=1:numel(begcixen); 
  for begcixix=numel(begcixen):-1:1;
    begcix=begcixen(begcixix); 
    allcixen = begcix-coverlap:begcix+cstrd-1+coverlap;
    cixen = allcixen(1 <= allcixen & allcixen <= size(LAT,2));
    disp({min(rixen),max(rixen),min(cixen),max(cixen)});
    dat(allrixen(roverlap+1:end-roverlap-1),allcixen(coverlap+1:end-coverlap-1)) = ...
        griddata(LON(rixen,cixen),LAT(rixen,cixen),DAT(rixen,cixen),MLON(rixen,cixen),MLAT(rixen,cixen));
  end; 
  clear begcixix begcix cixen

  toc,
end;
clear begrixix begrix rixen









for begrixix=1:numel(begrixen);
  fprintf(1,'%g of %g rows\n',begrixix,numel(begrixen)); 
  tic,

  begrix = begrixen(begrixix);
  allrixen = begrix-roverlap:begrix+rstrd-1+roverlap;
  rixen = allrixen(1 <= allrixen & allrixen <= size(LAT,1));
  for begcixix=1:numel(begcixen); 
    begcix=begcixen(begcixix); 
    allcixen = begcix-coverlap:begcix+cstrd-1+coverlap;
    cixen = allcixen(1 <= allcixen & allcixen <= size(LAT,2));
    disp({min(rixen),max(rixen),min(cixen),max(cixen)});
    d = griddata(LON(rixen,cixen),LAT(rixen,cixen),DAT(rixen,cixen),MLON(rixen,cixen),MLAT(rixen,cixen));
keyboard;
    %dat(rixen(roverlap:end-roverlap),cixen(coverlap:end-coverlap)) = d(roverlap:end-roverlap,coverlap:end-coverlap);
    dat(allrixen(roverlap:end-roverlap),allcixen(coverlap:end-coverlap)) = ...
        d(allrixen(roverlap:end-roverlap),allcixen(coverlap:end-coverlap));
    d=[]; clear d
  end; 
  clear begcixix begcix cixen

  toc,
end;
clear begrixix begrix rixen








for begrixix=1:numel(begrixen);
  begrix = begrixen(begrixix);
  rixen = begrix-roverlap:begrix+rstrd-1+roverlap;
  rixen(1 > rixen | rixen > size(LAT,1)) = [];
  fprintf(1,'%g of %g rows\n',begrixix,numel(begrixen)); 
  tic,
  for begcixix=1:numel(begcixen); 
    begcix=begcixen(begcixix); 
    cixen = begcix-coverlap:begcix+cstrd-1+coverlap;
    cixen(1 > cixen | cixen > size(LAT,2)) = [];
    disp({min(rixen),max(rixen),min(cixen),max(cixen)});
    d = griddata(LON(rixen,cixen),LAT(rixen,cixen),DAT(rixen,cixen),MLON(rixen,cixen),MLAT(rixen,cixen));
keyboard;
    dat(rixen(roverlap:end-roverlap),cixen(coverlap:end-coverlap)) = d(roverlap:end-roverlap,coverlap:end-coverlap);
    d=[]; clear d
  end; 
  clear begcixix begcix cixen
  toc,
end;
clear begrixix begrix rixen






% % rstrd = 1371;
% % cstrd =  551;
% rstrd = ceil(size(LAT,1)/40);
% cstrd = ceil(size(LAT,2)/40);





  %rixen = begrix:begrix+rstrd-1;
    %cixen = begcix:begcix+cstrd-1;






1;
%% SCRIPT process_F010.m
% Process "F010" USGS DEM high-res. (20 m, UTM grid) coastal bathymetry file
% for south  Florida: loads the ASCII file using the confusiingly-named
% USGS24KDEM (MATLAB Map Toolbox), interpolates onto a latitude-longitude
% grid, and saves the resulting variables (lon,lat,dat) in a .MAT file.
%
% Data is documented in file "F010_25081C1_BIG.dem.IS.FLORIDA.30m.bathy.txt"
%
% USGS DEM file format reference URL:
%  https://en.wikipedia.org/wiki/USGS_DEM
% Data source URL:
%  http://estuarinebathymetry.noaa.gov/bathy_htmls/F010.html
%
% Last Saved Time-stamp: <Thu 2016-07-07 17:44:54 Eastern Daylight Time gramer>

set_more off

coastpath = get_ecoforecasts_path('coast');
fname = 'F010_25081C1_BIG.dem';

if ( ~exist('LAT','var') ||~exist('LON','var') ||~exist('DAT','var') )
  DAT=[]; LAT=[]; LON=[]; clear DAT LAT LON
  intermed_matfname = fullfile(coastpath,[fname,'-intermed.mat']);
  if ( exist(intermed_matfname,'file') )
    disp(['Loading ',intermed_matfname]);
    load(intermed_matfname);
  else
    fpath = fullfile(coastpath,fname);
    disp(['USGS24KDEM(',fpath,',1)']);
    [LAT,LON,DAT] = usgs24kdem(fpath,1);
    disp(['Saving ',intermed_matfname]);
    save(intermed_matfname,'LAT','LON','DAT','fpath','-v7.3');
  end;
end;

LAT = LAT(

% Linearly interpolate back onto a latitude-longitude grid
lon = linspace(min(LON(:)),max(LON(:)),size(LON,2));
lat = linspace(min(LAT(:)),max(LAT(:)),size(LAT,1));

% Mirror the way LON and LAT are arranged by USGS24KDEM
[MLON,MLAT] = meshgrid(lon,lat(end:-1:1));


overlap = floor(size(LAT,1)/40);

% % rstrd = 1371;
% % cstrd =  551;
rstrd = ceil(size(LAT,1)/4);
cstrd = ceil(size(LAT,2)/4);
% rstrd = ceil(size(LAT,1)/40);
% cstrd = ceil(size(LAT,2)/40);
begrixen = 1:rstrd:size(LAT,1);
begcixen = 1:cstrd:size(LAT,2);

dat=[]; clear dat
dat = repmat(nan,size(DAT));

for begrixix=1:numel(begrixen);
  begrix = begrixen(begrixix);
  %rixen = begrix:begrix+rstrd-1;
  rixen = begrix-overlap:begrix+rstrd+overlap;
  rixen(1 > rixen | rixen > size(LAT,1)) = [];
  fprintf(1,'%g of %g rows\n',begrixix,numel(begrixen)); 
  tic,
  for begcixix=1:numel(begcixen); 
    begcix=begcixen(begcixix); 
    %cixen = begcix:begcix+cstrd-1;
    cixen = begcix-overlap:begcix+cstrd+overlap;
    cixen(1 > cixen | cixen > size(LAT,2)) = [];
    disp({min(rixen),max(rixen),min(cixen),max(cixen)});
    dat(rixen,cixen) = griddata(LON(rixen,cixen),LAT(rixen,cixen),DAT(rixen,cixen),MLON(rixen,cixen),MLAT(rixen,cixen));
  end; 
  clear begcixix begcix cixen
  toc,
end;
clear begrixix begrix rixen

MLAT=[]; MLON=[]; clear MLAT MLON

set_more

if 1;
  matfname = fullfile(coastpath,[fname,'.mat']);
  disp(['Saving ',matfname]);
  save(matfname,'lat','lon','dat','fname','-v7.3');
end;










for begrixix=numel(begrixen):-1:1;





    dy = distance_wgs84(bat.lat(latix(1)),bat.lon(lonix(1)),bat.lat(lonix(1)+N+1),bat.lon(lonix(1)));







function [bet,bat] = find_ngdc_slope(bat,lon,lat,N)
%function [bet,bat] = find_ngdc_slope(bat,lon,lat,N)
%
% Find seafloor slope BET from bathymetry struct BAT at location(s) LON and
% LAT.  If necessary, calls GRADIENTM to calculate 2-point finite difference
% and add fields .aspect_deg,.beta_deg,.beta_y,.beta_x to BAT. Calls INTERP2
% to linearly interpolate slope to coordinates LON,LAT; if N>0, then finite
% difference with an (N*2)+1 point template instead of a 2-point mean.
%

  if ( ~exist('N','var') || isempty(N) )
    N = 0;
  end;
  if ( ~isscalar(N) || ~isnumeric(N) || (floor(N) ~= N) )
    error('If specified, N must be an integer scalar!');
  end;
  if ( ~isfield(bat,'lon') || ~isfield(bat,'lat') || ~isfield(bat,'field') )
    error('BAT must be a bathymetry struct with fields .lon,.lat,.field');
  end;

  [LON,LAT] = meshgrid(bat.lon,bat.lat);

  if ( ~isfield(bat,'beta') )
    [bat.aspect_deg,bat.beta_deg,bat.beta_y,bat.beta_x] = gradientm(LAT,LON,bat.field,wgs84Ellipsoid);
    bat.beta = sqrt((bat.beta_y.^2)+(bat.beta_x.^2));
  end;

  if ( N < 1 )
    bet = interp2(LAT,LON,bat.beta,lat,lon,'linear',NaN);
    if ( any(isnan(bet)) )
      error('One or more points fall outside bathymetry region!');
    end;

  else
keyboard;
    lonix = interp2(LAT,LON,repmat(1:numel(bat.lon),[numel(bat.lat),1]),lat,lon,'nearest',NaN);
    latix = interp2(LAT,LON,repmat(1:numel(bat.lat),[1,numel(bat.lon)]),lat,lon,'nearest',NaN);
    if ( any(isnan(lonix)) || any(isnan(latix)) )
      error('One or more points fall outside bathymetry region!');
    end;

    dx = distance_wgs84(bat.lat(latix(1)),bat.lon(lonix(1)),bat.lat(latix(1)),bat.lon(lonix(1)+N+1));
    dy = distance_wgs84(bat.lat(latix(1)),bat.lon(lonix(1)),bat.lat(lonix(1)+N+1),bat.lon(lonix(1)));

    x.bath = bat;
    for ix=1:numel(lonix)
      ang = bat.beta_deg(latix(ix),lonix(ix));
      [lons,lats] = transect_wgs84(x,[dx,dy],ang,'bath')
      
    end;
    clear lons lats
    x=[]; clear x
  end;

  LON=[]; LAT=[]; clear LON LAT

return;












1;


if ( ~exist('doPrint','var') )
  doPrint = false;
end;
if ( doPrint )
  figpath = get_ecoforecasts_path('figs');
end;

if ( ~exist('doFKEYS','var') )
  doFKEYS = false;
end;
if ( ~exist('doGOM2','var') )
  doGOM2 = false;
end;

if ( doGOM2 )
  if ( ~exist('mdl','var') || ~isfield(mdl,'model') || ~strcmpi(mdl.model,'GOM2') )
    mdl=[]; clear mdl;
    angom2;
  end;
  dts = datenum(2014,1,[1:31],12,0,0);
  mdl.model = 'GOM2';
elseif ( doFKEYS )
  if ( ~exist('mdl','var') || ~isfield(mdl,'model') || ~strcmpi(mdl.model,'FKEYS') )
    mdl=[]; clear mdl;
    anfkeys;
  end;
  dts = datenum(2008,1,1,6,0,0):(6/24):datenum(2008,1,7,18,0,0);
  mdl.model = 'FKEYS';
else
  if ( ~exist('mdl','var') || ~isfield(mdl,'model') || ~strcmpi(mdl.model,'eFKEYS') )
    mdl=[]; clear mdl;
    anefkeys;
  end;
  dts = datenum(2012,1,1,6,0,0):(6/24):datenum(2012,1,7,18,0,0);
  mdl.model = 'eFKEYS';
end;

% Perform model-in situ comparisons

cstnms = {'fwyf1','mlrf1','lonf1','smkf1','looe1','sanf1'};
ncstnms = numel(cstnms);

for stix=1:ncstnms
  stnm = cstnms{stix};
  %mdl.(stnm)=[];
  %mdl = rmfield(mdl,stnm);
  if ( ~isfield(mdl,stnm) )

    mdl.(stnm) = get_station_from_station_name(stnm);
    mdl.(stnm) = station_optimal_isobath_orientation(mdl.(stnm));

    [mdl.(stnm).lonerr,mdl.(stnm).lonix] = min(abs(mdl.lon-mdl.(stnm).lon));
    [mdl.(stnm).laterr,mdl.(stnm).latix] = min(abs(mdl.lat-mdl.(stnm).lat));

    mdl.(stnm).t.z = mdl.z;
    mdl.(stnm).t.date = mdl.d;
    mdl.(stnm).t.prof = mdl.t(:,:,mdl.(stnm).latix,mdl.(stnm).lonix);
    mdl.(stnm).t.data = nanmean(mdl.(stnm).t.prof,2);
    for zix = 1:numel(mdl.z)
      fld = ['t',num2str(mdl.z(zix))];
      mdl.(stnm).(fld).z = mdl.z(1);
      mdl.(stnm).(fld).date = mdl.d;
      mdl.(stnm).(fld).data = mdl.t(:,zix,mdl.(stnm).latix,mdl.(stnm).lonix);
    end;
    clear fld zix

    if ( isfield(mdl,'s') )
      mdl.(stnm).s.z = mdl.z;
      mdl.(stnm).s.date = mdl.d;
      mdl.(stnm).s.prof = mdl.s(:,:,mdl.(stnm).latix,mdl.(stnm).lonix);
      mdl.(stnm).s.data = nanmean(mdl.(stnm).s.prof,2);
      for zix = 1:numel(mdl.z)
        fld = ['s',num2str(mdl.z(zix))];
        mdl.(stnm).(fld).z = mdl.z(1);
        mdl.(stnm).(fld).date = mdl.d;
        mdl.(stnm).(fld).data = mdl.s(:,zix,mdl.(stnm).latix,mdl.(stnm).lonix);
      end;
      clear fld zix
    end; %if ( isfield(mdl,'s') )

    mdl.(stnm).u.z = mdl.z;
    mdl.(stnm).u.date = mdl.d;
    mdl.(stnm).u.prof = mdl.u(:,:,mdl.(stnm).latix,mdl.(stnm).lonix);
    mdl.(stnm).u.data = nanmean(mdl.(stnm).u.prof,2);

    mdl.(stnm).v.z = mdl.z;
    mdl.(stnm).v.date = mdl.d;
    mdl.(stnm).v.prof = mdl.v(:,:,mdl.(stnm).latix,mdl.(stnm).lonix);
    mdl.(stnm).v.data = nanmean(mdl.(stnm).v.prof,2);

    % Reorient vector currents to local isobath: 'x'=cross-shore, 'l'=longshore
    mdl.(stnm).x.z = mdl.z;
    mdl.(stnm).x.date = mdl.d;
    mdl.(stnm).l.z = mdl.z;
    mdl.(stnm).l.date = mdl.d;
    [mdl.(stnm).x.prof,mdl.(stnm).l.prof] = ...
        reorient_vectors(mdl.(stnm).isobath_orientation,mdl.(stnm).u.prof,mdl.(stnm).v.prof);

    mdl.(stnm).x.data = nanmean(mdl.(stnm).x.prof,2);
    mdl.(stnm).l.data = nanmean(mdl.(stnm).l.prof,2);
  end; %if ( ~isfield(mdl,stnm) )

end; %for stix=1:ncstnms


min_t = +inf;
max_t = -inf;
min_s = +inf;
max_s = -inf;
min_curr = +inf;
max_curr = -inf;

dix = find(ismember(mdl.d,dts));

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


%stnixen = [1,2,4,5,6];
stnixen = [1,2,4,5,6];

if 1;
  fmg;
  ncols = numel(stnixen);
  colix = 1;
  %min_hist = floor(-max_spd); max_hist = ceil(+max_spd); d_hist = 0.025;
  min_hist = -0.2; max_hist = +0.2; d_hist = 0.01;
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
  if ( doPrint )
    print('-dpng',fullfile(figpath,[mdl.model,'.uv_hist.png']));
  end;
end;

%stnm = 'fwyf1';
%stnm = 'mlrf1';
%stnm = 'sanf1';
% if ( stnm )
%%for stix=1:ncstnms
%for stix=stnixen(:)'
for stix=stnixen([2,4])
  stnm = cstnms{stix};
  disp(stnm);

  % %fmg;
  % stn = []; clear stn
  % stn = get_station_from_station_name(stnm);
  % stn = station_optimal_isobath_orientation(stn);
  % stn = read_hires_bathymetry(stn,[12e3,12e3],[],false);

  % Get high-resolution bathymetry
  mdl.(stnm) = read_hires_bathymetry(mdl.(stnm),[12e3,12e3],[],false);

  % Draw bathymetry map
  mdl.(stnm) = plot_hires_bathymetry(mdl.(stnm),-[0:5:80],[12e3,12e3],true,[],false,[],false);
  fh = gcf;
  axis(axis);
  [c,h] = contour(mdl.lon,mdl.lat,squeeze(nanmean(mdl.t(:,1,:,:))),[floor(min_t):0.25:ceil(max_t)]);
  clabel(c,h);
  c=[]; clear c h


  %fmg;
  [res,mdl.(stnm)] = plot_bathy_transect(mdl.(stnm),10,mdl.(stnm).isobath_orientation+90,'ngdc_hires_bathy');
  [dx,az] = distance_wgs84(mdl.(stnm).lat,mdl.(stnm).lon,res.lat,res.lon);
  dx(round(az) ~= mdl.(stnm).isobath_orientation+90) = -dx(round(az) ~= mdl.(stnm).isobath_orientation+90);
  res.dx = dx;

  lonix = interp1(mdl.lon,1:mdl.nlon,res.lon,'nearest');
  latix = interp1(mdl.lat,1:mdl.nlat,res.lat,'nearest');
  ix = sub2ind([size(mdl.t,3),size(mdl.t,4)],latix,lonix);
  %t = squeeze(mdl.t(1,:,ix));
  t = squeeze(nanmean(mdl.t(:,:,ix)));
  disp([nanmin(t(:)),nanmax(t(:))]);
  u = squeeze(nanmean(mdl.u(:,:,ix)));
  v = squeeze(nanmean(mdl.v(:,:,ix)));
  % disp([nanmin(u(:)),nanmax(u(:))]);
  % disp([nanmin(v(:)),nanmax(v(:))]);
  [x,l] = reorient_vectors(mdl.(stnm).isobath_orientation,u,v);
  disp([nanmin(l(:)),nanmax(l(:))]);
  disp([nanmin(x(:)),nanmax(x(:))]);

  [c,h] = contourf(res.dx,-mdl.z,t,[11.0:0.5:24.5]);
  set(gca,'CLim',[11.0,24.5]);
  clabel(c,h);
  c=[]; clear c h
  plot(res.dx,res.field,'k-','LineWidth',2);
  xlim([-10,+10]); ylim([-140,0]);
  titlename([mdl.model,': ',upper(stnm),' bathymetry vs. T cross-shore profile (1 week mean)']);
  if ( doPrint )
    print('-dpng',fullfile(figpath,[mdl.model,'.',stnm,'.section_t_mean.png']));
  end;

  %fmg;
  plot_bathy_transect(mdl.(stnm),10,mdl.(stnm).isobath_orientation+90,'ngdc_hires_bathy');
  [c,h] = contourf(res.dx,-mdl.z,l,[-0.300:0.05:+1.200]);
  set(gca,'CLim',[-0.300,+1.200]);
  clabel(c,h);
  c=[]; clear c h
  plot(res.dx,res.field,'k-','LineWidth',2);
  xlim([-10,+10]); ylim([-140,0]);
  titlename([mdl.model,': ',upper(stnm),' bathymetry vs. alongshore V (1 week mean)']);
  if ( doPrint )
    print('-dpng',fullfile(figpath,[mdl.model,'.',stnm,'.section_ls_mean.png']));
  end;

  %fmg;
  plot_bathy_transect(mdl.(stnm),10,mdl.(stnm).isobath_orientation+90,'ngdc_hires_bathy');
  [c,h] = contourf(res.dx,-mdl.z,x,[-0.150:0.010:+0.150]);
  set(gca,'CLim',[-0.150,+0.150]);
  clabel(c,h);
  c=[]; clear c h
  plot(res.dx,res.field,'k-','LineWidth',2);
  xlim([-10,+10]); ylim([-140,0]);
  titlename([mdl.model,': ',upper(stnm),' bathymetry vs. cross-shore U (1 week mean)']);
  if ( doPrint )
    print('-dpng',fullfile(figpath,[mdl.model,'.',stnm,'.section_xs_mean.png']));
  end;

  % Draw cross-shore section on bathymetry map
  figure(fh);
  plot(res.lon,res.lat,'k-','LineWidth',2);
  clear fh
  if ( doPrint )
    print('-dpng',fullfile(figpath,[mdl.model,'.',stnm,'.surface_t_mean.png']));
  end;
end; %for stix=stnixen([2,4])'


clear ans cstnms stix stnm
clear az dx ig ix latix lonix s t
%clear datpath hycpath matfname ncols ncstnms ans
%stn=[]; clear res stn min_t max_t min_s max_s min_curr max_curr max_spd












  if ( ~exist('mdl','var') )






plot(snaf1.lon,snaf1.lat,'kp','MarkerFaceColor','r');





  relatedflds = {};


      elseif ( strncmpi(arg{1},'related',length('related')) )
        relatedflds(end+1) = { arg(2:end) };
        disp(['RELATED FIELDS: ']);
        disp(relatedflds{end});
        flds(ismember(flds,relatedflds{end}(2:end))) = [];







  if ( isempty(flds) )
  end;





function stn = qcstation(stn,varargin)
%function stn = qcstation(stn,[varname|{varnames...}|jump_initial_frame]...)
%
% Use BRUSH to allow caller to INTERACTIVELY  perform QC on all time-series
% (STRUCT fields with .date and .data subfields) in station STN. To limit
% work, User may specify a field name CHAR, or a CELLSTR of field names. Uses
% JUMPAX (v.) to jump through data graph of each variable; if caller passes
% in a numeric 2-vector JUMP_INITIAL_FRAME, initial XLIM of each data graph
% is set to that; if it is a numeric scalar < 1, initial XLIM is that pct. of
% the total record size (DEFAULT: 0.10); if >= 1, frame size is that many
% days; if +Inf, then JUMPAX is not called. Calls INPUT between fields.
%
% Last Saved Time-stamp: <Wed 2016-06-15 16:51:01 Eastern Daylight Time gramer>

  flds = {};

  % Default to one-tenth of total record size for JUMPAX (below)
  frame_size = 0.10;
  jump_initial_frame = [];

  args = varargin(:);
  nargs = numel(args);
  while ( nargs > 0 )
    if ( ischar(args{1}) )
      flds = cellstr(args{1});
    elseif ( iscellstr(args{1}) )
      flds = args{1};
    elseif ( isvector(args{1}) && isnumeric(args{1}) )
      if ( numel(args{1}) == 1 )
        frame_size = args{1};
      elseif ( numel(args{1}) == 2 )
        jump_initial_frame = args{1};
      end;
    else
      warning('ecoforecasts:qcstation:UnknownArg',...
              'Optional arg must be a CHAR, CELLSTR, or numeric vector');
    end;
    args(1) = [];
    nargs = nargs - 1;
  end;

  if ( isempty(flds) )
    flds = fieldnames(stn);
  end;

  for fldix = 1:numel(flds)
    fld = flds{fldix};
    if ( ~is_valid_ts(stn.(fld)) )
      disp(['Skipping STN.',fld]);

    else

      fldnm = strrep(fld,'_','\_');

      dts = stn.(fld).date;
      dat = stn.(fld).data;
      fh = fmg;
      %lh = plot_ts(stn.(fld));
      lh = plot(dts,dat); datetick3;
      titlename(['BRUSH points to remove from STN.',fldnm]);

      set(lh,'XDataSource','dts');
      set(lh,'YDataSource','dat');

      if ( isempty(jump_initial_frame) )
        if ( ~isinf(frame_size) )
          if ( frame_size >= 1 )
            jump_initial_frame = [dts(1),dts(1)+frame_size];
          elseif ( frame_size > 0 )
            jump_initial_frame = [dts(1),dts(floor(end*frame_size))];
          else
            error('Invalid FRAME SIZE %g',frame_size);
          end;
        end;
      end;
      if ( ~isempty(jump_initial_frame) )
        xlim(jump_initial_frame); 
        qcstation_jumpaction;
      end;

      brush(fh,'on');
      h = brush(fh);

      if ( ~isempty(jump_initial_frame) )
        [ax,fh] = jumpax([],[],@qcstation_jumpaction);
      end; %if ( ~isempty(jump_initial_frame) )
      if ( ishandle(fh) )
        close(fh);
      end;

      skip = false;
      inp = input('ENTER to save and move to next field, or "q"uit or "s"kip: ','s');
      if ( strncmpi(inp,'q',1) )
        break;
      elseif ( strncmpi(inp,'s',1) )
        skip = true;
      end;

keyboard;
      if ( ~skip )
        if ( numel(dts) ~= numel(stn.(fld).date) || any(dat ~= stn.(fld).data) )
          disp(['Modifying STN.',fld]);
          stn.(fld).date = dts;
          stn.(fld).data = dat;
        end;
      end;

    end;

  end;    

return;


function qcstation_jumpaction(ax)
  if ( exist('ax','var') && ishandle(ax) )
    ylim(ax,'default');
    datetick3(ax);
  else
    ylim('default');
    datetick3;
  end;
return;














function stn = qcstation(stn,varargin)
%function stn = qcstation(stn,[varname|{varnames...}|jump_initial_frame]...)
%
% Use BRUSH to allow caller to INTERACTIVELY  perform QC on all time-series
% (STRUCT fields with .date and .data subfields) in station STN. To limit
% work, User may specify a field name CHAR, or a CELLSTR of field names. Uses
% JUMPAX (v.) to jump through data graph of each variable; if caller passes
% in a numeric 2-vector JUMP_INITIAL_FRAME, initial XLIM of each data graph
% is set to that; if caller passes a numeric scalar < 1, initial XLIM is that
% percentage of the total record size (DEFAULT: 0.10); if >= 1, frame size is
% that many days; if +Inf, then JUMPAX is not called - calls PAUSE instead.
%
% Last Saved Time-stamp: <Wed 2016-06-08 15:27:02 Eastern Daylight Time gramer>

  flds = {};

  % Default to one-tenth of total record size for JUMPAX (below)
  frame_size = 0.10;
  jump_initial_frame = [];

  args = varargin(:);
  nargs = numel(args);
  while ( nargs > 0 )
    if ( ischar(args{1}) )
      flds = cellstr(args{1});
    elseif ( iscellstr(args{1}) )
      flds = args{1};
    elseif ( isvector(args{1}) && isnumeric(args{1}) )
      if ( numel(args{1}) == 1 )
        frame_size = args{1};
      elseif ( numel(args{1}) == 2 )
        jump_initial_frame = args{1};
      end;
    else
      warning('ecoforecasts:qcstation:UnknownArg',...
              'Optional arg must be a CHAR, CELLSTR, or numeric vector');
    end;
    args(1) = [];
    nargs = nargs - 1;
  end;

  if ( isempty(flds) )
    flds = fieldnames(stn);
  end;

  for fldix = 1:numel(flds)
    fld = flds{fldix};
    if ( ~is_valid_ts(stn.(fld)) )
      disp(['Skipping STN.',fld]);

    else

      fldnm = strrep(fld,'_','\_');

      dts = stn.(fld).date;
      dat = stn.(fld).date;
      fh = fmg;
      lh = plot_ts(stn.(fld));
      titlename(['BRUSH points to remove from STN.',fldnm]);

      set(lh,'XDataSource','dts');
      set(lh,'YDataSource','dat');

      if ( isempty(jump_initial_frame) )
        if ( ~isinf(frame_size) )
          if ( frame_size >= 1 )
            jump_initial_frame = [dts(1),dts(1)+frame_size];
          elseif ( frame_size > 0 )
            jump_initial_frame = [dts(1),dts(floor(end*frame_size))];
          else
            error('Invalid FRAME SIZE %g',frame_size);
          end;
        end;
      end;
      if ( ~isempty(jump_initial_frame) )
        xlim(jump_initial_frame); 
        qcstation_jumpaction;
      end;

      brush(fh,'on');
      h = brush(fh);

      skip = false;
      if ( ~isempty(jump_initial_frame) )
        [ax,fh] = jumpax([],[],@qcstation_jumpaction);
      else
        inp = input('Enter to save and move to next field, or ''q''uit or ''s''kip','s');
        if ( strncmpi(inp,'q',1) )
          done = true;
        elseif ( strncmpi(inp,'s',1) )
          skip = true;
        end;
      end;
      if ( ishandle(fh) )
        close(fh);
      end;

keyboard;
      if ( ~skip )
        if ( numel(dts) ~= numel(stn.(fld).date) || any(dat ~= stn.(fld).data) )
          disp(['Modifying STN.',fld]);
          stn.(fld).date = dts;
          stn.(fld).data = dat;
        end;
      end;

    end;

  end;    

return;


function qcstation_jumpaction(ax)
  if ( exist('ax','var') && ishandle(ax) )
    ylim(ax,'default');
    datetick3(ax);
  else
    ylim('default');
    datetick3;
  end;
return;











        initial_frame_dts = [dts(1),initial_frame_dts];







      if ( strncmpi(inp,'q',1) )
        break;
      end;









clear d ix m mx str



msix=find(ismember({stn.ps.station_name},{stn.ms.station_name}));
min_d = 5/111e3;
for ix=1:numel(stn.mlons);
  pix = find( abs(stn.plons(msix)-stn.mlons(ix))<min_d & ...
              abs(stn.plats(msix)-stn.mlats(ix))<min_d );
  if ( isempty(pix) )
    disp(stn.ps(ix).station_name);
  end;
end;







1;

if ( ~exist('stn','var') || ~isfield(stn,'ps') )
  stn=[]; clear stn
  stn = extract_xaymara_sites;
end;

%for stix = 1:stn.nps
for stix = 1:3
  disp([stn.ps(stix).station_name,':',stn.ps(stix).station_desc]);
  plot_hires_bathymetry(stn.ps(stix),-[0:1:50]);
  set(gca,'CLim',[-30,0]);
  axis(axis);

  plot(stn.lons,stn.lats,'ws','MarkerFaceColor','w');

  plot(stn.mlons,stn.mlats,'k.','MarkerSize',0.5);
  text(stn.mlons,stn.mlats,'\leftarrow M','Color','w','HorizontalAlignment','left','Rotation',[0:stn.nms-1]*7);

  plot(stn.plons,stn.plats,'r.','MarkerSize',0.5);
  text(stn.plons,stn.plats,'P \rightarrow ','Color','w','HorizontalAlignment','right','Rotation',[0:stn.nps-1]*7);
end;

clear ans stix




  %stn.ms = repmat(struct,[stn.nms,1]);
  %stn.ps = repmat(struct,[stn.nps,1]);




    % %x = plot_hires_bathymetry(x,-[0:1:50],[4e3,4e3],true,[],[],[],true);
    % x = plot_hires_bathymetry(x,-[0:1:50],[4e3,4e3],true,[],[],[],false);
    % set(gca,'CLim',[-50,0]);
    % axis(axis);
    % plot(stn.lons,stn.lats,'ws','MarkerFaceColor','w');
    % plot(stn.mlons,stn.mlats,'k.','MarkerSize',0.5);
    % plot(stn.plons,stn.plats,'r.','MarkerSize',0.5);







1;

datpath = get_ecoforecasts_path('data');
matfname = fullfile(datpath,'xaymara_sites.mat');

if ( exist(matfname,'file') )
  disp(['Loading ',matfname]);
  load(matfname);

else
tic,

  disp('Extracting Xaymara sites');
  stn.lons=[];
  stn.lats=[];

  stn.rawfname = fullfile(datpath,'Xaymara FL Keys sites for Lew.xlsx');
  rawin = importdata(stn.rawfname);

  stn.nms = size(rawin.data.FLKeysSitesMcav,1);
  %stn.ms = repmat(struct,[stn.nms,1]);

  stn.nps = size(rawin.data.FLKeysSitesPast,1);
  %stn.ps = repmat(struct,[stn.nps,1]);

  for stix = 1:stn.nms
    stn.ms(stix).station_name = rawin.textdata.FLKeysSitesMcav{stix+1,4};
    stn.ms(stix).station_desc = rawin.textdata.FLKeysSitesMcav{stix+1,3};
    stn.ms(stix).region = rawin.textdata.FLKeysSitesMcav{stix+1,1};
    stn.ms(stix).subregion = rawin.textdata.FLKeysSitesMcav{stix+1,2};
    stn.ms(stix).depth_category = rawin.textdata.FLKeysSitesMcav{stix+1,5};
    stn.ms(stix).lon = rawin.data.FLKeysSitesMcav(stix,3);
    stn.ms(stix).lat = rawin.data.FLKeysSitesMcav(stix,2);
    stn.ms(stix).depth = rawin.data.FLKeysSitesMcav(stix,1);
    stn.mlons(stix,1) = stn.ms(stix).lon;
    stn.mlats(stix,1) = stn.ms(stix).lat;
    stn.lons(end+1,1) = stn.ms(stix).lon;
    stn.lats(end+1,1) = stn.ms(stix).lat;
  end;

  for stix = 1:stn.nps
    stn.ps(stix).station_name = rawin.textdata.FLKeysSitesPast{stix+1,4};
    stn.ps(stix).station_desc = rawin.textdata.FLKeysSitesPast{stix+1,3};
    stn.ps(stix).region = rawin.textdata.FLKeysSitesPast{stix+1,1};
    stn.ps(stix).subregion = rawin.textdata.FLKeysSitesPast{stix+1,2};
    stn.ps(stix).depth_category = rawin.textdata.FLKeysSitesPast{stix+1,5};
    stn.ps(stix).lon = rawin.data.FLKeysSitesPast(stix,3);
    stn.ps(stix).lat = rawin.data.FLKeysSitesPast(stix,2);
    stn.ps(stix).depth = rawin.data.FLKeysSitesPast(stix,1);
    stn.plons(stix,1) = stn.ps(stix).lon;
    stn.plats(stix,1) = stn.ps(stix).lat;
    stn.lons(end+1,1) = stn.ps(stix).lon;
    stn.lats(end+1,1) = stn.ps(stix).lat;
  end;

  rawin=[]; clear rawin stix

  % disp(['NOT Saving ',matfname]);
  % %save(matfname,'nmstns','npstns','mstns','mlons','mlats','pstns','plons','plats','lons','lats','-v7.3');


  % Map wide-view bathymetry
  stn.lon = mean([min(stn.lons),max(stn.lons)]);
  stn.lat = mean([min(stn.lats),max(stn.lats)]);
  [ig,minix] = min(stn.lons);
  [ig,maxix] = max(stn.lons);
  stn.rady = ceil(max([distance_wgs84(stn.lats(minix),stn.lons(minix),stn.lat,stn.lons(minix)),...
                      distance_wgs84(stn.lats(maxix),stn.lons(maxix),stn.lat,stn.lons(maxix))]))*1e3;
  [ig,minix] = min(stn.lats);
  [ig,maxix] = max(stn.lats);
  stn.radx = ceil(max([distance_wgs84(stn.lats(minix),stn.lons(minix),stn.lats(minix),stn.lon),...
                      distance_wgs84(stn.lats(maxix),stn.lons(maxix),stn.lats(maxix),stn.lon)]))*1e3;

  clear ig minix maxix

  stn = plot_hires_bathymetry(stn,-[0:5:50,100:100:800],[stn.radx*1.05,stn.rady*1.10],true,[],[],[],false);
  axis(axis);
  plot(stn.lons,stn.lats,'ws','MarkerFaceColor','w');
  plot(stn.mlons,stn.mlats,'k.','MarkerSize',0.5);
  plot(stn.plons,stn.plats,'r.','MarkerSize',0.5);

  [LON,LAT] = meshgrid(stn.ngdc_hires_bathy.lon,stn.ngdc_hires_bathy.lat);
  [stn.aspect,stn.beta_deg,stn.beta_y,stn.beta_x] = gradientm(LAT,LON,stn.ngdc_hires_bathy.field,wgs84Ellipsoid);
  LON=[]; LAT=[]; clear LON LAT
  stn.beta = sqrt((stn.beta_y.^2)+(stn.beta_x.^2));


  for stix = 1:stn.nms
    x.lon = stn.ms(stix).lon;
    x.lat = stn.ms(stix).lat;
    x = read_hires_bathymetry(x,[4e3,4e3],[],false);
    %x = plot_hires_bathymetry(x,-[0:1:50],[4e3,4e3],true,[],[],[],true);
    x = plot_hires_bathymetry(x,-[0:1:50],[4e3,4e3],true,[],[],[],false);
    set(gca,'CLim',[-50,0]);
    axis(axis);
    plot(stn.lons,stn.lats,'ws','MarkerFaceColor','w');
    plot(stn.mlons,stn.mlats,'k.','MarkerSize',0.5);
    plot(stn.plons,stn.plats,'r.','MarkerSize',0.5);
    stn.mdeps(stix) = stn.ms(stix).depth;
    stn.mndeps(stix) = interp2(stn.ngdc_hires_bathy.lon,stn.ngdc_hires_bathy.lat,...
                               stn.ngdc_hires_bathy.field,...
                               stn.ms(stix).lon,stn.ms(stix).lat);
    stn.ms(stix).ngdc_hires_bathy = x.ngdc_hires_bathy;
    stn.ms(stix).ngdc_depth = stn.mndeps(stix);

    stn.mbetas(stix) = interp2(stn.ngdc_hires_bathy.lon,stn.ngdc_hires_bathy.lat,...
                               stn.beta,stn.ms(stix).lon,stn.ms(stix).lat);
    stn.ms(stix).ngdc_beta = stn.mbetas(stix);
    x=[]; clear x
  end;

  for stix = 1:stn.nps
    x.lon = stn.ps(stix).lon;
    x.lat = stn.ps(stix).lat;
    x = read_hires_bathymetry(x,[4e3,4e3],[],false);
    % %x = plot_hires_bathymetry(x,-[0:1:50],[4e3,4e3],true,[],[],[],true);
    % x = plot_hires_bathymetry(x,-[0:1:50],[4e3,4e3],true,[],[],[],false);
    % set(gca,'CLim',[-50,0]);
    % axis(axis);
    % plot(stn.lons,stn.lats,'ws','MarkerFaceColor','w');
    % plot(stn.mlons,stn.mlats,'k.','MarkerSize',0.5);
    % plot(stn.plons,stn.plats,'r.','MarkerSize',0.5);
    stn.pdeps(stix) = stn.ps(stix).depth;
    stn.pndeps(stix) = interp2(stn.ngdc_hires_bathy.lon,stn.ngdc_hires_bathy.lat,...
                               stn.ngdc_hires_bathy.field,...
                               stn.ps(stix).lon,stn.ps(stix).lat);

    stn.ps(stix).ngdc_hires_bathy = x.ngdc_hires_bathy;
    stn.ps(stix).ngdc_depth = stn.pndeps(stix);

    stn.pbetas(stix) = interp2(stn.ngdc_hires_bathy.lon,stn.ngdc_hires_bathy.lat,...
                               stn.beta,stn.ps(stix).lon,stn.ps(stix).lat);
    stn.ps(stix).ngdc_beta = stn.pbetas(stix);
    x=[]; clear x
  end;

  clear stix

  disp(['NOT Saving ',matfname]);
  %save(matfname,'stn','-v7.3');

toc,

end; %if ( exist(matfname,'file') )











  rawin = importdata(fullfile(datpath,'Xaymara FL Keys sites for Lew.xlsx'));







  plot(lon,lat,'kp','MarkerFaceColor','w');




  % plot(lon,lat,'kp','MarkerFaceColor','w');







for stix = 1:nmstns
  x.lon = mstns(stix).lon; x.lat = mstns(stix).lat;
  %x = plot_hires_bathymetry(x,-[0:1:30],[4e3,4e3],true,[],[],[],true);
  x = plot_hires_bathymetry(x,-[0:1:30],[4e3,4e3],true,[],[],[],false);
  set(gca,'CLim',[-30,0]);
  axis(axis);
  plot(lon,lat,'kp','MarkerFaceColor','w');
  ndeps(stix) = interp2(x.ngdc_hires_bathy.lon,x.ngdc_hires_bathy.lat,x.ngdc_hires_bathy.field,mstns(stix).lon,mstns(stix).lat);
  mstns(stix).ngdc_hires_bathy = x.ngdc_hires_bathy;
  mstns(stix).ngdc_depth = ndeps(stix);
  deps(stix) = mstns(stix).depth;
  x=[]; clear x
end;








  lon(nmstns+stix,1) = pstns(stix).lon;
  lat(nmstns+stix,1) = pstns(stix).lat;






  % Plot station location as a white pentagram
  if ( isfield(stn,'lon') && isfield(stn,'lat') )
    stnlon = stn.lon;
    stnlat = stn.lat;
    plot(stnlon,stnlat,'wp', 'MarkerEdgeColor','black', 'MarkerFaceColor','white');
  else
    stnlon = mean(stn.ngdc_hires_bathy.lon(:));
    stnlat = mean(stn.ngdc_hires_bathy.lat(:));
    plot(stnlon,stnlat,'wo', 'MarkerEdgeColor','black', 'MarkerFaceColor','white');
  end;








  [c,h] = contour(mdl.lon,mdl.lat,squeeze(mdl.t(1,1,:,:)),[floor(min_t):0.25:ceil(max_t)]);







  % % stn = plot_hires_bathymetry(stn,-[0:5:80],[12e3,12e3],true,[],false,[],false);
  % % fh = gcf;
  % % axis(axis);
  % % [c,h] = contour(mdl.lon,mdl.lat,squeeze(mdl.t(1,1,:,:)),[floor(min_t):0.25:ceil(max_t)]);
  % % clabel(c,h);
  % % c=[]; clear c h










%stnm = 'fwyf1';
%stnm = 'mlrf1';
%stnm = 'sanf1';
% if ( stnm )
for stix=1:ncstnms
  stnm = cstnms{stix};
  disp(stnm);

  %fmg;
  stn = []; clear stn
  stn = get_station_from_station_name(stnm);
  stn = station_optimal_isobath_orientation(stn);
  stn = read_hires_bathymetry(stn,[12e3,12e3],[],false);
  % stn = plot_hires_bathymetry(stn,-[0:5:80],[12e3,12e3],true,[],false,[],false);
  % fh = gcf;
  % axis(axis);
  % [c,h] = contour(mdl.lon,mdl.lat,squeeze(mdl.t(1,1,:,:)),[floor(min_t):0.25:ceil(max_t)]);
  % clabel(c,h);
  % c=[]; clear c h

  %fmg;
  [res,stn] = plot_bathy_transect(stn,10,stn.isobath_orientation+90,'ngdc_hires_bathy');
  [dx,az] = distance_wgs84(stn.lat,stn.lon,res.lat,res.lon);
  dx(round(az) ~= stn.isobath_orientation+90) = -dx(round(az) ~= stn.isobath_orientation+90);
  res.dx = dx;

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

  [c,h] = contourf(res.dx,-mdl.z,t,[11.0:0.5:24.5]);
  set(gca,'CLim',[11.0,24.5]);
  clabel(c,h);
  c=[]; clear c h
  plot(res.dx,res.field,'k-','LineWidth',2);
  xlim([-10,+10]); ylim([-140,0]);
  titlename([mdl.model,': ',upper(stnm),' bathymetry vs. cross-shore T (1 week mean)']);
  if ( doPrint )
    print('-dpng',fullfile(figpath,[mdl.model,'.',stnm,'.section_t_mean.png']));
  end;

  %fmg;
  plot_bathy_transect(stn,10,stn.isobath_orientation+90,'ngdc_hires_bathy');
  [c,h] = contourf(res.dx,-mdl.z,l,[0.0:0.01:+1.0]);
  set(gca,'CLim',[0.0,+1.0]);
  clabel(c,h);
  c=[]; clear c h
  plot(res.dx,res.field,'k-','LineWidth',2);
  xlim([-10,+10]); ylim([-140,0]);
  titlename([mdl.model,': ',upper(stnm),' bathymetry vs. cross-shore V (1 week mean)']);

  % % Draw cross-shore section on bathymetry map
  % figure(fh);
  % plot(res.lon,res.lat,'k-','LineWidth',2);
  % clear fh

  break;
end; %for stix=1:ncstnms

















    xlim([min_hist,max_hist]);



    xlim([-max_spd,+max_spd]);
    xlim([-max_spd,+max_spd]);



    hist(mdl.(stnm).l.data,min_hist:0.025:max_hist); 



    hist(mdl.(stnm).x.data,floor(-max_spd):0.025:ceil(+max_spd)); 





    %spt(2,floor(ncstnms/2),((stix-1)*ncols)+1); 




    hist(mdl.(stnm).x.data);





if 0;
  for stix=1:ncstnms
    fmg;
    stnm = cstnms{stix};
    spt(2,1,1); hist(mdl.(stnm).l.data); xlabel('Alongshore (m/s)');
    xlim([-max_spd,+max_spd]);
    titlename([upper(stnm),' currents']);
    spt(2,1,2); hist(mdl.(stnm).x.data); xlabel('Cross-shore (m/s)');
    xlim([-max_spd,+max_spd]);
    if ( doPrint )
      print('-dpng',fullfile(figpath,[mdl.model,'.',stnm,'.uv_hist.png']));
    end;
  end;
end;







  % For this analysis, restrict ourselves to just the FRT reefs
  mdl.lon_range = [-83.0,-79.5];
  mdl.lat_range = [+24.0,+27.5];
  % "Surface" (0 m) eFKEYS layer is almost identical to 1 m layer
  %mdl.z_range = [1,80];
  %mdl.z_range = [1,120];
  mdl.z_range = [0,140];






        % mdl.u(dtix,:,:,:) = squeeze(cast(nc{'u'}(:,:,:,:),'double'));
        % mdl.v(dtix,:,:,:) = squeeze(cast(nc{'v'}(:,:,:,:),'double'));
        % mdl.w(dtix,:,:,:) = squeeze(cast(nc{'w_velocity'}(:,:,:,:),'double'));
        % mdl.t(dtix,:,:,:) = squeeze(cast(nc{'temperature'}(:,:,:,:),'double'));
        % mdl.s(dtix,:,:,:) = squeeze(cast(nc{'salinity'}(:,:,:,:),'double'));





          % "Surface" (0 m) layer is almost identical to 1 m layer
          %mdl.nzix = 1:5;
          %mdl.nzix = 2:7;
          mdl.nzix = 2:9;
          mdl.z = mdl.true_z(mdl.nzix);







          mdl.true_lon = cast(nc{'Longitude'}(:),'double');
          mdl.true_lat = cast(nc{'Latitude'}(:),'double');
          mdl.true_nlon = numel(mdl.true_lon);
          mdl.true_nlat = numel(mdl.true_lat);

          %mdl.z = cast(nc{'Depth'}(:),'double');
          %mdl.nz = numel(mdl.z);
          mdl.true_z = cast(nc{'Depth'}(:),'double');
          mdl.true_nz = numel(mdl.true_z);
          % "Surface" (0 m) layer is almost identical to 1 m layer
          %mdl.nzix = 1:5;
          %mdl.nzix = 2:7;
          mdl.nzix = 2:9;
          mdl.z = mdl.true_z(mdl.nzix);





function idx = find_pct_good_dates(dts,grpfun,pct,goodfun)
%function idx = find_pct_good_dates(dts,grpfun,pct[,goodfun])
%
% Return all indices of DATENUM vector DTS, for which there are at least PCT
% values present in each time period: time periods are deterined by GRPFUN,
% as UNIQUE(GRPFUN(DTS)). Returns indices of all members of DTS which are in
% periods that have at least PCT good dates; also optionally returns

  if ( ~isnumeric(dts) || ~isvector(dts) )
    error('DTS must be a vector of DATENUM');
  end;
  if ( ~isa(grpfun,'function_handle') )
    error('GRPFUN must be a FUNCTION_HANDLE accepting a DATENUM vector');
  end;
  if ( 0 > pct || pct > 1 )
    error('PCT must be decimal between 0 and 1');
  end;

  d = grpfun(dts);
  ud = unique(d);

  dt = median(diff(dts));
  alldts = ud(1):dt:ud(end);
  alld = grpfun(alldts);
  ualld = unique(alld);

  % NOTE: The last period is unreliable - always include it for now
  idx = 1:numel(d);
  for uix = 1:numel(ud)-1
    dix = find(ismember(d,ud(uix)));
    allix = find(ismember(alld,ud(uix)));
    if ( (numel(dix)/numel(allix)) >= pct )
      idx(end+1:end+numel(dix)) = dix;
    end;
  end;
  dix = find(ismember(d,ud(end)));
  idx(end+1:end+numel(dix)) = dix;

return;






find_pct_good_dates


for stix=1:numel(cstnms)
  stnm = cstnms{stix};
  fmg;
  spt(2,1,1); hist(mdl.(stnm).l.data); xlabel('Alongshore');
  xlim([-max_spd,+max_spd]); %xlim([-0.4,+0.4]);
  titlename([upper(stnm),' currents']);
  spt(2,1,2); hist(mdl.(stnm).x.data); xlabel('Cross-shore');
  xlim([-max_spd,+max_spd]); %xlim([-0.4,+0.4]);
end;






  mdl.(stnm).t0.z = mdl.z(1);
  mdl.(stnm).t0.date = mdl.d;
  mdl.(stnm).t0.data = mdl.t(:,1,mdl.(stnm).latix,mdl.(stnm).lonix);
  mdl.(stnm).t1.z = mdl.z(2);
  mdl.(stnm).t1.date = mdl.d;
  mdl.(stnm).t1.data = mdl.t(:,2,mdl.(stnm).latix,mdl.(stnm).lonix);
  mdl.(stnm).t10.z = mdl.z(3);
  mdl.(stnm).t10.date = mdl.d;
  mdl.(stnm).t10.data = mdl.t(:,3,mdl.(stnm).latix,mdl.(stnm).lonix);
  mdl.(stnm).t20.z = mdl.z(4);
  mdl.(stnm).t20.date = mdl.d;
  mdl.(stnm).t20.data = mdl.t(:,4,mdl.(stnm).latix,mdl.(stnm).lonix);








dx = distance_wgs84(stn.lat,stn.lon,res.lat,res.lon);
dx(1:floor(numel(dx)/2)) = -dx(1:floor(numel(dx)/2));
res.dx = dx;



  min_t = nanmin([min_t,nanmin(mdl.(stnm).t0.data(:)),nanmin(mdl.(stnm).t20.data(:))]);
  max_t = nanmax([max_t,nanmax(mdl.(stnm).t0.data(:)),nanmax(mdl.(stnm).t20.data(:))]);
  min_curr = nanmin([min_curr,nanmin(mdl.(stnm).u.prof(:)),nanmin(mdl.(stnm).v.prof(:))]);
  max_curr = nanmax([max_curr,nanmax(mdl.(stnm).u.prof(:)),nanmax(mdl.(stnm).v.prof(:))]);



    disp(datestr(datenum(yr,1,0) + jd));




          disp(['PRE-Saving ',matfname]);
          save(matfname,'mdl');








clear lon lat z d u v w t s





   case 'pibhmc_bathy_5m_rota',
    corners = [ ...
        147,173 , 549,173 , 549,389 , 147,389 ; ...
        549,247 , 621,247 , 621,360 , 549,360 ; ...
        065,076 , 396,076 , 396,273 , 065,273 ; ...
        386,124 , 433,124 , 433,185 , 386,185 ; ...
        539,200 , 593,200 , 593,255 , 539,255 ; ...
              ];



figure(fh);
r = getrect;
b = rect2bbox(r);
ixes = bbox2ind(b,LON,LAT);
corners(end+1,1:8) = [ixes(1,1),ixes(1,2), ixes(2,1),ixes(2,2), ixes(3,1),ixes(3,2), ixes(4,1),ixes(4,2)];

rectangle('Position',r);
text(b(2),b(4),[' \leftarrow ',num2str(size(corners,1))]);

fill_dem_bathy;

if ( exist('fh2','var') && ishandle(fh2) )
  close(fh2);
end;
clear fh2;
figure(fh);

if ( startOver )
  startOver = false;
end;










1;

% Sample: c:/Users/gramer/Documents/rsmas/Coastal/thesis/data/hycom/eFKEYS_01-1-7-2012/eFKEYS_archv.2012_001_12_3zuvwts.nc

datpath = get_thesis_path('../data/hycom/eFKEYS_01-1-7-2012');

lon=[];
lat=[];
z=[];
d=[];
u=[];
v=[];
w=[];
t=[];
s=[];

dts = datenum(2012,1,1,[6:6:(7*24)-1],0,0);

for yr=2012;
  for jd=1:7;
    disp(datestr(datenum(yr,1,0) + jd));
    for hr = 0:6:18;
      ncname = sprintf('eFKEYS_archv.%04d_%03d_%02d_3zuvwts.nc',yr,jd,hr);
      ncpath = fullfile(datpath,ncname);
      if ( ~exist(ncpath,'file') )
        warning('Skipping %s',ncpath);
      else
        d(end+1,1) = datenum(yr,1,0,hr,0,0) + jd;
        try,
          nc = mDataset(ncpath);
          if ( isempty(lon) )
            lon = cast(nc{'Longitude'}(:),'double');
            lat = cast(nc{'Latitude'}(:),'double');
            z = cast(nc{'Depth'}(:),'double');
            nlon = numel(lon);
            nlat = numel(lat);
            nz = numel(z);
          end;
          u(end+1,1:nz,1:nlat,1:nlon) = squeeze(cast(nc{'u'}(:,:,:,:),'double'));
          v(end+1,1:nz,1:nlat,1:nlon) = squeeze(cast(nc{'v'}(:,:,:,:),'double'));
          w(end+1,1:nz,1:nlat,1:nlon) = squeeze(cast(nc{'w_velocity'}(:,:,:,:),'double'));
          t(end+1,1:nz,1:nlat,1:nlon) = squeeze(cast(nc{'temperature'}(:,:,:,:),'double'));
          s(end+1,1:nz,1:nlat,1:nlon) = squeeze(cast(nc{'salinity'}(:,:,:,:),'double'));
        catch,
          catchwarn;
        end;
        close(nc); clear nc;
      end;
    end;
  end;
end;







  % Tops(rix) = max(LAT(ixes(:,2)));
  % Rights(rix) = max(LON(ixes(:,1)));
  % Bottoms(rix) = min(LAT(ixes(:,2)));
  % Lefts(rix) = max(LON(ixes(:,1)));




  else
    maxerr = max(diff(unique(lon)));
    [err,ix(1)] = min(abs(lon-bbox(1)));
    if (err>maxerr+eps); warning('MINX is outside region'); end;
    [err,ix(2)] = min(abs(lon-bbox(2)));
    if (err>maxerr+eps); warning('MAXX is outside region'); end;

    maxerr = max(diff(unique(lat)));
    [err,ix(3)] = min(abs(lat-bbox(3)));
    if (err>maxerr+eps); warning('MINY is outside region'); end;
    [err,ix(4)] = min(abs(lat-bbox(4)));
    if (err>maxerr+eps); warning('MAXY is outside region'); end;
  end;







function ix = bbox2ind(bbox,lon,lat)
%function ix = bbox2ind(bbox,lon,lat)
% Return a 4-vector of indices IX in LON,LAT for Bounding BOX [MINX,MAXX,MINY,MAXY]
% If LON and LAT are not vectors, their non-singleton dimensions must be equal.

  if ( bbox(2)<=bbox(1) || bbox(4)<=bbox(3) )
    error('BBOX should have form [MINX,MAXX,MINY,MAXY]');
  end;
  if ( ~isnumeric(lon) || ~isnumeric(lat) )
    error('LON and LAT should be numeric matrices');
  end;

  if ( ~isvector(lon) || ~isvector(lat) )
    % Remove any singleton dimensions
    lon = squeeze(lon);
    lat = squeeze(lat);
    if ( ndims(lon) ~= ndims(lat) || any(size(lon) ~= size(lat)) )
      error('If matrices, LON and LAT should be equal-sized');
    end;

    warning('Unraveling LON and LAT: use IND2SUB(SIZE(LON),IX)');

    maxerr_ul = abs(lon(1,1)-lon(2,2)) + abs(lat(1,1)-lat(2,2));
    maxerr_lr = abs(lon(end-1,end-1)-lon(end,end)) + abs(lat(end-1,end-1)-lat(end,end));
    maxerr = max(maxerr_ul,maxerr_lr);

    [err,ix(1)] = min( abs(lon(:)-bbox(1)) + abs(lat(:)-bbox(3)) );
    if (err>maxerr+eps); warning('MINX,MINY may be outside region'); end;
    [err,ix(2)] = min( abs(lon(:)-bbox(2)) + abs(lat(:)-bbox(3)) );
    if (err>maxerr+eps); warning('MAXX,MINY may be outside region'); end;
    [err,ix(3)] = min( abs(lon(:)-bbox(1)) + abs(lat(:)-bbox(4)) );
    if (err>maxerr+eps); warning('MINX,MAXY may be outside region'); end;
    [err,ix(4)] = min( abs(lon(:)-bbox(2)) + abs(lat(:)-bbox(4)) );
    if (err>maxerr+eps); warning('MAXX,MAXY may be outside region'); end;


  else
    maxerr = max(diff(unique(lon)));
    [err,ix(1)] = min(abs(lon-bbox(1)));
    if (err>maxerr+eps); warning('MINX is outside region'); end;
    [err,ix(2)] = min(abs(lon-bbox(2)));
    if (err>maxerr+eps); warning('MAXX is outside region'); end;

    maxerr = max(diff(unique(lat)));
    [err,ix(3)] = min(abs(lat-bbox(3)));
    if (err>maxerr+eps); warning('MINY is outside region'); end;
    [err,ix(4)] = min(abs(lat-bbox(4)));
    if (err>maxerr+eps); warning('MAXY is outside region'); end;
  end;

return;








  if ( ~isnumeric(lon) || ~isnumeric(lat) || ~isvector(lon) || ~isvector(lat) )
    error('LON and LAT should be numeric vectors');
  end;





% % 'pibhmc_bathy_5m_saipan',
% Tops = [653,673,761,653,394,604,223,];
% Rights = [779,926,888,913,116,436,362,];
% Bottoms = [255,626,691,611,375,559,145,];
% Lefts = [480,872,816,869, 99,357,291,];




      280,336  ,  472,336  ,  472,494  ,  280,494 ; ...




switch ( batnm ),
 case 'pibhmc_bathy_5m_saipan',
  % % REGION CORNERS:
  % 480,255 / 779,255 / 779,653 / 480,653
  % 872,626 / 926,626 / 926,673 / 872,673
  % 816,691 / 888,691 / 888,761 / 816,761
  % 869,611 / 913,611 / 913,653 / 869,653
  %  99,375 / 116,375 / 116,394 /  99,394
  % 357,559 / 436,559 / 436,604 / 357,604
  % 291,145 / 362,145 / 362,223 / 291,223
  Tops = [653,673,761,653,394,604,223,];
  Rights = [779,926,888,913,116,436,362,];
  Bottoms = [255,626,691,611,375,559,145,];
  Lefts = [480,872,816,869, 99,357,291,];

 case 'pibhmc_bathy_5m_tinian.PARTIAL',
  % % REGION CORNERS:
  % 264,669 / 483,669 / 483,916 / 264,916
  % 170,581 / 513,581 / 513,680 / 170,680
  % 196,495 / 485,495 / 495,581 / 196,581

end;

for rix=1:size(corners,1)
  Tops(end+1) = corners(rix,6);
  Rights(end+1) = corners(rix,3);
  Bottoms(end+1) = corners(rix,2);
  Lefts(end+1) = corners(rix,2);
end;





    fmg; contourf(flon,flat,fdat,-1000:50:100); colorbar; titlename(['BEFORE ',num2str(ix)]);
    fdat(nanlandix(badlandix)),
    fmg; contourf(flon,flat,fdat,-1000:50:100); colorbar; titlename(['AFTER ',num2str(ix)]);
    fdat(nanlandix(badlandix)),







x=tic;
toc(x); clear x



tic,
disp('FIND'); toc,
tic,
disp('SUM'); toc,



tic,
disp('RESHAPE1'); toc,
tic,
disp('Trim'); toc,
tic,
disp('RESHAPE2'); toc,





nix=652:699; tix=260:294; n=lon(nix); t=lat(tix); d=dat(tix,nix); [N,T]=meshgrid(n,t); % THIS IS THE FULL BOX


nix=652:699; tix=240:294; n=lon(nix); t=lat(tix); d=dat(tix,nix); [N,T]=meshgrid(n,t); % *DEEP* BOX
nd=interp2(n,t,d,N(isnan(d)),T(isnan(d)),'spline',nan);


nix=652:699; tix=100:294; n=lon(nix); t=lat(tix); d=dat(tix,nix); [N,T]=meshgrid(n,t); % *DEEPER* BOX
nd=interp2(n,t,d,N(isnan(d)),T(isnan(d)),'spline',nan);

nix=652:699; tix=20:294; n=lon(nix); t=lat(tix); d=dat(tix,nix); [N,T]=meshgrid(n,t); % JUST OK
nd=interp2(n,t,d,N(isnan(d)),T(isnan(d)),'spline',nan);

nix=652:699; tix=19:294; n=lon(nix); t=lat(tix); d=dat(tix,nix); [N,T]=meshgrid(n,t); % BAD
nd=interp2(n,t,d,N(isnan(d)),T(isnan(d)),'spline',nan);


nix=300:699; tix=20:294; n=lon(nix); t=lat(tix); d=dat(tix,nix); [N,T]=meshgrid(n,t); % JUST OK
nd=interp2(n,t,d,N(isnan(d)),T(isnan(d)),'spline',nan);

nix=299:699; tix=20:294; n=lon(nix); t=lat(tix); d=dat(tix,nix); [N,T]=meshgrid(n,t); % BAD
nd=interp2(n,t,d,N(isnan(d)),T(isnan(d)),'spline',nan);






  %rawdat = interp2(dem.lon,dem.lat,dem.dat,RAWLON,RAWLAT,'cubic',nan);




basnm = 'cnmi-saipan_c3-ssl-REVERSE';
demnm = 'usgs_dem_10m_saipan';
batnm = 'pibhmc_bathy_5m_saipan';

cmbpath = fullfile(get_ecoforecasts_path('coast'),[batnm,'.dods.',num2str(target_res),'m.',demnm,'.dods.mat']);

if ( exist(cmbpath,'file') )
  disp(['Loading ',cmbpath]);
  load(cmbpath);
  
else

  dempath = fullfile(get_ecoforecasts_path('coast'),[demnm,'.dods.mat']);
  batpath = fullfile(get_ecoforecasts_path('coast'),[batnm,'.dods.mat']);

  if ( ~isempty(basnm) )
    baspath = fullfile(get_ecoforecasts_path('coast'),[basnm,'.mat']);
    if ( ~exist(baspath,'file') )
      error(['No base bathymetry map?? ',baspath]);
    end;
  end;




  bas = [];
  if ( ~isempty(baspath) )
    bas = load(baspath);
  end;




  % Up- (or down-)sample base and elevation model to bathymetry resolution and spatial coverage
  rawlon = bat.lon;
  rawlat = bat.lat;
  [RAWLON,RAWLAT]=meshgrid(rawlon,rawlat);

  basdat = interp2(bas.lon,bas.lat,bas.dat,RAWLON,RAWLAT,'linear',nan);














  newlon = lon(1:npts(1):end);
  newlat = lat(1:npts(2):end);
  [LON,LAT] = meshgrid(newlon,newlat);
tic,
  newdat = interp_field(lat,lon,dat,LAT,LON,{@nanmean,npts(1),npts(2)});
toc,
  LON=[]; LAT=[]; clear LON LAT




if ( min(diff(unique(bat.lon))) < min(diff(unique(bat.lon)))-eps )
  % If Bathymetry has finer resolution
  lon = dem.lon;
  lat = dem.lat;
  [LON,LAT]=meshgrid(lon,lat);
  %dat = interp2(bat.lon,bat.lat,bat.dat,LON,LAT,'nearest');
  dat = interp2(bat.lon,bat.lat,bat.dat,LON,LAT);
  LON=[]; LAT=[]; clear LON LAT
  dat(isnan(dat)) = dem.dat(isnan(dat));
else
  % If DEM has finer resolution
  lon = bat.lon;
  lat = bat.lat;
  [LON,LAT]=meshgrid(lon,lat);
  dat = interp2(dem.lon,dem.lat,dem.dat,LON,LAT);
  LON=[]; LAT=[]; clear LON LAT
  dat(isnan(dat)) = bat.dat(isnan(dat));
end;








% Control variable
cX = { ...
    'Ta', 'ndbc_air_t' ;...
    'ERAI Ta', 'erai_air_t' ;...
    'ERAI qa', 'erai_spechumid' ;...
    'ERAI QswI', 'erai_dsrf' ;...
    %'ERAI QlwI', 'erai_dlrf' ;...
    'ERAI Q0', 'erai_net_heat_flux' ;...
    %'NARR Ta', 'ncep_air_t' ;...
    %'NARR qa', 'ncep_spechumid' ;...
    %'NARR QswI', 'ncep_dsrf' ;...
    %'NARR QlwI', 'ncep_dlrf' ;...
    %'NARR Q0', 'ncep_net_heat_flux' ;...
     };



% Control variable
cX = { ...
    'Ta 1d', 'ndbc_air_t_1_d_avg' ;...
    'ERAI Ta 1d', 'erai_air_t_1_d_avg' ;...
    'ERAI qa 1d', 'erai_spechumid_1_d_avg' ;...
    'ERAI QswI 1d', 'erai_dsrf_1_d_avg' ;...
    %'ERAI QlwI 1d', 'erai_dlrf_1_d_avg' ;...
    'ERAI Q0 1d', 'erai_net_heat_flux_1_d_avg' ;...
    %'NARR Ta 1d', 'ncep_air_t_1_d_avg' ;...
    %'NARR qa 1d', 'ncep_spechumid_1_d_avg' ;...
    %'NARR QswI 1d', 'ncep_dsrf_1_d_avg' ;...
    %'NARR QlwI 1d', 'ncep_dlrf_1_d_avg' ;...
    %'NARR Q0 1d', 'ncep_net_heat_flux_1_d_avg' ;...
     };







cX = { ...
    'Ta', 'ndbc_air_t' ;...
    'ERAI Ta', 'erai_air_t' ;...
    'ERAI qa', 'erai_spechumid' ;...
    'ERAI QswI', 'erai_dsrf' ;...
    %'ERAI QlwI', 'erai_dlrf' ;...
    'ERAI Q0', 'erai_net_heat_flux' ;...
    %'NARR Ta', 'ncep_air_t' ;...
    %'NARR qa', 'ncep_spechumid' ;...
    %'NARR QswI', 'ncep_dsrf' ;...
    %'NARR QlwI', 'ncep_dlrf' ;...
    %'NARR Q0', 'ncep_net_heat_flux' ;...
     };



cX = { ...
    'Ta 1d', 'ndbc_air_t_1_d_avg' ;...
    'ERAI Ta 1d', 'erai_air_t_1_d_avg' ;...
    'ERAI qa 1d', 'erai_spechumid_1_d_avg' ;...
    'ERAI QswI 1d', 'erai_dsrf_1_d_avg' ;...
    %'ERAI QlwI 1d', 'erai_dlrf_1_d_avg' ;...
    'ERAI Q0 1d', 'erai_net_heat_flux_1_d_avg' ;...
    %'NARR Ta 1d', 'ncep_air_t_1_d_avg' ;...
    %'NARR qa 1d', 'ncep_spechumid_1_d_avg' ;...
    %'NARR QswI 1d', 'ncep_dsrf_1_d_avg' ;...
    %'NARR QlwI 1d', 'ncep_dlrf_1_d_avg' ;...
    %'NARR Q0 1d', 'ncep_net_heat_flux_1_d_avg' ;...
     };




%cX = {'Ta', 'ndbc_air_t'};
%cX = {'ERAI Ta', 'erai_air_t'};
%cX = {'NARR Ta', 'ncep_air_t'};
cX = {'ERAI qa', 'erai_spechumid'};
%cX = {'NARR qa', 'ncep_spechumid'};
%cX = {'ERAI QswI', 'erai_dsrf'};
%cX = {'NARR QswI', 'ncep_dsrf'};






clear ans cix mo moix n p r N P R yr yrs





for cix=1:numel(cstnms)
  stnm = cstnms{cix};







    if ( ~isempty(strfind(X,'ndbc')) || ~isempty(strfind(Y,'ndbc')) )
      stns{cix} = load_all_ndbc_data(stns{cix});
    end;
    if ( ~isempty(strfind(X,'erai')) || ~isempty(strfind(Y,'erai')) )
      stns{cix} = get_erai_station(stns{cix});
      stns{cix} = adjust_erai_station(stns{cix});
    end;
    if ( ~isempty(strfind(X,'ncep')) || ~isempty(strfind(Y,'ncep')) )
      stns{cix} = get_ncep_station(stns{cix},'narr');
    end;









% Plot monthly correlations
for cix=1:numel(cstnms)
  stnm = cstnms{cix};
  for mo=1:12;
    moix = find(get_month(taas{cix}.date)==mo);
    N(mo) = numel(moix);
    R(mo) = corr2(taas{cix}.data(moix),tsas{cix}.data(moix));
  end;
  N,
  R,

  fmg;
  bar(1:12,R.^2);
  axis([0.5,12.5,0,1]);
  titlename([upper(stnm),' Ta vs. Ts R^2 ',num2str(min(yrs)),'-',num2str(max(yrs))]);
end;


% Plot monthly correlations for individual years
if 0;

  yrs = 1992:2004;

  for yr=yrs(:)'
    skipYr=false;
    for cix=1:numel(cstnms)
      for mo=1:12;
        moix = find(get_month(taass{cix}.date)==mo & get_year(taass{cix}.date)==yr);
        N(mo) = numel(moix);
        R(mo) = corr2(taass{cix}.data(moix),tsass{cix}.data(moix));
      end;
      N,
      R,

      if ( min(N) < 18 )
        disp('********************');
        disp(['Removing ',num2str(yr),' at ',cstnms{cix}]);
        disp('********************');
        yrs(yrs==yr) = [];
        skipYr = true;
        continue;
      else
        fmg;
        bar(1:12,R.^2);
        axis([0.5,12.5,0,1]);
        titlename([upper(cstnms{cix}),' Ta - Ts R^2 ',num2str(yr)]);
      end;
    end;
    if ( skipYr )
      continue;
    end;
  end;

  for yr=yrs(:)'
    for cix=1:numel(cstnms)
      for mo=1:12;
        moix = find(get_month(taass{cix}.date)==mo & get_year(taass{cix}.date)==yr);
        N(mo) = numel(moix);
        R(mo) = corr2(taass{cix}.data(moix),tsass{cix}.data(moix));
      end;
      N,
      R,

      if ( min(N) < 18 )
        disp('********************');
        disp(['Oops?! At ',cstnms{cix}]);
        disp('********************');
        continue;
      else
        fmg;
        bar(1:12,R.^2);
        axis([0.5,12.5,0,1]);
        titlename([upper(cstnms{cix}),' Ta - Ts R^2 ',num2str(yr)]);
      end;
    end;
  end;

end;

clear ans cix mo moix N R stnm yr yrs

set_more











for cix=1:numel(cstnms)

  stnm = cstnms{cix};

  stns{cix} = get_station_from_station_name(stnm);
  stns{cix} = load_all_ndbc_data(stns{cix});
  % Coincident dates for Ta and Ts at this site
  [ta,ts] = intersect_tses(stns{cix}.ndbc_air_t,stns{cix}.ndbc_sea_t);

  [taa,taclim,tatid,tad] = anomalize_ts(ta,@get_jhour_no_leap);
  [tsa,tsclim,tstid,tsd] = anomalize_ts(ts,@get_jhour_no_leap);

  fmg;
  plot(tatid,taclim,'r-',tatid,taclim-tad,'r:',tatid,taclim+tad,'r:');
  plot(tstid,tsclim,'b-',tstid,tsclim-tsd,'b:',tstid,tsclim+tsd,'b:');
  axis([0,366,13,34]); datetick3('x','mmm');
  titlename([upper(stnm),' climatological Ta vs. Ts']);

  tas{cix} = ta;
  tss{cix} = ts;

  taas{cix} = taa;
  tsas{cix} = tsa;

end;






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
  corr2(ta.data(ondix),ts.data(ondix))
  % Anomaly
  tsa.date=ts.date; tsa.data=ts.data-scl(get_jday_no_leap(ts.date))';
  taa.date=ta.date; taa.data=ta.data-acl(get_jday_no_leap(ta.date))';
  corr2(taa.data,tsa.data)
  %scatter_fit_ts_seasons(taa,tsa,[],[],'Ta','Ts');










    asd = repmat(0,size(clim));
    for tix=1:numel(tid)
      ix = find(cumfun(anomts.date)==tid(tix));
      anomts.data(ix) = anomts.data(ix) - clim(tix);
    end;





function contour_field(fldstr,fld,ix,varargin)
%function contour_field(fldstr[,fld,ix,varargin])
%
% Use CONTOURF to plot field FLD with coordinates FLDSTR.lon and FLDSTR.lat.
% DEFAULT FLD is field FLDSTR.field. If IX is a non-empty index list, plot
% SQUEEZE(FLD(ix,:,:)), if a FUNCTION_HANDLE (v.), plot SQUEEZE(IX(FLD)).
% Uses DASPECT to set aspect ratio based on latitude. All other arguments are
% passed through to CONTOURF as args 4-N.
%
% Last Saved Time-stamp: <Mon 2016-05-02 12:57:01 Eastern Daylight Time gramer>

  if ( ~exist('fld','var') || isempty(fld) )
    fld = fldstr.field;
  elseif ( ischar(fld) )
    fld = fldstr.(fld);
  end;
  if ( exist('ix','var') && ~isempty(ix) )
    if ( isa(ix,'function_handle') )
      fld = ix(fld);
    else
      fld = fld(ix,:,:);
    end;
  end;

  fld = squeeze(fld);

  contourf(fldstr.lon,fldstr.lat,fld,varargin{:});
  daspect([1,cosd(fldstr.lat(1)),1]);
  if ( numel(varargin) > 0 && numel(varargin{1}) > 2 )
    set(gca,'CLim',[varargin{1}(1),varargin{1}(end)]);
  end;
  colorbar;

return;







  %[lat,lon] = utm2deg(949622.38,589331.52,'17 R')
  %[lat,lon] = minvtran(utmstruct,949622.38,589331.52)
  %[lat,lon] = minvtran(utmstruct,x,y);




    % But if user gave us vectors, pass back vectors in return
    lat = linspace(min(LAT(:)),max(LAT(:)),numel(unique(y(:))));
    lon = linspace(min(LON(:)),max(LON(:)),numel(unique(x(:))));







  %[x,y] = mfwdtran(utmstruct,latlim,lonlim)

  if ( ~isvector(x) && ~isvector(y) )
    %[lat,lon] = utm2deg(949622.38,589331.52,'17 R')
    %[lat,lon] = minvtran(utmstruct,949622.38,589331.52)
    [lat,lon] = minvtran(utmstruct,x,y);

  else
    % MINVTRAN relies on MESHGRID format
    [X,Y] = meshgrid(x,y);
    [LAT,LON] = minvtran(utmstruct,X,Y);

    % But if user gave us vectors, pass back vectors in return
    lat = linspace(min(LAT(:)),max(LAT(:)),numel(unique(y(:))));
    lon = linspace(min(LAT(:)),max(LAT(:)),numel(unique(y(:))));
    % But if user gave us vectors, pass back vectors in return
    latix = linspace(min(LAT(:)),max(LAT(:)),numel(unique(y(:))));
    lonix = linspace(min(LON(:)),max(LON(:)),numel(unique(x(:))));
    lat = interp1(unique(LAT),latix);
    lon = interp1(unique(LON),lonix);
keyboard;
  end;









% if AXLMS==[] (DEFAULT), set limits to include all points on all four plots.



  if ( ~exist('axLms','var') || isempty(axLms) )
    axLms = [];
  end;




  if ( isvector(axLms) )
    if ( numel(axLms) == 1 && isfinite(mn) && isfinite(mx) )
      if ( ischar(axLms) && strncmpi(axLms,'eq',2) )
        subplots_set(fh,'XLim',[xmn,xmx],'YLim',[ymn,ymx]);
      else
        subplots_set(fh,'XLim',[mn,mx],'YLim',[mn,mx]);
      end;
    elseif ( numel(axLms) == 2 )
      subplots_set(fh,'XLim',axLms,'YLim',axLms);
    elseif ( numel(axLms) == 4 )
      subplots_set(fh,'XLim',axLms(1:2),'YLim',axLms(3:4));
    else
      warning('Ignoring arg AXLMS...')
    end;
  end;






function catchwarn(msg)
%function catchwarn(msg)
%
% Display the last error which occured (see LASTERROR) as a warning, with
% complete stack trace. Display of this "warning" may also be disabled by
% calling WARNING('OFF',IDENT) with the identifier of the error. Normally
% called, e.g., in the CATCH of a TRY,CATCH (v.)
%
% Last Saved Time-stamp: <Fri 2016-03-18 14:25:41 Eastern Daylight Time gramer>

  if ( ~exist('msg','var') )
    msg = '';
  end;
  warnstr = ['CAUGHT ERROR: ',msg];
  e = lasterror;
  warnstr = sprintf('%s\n%s\n%s\n',warnstr,e.identifier,e.message);
  for ix=1:numel(e.stack)
    %warnstr = sprintf('%s\n%s:%g\n',warnstr,strrep(e.stack(ix).file,'\','\\'),e.stack(ix).line);
    warnstr = sprintf('%s\n%s:%g\n',warnstr,e.stack(ix).file,e.stack(ix).line);
  end;

  if ( isempty(e.identifier) )
    warning(warnstr);
  else
    s = warning('query',e.identifier);
    if ( ~strcmpi(s.state,'off') )
      warning(warnstr);
    end;
  end;

return;









      if ( STATUS ~= 0 )
        warning('Failed download! %s',url);
      end;





% Simple QA - assumes *cross-shore* component must generally be quite small
%badix = find(nanmax(abs(u_avg.prof),[],2)>0.5);
badix = find( ~isfinite(sp_sfc.data) | (nanmax(abs(x.prof),[],2)>1.0) );
if ( ~isempty(badix) )
  badtix = unique(find(ismember(sfc_t.date,sp_sfc.date(badix))));

  disp(['Clearing ',num2str(numel(badix)),' bad ADCP and ',num2str(numel(badtix)),' bad Temp. records']);
keyboard;
  sfc_t.date(badtix) = [];
  sfc_t.data(badtix) = [];

  sfc_bin(badix) = [];
  sfc_hgt(badix) = [];
  last_good_bin(badix) = [];

  sp_sfc.date(badix) = [];
  sp_sfc.data(badix) = [];

  dr_sfc.date(badix) = [];
  dr_sfc.data(badix) = [];

  x_sfc.date(badix) = [];
  x_sfc.data(badix) = [];
end;








function y = fixn(x,n) 
% 
% See also ROUNDN, FLOORN, CEILN, ROUND, FLOOR, CEIL, FIX.

  error(nargchk(2,2,nargin)) 
  error(nargoutchk(0,1,nargout)) 

  y = fix(x./n).*n; 

return;


function y = floorn(x,n) 
% 
% See also ROUNDN, ROUND, FLOOR, CEIL, FIX.

  error(nargchk(2,2,nargin)) 
  error(nargoutchk(0,1,nargout)) 

  y = floor(x./n).*n; 

return;


function y = ceiln(x,n) 
% 
% See also ROUNDN, ROUND, FLOOR, CEIL, FIX.

  error(nargchk(2,2,nargin)) 
  error(nargoutchk(0,1,nargout)) 

  y = ceil(x./n).*n; 

return;



  elseif ( numel(longform) == numel(seas) )
    strs = repmat(' ',size(season_names));
    ix = find(seas == 1);
    longix = find(longform(ix) == true); shortix = find(longform(ix) ~= true);
    strs(ix(longix),:) = long_season_names(seas,:);

    strs(longform,:) = long_season_names(seas(longform),:);
    strs(~longform,:) = season_names(seas(~longform),:);
  else
    error('LONGFORM must be scalar or have NUMEL == NUMEL(SEAS)');
  end;







function strs = get_season_name(seas,longform)
%function strs = get_season_name(seas,longform)

  if ( ~exist('longform','var') || isempty(longform) )
    longform = false;
  end;
  if ( numel(longform) ~= numel(seas) )
    longform = repmat(longform,size(seas));
  end;

  season_names      = ['JFM    ';'AMJ    ';'JAS    ';'OND    '];
  long_season_names = ['Jan-Mar';'Apr-Jun';'Jul-Sep';'Oct-Dec'];

  strs = repmat(' ',size(season_names));

  ix = find(seas == 1);
  longix = find(longform(ix) == true); shortix = find(longform(ix) ~= true);
  strs(ix(longix),:) = long_season_names(seas,:);

return;





  strsz = size(season_names,2);





    if ( samplehrs ~= floor(samplehrs) || nhrs ~= floor(nhrs) )

    % If reading frame doesn't evenly divide sample period, aliasing likely
    elseif ( gcd(samplehrs, nhrs) ~= samplehrs )




From scatter_fit_ts_seasons.m:
  % Note: SCATTER_FIT_TS calls TITLENAME (v.) - we must undo that
  set(fh,'Name',['Seasons: ',xlbl,' vs. ',ylbl]);
  % Return focus to upper-right panel in case caller wishes to assign a TITLE
  if ( ishandle(ax(1)) )
    axes(ax(1));
  elseif ( ishandle(ax(2)) )
    axes(ax(2));
  end;




function th = subplot_titlename(m,n,varargin)
%function th = subplot_titlename(m,n,varargin)
%
% Set Axes title, but also set Axes' Figure 'Name' property. Args are
% identical to TITLENAME (qv.). This version is for SUBPLOTS (v.): picks the
% middle subplot of the top row to apply the "TITLE" to. M and N are the args
% that were passed to SUBPLOT (or SUBPLOT_TIGHT).
% 
% Last Saved Time-stamp: <Fri 2015-08-14 16:13:47 Eastern Daylight Time gramer>

  
  th = titlename(varargin);

return;





function [bathfiles,stn_covered] = read_hires_bathymetry_check_file(stn,bbox,filebbox,file,bathfiles,stn_covered)
  if ( bboxint(bbox,filebbox) > eps )
    bathfiles{end+1} = file;
  end;
  if ( ~isempty(bboxinside(stn.lon,stn.lat,filebbox,false,0)) )
    stn_covered = true;
  end;
return;








  coastpath = get_ecoforecasts_path('coast');

  [lats,lons] = reckon_wgs84(stn.lat,stn.lon,rad([2,1,2,1])./1e3,[0,90,180,270]);
  bbox = [min(lons),max(lons),min(lats),max(lats)];

  if ( ~exist('bathfiles','var') || isempty(bathfiles) )
    bathfiles = {};

    % Do any of the files we've accumulated so far CONTAIN OUR LOCATION?
    stn_covered = false;

    % Use the finest available bathymetry - but don't mix resolutions
    if ( useHighestRes )

      % [bathfiles,stn_covered] = ...
      %     read_hires_bathymetry_check_file(stn,bbox,filebbox,bathfiles,stn_covered);

      % 1/3-arcsecond products
      if ( bboxint(bbox,[-88.30,-87.65,30.00,31.00]) > eps )
        bathfiles{end+1} = fullfile(coastpath,'mobile_al_mhw.grd');
      end;
      if ( bboxint(bbox,[-91.60,-88.80,28.60,29.70]) > eps )
        bathfiles{end+1} = fullfile(coastpath,'southern_louisiana_mhw.grd');
      end;
      % Palm Beach too big to fit in memory! Oh well...
      % if ( bboxint(bbox,[-?.30,-?.65,?.00,?.00]) > eps )
      %   bathfiles{end+1} = 'pb_mhw';
      % end;
      if ( bboxint(bbox,[-89.30,-88.30,29.70,30.60]) > eps )
        bathfiles{end+1} = 'biloxi_ms';
      end;
      if ( bboxint(bbox,[-86.10,-85.20,29.55,30.50]) > eps )
        bathfiles{end+1} = fullfile(coastpath,'panama_city_fl_1-3_arc-second_mhw_netcdf.grd');
      end;
      if ( bboxint(bbox,[-82.18,-81.27,23.98,25.02]) > eps )
        bathfiles{end+1} = 'key_west_fl_mhw';
      end;
      if ( bboxint(bbox,[-65.04,-64.40,17.61,17.88]) > eps )
        bathfiles{end+1} = 'st_croix_mhw_update';
      end;
      if ( bboxint(bbox,[-65.15,-64.64,18.17,18.48]) > eps )
        bathfiles{end+1} = 'sts_thoms_john_mhw_update';
      end;
      if ( bboxint(bbox,[-65.90,-65.35,18.05,18.60]) > eps )
        bathfiles{end+1} = 'fajardo';
      end;
      if ( bboxint(bbox,[-67.60,-67.10,17.90,18.60]) > eps )
        bathfiles{end+1} = 'mayaguez';
      end;
      if ( bboxint(bbox,[-67.10,-66.40,17.70,18.05]) > eps )
        bathfiles{end+1} = 'ponce';
      end;
      if ( bboxint(bbox,[144.30,145.20,13.00,13.90]) > eps )
        bathfiles{end+1} = 'guam_1_3s_20081014';
      end;

      % 1-arcsecond products
      if ( isempty(bathfiles) )
        if ( bboxint(bbox,[-65.15,-64.00,17.00,19.00]) > eps )
          %bathfiles{end+1} = 'us_virgin_islands_vi_1s_mhw';
          bathfiles{end+1} = 'us_virgin_islands_vi_1s_mhw_update';
        end;
        if ( bboxint(bbox,[-68.00,-65.00,17.00,19.00]) > eps )
          bathfiles{end+1} = 'pr_1s';
        end;
        if ( bboxint(bbox,[-90.75,-85.00,24.00,38.00]) > eps )
          bathfiles{end+1} = fullfile(coastpath,'northern_gulf_coast_mhw.grd');
        end;
        % How do I read that stupid south_florida 30 m file from NGDC?? 
      end;

    end;

    % 3-arcsecond products (6-arcsecond for Mariana region)
    if ( isempty(bathfiles) ) 
      if ( bboxint(bbox,[-87,-78,24,35]) > eps )
        bathfiles{end+1} = fullfile(coastpath,'fl_east_gom_crm_v1.nc');
      end;
      if ( bboxint(bbox,[-108,-94,24,38]) > eps )
        bathfiles{end+1} = fullfile(coastpath,'western_gom_crm_v1.nc');
      end;
      if ( bboxint(bbox,[-94,-87,24,36]) > eps )
        bathfiles{end+1} = fullfile(coastpath,'central_gom_crm_v1.nc');
      end;
      if ( bboxint(bbox,[-68,-64,16,20]) > eps )
        bathfiles{end+1} = fullfile(coastpath,'puerto_rico_crm_v1.nc');
      end;
      if ( bboxint(bbox,[-162,-152,18,24]) > eps )
        bathfiles{end+1} = fullfile(coastpath,'hawaii_crm_v1.nc');
      end;
      if ( bboxint(bbox,[-171.14,-169.30,-14.73,-13.83]) > eps )
        bathfiles{end+1} = 'pago_pago_3s';
      end;
      if ( bboxint(bbox,[139.00,150.00,10.00,22.00]) > eps )
        bathfiles{end+1} = 'MarianaTrench_6as_v1';
      end;
      if ( bboxint(bbox,[-85,-68,31,40]) > eps )
        bathfiles{end+1} = fullfile(coastpath,'se_atl_crm_v1.nc');
      end;
      if ( bboxint(bbox,[-128,-117,37,44]) > eps )
        bathfiles{end+1} = fullfile(coastpath,'central_pacific_crm_v1.nc');
      end;
      if ( bboxint(bbox,[-80,-64,40,48]) > eps )
        bathfiles{end+1} = fullfile(coastpath,'ne_atl_crm_v1.nc');
      end;
      if ( bboxint(bbox,[-128,-116,44,49]) > eps )
        bathfiles{end+1} = fullfile(coastpath,'nw_pacific_crm_v1.nc.nc');
      end;
    end;
  end;

  % Don't force caller to wrap a single filename string in '{}'
  if ( ischar(bathfiles) )
    bathfiles = {bathfiles};
  end;








  %[lats,lons] = reckon_wgs84(stn.lat,stn.lon,[rad,rad]./1e3,[0,180,90,270]);
  [lats,lons] = reckon_wgs84(stn.lat,stn.lon,rad./1e3,[0,180,90,270]);
  bbox = [min(lons),max(lons),min(lats),max(lats)];



          'TASP7', ...		% tasp7,14.122,145.159: Talakhaya/Sabana Watershed - Rota - CNMI



function [rng,az] = distance_wgs84(lat1,lon1,lat2,lon2)
%function [rng,az] = distance_wgs84(lat1,lon1,lat2,lon2)
%
% Calls DISTANCE (v.) with WGS-84 ellipticity and polar radius to calculate
% distance(s) [km] between pairs of lat/lon points. Args may be matrices.
% If LAT1 (LON1) scalar but LAT2 (LON2) is not, call REPMAT to match sizes.
% Otherwise, requires size(lat1) == size(lon1) == size(lat2) == size(lon2).
%
% Last Saved Time-stamp: <Fri 2015-07-17 16:31:31 Eastern Daylight Time gramer>

  if ( isscalar(lat1) );  lat1 = repmat(lat1,size(lat2));  end;
  if ( isscalar(lon1) );  lon1 = repmat(lon1,size(lon2));  end;

  if ( ndims(lat1)~=ndims(lon1) || any(size(lat1)~=size(lon1)) || ...
       ndims(lat1)~=ndims(lat2) || any(size(lat1)~=size(lat2)) || ...
       ndims(lon1)~=ndims(lon2) || any(size(lon1)~=size(lon2)) )
    error('LAT1,LON1 must be scalar or must match dimensions of LAT2,LON2!');
  end;

  [rng,az] = distance(lat1,lon1, lat2,lon2, wgs84_ellipsoid());

return;






      lonix = interp1(lon,1:numel(lon),x(bix).lon,'nearest','extrap');
      lonix(~isfinite(lonix)) = [];
      latix = interp1(lat,1:numel(lat),x(bix).lat,'nearest','extrap');
      latix(~isfinite(latix)) = [];







    [pth,fnm,fex]=fileparts(bathfile);
    doAsc = ( strcmpi(fex,'.asc') || strcmpi(fex,'.mat') );
    if ( doAsc )
      [all_lons,all_lats,all_zs] = asc2mat(bathfile);
    else
      nc = mDataset(bathfile);
      if ( isempty(nc) )
        error('MDATASET unable to open %s',bathfile);
      end;
      all_lats = cast(nc{'y'}(:),'double');
      all_lons = cast(nc{'x'}(:),'double');
    end;





    if ( ~isempty(latixen) && ~isempty(lonixen) )
      if ( doAsc )
        z = all_zs(latixen,lonixen)';
      else
        z = cast(nc{'z'}(latixen,lonixen),'double')';
        close(nc);
      end;
    end;

    x(bix).doAsc = doAsc;












  [lons,lats] = reckon_wgs84(stn.lat,stn.lon,[rad,rad]./1e3,[0,180,90,270]);
  bbox = [min(lons),max(lons),min(lats),max(lats)];

  if ( ~exist('bathfiles','var') || isempty(bathfiles) )
    bathfiles = {};

    % Use the finest available bathymetry - but don't mix resolutions
    if ( useHighestRes )

      % 1/3-arcsecond products
      if ( (30.00<=bbox(4) || bbox(3)<=31.00) && (-88.30<=bbox(2) || bbox(1)<=-87.65) )
        bathfiles{end+1} = fullfile(coastpath,'mobile_al_mhw.grd');
      end;
      % Palm Beach too big to fit in memory! Oh well...
      % if ( (?.??<=bbox(4) || bbox(3)<=?.??) && (-??.??<=bbox(2) || bbox(1)<=-?.??) )
      %   bathfiles{end+1} = 'pb_mhw';
      % end;
      if ( (23.98<=bbox(4) || bbox(3)<=25.02) && (-82.18<=bbox(2) || bbox(1)<=-81.27) )
        bathfiles{end+1} = 'key_west_fl_mhw';
      end;
      if ( (17.61<=bbox(4) || bbox(3)<=17.88) && (-65.04<=bbox(2) || bbox(1)<=-64.40) )
        bathfiles{end+1} = 'st_croix_mhw_update';
      end;
      if ( (18.17<=bbox(4) || bbox(3)<=18.48) && (-65.15<=bbox(2) || bbox(1)<=-64.64) )
        bathfiles{end+1} = 'sts_thoms_john_mhw_update';
      end;
      if ( (17.70<=bbox(4) || bbox(3)<=18.05) && (-67.10<=bbox(2) || bbox(1)<=-66.40) )
        bathfiles{end+1} = 'ponce';
      end;
      if ( (18.05<=bbox(4) || bbox(3)<=18.60) && (-65.90<=bbox(2) || bbox(1)<=-65.35) )
        bathfiles{end+1} = 'fajardo';
      end;
      if ( (17.90<=bbox(4) || bbox(3)<=18.60) && (-67.60<=bbox(2) || bbox(1)<=-67.10) )
        bathfiles{end+1} = 'mayaguez';
      end;

      % 1-arcsecond products
      if ( isempty(bathfiles) )
        if ( (17.00<=bbox(4) || bbox(3)<=19.00) && (-65.15<=bbox(2) || bbox(1)<=-64.00) )
          %bathfiles{end+1} = 'us_virgin_islands_vi_1s_mhw';
          bathfiles{end+1} = 'us_virgin_islands_vi_1s_mhw_update';
        end;
        if ( (17.00<=bbox(4) || bbox(3)<=19.00) && (-68.00<=bbox(2) || bbox(1)<=-65.00) )
          bathfiles{end+1} = 'pr_1s';
        end;
        if ( (24.00<=bbox(4) || bbox(3)<=38.00) && (-90.75<=bbox(2) || bbox(1)<=-85.0) )
          bathfiles{end+1} = fullfile(coastpath,'northern_gulf_coast_mhw.grd');
        end;
        % How do I read that stupid south_florida 30 m file from NGDC?? 
      end;

    end;

    % 3-arcsecond products
    if ( isempty(bathfiles) ) 
      if ( (24<=bbox(4) || bbox(3)<=35) && (-87<=bbox(2) || bbox(1) <=-78) )
        bathfiles{end+1} = fullfile(coastpath,'fl_east_gom_crm_v1.nc');
      end;
      if ( (24<=bbox(4) || bbox(3)<=38) && (-108<=bbox(2) || bbox(1)<=-94) )
        bathfiles{end+1} = fullfile(coastpath,'western_gom_crm_v1.nc');
      end;
      if ( (24<=bbox(4) || bbox(3)<=36) && (-94<=bbox(2) || bbox(1)<=-87) )
        bathfiles{end+1} = fullfile(coastpath,'central_gom_crm_v1.nc');
      end;
      if ( (16<=bbox(4) || bbox(3)<=20) && (-68<=bbox(2) || bbox(1)<=-64) )
        bathfiles{end+1} = fullfile(coastpath,'puerto_rico_crm_v1.nc');
      end;
      if ( (18<=bbox(4) || bbox(3)<=24) && (-162<=bbox(2) || bbox(1)<=-152) )
        bathfiles{end+1} = fullfile(coastpath,'hawaii_crm_v1.nc');
      end;
      if ( (31<=bbox(4) || bbox(3)<=40) && (-85<=bbox(2) || bbox(1)<=-68) )
        bathfiles{end+1} = fullfile(coastpath,'se_atl_crm_v1.nc');
      end;
      if ( (37<=bbox(4) || bbox(3)<=44) && (-128<=bbox(2) || bbox(1)<=-117) )
        bathfiles{end+1} = fullfile(coastpath,'central_pacific_crm_v1.nc');
      end;
      if ( (40<=bbox(4) || bbox(3)<=48) && (-80<=bbox(2) || bbox(1)<=-64) )
        bathfiles{end+1} = fullfile(coastpath,'ne_atl_crm_v1.nc');
      end;
      if ( (44<=bbox(4) || bbox(3)<=49) && (-128<=bbox(2) || bbox(1)<=-116) )
        bathfiles{end+1} = fullfile(coastpath,'nw_pacific_crm_v1.nc.nc');
      end;
    end;
  end;







  if ( numel(bathfiles) > 1 )
error('NEED TO PREALLOCATE STRUCT WITH PROPER-SIZED MATRICES!');
  end;






FROM read_hires_bathymetry.m:
  if ( ~exist('bathfiles','var') || isempty(bathfiles) )
    bathfiles = {};

    % Use the finest available bathymetry - but don't mix resolutions
    if ( useHighestRes && ...
         24<=stn.lat-(rad(2)/111e3) && stn.lat+(rad(2)/111e3)<=38 ...
         && -90.75<=stn.lon-(rad(1)/111e3) && stn.lon+(rad(1)/111e3)<=--85.0 )
      bathfiles{end+1} = fullfile(coastpath,'northern_gulf_coast_mhw.grd');
    else
if ( 24<=stn.lat && stn.lat<=35 && -87<=stn.lon && stn.lon<=-78 )
      bathfiles = fullfile(coastpath,'fl_east_gom_crm_v1.nc');
    elseif ( 24<=stn.lat && stn.lat<=38 && -108<=stn.lon && stn.lon<=-94 )
      bathfiles = fullfile(coastpath,'western_gom_crm_v1.nc');
    elseif ( 24<=stn.lat && stn.lat<=36 &&  -94<=stn.lon && stn.lon<=-87 )
      bathfiles = fullfile(coastpath,'central_gom_crm_v1.nc');
    elseif ( 16<=stn.lat && stn.lat<=20 &&  -68<=stn.lon && stn.lon<=-64 )
      bathfiles = fullfile(coastpath,'puerto_rico_crm_v1.nc');
    elseif ( 18<=stn.lat && stn.lat<=24 && -162<=stn.lon && stn.lon<=-152 )
      bathfiles = fullfile(coastpath,'hawaii_crm_v1.nc');
    else
      error('Please specify an NGDC CRM netCDF file including those coordinates??');
    end;
    disp(['Using ',bathfiles]);
  end;

  [lon,lat,dat]=asc2mat(bathfile);

  nc = mDataset(bathfile);








    xrad = floor(ncols/2);
    lon = xllc + ([-xrad:xrad].*cellsz);

    yrad = floor(nrows/2);
    lat = yllc + ([-yrad:yrad].*cellsz);









function [lon,lat,dat,matfname] = asc2mat(ascfname,matfname)
%function [lon,lat,dat,matfname] = asc2mat(ascfname,matfname)

  tic,

  matfname = strrep(ascfname,'.asc','.mat');








  cursegix = 1;
  while ( cursegix < size(cs,2) )
    if ( isempty(plotargs) )
      if ( isempty(zlevel) )
        plot(ax,cs(1,cursegix+1:cursegix+seglen),cs(2,cursegix+1:cursegix+seglen));
      elseif ( isscalar(zlevel) )
        if ( islogical(zlevel) && zlevel )
          zs=repmat(cs(1,cursegix),[1,seglen]);
          plot3(ax,cs(1,cursegix+1:cursegix+seglen),cs(2,cursegix+1:cursegix+seglen),zs);
          zs=[]; clear zs
        else
          zs=repmat(zlevel,[1,seglen]);
          plot3(ax,cs(1,cursegix+1:cursegix+seglen),cs(2,cursegix+1:cursegix+seglen),zs);
          zs=[]; clear zs
        end;
      end;
    else
      if ( isempty(zlevel) )
        plot(ax,cs(1,cursegix+1:cursegix+seglen),cs(2,cursegix+1:cursegix+seglen),plotargs{:});
      elseif ( islogical(zlevel) && isscalar(zlevel) && zlevel )
        zs=repmat(cs(1,cursegix),[1,seglen]);
        plot3(ax,cs(1,cursegix+1:cursegix+seglen),cs(2,cursegix+1:cursegix+seglen),zs,plotargs{:});
        zs=[]; clear zs
      else
        zs=repmat(zlevel,[1,seglen]);
        plot3(ax,cs(1,cursegix+1:cursegix+seglen),cs(2,cursegix+1:cursegix+seglen),zs,plotargs{:});
        zs=[]; clear zs
      end;
    end;

    cursegix = cursegix+seglen+1;
  end;





function [cs,h] = plot_contour_cs(varargin)
%function [cs,h] = plot_contour_cs([curaxes,zlevel,]cs[,plotargs])
%
% PLOT lines for all contours in contour matrix CS (v. CONTOURC) on the
% CURAXES (DEFAULT: GCA). Optional PLOTARGS is passed onto PLOT.
% (Hard to believe MATLAB doesn't have a way to do this, but... ?)
% If ZLEVEL is specified, use PLOT3; if ZLEVEL is the logical true
% scalar (LOGICAL(1)), use PLOT3 to plot the value levels from CS.
% 
% Last Saved Time-stamp: <Mon 2015-05-25 13:48:39 Eastern Daylight Time gramer>

  cs=[];
  h=[];

  args = varargin;
  if ( ishandle(args{1}) )
    ax = args{1};
    args(1) = [];
  elseif ( isempty(get(0,'Children')) )
    ax = gca;
    hold on;
  else
    ax = gca;
  end;

  zlevel = [];
  if ( isscalar(args{1}) )
    zlevel = args{1};
    args(1) = [];
  end;

  cs = args{1};
  args(1) = [];

  plotargs = args;

  if ( isempty(zlevel) )
    if ( isempty(plotargs) )
      plotfun = @(x,y)(plot(ax,x,y));
    else
      plotfun = @(x,y)(plot(ax,x,y,plotargs{:}));
    end;
  elseif ( isscalar(zlevel) )
    if ( isempty(plotargs) )
      plotfun = @(x,y,z)(plot3(ax,x,y,z));
    else
      plotfun = @(x,y,z)(plot3(ax,x,y,z,plotargs{:}));
    end;
  else
    error('If ZLEVEL is specified, it must be a double or logical scalar');
  end;

  cursegix = 1;
  while ( cursegix < size(cs,2) )
    seglen = floor(cs(2,cursegix));
    if ( isempty(plotargs) )
      if ( isempty(zlevel) )
        plot(ax,cs(1,cursegix+1:cursegix+seglen),cs(2,cursegix+1:cursegix+seglen));
      elseif ( isscalar(zlevel) )
        if ( islogical(zlevel) && zlevel )
          zs=repmat(cs(1,cursegix),[1,seglen]);
          plot3(ax,cs(1,cursegix+1:cursegix+seglen),cs(2,cursegix+1:cursegix+seglen),zs);
          zs=[]; clear zs
        else
          zs=repmat(zlevel,[1,seglen]);
          plot3(ax,cs(1,cursegix+1:cursegix+seglen),cs(2,cursegix+1:cursegix+seglen),zs);
          zs=[]; clear zs
        end;
      end;
    else
      if ( isempty(zlevel) )
        plot(ax,cs(1,cursegix+1:cursegix+seglen),cs(2,cursegix+1:cursegix+seglen),plotargs{:});
      elseif ( islogical(zlevel) && isscalar(zlevel) && zlevel )
        zs=repmat(cs(1,cursegix),[1,seglen]);
        plot3(ax,cs(1,cursegix+1:cursegix+seglen),cs(2,cursegix+1:cursegix+seglen),zs,plotargs{:});
        zs=[]; clear zs
      else
        zs=repmat(zlevel,[1,seglen]);
        plot3(ax,cs(1,cursegix+1:cursegix+seglen),cs(2,cursegix+1:cursegix+seglen),zs,plotargs{:});
        zs=[]; clear zs
      end;
    end;
    cursegix = cursegix+seglen+1;
  end;

return;










    elseif ( islogical(doLabels{2}) )
      inlineLabels = doLabels{2};
      doLabels = doLabels{1};







        if ( ~isempty(interpMinFin) )
          goodix = find(sum(isfinite(dat(:,:)')',2)>=interpMinFin);
        else
          goodix = 1:size(vals,1);
        end;
        vals(goodix,ix) = interpMethod(dat(:,:)')';
        % vals(:,ix) = interpMethod(dat(:,:)')';







  % Process all other data columns
  nGoodCols = 0;
  for col = begcol:endcol
    % Only process non-blank columns
    if ( ~isempty(hdr{col}) )
      nGoodCols = nGoodCols + 1;
      result.(hdr{col}).date = alltimestamps;
      % Return a vector of numbers where ever possible
      nnums = numel(find(cellfun(@ischar, rawcells(begrow:endrow,col))));
      if ( nnums == length(result.(hdr{col}).date) )
        result.(hdr{col}).data = [rawcells{begrow:endrow, col}]';
      else
        % Try a cellstr-to-number conversion if necessary
        result.(hdr{col}).data = str2double(rawcells(begrow:endrow, col));
        if ( numel(find(isnan(result.(hdr{col}).data))) > (0.5*length(result.(hdr{col}).data)) )
          % If most of the values in the column really are non-numeric,
          % return a cell array of strings instead.
          result.(hdr{col}).data = rawcells(begrow:endrow, col);
        end;
      end;
    end;
  end;







  for col = begcol:endcol
    % Only process non-blank columns
    if ( ~isempty(hdr{col}) )
      nGoodCols = nGoodCols + 1;
      result.(hdr{col}).date = alltimestamps;
      result.(hdr{col}).data = repmat(nan,size(result.(hdr{col}).date));
      % Return a vector of numbers where ever possible
      if ( ~iscellstr(rawcells{begrow:endrow, col}) )
      try,
        warning('OFF','MATLAB:nonIntegerTruncatedInConversionToChar');
        result.(hdr{col}).data = [rawcells{begrow:endrow, col}]';
        warning('ON','MATLAB:nonIntegerTruncatedInConversionToChar');
      catch,
      end;
      if ( numel(find(isnan(result.(hdr{col}).data))) > (0.5*length(result.(hdr{col}).data)) )
        % Try a cellstr-to-number conversion if necessary
        result.(hdr{col}).data = str2double(rawcells(begrow:endrow, col));
        if ( numel(find(isnan(result.(hdr{col}).data))) > (0.5*length(result.(hdr{col}).data)) )
          % If most of the values in the column really are non-numeric,
          % return a cell array of strings instead.
          result.(hdr{col}).data = rawcells(begrow:endrow, col);
        end;
      end;
    end;
  end;







From INTERSECT_DATES.m:
  % Make sure numel(ix1) ~= numel(ix2) after unique'ing
  [ig,newix] = unique(ix1);
  ix1 = ix1(newix);
  ix2 = ix2(newix);
  [ig,newix] = unique(ix2);
  ix1 = ix1(newix);
  ix2 = ix2(newix);

  % Assume dates came in sorted - ix2 was just sorted by UNIQUE 
  ix1 = sort(ix1);





function stn = read_sbe(fnames,stn)
%function stn = read_sbe(fnames,stn)
%
% Read one or more .ASC files, or a directory of .ASC files, containing SBE37
% or SBE39 ASC-formatted data. This supplants the seemingly FUBAR functions
% named SBE* (SBE37PARSE, etc.) contained in the CSIRO's IMOS-Toolbox.
%
% Last Saved Time-stamp: <Fri 2015-04-24 20:19:34 Eastern Daylight Time gramer>

  % Did user forget the first arg?
  if ( ~exist('fnames','var') )
    fnames = [];
  elseif ( ischar(fnames) )
    % Is it a single filename?
    if ( exist(fnames,'file') )
      fnames = {fnames};
    % Is it a directory?
    elseif ( exist(fnames,'dir') )
      fs = dir(fullfile(fnames,'*.asc'));
      fnames = {fs.name};
    % Oops
    else
      fnames = [];
    end;
  % Is it a cell array of filename strings?
  elseif ( ~iscellstr(fnames) )
    fnames = [];
  end;

  % No matter what FNAMES was passed in as, it is now either empty, ,or a
  % cell array of (one or more) pathnames for data files.
  if ( isempty(fnames) )
      error('First arg must be CELLSTR, or  individual directory or filename string');
  end;

  if ( ~exist('stn','var') )
    stn = [];
  end;

tic,
  fid = fopen('NSUOC\SFOMC\data\c-buoy\mc\c0842.asc','r');
  ccstr = textscan(fid,'%[^\n]\n');  fclose(fid);

  cstr=ccstr{:}; clear ccstr;

  begix = strmatch('start sample number',cstr); begix = begix + 1;
  while ( isempty(cstr{begix}) ); begix = begix + 1; end;

  ncoms = numel(strfind(cstr{begix},','));

  str = sprintf('%s\n',cstr{begix:end});

  % T,P,C,date,time
  if ( ncoms == 4 )
    dat = textscan(str,'%f,%f,%f %[^ ] %f,%f:%f:%f\n');
    tmp = dat{1};
    prs = dat{2};
    cnd = dat{3};
    dts = datenum(dat{6},month_str2num(dat{5}),dat{4},dat{7},dat{8},dat{9});
  % T,P,date,time
  elseif ( ncoms == 3 )
    dat = textscan(str,'%f,%f,%f %[^ ] %f,%f:%f:%f\n');
    tmp = dat{1};
    prs = dat{2};
    dts = datenum(dat{5},month_str2num(dat{4}),dat{3},dat{6},dat{7},dat{8});
  % T,date,time
  elseif ( ncoms == 2 )
    dat = textscan(str,'%f,%f,%f %[^ ] %f,%f:%f:%f\n');
    tmp = dat{1};
    dts = datenum(dat{4},month_str2num(dat{3}),dat{2},dat{5},dat{6},dat{7});
  else
    error('Unknown format! First data line: "%s"',cstr{begix});
  end;
toc,
return;







  str = strcat(char(cstr{begix:end-1}),sprintf('\n'));





function yrs = get_yearfrac(dts)
%function yrs = get_yearfrac(dts)
%
% Return "year-fractions" for each date in the vector of DATENUM (qv.) DTS.
% E.g., the year-fraction for DATENUM(2003,1,1,0,0,0) is 2003.0000, while
% that for DATENUM(2004,06,01,23,59,59) is 2004.0000
%
% Last Saved Time-stamp: <Thu 2014-10-23 14:32:40 Eastern Daylight Time gramer>

  [yrs,ig,ig] = datevec(dts);
  yds = dts - datenum(yrs,1,1);

return;




  if ( ~exist('doNorm', 'var') || isempty(doNorm) )
    %doNorm = false;
    doNorm = true;
  end;


  if ( doNorm )
    normFac = nanstd(vecs,0,2);
    normAdd = nanmean(vecs,2);
    vecs = (vecs - repmat(normAdd,[1,size(vecs,2)])) .* repmat(normFac,[1,size(vecs,2)]);
  end;


  if ( doNorm )
    nds = (sm.codebook .* normFac) + normAdd;
  else
    nds = sm.codebook;
  end;


  cmins = min(nds); cmn = min(reshape(cmins, [npts nvars]));
  cmaxs = max(nds); cmx = max(reshape(cmaxs, [npts nvars]));
  cmin = min([-abs(cmn) ; -abs(cmx)]); cmax = max([abs(cmn) ; abs(cmx)]);
  cmin = cmn; cmax = cmx;
  %%%% ???  cmin = repmat(0, size(cmn)); cmax = cmx;
  for orderix = 1:length(disp_order)
    ix = disp_order(orderix);
    subplot(mapdims(2), mapdims(1), orderix);
    hold on;
    x = [1:npts]/24;
    cdbk = reshape(nds(ix,:), [npts nvars]);










    % Catenate criterion time-series at end of data vector
    evtix = find(cfdts(ix)-bakdys<=cts.date & cts.date<=cfdts(ix)+fwddys);
    [ig,hrix] = unique(get_yearhour(cts.date(evtix)));
    datix = unique(round(interp1(cfdts(ix)+rdts,1:numel(rdts),get_yearhour(cts.date(evtix)))));
    evtix(~isfinite(datix)) = [];
    datix(~isfinite(datix)) = [];
    vecs(ix,npts+datix) = cts.data(evtix);







  disp('Training Self-Organizing Map...');

  %initStr = 'randinit';
  initStr = 'lininit';

  % % % sm = som_make(vecs, 'msize',mapdims, 'lininit', 'batch', 'tracking',0, ...
  % % %               'shape','sheet', 'training','long', 'neigh','ep');
  % % sm = som_make(vecs, 'msize',mapdims, 'randinit', 'seq', 'tracking',0, ...
  % %               'shape','sheet', 'training','long', 'neigh','ep');
  % sm = som_make(vecs, 'msize',mapdims, 'randinit', 'batch', 'tracking',0, ...
  %               'shape','sheet', 'training','long', 'neigh','ep');
  sm = som_make(vecs, 'msize',mapdims, initStr, 'batch', 'tracking',0, ...
                'shape','sheet', 'training','long', 'neigh','ep');








  if ( any(~isfinite(vecs(:))) )
    disp('NaNs force random initialization');
    initStr = 'randinit';
  else
    initStr = 'lininit';
  end;





crit_ts = stn.(cfld);
%crit_fn = str2func(crit);



  evtix = find(cfdts(ix)-bakdys<=cent_ts.date & cent_ts.date<=cfdts(ix)+fwddys);
  % Glean from non-hourly time series
  evtix = unique(round(interp1(cfdts(ix)+rdts,1:numel(cent_ts.date),cent_ts.date)));




    [ig,startix] = cent_ts.date(evtix(1))
    datix = evtix-startix+1;
    vecs(ix,datix) = cent_ts.data(evtix);




%center_ts = stn.ndbc_barom;
center_ts = crit_ts;

% Center each event based on a local extremum of the centering time-series
for ix = 1:numel(cfdts)
  evtix = find(cfdts(ix)-3<=center_ts.date & center_ts.date<cfdts(ix)+3);
  % If we even had centering TS data during this event!
  if ( ~isempty(evtix) )
    % Peak value
    [ig,peakix] = max(center_ts.data(evtix));
    % % Lowest value
    % [ig,peakix] = min(center_ts.data(evtix));
    cfdts(ix) = center_ts.date(evtix(peakix));
  end;
end;




% Center each event based on a local extremum of barometric pressure
for ix = 1:numel(cfdts)
  evtix = find(cfdts(ix)-3<=stn.ndbc_barom.date & stn.ndbc_barom.date<cfdts(ix)+3);
  % If we even had pressure data during this event!
  if ( ~isempty(evtix) )
    % % Peak barometric pressure
    % [ig,peakix] = max(stn.ndbc_barom.data(evtix));
    % Lowest barometric pressure
    [ig,peakix] = min(stn.ndbc_barom.data(evtix));
    cfdts(ix) = stn.ndbc_barom.date(evtix(peakix));
  end;
end;

% % Center each event based on a local extremum of the criterion time-series
% for ix = 1:numel(cfdts)
%   evtix = find(cfdts(ix)-3<=crit_ts.date & crit_ts.date<cfdts(ix)+3);
%   % Peak value
%   [ig,peakix] = max(crit_ts.data(evtix));
%   % % Lowest value
%   % [ig,peakix] = min(crit_ts.data(evtix));
%   cfdts(ix) = crit_ts.date(evtix(peakix));
% end;






  %evtix = find(cfdts(ix)<=stn.ndbc_barom.date & stn.ndbc_barom.date<cfdts(ix)+3);
  %evtix = find(cfdts(ix)-3<=stn.ndbc_barom.date & stn.ndbc_barom.date<cfdts(ix)+1);




% Center each event on the peak barometric pressure
for ix = 1:numel(cfdts)
  %evtix = find(cfdts(ix)<=stn.ndbc_barom.date & stn.ndbc_barom.date<cfdts(ix)+3);
  evtix = find(cfdts(ix)-3<=stn.ndbc_barom.date & stn.ndbc_barom.date<cfdts(ix)+3);
  [ig,peakix] = max(stn.ndbc_barom.data(evtix));
  cfdts(ix) = stn.ndbc_barom.date(evtix(peakix));
end;

% Center each event on the lowest barometric pressure
for ix = 1:numel(cfdts)
  %evtix = find(cfdts(ix)-3<=stn.ndbc_barom.date & stn.ndbc_barom.date<cfdts(ix)+1);
  evtix = find(cfdts(ix)-3<=stn.ndbc_barom.date & stn.ndbc_barom.date<cfdts(ix)+3);
  [ig,minix] = min(stn.ndbc_barom.data(evtix));
  % If we had pressure data during this event!
  if ( ~isempty(minix) )
    cfdts(ix) = stn.ndbc_barom.date(evtix(minix));
  end;
end;




% % % % % Find cold fronts: pressure >= 99th percentile (around 1026 hPa)
% % % % cfdts = unique(floor(stn.ndbc_barom.date(stn.ndbc_barom.data>=prctile(stn.ndbc_barom.data,99))));

% % % % Alternate method: pressure >= 96th percentile for that part of year (warm or cool season)
% % % clear jdlim
% % % jdlim(ts_boreal_cool(stn.ndbc_barom),1) = prctile(stn.ndbc_barom.data(ts_boreal_cool(stn.ndbc_barom)),98);
% % % jdlim(ts_boreal_warm(stn.ndbc_barom),1) = prctile(stn.ndbc_barom.data(ts_boreal_warm(stn.ndbc_barom)),98);
% % % cfdts = unique(floor(stn.ndbc_barom.date(stn.ndbc_barom.data>=jdlim)));

% % % Alternate method: 5-day pressure range > 10 hPa
% % cfdts = unique(floor(stn.ndbc_barom_5_d_range.date(stn.ndbc_barom_5_d_range.data>=10)));

% % Alternate method: 5-day pressure variance > 96th percentile
% cfdts = unique(floor(stn.ndbc_barom_5_d_var.date(stn.ndbc_barom_5_d_var.data>=prctile(stn.ndbc_barom_5_d_var.data,96))));

% Alternate method: 3-day vectorial wind variance > 96th percentile
cfdts = unique(floor(stn.ndbc_wind1_u_3_d_var_0_d_asof_sum_ndbc_wind1_v_3_d_var.date(stn.ndbc_wind1_u_3_d_var_0_d_asof_sum_ndbc_wind1_v_3_d_var.data>=prctile(stn.ndbc_wind1_u_3_d_var_0_d_asof_sum_ndbc_wind1_v_3_d_var.data,96))));










1;

% if ( ~exist('stn','var') )
%     stn = get_station_from_station_name('lonf1'); stn = load_all_ndbc_data(stn);
% end;
% if ( ~isfield(stn,'erai_spechumid') )
%   stn = get_erai_station(stn);
% end;
% stn = station_spddir_to_uv(stn,'ndbc_wind1_speed','ndbc_wind1_dir');
if ( ~exist('stn','var') )
    stn = optimize_station_heat_budget('lonf1','erai','avhrr_weekly','ndbc','tpxo_tide','erai');
    stn = verify_variable(stn,'erai_ndbc_arf_1_d_sum');
    stn = verify_variable(stn,'erai_ndbc_srf_1_d_sum');
    stn = verify_variable(stn,{'ndbc_erai_erai_30a_latent_flux_term_1_d_sum','ndbc_erai_erai_30a_sensible_flux_term_1_d_sum','erai_ndbc_arf_term_1_d_sum','ndbc_erai_erai_30a_net_flux_term_1_d_sum'});
end;
stn = verify_variable(stn,'ndbc_wind1_speed_1_d_median');
stn = verify_variable(stn,'ndbc_wind1_u_1_d_median');
stn = verify_variable(stn,'ndbc_wind1_v_1_d_median');

% if ( ~exist('smkf1','var') )
%     smkf1 = get_station_from_station_name('smkf1'); smkf1 = load_all_ndbc_data(smkf1);
% end;
% smkf1 = station_spddir_to_uv(smkf1,'ndbc_wind1_speed','ndbc_wind1_dir');

% if ( ~exist('mlrf1','var') )
%     mlrf1 = get_station_from_station_name('mlrf1'); mlrf1 = load_all_ndbc_data(mlrf1);
% end;
% mlrf1 = station_spddir_to_uv(mlrf1,'ndbc_wind1_speed','ndbc_wind1_dir');

% if ( ~exist('fwyf1','var') )
%     fwyf1 = get_station_from_station_name('fwyf1'); fwyf1 = load_all_ndbc_data(fwyf1);
% end;
% fwyf1 = station_spddir_to_uv(fwyf1,'ndbc_wind1_speed','ndbc_wind1_dir');

% if ( ~exist('sanf1','var') )
%     sanf1 = get_station_from_station_name('sanf1'); sanf1 = load_all_ndbc_data(sanf1);
% end;
% sanf1 = station_spddir_to_uv(sanf1,'ndbc_wind1_speed','ndbc_wind1_dir');


%% ISOLATE COLD FRONT "EVENTS"

% % % Find cold fronts: pressure >= 1026 hPa
% % cfdts = unique(floor(stn.ndbc_barom.date(stn.ndbc_barom.data>=1026)));

% % Alternate method: pressure >= 99th percentile
% cfdts = unique(floor(stn.ndbc_barom.date(stn.ndbc_barom.data>=prctile(stn.ndbc_barom.data,99))));

% % % "Alternate alternate" method: 1 d median wind <= 10 then >= 20 kts in 48 hours
% % cfdts = unique(floor(stn.ndbc_wind1_speed_1_d_median.date(stn.ndbc_wind1_speed_1_d_median.data<=10)));
% % cfdts = intersect(cfdts,floor(stn.ndbc_wind1_speed_1_d_median.date(stn.ndbc_wind1_speed_1_d_median.data>=20))-2);


% Alternate method: pressure >= 96th percentile for that season (JFM,AMJ,JAS,OND)
for seas=1:4
  seaslim(seas) = prctile(stn.ndbc_barom.data(get_season(stn.ndbc_barom.date)==seas),96);
end;
jdlim = seaslim(get_season(stn.ndbc_barom.date))';
cfdts = unique(floor(stn.ndbc_barom.date(stn.ndbc_barom.data>=jdlim)));

% Alternate method: pressure >= 96th percentile for that part of year (warm or cool season)
jdlim(ts_boreal_cool(stn.ndbc_barom)) = prctile(stn.ndbc_barom.data(ts_boreal_cool(stn.ndbc_barom)),96);
jdlim(ts_boreal_warm(stn.ndbc_barom)) = prctile(stn.ndbc_barom.data(ts_boreal_warm(stn.ndbc_barom)),96);
cfdts = unique(floor(stn.ndbc_barom.date(stn.ndbc_barom.data>=jdlim)));

% Isolate them to individual non-contiguous days of high pressure
cfdts(find(diff(cfdts)<3)+1) = [];


% % % Center each event on the peak barometric pressure
% % for ix = 1:numel(cfdts)
% %   evtix = find(cfdts(ix)<=stn.ndbc_barom.date & stn.ndbc_barom.date<cfdts(ix)+3);
% %   [ig,peakix] = max(stn.ndbc_barom.data(evtix));
% %   cfdts(ix) = stn.ndbc_barom.date(evtix(peakix));
% % end;

% % Center each event on the lowest barometric pressure
% for ix = 1:numel(cfdts)
%   evtix = find(cfdts(ix)-3<=stn.ndbc_barom.date & stn.ndbc_barom.date<cfdts(ix)+1);
%   [ig,minix] = min(stn.ndbc_barom.data(evtix));
%   cfdts(ix) = stn.ndbc_barom.date(evtix(minix));
% end;


% for dt = cfdts(:)'
%     fh=fmg;
%     fn = @(x)(find(abs(x.date-dt)<7));
%     subplot(3,1,1); plot_ts(subset_ts(stn.ndbc_wind1_speed,fn)); xlabel('U10'); axis([dt-7,dt+6.9,0,40]); datetick3;
%     subplot(3,1,2); plot_ts(subset_ts(stn.ndbc_air_t,fn),subset_ts(stn.ndbc_sea_t,fn)); xlabel('Ta'); axis([dt-7,dt+6.9,7,27]); datetick3;
%     subplot(3,1,3); plot_ts(subset_ts(stn.ndbc_barom,fn)); xlabel('Pr'); axis([dt-7,dt+6.9,1000,1033]); datetick3;
%     pause;
%     close(fh); clear fh;
% end;

if ( ~exist('fh','var') || ~ishandle(fh) )
    fh = [];
end;
if ( isempty(fh) )
    fh=fmg;
    %DEBUG:    disp('New fig');

    % %subplot(3,1,1); plot_ts(stn.ndbc_wind1_speed); xlabel('U10'); grid on;
    % subplot(3,1,1); plot_ts(stn.ndbc_wind1_u,stn.ndbc_wind1_v); xlabel('U10, V10'); grid on;
    % hold on; plot_ts(sanf1.ndbc_wind1_u,'-',sanf1.ndbc_wind1_v,'-');     
    % hold on; plot_ts(stn.ndbc_wind1_speed_1_d_median,'r-');     

    % subplot(3,1,2); plot_ts(stn.ndbc_air_t,stn.ndbc_sea_t); xlabel('Ta, Ts'); grid on;
    % hold on; plot_ts(sanf1.ndbc_air_t,'-',sanf1.ndbc_sea_t,'-');
    % hold on; plot_ts(stn.erai_air_t,'m-');
    % hold on; plot_ts(ts_op(stn.erai_spechumid,1e3,'*'),'c-');

    % subplot(3,1,3); plot_ts(stn.ndbc_barom); xlabel('Pr'); grid on;
    % hold on; plot_ts(sanf1.ndbc_barom,'-'); 
    % hold on; plot_ts(stn.erai_barom,'m-'); 

    subplot(3,1,1); plot_ts(stn.ndbc_air_t,stn.ndbc_sea_t); xlabel('Ta, Ts'); grid on;
    % subplot(3,1,2); plot_ts(stn.erai_ndbc_arf_1_d_sum,stn.erai_ndbc_srf_1_d_sum); xlabel(',Q_S_W'); grid on;
    % subplot(3,1,2); plot_ts(stn.ndbc_erai_erai_30a_latent_flux,stn.ndbc_erai_erai_30a_sensible_flux,stn.ndbc_erai_erai_30a_net_flux); xlabel('Q_L_H,Q_S_H,Q_0'); grid on;
    % stn.erai_ndbc_arf_term_dly \Sigma_1_d(Q_S_W(\gamma)+Q_L_W)
    subplot(3,1,2); plot_ts(stn.ndbc_erai_erai_30a_latent_flux_term_1_d_sum,stn.ndbc_erai_erai_30a_sensible_flux_term_1_d_sum,stn.erai_ndbc_arf_term_1_d_sum,stn.ndbc_erai_erai_30a_net_flux_term_1_d_sum); xlabel('\Sigma_1_dQ_L_H,Q_S_H,Q_S_W(\gamma)+Q_L_W,Q_0'); grid on;
    legend('Q_L_H','Q_S_H','Q_S_W(\gamma)+Q_L_W','Q_0');
    subplot(3,1,3); plot_ts(stn.ndbc_barom); xlabel('Pr'); grid on;

end;

if ( ~exist('dtix','var') )
    dtix = 1;
end;
for dtix = dtix:numel(cfdts)
    dt = cfdts(dtix);
    % % %subplot(3,1,1); axis([dt-7,dt+6.9,0,40]); datetick3;
    % % subplot(3,1,1); axis([dt-7,dt+6.9,-20,25]); datetick3;
    % % titlename([num2str(dtix),' of ',num2str(numel(cfdts))]);
    % % subplot(3,1,2); axis([dt-7,dt+6.9,5,29]); datetick3;
    % % subplot(3,1,3); axis([dt-7,dt+6.9,1005,1033]); datetick3;
    % subplot(3,1,1); axis([dt-7,dt+6.9,-600,600]); datetick3;
    % subplot(3,1,2); axis([dt-7,dt+6.9,0,5500]); datetick3;
    subplot(3,1,1); axis([dt-7,dt+6.9,5,29]); datetick3;
    titlename([num2str(dtix),' of ',num2str(numel(cfdts))]);
    subplot(3,1,2); axis([dt-7,dt+6.9,-8,8]); datetick3;
    subplot(3,1,3); axis([dt-7,dt+6.9,1002,1033]); datetick3;
    pause;
    if ( ~ishandle(fh) )
        disp('Quit');
        break;
    end;
    % r = input('"Q"uit, "b"ack, or forward (any) ','s');
    % if ( strncmpi(r,'q',1) || ~ishandle(fh) )
    %     disp('Quit');
    %     break;
    % elseif ( strncmpi(r,'b',1) )
    %     dtix = max(1,dtix - 1);
    % end;
end;









    [ix,jx] = ind2sub(size(stn.(['ndbc_adcp_dir']).prof),hix);





      % Invalid ("999") speeds and dirs may correspond with "valid" depths
      [stn.(['ndbc_adcp_spd_',bns]),stn.(['ndbc_adcp_spd_',bns]),stn.(['ndbc_adcp_spd_',bns])] = ...
          intersect_tses(







    [h,stat] = boxplot_ts(tses{tix},grpfun,'allcol',clr,'widths',wid,varargin{:});
    if ( nargout > 0 )
      hs{end+1} = h;
      if ( nargout > 1 )
        stats{end+1} = stat;
      end;
    end;
  end;




    %fill3(sofla_coast(:,1), sofla_coast(:,2), repmat(0,size(sofla_coast(:,1))), [.4 .3 .2]);




function [ns,xs,fh] = hist_ts_seasons(ts,varargin)
%function [ns,xs,fh] = hist_ts_seasons(ts,varargin)
%
% Make a figure with four subplots, each showing the distribution of time
% series TS for the corresponding season (v. GET_SEASON), 1 through 4.
%
% Last Saved Time-stamp: <Fri 2013-10-25 12:22:33 Eastern Daylight Time gramer>

  fh = fmg;
  for seas = 1:4
    subplot(2,2,seas);
    dat = ts.data(get_season(ts.date)==seas);
    [ns{seas},xs{seas}] = hist(dat,varargin{:});
    xlabel(['Season ',num2str(seas),'  (N=',num2str(numel(dat)),')']);
    dat=[]; clear dat
  end;
  linkaxes;

  if ( nargout < 3 )
    fh=[]; clear fh;
    if ( nargout < 2 )
      xs=[]; clear xs;
      if ( nargout < 1 )
        ns=[]; clear ns;
      end;
    end;
  end;

return;





    if ( isempty(boxclr) )
      plot(cum,['b',meanMarker],'MarkerSize',8);
    elseif ( ischar(boxclr) )
      plot(cum,[boxclr,meanMarker],'MarkerSize',8);
    elseif ( isnumeric(boxclr) && numel(boxclr) == 3 )
      plot(cum,meanMarker,'Color',boxclr,'MarkerSize',8);
    else
      error('ALLCOLOR arg must be a color string or numeric 3-vector!');
    end;







      if ( ~doWarm )

        % Perform Fairall et al. (TOGA-COARE 3.0a) calculation
        disp('COR30 (no warm-layer)');
        disp([num2str(length(w)) ' data points']);
        nby10 = floor(length(w)/10);
        for ix = 1:length(w)
          res = cor30a([w(ix),ou(ix),t(ix),a(ix),ss(ix),sa(ix),dsrf(ix),dlrf(ix),prcp(ix),pblz(ix),p(ix),wz,az,az,stn.lat,1,1,wvper(ix),wvhgt(ix)]);
          % Preallocate for speed
          if ( isempty(result) )
            result = repmat(nan,[length(w) length(res)]);
          end;
          % Algorithm sometimes returns complex numbers??
          result(ix,:) = real( res );
          if (mod(ix,nby10)==0); disp(ix); end;
        end;

      else

        % Perform Fairall et al. (TOGA-COARE 3.0a + WARM CORRECTION) calculation
        disp('COR30A (warm-layer)');
        % Algorithm sometimes returns complex numbers??
        % result = real( cor30a_warm(dts,w,ou,t,a,ss,sa,dsrf,dlrf,prcp,pblz,p,stn.lon,stn.lat,wz,az,az,stz,1,1,1,wvper,wvhgt) );
        % Lew.Gramer@noaa.gov, 2013 Feb 06:
        % JWave==2: Use Taylor & Yelland z0
        result = real( cor30a_warm(dts,w,ou,t,a,ss,sa,dsrf,dlrf,prcp,pblz,p,stn.lon,stn.lat,wz,az,az,stz,1,1,2,wvper,wvhgt,max_wl) );
      end; %if ~doWarm else






  x.ts = ts;
  x.result.date = ts.date(1:dt:end-1);
  x.result.data = diff(ts.data(1:dt:end));




  if ( doMean )
    [cum,tid]=grp_ts(ts.data,ts.date,grpfun,@nanmean,[]);
    x = get(get(h(1),'Parent'),'XTick');
    xs = get(get(h(1),'Parent'),'XTickLabel');
    if ( isempty(boxclr) )
      plot(x,cum,'bs','MarkerSize',6);
    else
      plot(x,cum,[boxclr,'s'],'MarkerSize',6);
    end;
    set(get(h(1),'Parent'),'XTickLabel',xs);
  end;





%{
  if ( isempty(boxclr) )
    h = boxplot(ts.data,dt, 'notch','on', 'whisker',2, varargin{:});
  else
    h = boxplot(ts.data,dt, 'notch','on', 'whisker',2, ...
                'symbol',[boxclr,'x'],'colors',[boxclr;boxclr;boxclr], varargin{:});
  end;

  if ( doMean )
    [cum,tid]=grp_ts(ts.data,ts.date,grpfun,@nanmean,[]);
    x = get(get(h(1),'Parent'),'XTick');
    if ( isempty(boxclr) )
      plot(x,cum,'bs','MarkerSize',6);
    else
      plot(x,cum,[boxclr,'s'],'MarkerSize',6);
    end;
  end;
%}

  if ( doMean )
    [cum,tid]=grp_ts(ts.data,ts.date,grpfun,@nanmean,[]);
    if ( isempty(boxclr) )
      plot(cum,'bs','MarkerSize',6);
    else
      plot(cum,[boxclr,'s'],'MarkerSize',6);
    end;
  end;
  if ( isempty(boxclr) )
    h = boxplot(ts.data,dt, 'notch','on', 'whisker',2, varargin{:});
  else
    h = boxplot(ts.data,dt, 'notch','on', 'whisker',2, ...
                'symbol',[boxclr,'x'],'colors',[boxclr;boxclr;boxclr], varargin{:});
  end;






    for curs = 1:10:length(line)
        idx = idx + 1;
        % hdr{idx} = strtrim(line(curs:curs+9));
        len = min(9,length(line)-curs);
        hdr{idx} = strtrim(line(curs:curs+len));
        % Try to make sure we have what looks like a valid header line
        if ( ~isnan(str2double(hdr{idx})) )
            error('Invalid header line?? Field %d is the number "%s"!', ...
                  idx, hdr{idx});
        end;




From ROSEWGTD.m:
  % Could use ceil(nansum(sts.data)) to preallocate here
  wd=[];

    % Slow (growing a vector) but simple
    wd(end+1:end+n(dr))=repmat(drs(dr),[1,n(dr)]);

  wd(endix+1:end)=[];





            % Otherwise, merge all values - overwriting old data with new
            else
                % First do the delicate merge of datenums
                % Use a ~36-second (0.01-day) tolerance for matching datenums
                % (Narrow window ensures we only remove LIKELY duplicates.)
                newdatenum = union(roundn(stn.(fld).date), roundn(result.(fld).date));
                newdatenum = union(stn.(fld).date, result.(fld).date);
                si = find(ismember(newdatenum, stn.(fld).date));
                ri = find(ismember(newdatenum, result.(fld).date));
                stn.(fld).date = newdatenum;
                stn.(fld).data(si,1) = stn.(fld).data;
                % If the input field was a cell array, keep it that way
                if ( iscell(stn.(fld).data) == iscell(result.(fld).data) )
                    stn.(fld).data(ri,1) = result.(fld).data;
                elseif ( iscell(stn.(fld).data) )
                    stn.(fld).data(ri,1) = cellstr(num2str(result.(fld).data(:)));
                else
                    stn.(fld).data(ri,1) = str2double(result.(fld).data);
                end;
                % Frequently triggered by loading .date as row-, and .data as
                % column-vectors, for example, or vice versa. Stupid MATLAB. 
                if ( any(size(stn.(fld).date) ~= size(stn.(fld).data)) )
                    warning('Fields %s.date and .data do not match!', fld);
                    size(stn.(fld).date), size(stn.(fld).data),
                end;

            end;






function [anomts,clim,tms] = anomalize_daily_mean_ts(ts,varargin)
%function [anomts,clim,tms] = anomalize_daily_mean_ts(ts[,GRP_TS args])
%
% Calculate a year-day climatology CLIM for time series struct TS. Subtract
% the climatology from each year-day of TS to form new time series ANOMTS.
%
% Last Saved Time-stamp: <Wed 2012-10-17 14:04:49 Eastern Daylight Time lew.gramer>

  [clim,tms] = grp_ts(ts.data,ts.date,varargin{:});
  anomts = ts;
  for jdix=1:numel(tms)
    ix = find(get_jday(anomts.date)==tms(jdix));
    anomts.data(ix) = anomts.data(ix) - clim(jdix);
  end;

return;





get_ngdc_bathy_station.m:
  if ( isfield(stn,'station_name') )
    stnm = lower(stn.station_name);
  elseif ( ischar(stn) )
    stnm = lower(stn);
    clear stn;
    stn.station_name = upper(stnm);
  else
    error('First arg should be either station name or struct with .station_name field!');
  end;








  if ( ischar() )
    station = [];
    stnm = station_or_stnm;
  else
    station = station_or_stnm;
    stnm = station.station_name;
  end;




  if ( ~isempty(fh) )
    % Just make sure: calling LOGLOG/SEMILOGX multiple times can be squirrelly
    set(gca,'XScale','log');
    if ( ~presVar )
      set(gca,'YScale','log');
    end;





dat = 24*pos(2)/(2*pi);
dat = pos(2);



error('Wow - I have no idea what DataCursorMode UpdateFcn to use for CONTOURF (which is what WT calls)!');





  % set_hovmuller_cursor;
  h = datacursormode(gcf);
  set(h, 'UpdateFcn', @wt_select_cb);




  % dts = [ts.date(ix(1)):(1/24):ts.date(ix(end))]';
  % dat = interp1(ts.date,ts.data,dts);
  dts = ts.date;
  dat = ts.data;



    min_W_day = min([min_W_day(:);W_day(:)]);
    max_W_day = max([max_W_day(:);W_day(:)]);
    min_Pxx = min([min_Pxx(:);Pxx(:)]);
    max_Pxx = max([max_Pxx(:);Pxx(:)]);




    if ( presVar )
      disp(['Spectrum of ',var,' preserves variance']); 
      %Pxx = Pxx./W_day;
      Pxx = W_day.*Pxx;
      %Pxx = W_hour.*Pxx;
    end;





    if ( presVar )
      disp(['Spectrum of ',var,' preserves variance']); 
      Pxx = W.*Pxx;
      %Pxx = W_day.*Pxx;
      %Pxx = W_hour.*Pxx;
    end;







    [Pxx, W] = feval(specfunc, dat, specargs{:});

    W_hour = (2*pi/per_hr) ./ W;
    W_day = (2*pi/per_day) ./ W;

    min_W_day = min(min_W_day,W_day);
    max_W_day = max(max_W_day,W_day);
    min_Pxx = min(min_Pxx,Pxx);
    max_Pxx = max(max_Pxx,Pxx);

    if ( nvars > 1 )
        Pxxes{end+1} = Pxx;
        W_hours{end+1} = W_hour;
        W_days{end+1} = W_day;
    else
        Pxxes = Pxx;
        W_hours = W_hour;
        W_days = W_day;
    end;

    if ( ~isempty(fh) )
      if ( presVar )
        lhs(end+1) = feval(plotfunc, W_day,W_day.*Pxx, 'Color',co(mod(varix-1,ncos)+1,:));
        %lhs(end+1) = feval(plotfunc, W_hour,W_hour.*Pxx, 'Color',co(mod(varix-1,ncos)+1,:));
      else
        lhs(end+1) = feval(plotfunc, W_day,Pxx, 'Color',co(mod(varix-1,ncos)+1,:));
        %lhs(end+1) = feval(plotfunc, W_hour,Pxx, 'Color',co(mod(varix-1,ncos)+1,:));
      end;
    end;
  end;

  if ( ~isempty(fh) )
    % Just make sure: calling LOGLOG/SEMILOGX multiple times can be squirrelly
    set(gca,'xscale','log');
    if ( ~presVar )
      set(gca,'yscale','log');
    end;




GET_AVHRR_FIELD:
  if ( calcFieldTerms )
    % Do NOT try to construct hourly time series from synoptic data!
    if ( isfield(stn,'lat') && isfield(stn,'lon') )
      stn = calc_field_terms(stn,'avhrr_sst_field','avhrr_sst',interpMethod,stn.lat,stn.lon,npts,false);
    else
      stn = calc_field_terms(stn,'avhrr_sst_field',[],[],[],[],npts,false);
    end;
  end;






  results = repmat(struct('lon',[],'lat',[],'field',[]),size(azs));





function result = station_field_transect(stn,fldnm,rng,az,dts)
%function result = station_field_transect(stn,fldnm,rng,az,dts)
%
% Extract a transect from field STN.(FLDNM) along azimuth AZ out to distance
% RNG (in [km]), beginning at the location of STN (STN.lon,STN.lat). If
% STN.(FLDNM) is a time series of fields, extract transect for each datenum
% in DTS (DEFAULT: all times in time series) using INTERP3 (v.). Otherwise,
% use INTERP2 (v.). In either case, use 'spline'. CALLS: TRANSECT_WGS84.
%
% Last Saved Time-stamp: <Tue 2012-05-08 15:20:13  Lew.Gramer>

  if ( any(~isfield(stn,{'lon','lat'})) )
    error('First arg STN must be a struct with .lon and .lat fields');
  end;

  if ( ~exist('dts','var') || isempty(dts) )
    if ( isfield(stn.(fldnm),'date') )
      dts = stn.(fldnm).date;
    else
      dts = [];
    end;
  end;

  [result.lon,result.lat] = transect_wgs84(stn.lon,stn.lat,rng,az);

  if ( isfield(stn.(fldnm),'date') )
    [X,Y,Z] = meshgrid(stn.(fldnm).lon,stn.(fldnm).lat,stn.(fldnm).date);
    result.date = dts;
    result.field = interp3(X,Y,Z,...
                           permute(stn.(fldnm).field,[2 3 1]), ...
                           result.lon,result.lat,result.date,'spline');
  else
    if ( isvector(stn.(fldnm).lon) || isvector(stn.(fldnm).lat) )
      [LONS,LATS] = meshgrid(unique(stn.(fldnm).lon),unique(stn.(fldnm).lat));
      % result.field = interp2(stn.(fldnm).lon,stn.(fldnm).lat,stn.(fldnm).field, ...
      result.field = interp2(LONS,LATS,stn.(fldnm).field', ...
                             result.lon,result.lat,'nearest');
    else
      LONS = stn.(fldnm).lon;
      LATS = stn.(fldnm).lat;
      result.field = interp2(LONS,LATS,stn.(fldnm).field', ...
                             result.lon,result.lat,'nearest');
    end;
  end;

return;







function result = station_field_transect(stn,fldnm,rng,az,dts)
%function result = station_field_transect(stn,fldnm,rng,az,dts)
%
% Extract a transect from field STN.(FLDNM) along azimuth AZ out to distance
% RNG (in [km]), beginning at the location of STN (STN.lon,STN.lat). If
% STN.(FLDNM) is a time series of fields, extract transect for each datenum
% in DTS (DEFAULT: all times in time series) using INTERP3 (v.). Otherwise,
% use INTERP2 (v.). In either case, use 'spline'. CALLS: TRANSECT_WGS84.
%
% Last Saved Time-stamp: <Fri 2012-04-13 12:19:00  lew.gramer>

  if ( any(~isfield(stn,{'lon','lat'})) )
    error('First arg STN must be a struct with .lon and .lat fields');
  end;

  if ( ~exist('dts','var') || isempty(dts) )
    if ( isfield(stn.(fldnm),'date') )
      dts = stn.(fldnm).date;
    else
      dts = [];
    end;
  end;

  [result.lon,result.lat] = transect_wgs84(stn.lon,stn.lat,rng,az);

  if ( isfield(stn.(fldnm),'date') )
    [X,Y,Z] = meshgrid(stn.(fldnm).lon,stn.(fldnm).lat,stn.(fldnm).date);
    result.date = dts;
    result.field = interp3(X,Y,Z,...
                           permute(stn.(fldnm).field,[2 3 1]), ...
                           result.lon,result.lat,result.date,'spline');
  else
    [LONS,LATS] = meshgrid(unique(stn.(fldnm).lon),unique(stn.(fldnm).lat));
    % result.field = interp2(stn.(fldnm).lon,stn.(fldnm).lat,stn.(fldnm).field, ...
    result.field = interp2(LONS,LATS,stn.(fldnm).field', ...
                           result.lon,result.lat,'nearest');
  end;

return;








    if ( ~isempty(fh) )
      % if ( numel(vars) == 1 )
      %   feval(plotfunc, W_day, Pxx);
      %   %feval(plotfunc, W_hour, Pxx);
      % else
        feval(plotfunc, W_day,Pxx, 'Color',co(mod(varix-1,ncos)+1,:));
      % end;
    end;






/* EXPERIMENT */
/*
raw_avhrr_weekly_sst_field,1993,03,09,1993,05,03
raw_avhrr_weekly_sst_field,1994,03,09,1994,05,03
raw_avhrr_weekly_sst_field,1995,03,09,1995,05,03
raw_avhrr_weekly_sst_field,1996,03,09,1996,05,03
raw_avhrr_weekly_sst_field,1997,03,09,1997,05,03
raw_avhrr_weekly_sst_field,1998,03,09,1998,05,03
raw_avhrr_weekly_sst_field,1999,03,09,1999,05,03
raw_avhrr_weekly_sst_field,2000,03,09,2000,05,03
raw_avhrr_weekly_sst_field,2001,03,09,2001,05,03
raw_avhrr_weekly_sst_field,2002,03,09,2002,05,03
raw_avhrr_weekly_sst_field,2003,03,09,2003,05,03
raw_avhrr_weekly_sst_field,2004,03,09,2004,05,03
raw_avhrr_weekly_sst_field,2005,03,09,2005,05,03
raw_avhrr_weekly_sst_field,2006,03,09,2006,05,03
raw_avhrr_weekly_sst_field,2007,03,09,2007,05,03
raw_avhrr_weekly_sst_field,2008,03,09,2008,05,03
raw_avhrr_weekly_sst_field,2009,03,09,2009,05,03
raw_avhrr_weekly_sst_field,2010,03,09,2010,05,03
raw_avhrr_weekly_sst_field,2011,03,09,2011,05,03
raw_avhrr_weekly_sst_field,2012,03,09,2012,05,03
*/



  if ( ischar(sal) )
    if ( ~isfield(stn,sal) )
      error('No salinity field STN.%s found!',sal);
    else
      mindt = max(mindt,stn.(sal).date(1));
      maxdt = min(maxdt,stn.(sal).date(end));
    end;
  end;






  % Theta_a - Theta_t
  [ix1,ix2] = intersect_dates(stn.ndbc_air_t.date, stn.ndbc_sea_t.date);
  stn.ndbc_air_sea_t.date = stn.ndbc_air_t.date(ix1);
  stn.ndbc_air_sea_t.data = stn.ndbc_air_t.data(ix1) - ...
      stn.ndbc_sea_t.data(ix2);




  % Convert data to standard 10m height assuming a simple log-layer
  % Uses simple algorithms from COR30a by Fairall et al. (2003)
  if ( ~isfield(stn,'ndbc_air_sea_t10m') )
    stn.ndbc_air_sea_t10m.date = stn.ndbc_air_sea_t.date;
    stn.ndbc_air_sea_t10m.data = stn.ndbc_air_sea_t.data - (0.0098*az);
  end;
  if ( ~isfield(stn,'ndbc_wind1_u10m') || ~isfield(stn,'ndbc_wind1_v10m') )
    stn.ndbc_wind1_u10m.date = stn.ndbc_wind1_u.date;
    stn.ndbc_wind1_u10m.data = stn.ndbc_wind1_u.data .* log(10/1e-4)/log(wz/1e-4);
    stn.ndbc_wind1_v10m.date = stn.ndbc_wind1_v.date;
    stn.ndbc_wind1_v10m.data = stn.ndbc_wind1_v.data .* log(10/1e-4)/log(wz/1e-4);
  end;













  holdState = get(ax,'NextPlot');
  hold on;


  % Restore caller's HOLD state
  % set(gca,'NextPlot',holdState);
  if ( strcmpi(holdState,'replace') )
    hold off;
  end;





  % Was 0.96 - why?
  epsilon = 0.97;   % Ocean surface emissivity (see Anderson, 1952)
                    % Xue et al (1998) say: "0.97 after Anderson (1952) and Reed (1976)"
                    % Kraus and Businger suggest 0.98





  epsilon = 0.96;   % Ocean surface emissivity (see Anderson, 1952)







%%%%???DEBUG
  % stn = station_filter_bad_dates(stn);

  % % Offset raw data to time midway between start and end of week (+3.5)
  % stn.raw_avhrr_weekly_sst_field.date = stn.raw_avhrr_weekly_sst_field.date + 3.5;
%%%%???DEBUG





  %% Let user specify BAD weekly composites in each station's AVHRR data
  % (See STATION_FILTER_BAD_DATES.m)

  % stn.avhrr_weekly_sst_field = station.avhrr_weekly_sst;
  stn.raw_avhrr_weekly_sst_field = station.avhrr_weekly_sst;
  % First week of the year always begins on 01 Jan
  wkix = find(ismember(get_jday(station.avhrr_weekly_sst.date),[1:7:358]));
  stn.raw_avhrr_weekly_sst_field.date = station.avhrr_weekly_sst.date(wkix);
  stn.raw_avhrr_weekly_sst_field.field = station.avhrr_weekly_sst.field(wkix,:,:);
  % Field .date should be a column vector (Nx1)
  stn.raw_avhrr_weekly_sst_field.date = stn.raw_avhrr_weekly_sst_field.date(:);
  station = []; clear station;

  stn = station_filter_bad_dates(stn);


  %% Spatially interpolate time series from weekly field

  if ( isfield(stn,'lat') && isfield(stn,'lon') )
    stn.raw_avhrr_weekly_sst.date = stn.raw_avhrr_weekly_sst_field.date;
    stn.raw_avhrr_weekly_sst.data = ...
        interp_field(stn.raw_avhrr_weekly_sst_field.lat,stn.raw_avhrr_weekly_sst_field.lon,...
                     stn.raw_avhrr_weekly_sst_field.field,stn.lat,stn.lon,interpMethod);
  end;


  %% Interpolate daily field from weekly field

  % Redo: EXTRACT_AVHRR_WEEKLY_FIELD did not call STATION_FILTER_BAD_DATES
  dts = stn.raw_avhrr_weekly_sst_field.date;
  sst = stn.raw_avhrr_weekly_sst_field.field;
  % Interpolate to time midway between start and end of week (+3.5)
  alldts = [dts(1):dts(end)]+3.5;
  % Quick and dirty way of finding all non-land pixels
  meansst = squeeze(nanmean(sst,1));
  seaix = find(isfinite(meansst));
  nwks = size(sst,1);

  stn.avhrr_weekly_sst_field.date = alldts(:);
  stn.avhrr_weekly_sst_field.lon = stn.raw_avhrr_weekly_sst_field.lon;
  stn.avhrr_weekly_sst_field.lat = stn.raw_avhrr_weekly_sst_field.lat;
  stn.avhrr_weekly_sst_field.N = repmat(nan, [size(sst,2) size(sst,3)]);
  stn.avhrr_weekly_sst_field.pctused = repmat(nan, [size(sst,2) size(sst,3)]);
  stn.avhrr_weekly_sst_field.field = repmat(nan, [numel(alldts) size(sst,2) size(sst,3)]);
  for ix = seaix(:)'
    usedix = find(isfinite(sst(:,ix)));
    stn.avhrr_weekly_sst_field.N(ix) = numel(usedix);
    stn.avhrr_weekly_sst_field.pctused(ix) = (numel(usedix)/nwks);
    % Interpolate to time midway between start and end of week (+3.5)
    stn.avhrr_weekly_sst_field.field(:,ix) = ...
        interp1(dts(usedix)+3.5,sst(usedix,ix),alldts,'pchip',nan);
  end;
  dts=[]; sst=[]; alldts=[]; clear dts sst alldts;










  % stn.avhrr_weekly_sst_field = station.avhrr_weekly_sst;
  stn.raw_avhrr_weekly_sst_field.lon = station.avhrr_weekly_sst.lon;
  stn.raw_avhrr_weekly_sst_field.lat = station.avhrr_weekly_sst.lat;
  stn.raw_avhrr_weekly_sst_field.date = station.avhrr_weekly_sst.date(1:7:end);
  stn.raw_avhrr_weekly_sst_field.field = station.avhrr_weekly_sst.field(1:7:end,:,:);
  stn.raw_avhrr_weekly_sst_field.N = station.avhrr_weekly_sst.N(1:7:end,:,:);
  stn.raw_avhrr_weekly_sst_field.pctused = station.avhrr_weekly_sst.pctused(1:7:end,:,:);






          elseif ( ~isfield(stn,fld) )
            % Allow this without annoying warnings for now
            %warning('Ecoforecasts:BadDates:BadField','No field "%s" in STN',fld);





                if ( ndims(stn.(fld).field)==2 )
                  stn.(fld).field(dtix,:) = [];
                elseif ( ndims(stn.(fld).field)==3 )
                  stn.(fld).field(dtix,:,:) = [];
                elseif ( ndims(stn.(fld).field)==4 )
                  stn.(fld).field(dtix,:,:,:) = [];
                end;






function newts = subset_ts(ts,idx)
%function newts = subset_ts(ts,idx)
%
% Return a new time series struct NEWTS containing only those indices of time
% series TS selected by IDX, a logical or index vector, or a function_handle
% for a function accepting a time series and returning a vector of indices.
% TS may also be an array or a cell array of time series structs, in which
% case IDX is applied
%
% Last Saved Time-stamp: <Fri 2012-03-16 12:31:08  Lew.Gramer>

  if ( iscell(ts) )

  else

    if ( isa(idx,'function_handle') )
      idx = idx(ts);
    elseif ( (~isnumeric(idx) && ~islogical(idx)) || ~isvector(idx) )
      error('Second arg IDX must be a function handle, numeric- or logical-vector!');
    end;

  newts.date = ts.date(idx);
  newts.data = ts.data(idx);

return;






    stn.ndbc_stats{fldix,1} = fld;
    stn.ndbc_stats{fldix,2} = mea_tid;
    stn.ndbc_stats{fldix,3} = mea;
    stn.ndbc_stats{fldix,4} =  sd;
    stn.ndbc_stats{fldix,5} = med;
    stn.ndbc_stats{fldix,6} = iqr;




  % ok(end+1) = input(['Week #',num2str(wk),' (',datestr(fld.date(ixen(wk))),') OK? ']);
  yorn = input(['Week #',num2str(wk),' (',datestr(fld.date(ixen(wk))),') OK? '],'s');
  if ( islogical(yorn) )
  elseif ( ischar(yorn) )
    yorn = strcmpi(yorn(1),'y');
  elseif ( isnumeric(yorn) )
    yorn = (yorn ~= 0);
  else
    error('User response should be logical, numeric, or "Y" or "N"');
  end;
  ok(end+1) = yorn;








    fdir_clear_sky = 0.75; % Approximation for Keys - should depend on zenith angle, etc.





    % Fraction of direct vs. total (direct+diffuse) insolation
    % Clear sky
    % fdir_clear_sky = 0.85; % From Jin et al. (2011) Fig. 8, ignoring zenith angle
    fdir_clear_sky = 0.95; % Approximation for Keys - should depend on zenith angle, etc.

    % Varying cloud conditions (e.g., from airport reports or reanalysis)
    fdir = fdir_clear_sky .* (1 - C');





% badix = find( ~isfinite(x) | ~isfinite(y) );
% x(badix) = [];
% y(badix) = [];





  % Extract time series structs to intersect dates for
  ts = {};
  nts = 0;
  while ( argix <= nargin )
    arg = varargin{argix};
    if ( isstruct(arg) && isfield(arg,'date') && isfield(arg,'data') )
      % Gracefully deal with vectors OR MATRICES of structs passed as a single argument
      for strcix = 1:numel(arg)
        strc = arg(strcix);
        if ( isempty(strc.date) || any(size(strc.date) ~= size(strc.date)) )
          error('Struct arguments must be non-empty time series!');
        end;
        nts = nts + 1;
        ts{nts} = strc;
      end;
    else
      error('All args after optional TOL must be time series structs or arrays of TS structs!');
    end;
    argix = argix + 1;
  end;







    % ds = dir('../Xiaofang/AIMS_*.nc');




  [climts.data,climts.date,nPerYr,dt] = grp_ts(ts.data,ts.date,grpargs{:});

  dtstr=[];
  if ( 8785 <= max(climts.date) && max(climts.date) < 17569 )
    % Year-half-hour
    climts.date=datenum(1,1,1)+(climts.date.*48);
    dtstr=6;
  elseif ( 367 <= max(climts.date) && max(climts.date) < 8785 )
    % Year-hour
    climts.date=datenum(1,1,1)+(climts.date.*24);
    dtstr=6;
  elseif ( 74 <= max(climts.date) && max(climts.date) < 367 )
    % Year-day or Julian day
    climts.date=datenum(1,1,1)+(climts.date-1);
    dtstr=6;
  elseif ( 53 <= max(climts.date) && max(climts.date) < 74 )
    % Pentad
    climts.date=datenum(1,1,1)+((climts.date-1).*5);
    dtstr=6;
  elseif ( 25 <= max(climts.date) && max(climts.date) < 53 )
    % Week
    climts.date=datenum(1,1,1)+((climts.date-1).*7);
    dtstr=6;
  elseif ( 13 <= max(climts.date) && max(climts.date) < 25 )
    % Hour of day
    climts.date=datenum(1,1,1)+(climts.date./24);
    dtstr=15;
  elseif ( 5 <= max(climts.date) && max(climts.date) < 13 )
    % Month
    climts.date=datenum(1,climts.date,1);
    dtstr=3;
  elseif ( max(climts.date) < 5 )
    % Season / Quarter
    climts.date=datenum(1,(((climts.date-1).*3)+1),1);
    dtstr=18;
  end;









  [climts.data,climts.date,nPerYr] = grp_ts(ts.data,ts.date,grpargs{:});

  dtstr=[];
  if ( 8785 <= max(climts.date) && max(climts.date) < 17569 )
    % Year-half-hour
    climts.date=datenum(1,1,1)+(climts.date.*48);
    dtstr=6;
  elseif ( 367 <= max(climts.date) && max(climts.date) < 8785 )
    % Year-hour
    climts.date=datenum(1,1,1)+(climts.date.*24);
    dtstr=6;
  elseif ( 74 <= max(climts.date) && max(climts.date) < 367 )
    % Year-day or Julian day
    climts.date=datenum(1,1,1)+(climts.date-1);
    dtstr=6;
  elseif ( 53 <= max(climts.date) && max(climts.date) < 74 )
    % Pentad
    climts.date=datenum(1,1,1)+((climts.date-1).*5);
    dtstr=6;
  elseif ( 25 <= max(climts.date) && max(climts.date) < 53 )
    % Week
    climts.date=datenum(1,1,1)+((climts.date-1).*7);
    dtstr=6;
  elseif ( 13 <= max(climts.date) && max(climts.date) < 25 )
    % Hour of day
    climts.date=datenum(1,1,1)+(climts.date./24);
    dtstr=15;
  elseif ( 5 <= max(climts.date) && max(climts.date) < 13 )
    % Month
    climts.date=datenum(1,climts.date,1);
    dtstr=3;
  elseif ( max(climts.date) < 5 )
    % Season / Quarter
    climts.date=datenum(1,(((climts.date-1).*3)+1),1);
    dtstr=18;
  end;










    % Season / Quarter
    climts.date=datenum(1,1,1)+((climts.date-1)*92);






function h = boxplot_ts(ts,grpfun,varargin)
%function h = boxplot_ts(ts,grpfun,varargin)
%
% Create a BOXPLOT of the time series struct TS, grouping elements of TS by
% GRPFUN(TS.date). GRPFUN may be a function handle (DEFAULT: @GET_MONTH) or a
% string specifying grouping: 'hour', 'day', 'jday', 'annual' (uses grouping
% function GET_YEARDAY), 'pentad', 'week', 'month', 'year', or the name of
% any MATLAB function. All other args are name-value pairs: see BOXPLOT for a
% list, except 'fh' which specifies a figure handle; 'title' a TITLENAME (v.)
% or title for the figure; 'allcolors' a single color to use for both box and
% outliers; and 'indices' a vector of indices to include in BOXPLOT, or else
% a function-handle which accepts a time series and returns an index vector.
%
% DEFAULTS for BOXPLOT (v.): 'notch','on', 'whisker',2
%
% Last Saved Time-stamp: <Tue 2011-11-08 11:16:16  Lew.Gramer>

  h = [];

  if ( ~exist('ts','var') || ~is_valid_ts(ts) )
    error('First arg TS must be a valid time series struct');
  end;

  if ( ~exist('grpfun','var') || isempty(grpfun) )
    grpfun = @get_month;
  elseif ( ischar(grpfun) && isvector(grpfun) )
    f = grpfun(1:min(3,length(grpfun)));
    switch (lower(f)),
     case 'hou', grpfun=@get_hour;
     case 'day', grpfun=@floor;
     case 'jda', grpfun=@get_jday;
     case 'ann', grpfun=@get_yearday;
     case 'pen', grpfun=@get_pentad;
     case 'wee', grpfun=@get_week;
     case 'mon', grpfun=@get_month;
     case 'yea', grpfun=@get_year;
     otherwise,  grpfun=str2func(grpfun);
    end;
  end;
  if ( ~isa(grpfun,'function_handle') )
    error('Optional 2nd arg GRPFUN must be a grouping-period name or function name/handle');
  end;

  [fh,varargin] = getarg(varargin,'FH');
  if ( isempty(fh) )
    fh = fmg;
  end;

  [ttl,varargin] = getarg(varargin,'titl');
  if ( ~isempty(ttl) )
    if ( exist('titlename') > 1 )
      titlename(ttl);
    else
      title(ttl);
    end;
  end;

  [idxs,varargin] = getarg(varargin,'IND','default',1:numel(ts.date));
  if ( ischar(idxs) || isa(idxs,'function_handle') )
    idxs = feval(idxs,ts);
  end;

  [boxclr,varargin] = getarg(varargin,'ALLCOL');
  if ( isempty(boxclr) )
    h = boxplot(ts.data,grpfun(ts.date), 'notch','on', 'whisker',2, varargin{:});
  else
    h = boxplot(ts.data,grpfun(ts.date), 'notch','on', 'whisker',2, ...
                'symbol',[boxclr,'x'],'colors',[boxclr,boxclr,boxclr], varargin{:});
  end;

  if ( nargout < 1 )
    clear h;
  end;

return;













  %[idxs,varargin] = getarg(varargin,'ind','default',1:numel(ts.date));
  ix = strmatch('ind',lower(varargin));
  if ( isempty(ix) )
    idxs = 1:numel(ts.date);
  else
    idxs = varargin{ix+1};
    varargin = varargin([1:ix-1,ix+2:end]);
    if ( ischar(idxs) || isa(idxs,'function_handle') )
      idxs = feval(idxs,ts);
    end;
  end;



  % [fh,varargin] = getarg(varargin,'fh');
  ix = strmatch('fh',lower(varargin));
  if ( isempty(ix) )
    fh = fmg;
  else
    fh = varargin{ix+1};
    varargin = varargin([1:ix-1,ix+2:end]);
  end;


  % [boxclr,varargin] = getarg(varargin,'allcol');
  ix = strmatch('allcol',lower(varargin));
  if ( isempty(ix) )
    boxplot(ts.data,grpfun(ts.date), 'notch','on', 'whisker',2, varargin{:});
  else
    boxclr = varargin{ix+1};
    varargin = varargin([1:ix-1,ix+2:end]);
    boxplot(ts.data,grpfun(ts.date), 'notch','on', 'whisker',2, ...
            'symbol',[boxclr,'x'],'colors',[boxclr,boxclr,boxclr], varargin{:});
  end;









  if ( ~exist('boxclr','var') || isempty(boxclr) )
  else
    boxplot(ts.data,grpfun(ts.date), 'notch','on', 'whisker',2, ...
            'symbol',[boxclr,'x'],'colors',[boxclr,boxclr,boxclr]);
  end;





boxclr,fh


    stn = calc_field_terms(stn,rotfld,[],[],[],[],grdtmplt);








    x.rf.field = stn.(rotfld).field







  ulon = unique(stn.(fld).lon(:));
  nlon = numel(ulon);
  ulat = unique(stn.(fld).lat(:));
  nlon = numel(nlon);

  if ( ndims(stn.(fld).field) == 2 )
    % Special handling for NxM (non-time series) fields
    dim1 = size(stn.(fld).field,1);
    dim2 = size(stn.(fld).field,2);
  else
    dim1 = size(stn.(fld).field,2);
    dim2 = size(stn.(fld).field,3);
  end;

  if ( dim1 ~= nlon && dim2 ~= nlat )
    if ( dim2 ~= nlon && dim1 ~= nlat )
      error('Size mismatch in STN.(%s)',fld);
    else
    end;
  end;

  [vlat,vlon]=reckon(stn.lat,stn.lon,[-6:6]/111,ori);
  [hlat,hlon]=reckon(stn.lat,stn.lon,[-6:6]/111,ori-90);

  if ( ~exist('ufld','var') || isempty(ufld) )
    ufld = 'gradient_x';
  end;
  if ( ~exist('vfld','var') || isempty(vfld) )
    vfld = 'gradient_y';
  end;


  % If caller gave us none, try to munge reasonable fieldnames. But NEVER
  % overwrite U and V fields, unless caller told us to explicitly. NOTE: The
  % UFLD~'_x' and VFLD~'_y' cases are useful for field GRADIENT time series.

  if ( ~exist('xfld','var') || isempty(xfld) )
    xfld = regexprep(ufld,'_u$','_xshore');
    if ( strcmpi(ufld,xfld) || strcmpi(vfld,xfld) )
      xfld = regexprep(ufld,'_u_','_xshore_');
      if ( strcmpi(ufld,xfld) || strcmpi(vfld,xfld) )
        xfld = regexprep(ufld,'_x$','_xshore');
        if ( strcmpi(ufld,xfld) || strcmpi(vfld,xfld) )
          error('Specify a unique XFLD name, or use UFLD again to overwrite!');
        end;
      end;
    end;
  end;
  if ( ~exist('lfld','var') || isempty(lfld) )
    lfld = regexprep(vfld,'_v$','_lshore');
    if ( strcmpi(ufld,lfld) || strcmpi(vfld,lfld) )
      lfld = regexprep(vfld,'_v_','_lshore_');
      if ( strcmpi(ufld,lfld) || strcmpi(vfld,lfld) )
        lfld = regexprep(vfld,'_y$','_lshore');
        if ( strcmpi(ufld,lfld) || strcmpi(vfld,lfld) )
          error('Specify a unique LFLD name, or use VFLD again to overwrite!');
        end;
      end;
    end;
  end;

  u = stn.(fld).(ufld);
  v = stn.(fld).(vfld);

  % NOTE: REORIENT_VECTORS preserves input matrix shape, incl. time dimension
  [x,l] = reorient_vectors(ori,u,v);

  stn.(fldnm).(xfld) = x;
  stn.(fldnm).(lfld) = l;

  if ( nargout < 2 )
    clear xfld lfld;
  end;












,ufld,vfld,xfld,lfld







    topdlon = min(distance_wgs84(lats(1),lons(1),lats(1),lons(2)));
    btmdlon = min(distance_wgs84(lats(end),lons(1),lats(end),lons(2)));
    topdlat = min(distance_wgs84(lats(1:end-1),lons(1),lats(2:end),lons(1)));
    btmdlat = min(distance_wgs84(lats(1:end-1),lons(end),lats(2:end),lons(end)));
    dx = min([topdlon,btmdlon,topdlat,btmdlat]);






    btmdlon = min(distance_wgs84(lats(end),lons(1:end-1),lats(end),lons(2:end)));
    topdlat = min(distance_wgs84(lats,lons(1),lats,lons(2:end)));
    btmdlat = min(distance_wgs84(lats,lons(end),lats,lons(1:end-1)));




  [llats,llons]=reckon_wgs84(sitelat,sitelon,rngs,ori);
  [xlats,xlons]=reckon_wgs84(sitelat,sitelon,rngs,ori-90);

  newlats = llats;
  newlons = xlons;
  [NEWLAT,NEWLON] = meshgrid(newlats,newlons);






, 'positions',[1:52]



    xlbl = [xlbl ' (' upper(strrep(func2str(fix1),'_','\_')) ')'];
    ylbl = [ylbl ' (' upper(strrep(func2str(fix2),'_','\_')) ')'];





  if ( isstruct() )
    stn = stnm_or_stn;
  elseif ( ischar(stnm_or_stn) )
    stn.station_name = stnm_or_stn;
  else
    error('STNM_OR_STN must either be station struct or station name string!');
  end;




  % Conversion factor from Insolation to PAR
  % See Papaioannou, Papanikolaou and Retalis, 1993
  PAR_PER_INSOL = 0.473;

  % Conversion factor from W/m^2 to micromole quanta/m^2.s
  % See Morel and Smith, 1974
  UMOL_PER_WATT = 4.1513469579;
  % NOTE: Dye (JGR-A, 2004) recommends 4.56 mumol/J instead!






function g = par_to_insol(p)
%function g = par_to_insol(p)
%
% P is a vector or scalar of PAR in micro-mole quanta / m^2 / s. Returned
% vector or scalar G is insolation in Watts/m^2. Conversion factor from
% Insolation to PAR is derived from Papaioannou, Papanikolaou and Retalis,
% 1993. Conversion factor from W/m^2 to micromole quanta/m^2.s is derived
% from Morel and Smith, 1974.
%
% SEE ALSO: INSOL_TO_PAR
%
% Last Saved Time-stamp: <Tue 2011-11-01 15:21:10  lew.gramer>





function [p,w] = insol_to_par(g)
%function [p,w] = insol_to_par(g)
%
% G is a vector or scalar of insolation in Watts/m^2. Returned vector or
% scalar P is PAR in micro-mole quanta / m^2 / s, W is PAR in
% Watts/m^2. Conversion factor from Insolation to PAR is derived from
% Papaioannou, Papanikolaou and Retalis, 1993. Conversion factor from W/m^2
% to micromole quanta/m^2.s is derived from Morel and Smith, 1974.
%
% SEE ALSO: PAR_TO_INSOL
%
% Last Saved Time-stamp: <Sat 2011-10-22 17:54:22  lew.gramer>

  % Conversion factor from Insolation to PAR
  % See Papaioannou, Papanikolaou and Retalis, 1993
  PAR_PER_INSOL = 0.473;

  % Conversion factor from W/m^2 to micromole quanta/m^2.s
  % See Morel and Smith, 1974
  UMOL_PER_WATT = 4.1513469579;
  % NOTE: Dye (JGR-A, 2004) recommends 4.56 mumol/J instead!

  w = g .* PAR_PER_INSOL;
  p = w .* UMOL_PER_WATT;

  if ( nargout < 2 )
    clear w;
  end;

return;







    set(gca,'XMinorGrid','off');
    set(gca,'YMinorGrid','on');



    pwrrng = log10(max(pwrlims))-log10(min(pwrlims));
    topsz = power(10,pwrrng);
    btmsz = power(10,-pwrrng);




    xmin = min([W_day(:);(2/24)]);
    xmax = min([max(W_day),(365*25)]);
    ymin = min([Pxx(:);1e-7]);
    ymax = min([max(Pxx),1e+7]);





%DEBUG
save('oq0_in.mat','w','ou','t','a','ss','sa','dsrf','dlrf','prcp','pblz','p','wz','az','az','wvper','wvhgt');


%DEBUG
save('shb_in.mat','w','ou','t','a','ss','sa','dsrf','dlrf','prcp','pblz','p','wz','az','az','wvper','wvhgt');



    C = textscan(fid,'%d,%[^,],%s','Headerlines',2,'Delimiter','');
    C = textscan(fid,'%d,%[^,],%s','Headerlines',2,'Delimiter','');




  figure;
  maxigraph;
  hold on;
  grid on;





% Simple QA - assumes *cross-shore* component must generally be quite small
%badix = find(nanmax(abs(u.prof),[],2)>0.5);
badix = find( ~isfinite(u_sfc.data) | (nanmax(abs(x.prof),[],2)>1.0) );
if ( ~isempty(badix) )
  disp(['Clearing ' num2str(numel(badix)) ' bad records']);
  sfc_bin(badix) = nan;
  sfc_hgt(badix) = nan;
  last_top_bin(badix) = nan;
  u_avg.data(badix) = nan;
  u_sfc.data(badix) = nan;
  u_btm.data(badix) = nan;
  v_avg.data(badix) = nan;
  v_sfc.data(badix) = nan;
  v_btm.data(badix) = nan;
end;





  pngbytes = imread(fpath);

  % Calculate SST from 1-byte PNG (colormap) values
  sstbytes = cast(pngbytes,'double');
  sst = (22.0 * sstbytes / 235.0) + 10.0;

  % Per Chuanmin - values above 235 (32oC) are mask values(?)
  sst(sstbytes >= 236) = nan;

  % % Just in case our assumptions about cloud mask are wrong!
  % sst(minval > sst | sst > maxval) = nan;

  pngbytes = []; clear pngbytes
  sstbytes = []; clear sstbytes






          %DEBUG:
          if ( length(ts.date) < 0.8*length(oldts.date) )
            fmg; plot_ts(oldts,ts,'r'); titlename(strrep([stnm ' ' fname],'_','\_'));
          end;









    endix = find(get_hour(ts.date(1:spikeix(endspikeix)-1))>=4,1,'last');


    begix = spikeix(begspikeix)+1;
    begix = find(get_hour(ts.date(spikeix(begspikeix)+1:end))==5,1,'first');







          fmg; plot_ts(ts); titlename(strrep([stnm ' ' fname],'_','\_'));






  linspec = cell(size(u(:))); linspec(:) = {'k'}; linspec{end} = 'r';
  qh = quiver( [LON(:);legendlon],...
               [LAT(:);legendlat],...
               [u(:);1.0],[v(:);0.0], 0.9, linspec );






get(qh),


qh2=quiver(legendlon,legendlat,1.0,0.0,0.9,'r');
get(qh2),






  text(legendlon+dlon/2,legendlat-dlat/2,'  1 m/s  ', ...
       'HorizontalAlignment','Left', ...
       'VerticalAlignment','Bottom', ...
       'BackgroundColor',[1,1,1]);



  legend(qh(end),'1 m/s', 'Location','NorthWest');







  if ( isscalar(dtrng) )
    [ig,dtix]=min(abs(dtrng-stn.(ufld).date));
    dtstr = datestr(stn.(ufld).date(dtix));
  else
    dtix=find(min(dtrng)<=stn.(ufld).date & stn.(ufld).date<=max(dtrng));
    dtstr = [datestr(stn.(ufld).date(dtix(1))) ' - ' datestr(stn.(ufld).date(dtix(end)))];
  end;







  if ( isscalar(dtrng) )
    [ig,dtix]=min(abs(dtrng-stn.(ufld).date));
    dtstr = datestr(stn.(ufld).date(dtix));
  else
    dtix=find(min(dtrng)<=stn.(ufld).date & stn.(ufld).date<=max(dtrng));
    dtstr = [datestr(stn.(ufld).date(dtix(1))) ' - ' datestr(stn.(ufld).date(dtix(end)))];
  end;






        %dat = str2num(x(valix})';





        x{openix} = '-9.00';
        valcelix = regexp(x(openix(1):end),'^[0-9]+[.][0-9]+$');
        valix = find(~cellfun(@isempty,valcelix));
        if ( isempty(valix) )
          warning('No sea temperature data found in "%s"',fpath);
          continue;
        end;
        dts = begdt + ([0:numel(valix)-1]*(2/24))';
        %dat = str2num(x(valix})';
        dat = cellfun(@str2num,x(valix),'ErrorHandler',@get_fknms_thermistors_naner)';
        dts(0>=dat | dat>=40) = [];
        dat(0>=dat | dat>=40) = [];
        if ( isempty(dts) )
          warning('No valid sea temperature data found in "%s"',datfpath);
          dts=[]; dat=[]; clear dts dat;
          continue;
        end;






        %dat = str2num(x(valix})';
        datcel = cellfun(@str2num,x(valix),'ErrorHandler',@get_fknms_thermistors_naner);
        dat = [datcel{:}]';





        datcel = cellfun(@str2num,x(valix),'UniformOutput',false);




        fpath = fullfile(fknmspath,stns.(stnm).dir_name,fname);
        x = textread(fpath,'%s');
        dtix = strmatch('Deployed:',x);
        if ( isempty(dtix) )
          warning('No "Deployed" line found in "%s"',fpath);
          continue;
        end;
        begdt = datenum([x{dtix+1} ' ' x{dtix+2}]);
        openix = strmatch('Open',x);
        if ( isempty(openix) )
          warning('No "Open" tags found in "%s"',fpath);
          continue;
        end;
        valcelix = regexp(x(openix(end)+1:end),'^[0-9]+[.][0-9]*$');
        valix = find(~cellfun(@isempty,valcelix));
        if ( isempty(valix) )
          warning('No sea temperature data found in "%s"',fpath);
          continue;
        end;
        dts = begdt + ( ([1:numel(valix)]+numel(openix)-1)*(2/24) )';
        %dat = str2num(x(valix})';
        datcel = cellfun(@str2num,x(valix),'UniformOutput',false);
        dat = [datcel{:}]';
        result.fknms_seatemp.date = dts;
        result.fknms_seatemp.data = dat;
        stns.(stnm) = merge_station_data(stns.(stnm),result);
        result=[]; clear result;
        %DEBUG:
        disp(['Processed ' fpath]);









    stns.(stnm).fknms_seatemp = struct('date',[],'data',[]);





  % Add datestrings to X labels
  if ( exist('datetick3','file') )
    datetick3;
  elseif ( exist('datetick2','file') )
    datetick2('x',2,'keeplimits','keepticks');
  else
    datetick('x',2,'keeplimits','keepticks');
  end;






          % 145.7642780, ...      % LLBP7
          % 15.161, ...           % LLBP7




% h3 = ...
%     plot(ax(2),...
%          lciy2.bic_deep_330nm_daily_dose_3_day_average.date,...
%          lciy2.bic_deep_330nm_daily_dose_3_day_average.data);




  for ix = 1:size(x.data)
    stns(ix).station_name = regexprep(x.textdata{ix+1,1},'[^A-Z0-9_]','_');
    stns(ix).long_name = x.textdata{ix+1,2};

    stns(ix).lat = anfknms_dm2degree(x.textdata{ix+1,8});
    stns(ix).lon = -anfknms_dm2degree(x.textdata{ix+1,9});

    stns(ix).depth = x.data(ix) ./ 3.2808399;
  end;






          'DewPt',			'air_t_dewp' ; ...
          'DewPt_60',			'air_t_dewp' ; ...
          'Dew_Point',			'air_t_dewp' ; ...





function stn = station_relhumid_to_dewp(stn, afld, qfld, dfld)
%function stn = station_relhumid_to_dewp(stn, afld, qfld, dfld)
%
% Call RELHUMID_TO_DEWP (v.) on data in fields name AFLD (air temperature)
% and QFLD (Relative Humidity 0-100%, assumed to be measured at the same
% elevation), to create a new field DFLD containing dewpoint temperature.
%
% Last Saved Time-stamp: <Sat 2011-08-06 23:42:44  Lew.Gramer>

  ix = find(~isnan(stn.(qfld).data));

  [aix,qix] = intersect_dates(stn.(afld).date, stn.(qfld).date(ix));

  if ( isfield(stn, dfld) )
    stn = rmfield(stn, dfld);
  end;
  stn.(dfld).date = stn.(qfld).date(ix(qix));
  stn.(dfld).data = relhumid_to_dewp(stn.(afld).data(aix), stn.(qfld).data(ix(qix)));

  stn.(dfld).date(isnan(stn.(dfld).data)) = [];
  stn.(dfld).data(isnan(stn.(dfld).data)) = [];

return;






function stn = station_relhumid_to_spechumid(stn,afld,qfld,sfld)
%function stn = station_relhumid_to_spechumid(stn,afld,qfld,sfld)
%
% Call RELHUMID_TO_SPECHUMID (v.) on data in the two fields STN.(AFLD) (air
% temperature) and STN.(QFLD) (Relative Humidity - assumed to be measured or
% calculated for the same elevation), to create a new field STN.(SFLD)
% containing Specific Humidity [kg/kg]. (NOTE: STN.SFLD is removed first!)
%
% Last Saved Time-stamp: <Fri 2011-07-01 08:23:27  Lew.Gramer>

  ix = find(~isnan(stn.(qfld).data));

  [aix,qix] = intersect_dates(stn.(afld).date, stn.(qfld).date(ix));

  if ( isfield(stn, sfld) )
    stn = rmfield(stn, sfld);
  end;
  stn.(sfld).date = stn.(qfld).date(ix(qix));
  stn.(sfld).data = relhumid_to_spechumid(stn.(afld).data(aix), stn.(qfld).data(ix(qix)));

  stn.(sfld).date(isnan(stn.(sfld).data)) = [];
  stn.(sfld).data(isnan(stn.(sfld).data)) = [];

return;







    %%%% ??? HACKOLA - need to figure this out on the fly...
    region = 'freef';

    % Get region boundaries for later reference
    [lons,lats,ig,dx] = read_misst_region(region, 0, 0);

    for yr = 2002 : 2010
      % Ignore leap years until we get more data
      for jd = 1 : 365
        ix = ix + 1;
        dt = datenum(yr,1,1, 12,0,0) + jd - 1;
        [ig,mn,dy] = datevec(dt);
        wk = floor((jd-1)/7) + 1; wk(wk > 52) = 52;

        [lon,lat,sst] = read_misst_region(region, yr, jd);
        if ( ~isempty(sst) )












  % For most coral reef areas, "bathymetry" above MHHW is just a nuisance
  contourf(stn.ngdc_92m_bathy.lon,stn.ngdc_92m_bathy.lat,stn.ngdc_92m_bathy.field,my_contours);
  maxigraph;
  colorbar;
  % Plot station location as a white pentagram
  if ( isfield(stn,'lon') && isfield(stn,'lat') )
    stnlon = stn.lon;
    stnlat = stn.lat;
  else
    stnlon = mean(stn.ngdc_92m_bathy.lon(:));
    stnlat = mean(stn.ngdc_92m_bathy.lat(:));
  end;
  plot(stnlon,stnlat,'wp', 'MarkerEdgeColor','black', 'MarkerFaceColor','white');









          -82.862, ...          % DRYF1
          24.638, ...           % DRYF1
          1.0, ...              % DRYF1


          -82.770, ...          % PLSF1
          24.690, ...           % PLSF1
          1.0, ...              % PLSF1







  if ( iscellstr(flds) )

  else
    s = rmfield(s,flds(~badix,:));
    for ix=find(badix)'
      warning('rmfield:NoField','Missing field "%s"',strtrim(flds(ix,:)));
    end;    
  end;


function s = safe_rmfield(s,flds)
%function s = safe_rmfield(s,flds)
%
% Identical to RMFIELD (v.), except that when any field in FLDS does not
% exist in struct (or matrix of structs) S, no error is returned: instead, a
% WARNING with ID 'rmfield:NoField' (which may be toggled OFF) is displayed.
%
% Last Saved Time-stamp: <Thu 2011-06-02 10:30:04  lew.gramer>

  if ( ~iscellstr(flds) && ~ischar(flds) )
    error('Second arg must be a string array or cell array of strings!');
  end;

  all_flds = fieldnames(s);
  badix = (~ismember(flds,all_flds));
  if ( iscellstr(flds) )
    s = rmfield(s,flds(~badix));
    for ix=find(badix)'
      warning('rmfield:NoField','Missing field "%s"',flds{ix});
    end;    
  else
    s = rmfield(s,flds(~badix,:));
    for ix=find(badix)'
      warning('rmfield:NoField','Missing field "%s"',strtrim(flds(ix,:)));
    end;    
  end;

return;









  % Remove every whole day missing more than 3 hours of data
  badix = find(n<21);
  for ix=badix(:)'
    newts.data(ismember(dys,dys(badix)) = [];
    newts.date(dys==dys(ix)) = [];
  end;  




FROM function output_txt = tyz_select_cb(obj,event_obj):
output_txt = {['X: ',datestr(pos(1)),' (',jday,')'],...
    ['Y: ',sprintf('%f', pos(2)),pointno]};





  % % Project wind onto peak wave source direction - ??
  % w = w .* abs(cosd(wd - td));

  stn.(ssfld).date = stn.(wsfld).date(wsix);
  stn.(ssfld).data = stokes_drift(w,hs,tp);

  % Convert wave/wind SOURCE direction into Stokes TARGET direction
  td = td - 180;
  td(td<0) = td(td<0) + 360;

  stn.(sdfld).date = tds.date(tdix);
  % stn.(sdfld).data = td;
  % % Assume quasi-Eulerian flow directed 35o to right of wind direct?!
  % wd = wd + 35;
  % wd(wd>=360) = wd(wd>=360) - 360;
  stn.(sdfld).data = wd;


  [ssix,sdix] = intersect_dates(stn.(ssfld).date,stn.(sdfld).date);

  stn.(sufld).date = stn.(ssfld).date(ssix);
  stn.(svfld).date = stn.(ssfld).date(ssix);
  [stn.(sufld).data,stn.(svfld).data] = ...
      spddir_to_uv(stn.(ssfld).data(ssix),stn.(sdfld).data(sdix));









  cm = [ repmat([0 0 1],[n,1]) ; ...
         1 1 1 ; ...
         repmat([1 0 0],[n,1]) ];         

  cm = cm .* log(







  % [rng,az] = distance(lat1,lon1, lat2,lon2, [6356.752 (1/298.25722356)]);

  [rng,az] = distance(lat1,lon1, lat2,lon2, wgs84_ellipsoid());





      % If input matrices had the same SIZE, RESHAPE output matrices to it 
      if ( ndims(u)==ndims(v) && all(size(u)==size(v)) )
        x = reshape(x,size(u));
        l = reshape(l,size(u));
      end;





  if ( ~isempty(xlbl) )
    xlabel(xlbl);
  end;
  if ( ~isempty(ylbl) )
    ylabel(ylbl);
  end;








,doPrint

  figspath = get_ecoforecasts_path('figs');

  if ( ~exist('doPrint','var') || isempty(doPrint) )
    doPrint = false;
  end;

  if ( doPrint )
    print('-dtiff',fullfile(figspath,sprintf('%s-%s-seasonal-boxplot.tiff',stnm,fld)));
  end;




FROM read_oisst2.m:
                    % dat(dat<0) = nanmean(dat(dat>0));




FROM read_pathfinder_pentad.m:
FROM read_pathfinder_cortad.m:
  if ( ~exist('interpMethod','var') || isempty(interpMethod) )
    interpMethod = 'linear';
  end;
  if ( ischar(interpMethod) && interpMethod(1) == '*' )
    interpMethod = interpMethod(2:end);
  end;




FROM read_pathfinder_pentad.m:
  %DEBUG:  yrs = 1993:2009;
  %DEBUG:  yrs = 2009;




FROM read_pathfinder_pentad.m:
                    % stations(ix).(fld).data(dtix,1) = ...
                    %     interp2(lon{ix},lat{ix},dat,stations(ix).lon,stations(ix).lat,interpMethod);





                    % dat(dat<0) = nanmean(dat(dat>0));


            % dat(dat<0) = nanmean(dat(dat>0));





  % STRREP calls handle special cases to avoid "drooping letter syndrome".
  % Underscore in MATLAB LaTex means "subscript", but not in Window Names!

  th = title(varargin{:});
  txt = get(th, 'String');
  txt = strrep(txt,'\_','_');
  txt = strrep(txt,'_','\_');
  set(th, 'String',txt);



  % Just in case this TITLE was preprocessed to avoid "drooping letter
  % syndrome".  The underscore in MATLAB LaTex text means "subscript".





  % STRREP calls handle special cases to avoid "drooping letter syndrome".
  % Underscore in MATLAB LaTex means "subscript", but not in Window Names!
  txt = strrep(txt,'\_','_');
  txt = strrep(txt,'_','\_');






%%%% ??? DEBUG
  lon = ( ([lonoff:(lonoff+lonlen-1)]) * dx ) - (dx/2);
  lat = ( ([latoff:(latoff+latlen-1)]) * dx ) - (dx/2) - 90;



%%%% ??? DEBUG
  lon = ( ([lonoff:(lonoff+lonlen-1)]) * dx );
  lat = ( ([latoff:(latoff+latlen-1)]) * dx ) - 90;




%   lon = (([lonoff:(lonoff+lonlen-1)] - 1) * dx) + (dx/2);
%   lat = (([latoff:(latoff+latlen-1)] - 1) * dx) - 90 + (dx/2);
%   % lon = (([lonoff:(lonoff+lonlen-1)] - 0) * dx) + (dx/2);
%   % lat = (([latoff:(latoff+latlen-1)] + 1) * dx) - 90 + (dx/2);

% FOR SMKF1, LONF1, MLRF1, FWYF1 - NO! This EXACTLY MATCHES netCDF version! 
  lon = (([lonoff:(lonoff+lonlen-1)] - 0) * dx) + (dx/2);
  lat = (([latoff:(latoff+latlen-1)] + 0) * dx) - 90 + (dx/2);

% % FOR DBJM1, LCIY2
%   lon = (([lonoff:(lonoff+lonlen-1)] - 0) * dx) + (dx/2);
%   lat = (([latoff:(latoff+latlen-1)] - 1) * dx) - 90 + (dx/2);






%%%% ??? HACK
      % needix(needix == ix) = [];
%%%% ??? HACK

    else
%%%% ??? HACK
warning('WHOOPS! No MAT file for "%s"?!',stns(ix).station_name);
%%%% ??? HACK





  yrs = 1994:get_year(now);




function stns = load_misst_region_stations(rgn_or_stns,region)
%function stns = load_misst_region_stations(rgn_or_stns,region)
%
% Calculate MISST indices and load MISST time series for multiple sites
%
% Last Saved Time-stamp: <Sun 2011-04-10 19:33:08 Eastern Daylight Time gramer>

  set_more off;
  %DEBUG:
  tic,

  datapath = get_ecoforecasts_path('data');

  region = '';
  if ( ischar(rgn_or_stns) )
    region = rgn_or_stns;
    stns = load_region_metadata(region);
  elseif ( isstruct(rgn_or_stns) )
    stns = rgn_or_stns;
    if ( isfield(stns,'misst_region') )
      region = stns(1).misst_region;
    end;
  else
    error('First arg must be MISST region name string or vector of STN structs!');
  end;

  region = lower(region);

  worldOnly = ( strcmpi(region,'world') );

  disp(worldOnly);

  needix = 1:numel(stns);
  for ix = 1:numel(stns)
    matfname = fullfile(datapath,[stns(ix).station_name '_misst_' region '.mat']);
    % If we already loaded data one MISST file at a time, don't do it again
    if ( exist(matfname, 'file') )
      disp(['Reloading ' matfname]);
      load(matfname,'station');
      stations(ix) = station;
      station = []; clear station
      needix(needix == ix) = [];
    else
      flds = {'station_name','name','fname','lon','lat','depth','misst_region','misst_lonix','misst_latix'};
      for fldix = 1:length(flds)
        fld = flds{fldix};
        if ( isfield(stns,fld) )
          stations(ix).(fld) = stns(ix).(fld);
        end;
      end;
    end;
  end;


  % Extract MISST SSTs as needed
  if ( ~isempty(needix) )

    disp('Extracting SST from original MISST binary files...');

    % Find lat/lon indices for each station in MISST dataset for REGION
    if ( worldOnly )
      lonfld = 'misst_lonix';
      latfld = 'misst_latix';
    else
      lonfld = ['misst_' region '_lonix'];
      latfld = ['misst_' region '_latix'];
    end;
    if ( ~isfield(stations,lonfld) )
      stations(1).(lonfld) = [];
    end;
    if ( ~isfield(stations,latfld) )
      stations(1).(latfld) = [];
    end;
    [lon,lat,ig,dx] = read_misst_region(region);
    stations(needix) = load_misst_cfg(stations(needix),region);

    % If region is not 'world', find corresponding indices for global dataset also
    if ( ~worldOnly )
      [wlon,wlat,ig,wdx] = read_misst_region('world');
      rgnlonix = interp1(wlon,1:length(wlon),lon,'nearest');
      rgnlatix = interp1(wlat,1:length(wlat),lat,'nearest');

      % stations(needix) = load_misst_cfg(stations(needix),'world');
      for ix=needix(:)'
        stations(ix).misst_lonix = interp1(wlon,1:length(wlon),lon(stations(ix).(lonfld)),'nearest');
        stations(ix).misst_latix = interp1(wlat,1:length(wlat),lat(stations(ix).(latfld)),'nearest');
      end;
    end;


    % Once per day, Julian Day 2002-185 to "the present"
    begyear = 2002;  begjday = 185;
    endyear = 2011;  endjday = 61;
    alldts = datenum(begyear,1,begjday):datenum(endyear,1,endjday);

    % Make subsets rectangular, to avoid dyslexic errors
    xrad = 4;
    yrad = 3;

    for ix = needix(:)'
      stations(ix).misst_sst.date = repmat(nan,[length(alldts) 1]);
      stations(ix).misst_sst.data = repmat(nan,[length(alldts) 1]);

      stations(ix).misst_sst_field.date = repmat(nan,[length(alldts) 1]);

      xix = min(stations(ix).(lonfld))-xrad:max(stations(ix).(lonfld))+xrad;
      yix = min(stations(ix).(latfld))-yrad:max(stations(ix).(latfld))+yrad;
      stations(ix).misst_sst_field.lon = lon(xix);
      stations(ix).misst_sst_field.lat = lat(xix);

      stations(ix).misst_sst_field.field = repmat(nan,[length(alldts),length(lat(yix)),length(lon(xix))]);
    end;


    warning('off','MISST:NoRawFile');

    yrs = get_year(alldts(1)):get_year(alldts(end));

    for yr = yrs(:)'

      %DEBUG:
      disp(yr);

      switch (yr),
       case begyear,
        jds = begjday:365;
       case endyear,
        jds = 1:endjday;
       otherwise,
        jds = 1:365;
        if ( mod(yr,4) == 0 )
          jds = 1:366;
        end;
      end;

      for jd = jds(:)'
        dt = datenum(yr,1,1) + jd - 1;
        %DEBUG:        disp(datestr(dt));

        dtix = find(alldts == dt);
        if ( isempty(dtix) )
          error('Found no matching date for %g ("%s")?!',dt,datestr(dt));
        end;

        [ig,ig,sst] = read_misst_region(region,yr,jd);
        if ( isempty(sst) && ~worldOnly )
          [ig,ig,wsst] = read_misst_region('world',yr,jd);
          % Subset so each site's LATFLD and LONFLD will still index properly!
          if ( ~isempty(wsst) )
            if ( yr < 2010 )
              warning('Using WORLD file "%s"',datestr(dt));
            end;
            sst = wsst(rgnlatix,rgnlonix);
          end;
        end;

        if ( isempty(sst) )
          warning('MISST:MissingDay','No MISST data for year %g day %g',yr,jd);

        else
          for ix = needix(:)'
            dat = sst(stations(ix).(latfld),stations(ix).(lonfld));

            stations(ix).misst_sst.date(dtix,1) = dt;
            stations(ix).misst_sst.data(dtix,1) = nanmedian(dat(:));

            xix = min(stations(ix).(lonfld))-xrad:max(stations(ix).(lonfld))+xrad;
            yix = min(stations(ix).(latfld))-yrad:max(stations(ix).(latfld))+yrad;
            stations(ix).misst_sst_field.date(dtix,1) = dt;
            stations(ix).misst_sst_field.field(dtix,:,:) = sst(yix,xix);
          end; %for ix = needix(:)'

        end; %if ( isempty(sst) ) else

      end; %for jd = jds(:)'

    end; %for yr = yrs(:)'

    warning('on','MISST:NoRawFile');


    for ix = needix(:)'
      % Basic QA - remove days with missing or bad files
      badix = find( ~isfinite(stations(ix).misst_sst.date) | ...
                    ~isfinite(stations(ix).misst_sst_field.date) );
      %             ~isfinite(stations(ix).misst_sst.data) | ...
      %             all(~isfinite(stations(ix).misst_sst_field.field(:))) | ...
      stations(ix).misst_sst.date(badix) = [];
      stations(ix).misst_sst.data(badix) = [];
      stations(ix).misst_sst_field.date(badix) = [];
      stations(ix).misst_sst_field.field(badix,:,:) = [];

      matfname = fullfile(datapath,[stations(ix).station_name '_misst_' region '.mat']);
      disp(['Saving ' matfname]);
      station = stations(ix);
      save(matfname,'station');
      station = []; clear station
    end; %for ix = needix(:)'

  end; %if ( ~isempty(needix) )


  flds = fieldnames(stations);
  for fldix = 1:length(flds)
    fld = flds{fldix};
    for ix = 1:length(stations)
      stns(ix).(fld) = stations(ix).(fld);
      stations(ix).(fld) = [];
    end;
  end;
  stations = []; clear stations

  %DEBUG:
  toc,
  set_more;

return;










From _cortad.m:
        for ix=needix(:)'
          stations(ix).(fld).date = dts(:);
          if ( ~isempty(strfind(fld,'_field')) )
            dat = cast(nc{var}(yix(ix)-yrad:yix(ix)+yrad,xix(ix)-xrad:xix(ix)+xrad,:),'double');
            stations(ix).(fld).field = permute(dat,[3 2 1]);
          else
            dat = cast(nc{var}(yix(ix)-2:yix(ix)+2,xix(ix)-3:xix(ix)+3,:),'double');
            % dat(dat<0) = nanmean(dat(dat>0));
            stations(ix).(fld).data = ...
                interp2(lon{ix},lat{ix},dat,stations(ix).lon,stations(ix).lat,interpMethod);
          end; %if ( ~isempty(strfind(fld,'_field')) )
        end; %for ix=needix(:)'






  warnstat = warning('query','BuildDerivedVar:EmptyField');

  warning(warnstat);




function stns = verify_variable_multi(stns,varname,stnix)
%function stns = verify_variable_multi(stns,varname,stnix)
%
% Calls build_combo_var and build_derived_var (q.v.)
%
% Last Saved Time-stamp: <Mon 2011-04-04 12:31:54  Lew.Gramer>

  warningOff = warning('query','BuildDerivedVar:EmptyField');

  if ( ~exist('stnix','var') )
    stnix = 1:numel(stns);
  elseif ( numel(stns(stnix)) < numel(stns) )
    if ( ~isfield(stns,varname) )
      stns.(varname) = [];
      warning('off','BuildDerivedVar:EmptyField');
      warningOff = true;
    end;
  end;

  for ix = stnix(:)'
    % First ensure overarching 'combo' variable is built, if any
    stns(ix) = build_combo_var(stns(ix), varname);

    % Then make sure all other 'derived' variables get built too
    stns(ix) = build_derived_var(stns(ix), varname);
  end;

  % (Note: Building a Combo Var automatically forces its second input
  % variable to be built (if it is not yet). Separate build_derived_var
  % call above ensures any vars derived FROM a combo var also get built.)

  if ( warningOff )
    warning('on','BuildDerivedVar:EmptyField');
  end;

return;






%%%% DEBUG ???
    endyear = 2003;  endjday = 61;




  if ( strcmpi(region,'world') )





% %%%% HACK??? Already did 1993:2004
% yrs = 2005:2009;


% %%%% HACK??? Already did 1993 and 1994
      % needix(needix == ix) = [];


% %%%% HACK??? Already did 1993 and 1994
% warning('STATION STARTING FROM SCRATCH! "%s"',stations(ix).station_name);




% %%%% HACK??? Already did 1993 and 1994
% % SAVE MAT FILES AT THE END OF EVERY YEAR!
% % This extremely slow, annoying workaround is necessitated by the
% % frequent failures and time outs while accessing remote HDF files.
% for ix=needix(:)'
%   matfname = fullfile(datapath,[lower(stations(ix).station_name) '_pathfinder_pentad.mat']);
%   disp(['Saving MAT file ' matfname]);
%   station = stations(ix);
%   save(matfname,'station');
%   station=[]; clear station;
% end; %for ix=needix(:)'









        if ( isempty(dtix) )
          error('Unmatched date "%s"!',datestr(dt));
        end;





        % lerr = lasterror;
        % msg = 'CAUGHT ERROR!';
        % if ( isfield(lerr,'identifier') )
        %   msg = [ msg ' ' lerr.identifier ];
        % end;
        % if ( isfield(lerr,'message') )
        %   msg = [ msg ' : ' lerr.message ];
        % end;
        % warning(msg);




  if ( ~exist('flds','var') || isempty(flds) )
    flds = lower(vars);
  end;



  if ( ~exist('vars','var') || isempty(var) )
    vars = { 'SST4_daynitavg', 'STDV4_daynitavg', 'COUNTS4_daynitavg', ...
             'SST4_nighttime', 'STDV4_nighttime', 'COUNTS4_nighttime', ...
             'SST4_daynitavg', ...
           };
  end;
  if ( ~exist('flds','var') || isempty(flds) )
    flds = { 'sst',            'sst_std',         'sst_n', ...
             'night_sst',      'night_sst_std',   'night_sst_n', ...
             'sst_field', ...
           };
  end;






FROM load_10col_data.m:
    % Finally, remove any obviously INVALID values from each time series
    sens = fieldnames(result);
    for isen = 1:length(sens)
        if ( ~iscell(result.(sens{isen}).data) )
            rng = valid_var_range(sens{isen});
            goodidx = find( rng(1) <= result.(sens{isen}).data & ...
                            result.(sens{isen}).data <= rng(2) );
% %%%% DEBUG???
% badix = setdiff(1:length(result.(sens{isen}).data),goodidx);
% if ( ~isempty(badix) ); fprintf(2,'Removing %d points from %s\n',length(badix),sens{isen}); end;
% %%%% DEBUG???
            result.(sens{isen}).date = result.(sens{isen}).date(goodidx);
            result.(sens{isen}).data = result.(sens{isen}).data(goodidx);
        end;
    end;










    dy = (1-titlebuf)/n;
    dx = 1/m;
    % Tom and Dick
    [mix,nix] = ind2sub([m,n],ix(:));
    if ( max(nix) > n )
      error('Ecoforecasts:subplot_tight:SubplotIndexTooLarge',...
            'Index exceeds number of subplots.');
    end;
    nm = length(unique(mix));
    nn = length(unique(nix));
    ax = subplot('position',...
                 [((min(mix)-1)*dx)+buf+buf,((n-max(nix))*dy)+buf+titlebuf-0.01,(dx*nm)-buf4,(dy*nn)-buf4]);
    if ( nargs > 3 )
      set(ax,args{:});
    end;









    ax = subplot('position',...
                 [((min(mix)-1)*dx)+buf+0.03,((n-max(nix))*dy)+buf+titlebuf-0.01,(dx*nm)-buf4,(dy*nn)-buf4]);







disp(['DATETICK* ' num2str(gca)]);
disp(['SET_DATETICK_CURSOR ' num2str(gca)]);
disp(['DATETICK3 ' num2str(gca)]);





  persistent dtfn;





disp(['DATETICK3 ' num2str(gca)]);
disp(['DATETICK3 ' num2str(gca)]);




function vals = interp_field(lats,lons,fld,sitelat,sitelon,dlat,dlon,skipNaNs)
%function vals = interp_field(lats,lons,fld,sitelat,sitelon,[dlat],[dlon],[skipNaNs])
%
% Return bilinear interpolation on a *time series field*: LATS and LONS are
% monotonic vectors of length M and N, resp.; FLD is a matrix of size DxNxM;
% SITELAT,SITELON are both Px1 vectors of site coordinates. *Optional* args DLAT,
% DLON are inter-grid spacing in latitude, longitude, and are calculated if
% needed. VALS is a DxP vector (time series) of values interpolated on FLD.
% If optional SKIPNANS==True, treat NaNs as missing values in interpolation.
%
% This function is offered as a kindness to fellow dyslexics: the orgy of row
% flipping, dimension permuting, and plaid'ing required by INTERP3, and the
% glacial slowness of interpreted looping on INTERP2, all make it seemingly
% impossible in MATLAB to do cleanly and quickly what this function does.
%
% Last Saved Time-stamp: <Wed 2011-03-23 09:36:38  lew.gramer>

  ulats = unique(lats(:));
  ulons = unique(lons(:));
  nlats = length(ulats);
  nlons = length(ulons);

  if ( ndims(fld) ~= 3 || size(fld,2) ~= nlats || size(fld,3) ~= nlons )
    error('LATS,LONS,FLD should be Nx1, Mx1, and DxNxM, resp.!');
  end;
  if ( ndims(sitelat) ~= 2 || length(sitelat) ~= length(sitelon) )
    error('SITELAT,SITELON should be numerical vectors of identical length!');
  end;
  if ( ~exist('dlat','var') || isempty(dlat) )
    dlat = mean(diff(ulats));
  end;
  if ( ~exist('dlon','var') || isempty(dlon) )
    dlon = mean(diff(ulons));
  end;
  if ( ~exist('skipNaNs','var') || isempty(skipNaNs) )
    skipNaNs = false;
  end;

  vals = repmat(nan,[size(fld,1) length(sitelat)]);


  % We loop on sites rather than dates - which is normally MUCH faster.
  % If however we are called with thousands of sites, this will bog down.
  for ix = 1:length(sitelat)

    lat = sitelat(ix);
    lon = sitelon(ix);

    if ( lat<min(ulats) || lat>max(ulats) || lon<min(ulons) || lon>max(ulons) )
      % This could be considered a caller ERROR, but be kind
      warning('interp_field:OutOfBounds',...
              'Site (%g;%g) outside bounding box (%g,%g;%g,%g)',...
              lat,lon,min(ulats),max(ulats),min(ulons),max(ulons));
      %%%%%%%%%%%
      % CONTINUE
      %%%%%%%%%%%
      continue;
    end;

    % If dataset is column- or (as is usual) row-inverted, undo that for us!
    if ( lats(1) > lats(2) )
      lats = lats(end:-1:1);
      fld = fld(:,end:-1:1,:);
    end;
    if ( lons(1) > lons(2) )
      lons = lons(end:-1:1);
      fld = fld(:,:,end:-1:1);
    end;

    % Find indices of the nearest dataset gridpoint for our site
    [ig,yix] = min(abs(lat-lats));
    [ig,xix] = min(abs(lon-lons));

    % Find interpolation scales and gridpoint vertices
    yerr = (lat-lats(yix))/dlat;
    y = yix + yerr;
    if ( yerr >= 0 )
      y1=yix;
      y2=yix+1;
    else
      y1=yix-1;
      y2=yix;
    end;

    xerr = (lon-lons(xix))/dlon;
    x = xix + xerr;
    if ( xerr >= 0 )
      x1=xix;
      x2=xix+1;
    else
      x1=xix-1;
      x2=xix;
    end;

    % Brute force corner- and edge-handling. NOTE: We will make no attempt to
    % skip NaNs when the caller requests a point on an edge or a vertex.
    if (y==1)
      if (x==1)		vals(:,ix) = fld(:,1,1);
      elseif (x==nlons)	vals(:,ix) = fld(:,1,end);
      else		vals(:,ix) = (fld(:,1,x1).*(x2-x)) + (fld(:,1,x2).*(x-x1));
      end;
    elseif (y==nlats)
      if (x==1)		vals(:,ix) = fld(:,end,1);
      elseif (x==nlons)	vals(:,ix) = fld(:,end,end);
      else		vals(:,ix) = (fld(:,end,x1).*(x2-x)) + (fld(:,end,x2).*(x-x1));
      end;
    elseif (x==1)
      vals(:,ix) = (fld(:,y1,1).*(y2-y)) + (fld(:,y2,1).*(y-y1));
    elseif (x==nlons)
      vals(:,ix) = (fld(:,y1,end).*(y2-y)) + (fld(:,y2,end).*(y-y1));
    else
      f11 = fld(:,y1,x1);
      f12 = fld(:,y2,x1);
      f21 = fld(:,y1,x2);
      f22 = fld(:,y2,x2);
      % If caller requested NaN-skipping, do the best we can
      if ( skipNaNs )
        finiteF = nanmean([f11,f12,f21,f22]')';
        f11(isnan(f11)) = finiteF(isnan(f11));
        f12(isnan(f12)) = finiteF(isnan(f12));
        f21(isnan(f21)) = finiteF(isnan(f21));
        f22(isnan(f22)) = finiteF(isnan(f22));
      end;
      vals(:,ix) = (f11.*(x2-x).*(y2-y)) + (f21.*(x-x1).*(y2-y)) + ...
          (f12.*(x2-x).*(y-y1)) + (f22.*(x-x1).*(y-y1));
    end;

  end; %for ix=1:length(sitelat)

return;










          % -80.097, ...          % FWYF1
          % -80.380, ...          % MLRF1
          % -80.7833, ...         % TNRF1
          % -81.1100, ...         % SMKF1

          % 25.590, ...           % FWYF1
          % 25.010, ...           % MLRF1
          % 24.750, ...           % TNRF1
          % 24.62666, ...         % SMKF1









    error('See HELP SUBPLOT and HELP SUBPLOT_TIGHT for available arguments');




    % ax = subplot('position',...
    %              [((m-mix)*dx)+buf+0.01,((n-nix)*dy)+buf+0.01,dx-buf2,dy-buf2],...
    %              varargin{:});



%%%% DEBUG ??? MAYBE NDBC HEIGHTS ARE ALL NORMALIZED? 2011 Mar 15
    d = [10,10,10,2.7,nan,nan,nan];






% Identical to SUBPLOT (v., only certain calling forms), except that the gap
% or buffer space between subplots is significantly reduced. The following
% three calling forms for SUBPLOT are *not* supported by SUBPLOT_TIGHT:
%     SUBPLOT(m,n,P), where P is a vector
%     SUBPLOT('position',[left bottom width height])
%     SUBPLOT(111)




%   if ( ischar(zfld) )
%     stn = verify_variable(stn, zfld);
%   end;




  % isobathColor = [.6 .5 .4];




function [isobath,transects,cs,ch] = map_freef(boundingbox,isodepth,offsetLabels,axesOpts)
%function [isobath,transects,cs,ch] = map_freef(boundingbox,isodepth,offsetLabels,axesOpts)
%
% Draw a map of the Florida Reef Tract coastline, together with one or more
% isobaths at depth(s) ISODEPTH (DEFAULT [-350m]), suitable for plotting
% current vectors or property contour fields. BOUNDINGBOX defines the corners
% of the map drawn. If OFFSETLABELS is given, non-empty, not false (DEFAULT
% False), isobath contour labels are offset from contour lines (v. CONTOURF).
% AXESOPTS (DEFAULT {'same','nohit'}) controls AXES where contours are drawn:
% 'top'/'bot[tom]'=new AXES on top/bottom, 'same'=draw contours in current
% AXES; 'hit'/'nohit' = AXES with contours is/is not selectable (the property
% 'HitTest' in AXES and ContourGroup objects is set to 'on' or 'off', resp.)
%
% Returns coordinates of inner (Florida near-shore) isobath, coordinates of
% periodic TRANSECT lines perpendicular to that isobath, and the matrix and
% CONTOURGROUP handle returned by the CONTOUR of isobaths. If ISODEPTH is all
% NaN values, no isobaths are plotted or returned.
%
% DEFAULT: BOUNDINGBOX = [-80.50 -79.10 +24.80 +25.90];
%
% Last Saved Time-stamp: <Wed 2011-03-16 13:06:49  lew.gramer>

  global sofla_coast;
  global freef_topo;

  if ( ~exist('boundingbox', 'var') || isempty(boundingbox) )
    boundingbox = [-80.50 -79.10 +24.80 +25.90];
  end;
  if ( ~exist('isodepth', 'var') || isempty(isodepth) )
    isodepth = -350;
  end;
  if ( ~isnumeric(isodepth) )
    if ( strcmpi(isodepth,'none') )
      isodepth = [nan nan];
    else
      error('Optional 2nd arg ISODEPTH should be numeric or the string "none"!');
    end;
  elseif ( isscalar(isodepth) && ~isnan(isodepth) )
    isodepth = [ isodepth isodepth ];
  end;
  if ( ~exist('axesOpts','var') || isempty(axesOpts) )
    axesOpts = {'same','nohit'};
  end;

  % AXES Options - see above. Allow caller to specify simple string
  if ( ischar(axesOpts) )
    % But don't be picky - let caller specify comma-separated string as well
    axesOpts = strread(axesOpts,'%s','delimiter',',; ')
  end;
  axesOpts = lower(axesOpts);


  isobath = [];
  transects = [];
  cs = [];
  ch = [];

  cax = gca;
  if ( ismember('same',axesOpts) )
    hold on;
    ax = cax;
    nohit = false;

  else  % 'top' or 'bottom'
    pos = get(cax,'Position');
    ax = axes('position',pos);
    set(ax,'Color','none');
    nohit = true;

    if ( ismember('bottom',axesOpts) || ismember('bot',axesOpts) )
      %%%% ??? What does this do exactly?!
    end;
  end;

  axes(ax);
  %DEBUG:  cax,ax,gca,

  if ( ismember('hit',axesOpts) )
    nohit = false;
    set(ax,'HitTest','on');
  elseif ( ismember('nohit',axesOpts) )
    nohit = true;
    set(ax,'HitTest','off');
  end;

  % Load and draw (low-) medium- (full-)resolution coastline

  if ( ~exist('sofla_coast', 'var') || isempty(sofla_coast) )
    disp('Reloading coastline');
    % sofla_coast = load('sofla_coast_low.dat');
    sofla_coast = load('sofla_coast_medium.dat');
    % sofla_coast = load('sofla_coast.dat');
  end;

  %line(sofla_coast(:,1), sofla_coast(:,2), ...
  %     'Color', [.4 .3 .2]);
  cobj = fill(sofla_coast(:,1), sofla_coast(:,2), [.4 .3 .2]);
  if ( ~nohit )
    set(cobj,'HitTest','on');
  else
    set(cobj,'HitTest','off');
  end;


  % Load and draw requested isobath from ETOPO1 1-arcmin global relief model:
  % http://www.ngdc.noaa.gov/mgg/global/relief/ 
  % Amante, C. and B. W. Eakins, ETOPO1 1 Arc-Minute Global Relief Model: Procedures, Data Sources and Analysis. NOAA Technical Memorandum NESDIS NGDC-24, 19 pp, March 2009.

  if ( ~exist('freef_topo', 'var') || isempty(freef_topo) )
    disp('Reloading freef_topo...');

    freef_topo = load('freef.topo.mat');

    % This mesh would be HUGE: We don't really seem to need it??
    % [freef_topo.lons,freef_topo.lats] = meshgrid(freef_topo.tlon,freef_topo.tlat);
  end;

  latix = find(boundingbox(3) <= freef_topo.tlat & freef_topo.tlat <= boundingbox(4));
  lonix = find(boundingbox(1) <= freef_topo.tlon & freef_topo.tlon <= boundingbox(2));
  lon = freef_topo.tlon(lonix);
  lat = freef_topo.tlat(latix);
  topo = freef_topo.topo(latix,lonix);

  if ( numel(topo) < 4 )
    warning('map_freef:NotEnoughTopo', ...
            'Too few topo points in bounding box: no isobaths rendered');
    map_freef_return(cax,ax,boundingbox);
    return;
  end;

  if ( all(isnan(isodepth)) )
    disp('Isobars not plotted per user request...');
    map_freef_return(cax,ax,boundingbox);
    return;
  end;


  [cs, ch] = contour(ax, lon, lat, topo, ...
                     isodepth, 'LineColor', [.6 .5 .4]);
  if ( ~nohit )
    set(ch,'HitTest','on');
  else
    set(ch,'HitTest','off');
  end;
  if ( isempty(cs) )
    warning('map_freef:NoContourLines', ...
            'No isobath contours found in bounding box!');
    map_freef_return(cax,ax,boundingbox);
    return;
  end;

  % Put "+" within contours, and place labels nearby instead...
  if ( exist('offsetLabels','var') && ~isempty(offsetLabels) && (offsetLabels ~= false) )
    clhs = clabel(cs, 'Color',[.6 .5 .4]);
  else
    clhs = clabel(cs, ch, 'Color',[.6 .5 .4]);
  end;
  for clh=clhs(:)'
    if ( ~nohit )
      set(clh,'HitTest','on');
    else
      set(clh,'HitTest','off');
    end;
  end;


  set(ax, 'xlim', [boundingbox(1) boundingbox(2)]);
  set(ax, 'ylim', [boundingbox(3) boundingbox(4)]);

  if ( nargout > 0 )
    % Find lat/lon positions (nearest Florida) for desired isobath
    isobath = cs(1:2, 2:(cs(2,1) + 1));

    if ( nargout > 1 )
      % Subset gradient vectors at each interior point of isobath
      %[topx, topy] = gradient(sofla_topo);
      %lnix = ismember(sofla_topo_lons, isobath(1,:));
      %ltix = ismember(sofla_topo_lats, isobath(2,:));
      %topx = topx(lnix & ltix);
      %topy = topy(lnix & ltix);
      %transects(1,:,:) = isobath(:,2:end-1) - [topx ; topy];
      %transects(2,:,:) = isobath(:,2:end-1) + [topx ; topy];

      % Calculate secant and normal slopes at each interior point
      sct = isobath(:,3:end) - isobath(:,1:end-2);
      nrm = [sct(2,:) ; -sct(1,:)];
      % Calculate straight transects orthogonal to our isobath
      transects = repmat(nan, [2 2 length(isobath(1,1:end-2))]);
      transects(1,:,:) = isobath(:,2:end-1) - (2.*nrm);
      transects(2,:,:) = isobath(:,2:end-1) + (2.*nrm);

      %scts = [];
      %line(scts(:,1,:), scts(:,2,:));
    end;

  end;

  map_freef_return(cax,ax,boundingbox);

return;



%%%%%%%%%%%%%%%%%%%
% PRIVATE FUNCTIONS
%%%%%%%%%%%%%%%%%%%

function map_freef_return(cax,ax,boundingbox)
  set(ax, 'XLim', [boundingbox(1) boundingbox(2)]);
  set(ax, 'YLim', [boundingbox(3) boundingbox(4)]);
  % linkaxes([ax,cax],'xy');
  hlink = linkprop([cax,ax],{'DataAspectRatio','Position','View','XLim','YLim',});
  set(ax,'UserData',hlink);
  axes(cax);
  %DEBUG:  cax,ax,gca,
return;














function [isobath,transects,cs,ch] = map_freef(boundingbox,isodepth,offsetLabels,axesOpts)
%function [isobath,transects,cs,ch] = map_freef(boundingbox,isodepth,offsetLabels,axesOpts)
%
% Draw a map of the Florida Reef Tract coastline, together with one or more
% isobaths at depth(s) ISODEPTH (DEFAULT [-350m]), suitable for plotting
% current vectors or property contour fields. BOUNDINGBOX defines the corners
% of the map drawn. If OFFSETLABELS is given, non-empty, not false (DEFAULT
% False), isobath contour labels are offset from contour lines (v. CONTOURF).
% AXESOPTS (DEFAULT {'same','nohit'}) controls AXES where contours are drawn:
% 'top'/'bot[tom]'=new AXES on top/bottom, 'same'=draw contours in current
% AXES; 'hit'/'nohit' = AXES with contours is/is not selectable (the property
% 'HitTest' in AXES and ContourGroup objects is set to 'on' or 'off', resp.)
%
% Returns coordinates of inner (Florida near-shore) isobath, coordinates of
% periodic TRANSECT lines perpendicular to that isobath, and the matrix and
% CONTOURGROUP handle returned by the CONTOUR of isobaths. If ISODEPTH is all
% NaN values, no isobaths are plotted or returned.
%
% DEFAULT: BOUNDINGBOX = [-80.50 -79.10 +24.80 +25.90];
%
% Last Saved Time-stamp: <Thu 2011-03-10 12:25:09  lew.gramer>

  global sofla_coast;
  global freef_topo;

  if ( ~exist('boundingbox', 'var') || isempty(boundingbox) )
    boundingbox = [-80.50 -79.10 +24.80 +25.90];
  end;
  if ( ~exist('isodepth', 'var') || isempty(isodepth) )
    isodepth = -350;
  end;
  if ( ~isnumeric(isodepth) )
    if ( strcmpi(isodepth,'none') )
      isodepth = [nan nan];
    else
      error('Optional 2nd arg ISODEPTH should be numeric or the string "none"!');
    end;
  elseif ( isscalar(isodepth) && ~isnan(isodepth) )
    isodepth = [ isodepth isodepth ];
  end;
  if ( ~exist('axesOpts','var') || isempty(axesOpts) )
    axesOpts = {'same','nohit'};
  end;

  % AXES Options - see above. Allow caller to specify simple string
  if ( ischar(axesOpts) )
    % But don't be picky - let caller specify comma-separated string as well
    axesOpts = strread(axesOpts,'%s','delimiter',',; ')
  end;
  axesOpts = lower(axesOpts);


  isobath = [];
  transects = [];
  cs = [];
  ch = [];

  cax = gca;
  if ( ismember('same',axesOpts) )
    hold on;
    ax = cax;
    allowHits = true;
  else  % 'top' or 'bottom'
    pos = get(cax,'Position');
    ax = axes('position',pos);
    set(ax,'Color','none');
    allowHits = false;

    if ( ismember('bottom',axesOpts) || ismember('bot',axesOpts) )
      %%%% ??? What does this do exactly?!
    end;
  end;

  axes(ax);
  %DEBUG:  cax,ax,gca,

  if ( ismember('hit',axesOpts) )
    allowHits = true;
    set(ax,'HitTest','on');
  elseif ( ismember('nohit',axesOpts) )
    allowHits = false;
    set(ax,'HitTest','off');
  end;

  % Load and draw (low-) medium- (full-)resolution coastline

  if ( ~exist('sofla_coast', 'var') || isempty(sofla_coast) )
    disp('Reloading coastline');
    % sofla_coast = load('sofla_coast_low.dat');
    sofla_coast = load('sofla_coast_medium.dat');
    % sofla_coast = load('sofla_coast.dat');
  end;

  %line(sofla_coast(:,1), sofla_coast(:,2), ...
  %     'Color', [.4 .3 .2]);
  cobj = fill(sofla_coast(:,1), sofla_coast(:,2), [.4 .3 .2]);
  if ( allowHits )
    set(cobj,'HitTest','on');
  else
    set(cobj,'HitTest','off');
  end;


  % Load and draw requested isobath from ETOPO1 1-arcmin global relief model:
  % http://www.ngdc.noaa.gov/mgg/global/relief/ 
  % Amante, C. and B. W. Eakins, ETOPO1 1 Arc-Minute Global Relief Model: Procedures, Data Sources and Analysis. NOAA Technical Memorandum NESDIS NGDC-24, 19 pp, March 2009.

  if ( ~exist('freef_topo', 'var') || isempty(freef_topo) )
    disp('Reloading freef_topo...');

    freef_topo = load('freef.topo.mat');

    % This mesh would be HUGE: We don't really seem to need it??
    % [freef_topo.lons,freef_topo.lats] = meshgrid(freef_topo.tlon,freef_topo.tlat);
  end;

  latix = find(boundingbox(3) <= freef_topo.tlat & freef_topo.tlat <= boundingbox(4));
  lonix = find(boundingbox(1) <= freef_topo.tlon & freef_topo.tlon <= boundingbox(2));
  lon = freef_topo.tlon(lonix);
  lat = freef_topo.tlat(latix);
  topo = freef_topo.topo(latix,lonix);

  if ( numel(topo) < 4 )
    warning('map_freef:NotEnoughTopo', ...
            'Too few topo points in bounding box: no isobaths rendered');
    map_freef_return(cax,ax,boundingbox);
    return;
  end;

  if ( all(isnan(isodepth)) )
    disp('Isobars not plotted per user request...');
    map_freef_return(cax,ax,boundingbox);
    return;
  end;


  [cs, ch] = contour(ax, lon, lat, topo, ...
                     isodepth, 'LineColor', [.6 .5 .4]);
  if ( allowHits )
    set(ch,'HitTest','on');
  else
    set(ch,'HitTest','off');
  end;
  if ( isempty(cs) )
    warning('map_freef:NoContourLines', ...
            'No isobath contours found in bounding box!');
    map_freef_return(cax,ax,boundingbox);
    return;
  end;

  % Put "+" within contours, and place labels nearby instead...
  if ( exist('offsetLabels','var') && ~isempty(offsetLabels) && (offsetLabels ~= false) )
    clhs = clabel(cs, 'Color',[.6 .5 .4]);
  else
    clhs = clabel(cs, ch, 'Color',[.6 .5 .4]);
  end;
  for clh=clhs(:)'
    if ( allowHits )
      set(clh,'HitTest','on');
    else
      set(clh,'HitTest','off');
    end;
  end;


  set(ax, 'xlim', [boundingbox(1) boundingbox(2)]);
  set(ax, 'ylim', [boundingbox(3) boundingbox(4)]);

  if ( nargout > 0 )
    % Find lat/lon positions (nearest Florida) for desired isobath
    isobath = cs(1:2, 2:(cs(2,1) + 1));

    if ( nargout > 1 )
      % Subset gradient vectors at each interior point of isobath
      %[topx, topy] = gradient(sofla_topo);
      %lnix = ismember(sofla_topo_lons, isobath(1,:));
      %ltix = ismember(sofla_topo_lats, isobath(2,:));
      %topx = topx(lnix & ltix);
      %topy = topy(lnix & ltix);
      %transects(1,:,:) = isobath(:,2:end-1) - [topx ; topy];
      %transects(2,:,:) = isobath(:,2:end-1) + [topx ; topy];

      % Calculate secant and normal slopes at each interior point
      sct = isobath(:,3:end) - isobath(:,1:end-2);
      nrm = [sct(2,:) ; -sct(1,:)];
      % Calculate straight transects orthogonal to our isobath
      transects = repmat(nan, [2 2 length(isobath(1,1:end-2))]);
      transects(1,:,:) = isobath(:,2:end-1) - (2.*nrm);
      transects(2,:,:) = isobath(:,2:end-1) + (2.*nrm);

      %scts = [];
      %line(scts(:,1,:), scts(:,2,:));
    end;

  end;

  map_freef_return(cax,ax,boundingbox);

return;



%%%%%%%%%%%%%%%%%%%
% PRIVATE FUNCTIONS
%%%%%%%%%%%%%%%%%%%

function map_freef_return(cax,ax,boundingbox)
  set(ax, 'XLim', [boundingbox(1) boundingbox(2)]);
  set(ax, 'YLim', [boundingbox(3) boundingbox(4)]);
  % linkaxes([ax,cax],'xy');
  hlink = linkprop([cax,ax],{'DataAspectRatio','Position','View','XLim','YLim',});
  set(ax,'UserData',hlink);
  axes(cax);
  %DEBUG:  cax,ax,gca,
return;
























function [isobath,transects,cs,ch] = map_freef(boundingbox,isodepth,offsetLabels,axesOpts)
%function [isobath,transects,cs,ch] = map_freef(boundingbox,isodepth,offsetLabels,axesOpts)
%
% Draw a map of the Florida Reef Tract coastline, together with one or more
% isobaths at depth(s) ISODEPTH (DEFAULT [-350m]), suitable for plotting
% current vectors or property contour fields. BOUNDINGBOX defines the corners
% of the map drawn. If OFFSETLABELS is given, non-empty, not false (DEFAULT
% False), isobath contour labels are offset from contour lines (v. CONTOURF).
% AXESOPTS (DEFAULT {'same','nohit'}) controls AXES where contours are drawn:
% 'top'/'bot[tom]'=new AXES on top/bottom, 'same'=draw contours in current
% AXES; 'hit'/'nohit' = AXES with contours is/is not selectable (the property
% 'HitTest' in AXES and ContourGroup objects is set to 'on' or 'off', resp.)
%
% Returns coordinates of inner (Florida near-shore) isobath, coordinates of
% periodic TRANSECT lines perpendicular to that isobath, and the matrix and
% CONTOURGROUP handle returned by the CONTOUR of isobaths. If ISODEPTH is all
% NaN values, no isobaths are plotted or returned.
%
% DEFAULT: BOUNDINGBOX = [-80.50 -79.10 +24.80 +25.90];
%
% Last Saved Time-stamp: <Thu 2011-03-10 12:13:25  lew.gramer>

  global sofla_coast;
  global freef_topo;

  if ( ~exist('boundingbox', 'var') || isempty(boundingbox) )
    boundingbox = [-80.50 -79.10 +24.80 +25.90];
  end;
  if ( ~exist('isodepth', 'var') || isempty(isodepth) )
    isodepth = -350;
  end;
  if ( ~isnumeric(isodepth) )
    if ( strcmpi(isodepth,'none') )
      isodepth = [nan nan];
    else
      error('Optional 2nd arg ISODEPTH should be numeric or the string "none"!');
    end;
  elseif ( isscalar(isodepth) && ~isnan(isodepth) )
    isodepth = [ isodepth isodepth ];
  end;
  if ( ~exist('axesOpts','var') || isempty(axesOpts) )
    axesOpts = {'same','nohit'};
  end;

  % AXES Options - see above. Allow caller to specify simple string
  if ( ischar(axesOpts) )
    % But don't be picky - let caller specify comma-separated string as well
    axesOpts = strread(axesOpts,'%s','delimiter',',; ')
  end;
  axesOpts = lower(axesOpts);


  isobath = [];
  transects = [];
  cs = [];
  ch = [];

  cax = gca;
  if ( ismember('same',axesOpts) )
    hold on;
    ax = cax;
    allowHits = true;
  else  % 'top' or 'bottom'
    pos = get(cax,'Position');
    ax = axes('position',pos);
    set(ax,'Color','none');
    allowHits = false;

    if ( ismember('bottom',axesOpts) || ismember('bot',axesOpts) )
      %%%% ??? What does this do exactly?!
    end;
  end;

  axes(ax);
  %DEBUG:  cax,ax,gca,

  if ( ismember('hit',axesOpts) )
    allowHits = true;
    set(ax,'HitTest','on');
  elseif ( ismember('nohit',axesOpts) )
    allowHits = false;
    set(ax,'HitTest','off');
  end;

  % Load and draw (low-) medium- (full-)resolution coastline

  if ( ~exist('sofla_coast', 'var') || isempty(sofla_coast) )
    disp('Reloading coastline');
    % sofla_coast = load('sofla_coast_low.dat');
    sofla_coast = load('sofla_coast_medium.dat');
    % sofla_coast = load('sofla_coast.dat');
  end;

  %line(sofla_coast(:,1), sofla_coast(:,2), ...
  %     'Color', [.4 .3 .2]);
  cobj = fill(sofla_coast(:,1), sofla_coast(:,2), [.4 .3 .2]);
  if ( allowHits )
    set(cobj,'HitTest','on');
  else
    set(cobj,'HitTest','off');
  end;


  % Load and draw requested isobath from ETOPO1 1-arcmin global relief model:
  % http://www.ngdc.noaa.gov/mgg/global/relief/ 
  % Amante, C. and B. W. Eakins, ETOPO1 1 Arc-Minute Global Relief Model: Procedures, Data Sources and Analysis. NOAA Technical Memorandum NESDIS NGDC-24, 19 pp, March 2009.

  if ( ~exist('freef_topo', 'var') || isempty(freef_topo) )
    disp('Reloading freef_topo...');

    freef_topo = load('freef.topo.mat');

    % This mesh would be HUGE: We don't really seem to need it??
    % [freef_topo.lons,freef_topo.lats] = meshgrid(freef_topo.tlon,freef_topo.tlat);
  end;

  latix = find(boundingbox(3) <= freef_topo.tlat & freef_topo.tlat <= boundingbox(4));
  lonix = find(boundingbox(1) <= freef_topo.tlon & freef_topo.tlon <= boundingbox(2));
  lon = freef_topo.tlon(lonix);
  lat = freef_topo.tlat(latix);
  topo = freef_topo.topo(latix,lonix);

  if ( numel(topo) < 4 )
    warning('map_freef:NotEnoughTopo', ...
            'Too few topo points in bounding box: no isobaths rendered');
    map_freef_return(cax,ax,boundingbox);
    return;
  end;

  if ( all(isnan(isodepth)) )
    disp('Isobars not plotted per user request...');
    map_freef_return(cax,ax,boundingbox);
    return;
  end;


  [cs, ch] = contour(ax, lon, lat, topo, ...
                     isodepth, 'LineColor', [.6 .5 .4]);
  if ( allowHits )
    set(ch,'HitTest','on');
  else
    set(ch,'HitTest','off');
  end;
  if ( isempty(cs) )
    warning('map_freef:NoContourLines', ...
            'No isobath contours found in bounding box!');
    map_freef_return(cax,ax,boundingbox);
    return;
  end;

  % Put "+" within contours, and place labels nearby instead...
  if ( exist('offsetLabels','var') && ~isempty(offsetLabels) && (offsetLabels ~= false) )
    clhs = clabel(cs, 'Color',[.6 .5 .4]);
  else
    clhs = clabel(cs, ch, 'Color',[.6 .5 .4]);
  end;
  for clh=clhs(:)'
    if ( allowHits )
      set(clh,'HitTest','on');
    else
      set(clh,'HitTest','off');
    end;
  end;


  set(ax, 'xlim', [boundingbox(1) boundingbox(2)]);
  set(ax, 'ylim', [boundingbox(3) boundingbox(4)]);

  if ( nargout > 0 )
    % Find lat/lon positions (nearest Florida) for desired isobath
    isobath = cs(1:2, 2:(cs(2,1) + 1));

    if ( nargout > 1 )
      % Subset gradient vectors at each interior point of isobath
      %[topx, topy] = gradient(sofla_topo);
      %lnix = ismember(sofla_topo_lons, isobath(1,:));
      %ltix = ismember(sofla_topo_lats, isobath(2,:));
      %topx = topx(lnix & ltix);
      %topy = topy(lnix & ltix);
      %transects(1,:,:) = isobath(:,2:end-1) - [topx ; topy];
      %transects(2,:,:) = isobath(:,2:end-1) + [topx ; topy];

      % Calculate secant and normal slopes at each interior point
      sct = isobath(:,3:end) - isobath(:,1:end-2);
      nrm = [sct(2,:) ; -sct(1,:)];
      % Calculate straight transects orthogonal to our isobath
      transects = repmat(nan, [2 2 length(isobath(1,1:end-2))]);
      transects(1,:,:) = isobath(:,2:end-1) - (2.*nrm);
      transects(2,:,:) = isobath(:,2:end-1) + (2.*nrm);

      %scts = [];
      %line(scts(:,1,:), scts(:,2,:));
    end;

  end;

  map_freef_return(cax,ax,boundingbox);

return;



%%%%%%%%%%%%%%%%%%%
% PRIVATE FUNCTIONS
%%%%%%%%%%%%%%%%%%%

function map_freef_return(cax,ax,boundingbox)
  set(ax, 'XLim', [boundingbox(1) boundingbox(2)]);
  set(ax, 'YLim', [boundingbox(3) boundingbox(4)]);
  % linkaxes([ax,cax],'xy');
  hlink = linkprop([cax,ax],{'DataAspectRatio','Position','View','XLim','YLim',});
  set(ax,'UserData',hlink);
  axes(cax);
  %DEBUG:  cax,ax,gca,
return;


















    % axesOpts = {'top','nohit'};





    linkaxes([ax,cax],'xy');
    % linkaxes([ax,cax]);






  if ( all(isnan(isodepth)) )
    set(ax, 'xlim', [boundingbox(1) boundingbox(2)]);
    set(ax, 'ylim', [boundingbox(3) boundingbox(4)]);

    disp('Isobars not plotted per user request...');
    return;
  end;








  if ( exist('fh','var') && ~isempty(fh) )
    % Useful to pass in fh for, e.g., use in subplots
    figure(fh);
    hold on;
  else
    fh = figure;
    set(fh, 'units','normalized', 'outerposition',[0 0 1 1]);
    hold on;
  end;







  %% psst(mask) = nan;

  % osst = interp2(psst,[1:subs:size(psst,2)]',[1:subs:size(psst,1)],'*cubic');
  osst = interp2(psst,[1:subs:size(psst,2)]',[1:subs:size(psst,1)],'*linear');

  % Propagate land mask to oversampled grid in some reasonable way
  omask = logical(interp2(mask,[1:subs:size(mask,2)]',[1:subs:size(mask,1)],'*linear'));
  osst(omask) = nan;





  % [nanix,nanjx]=find(~isfinite(psst) & ~mask);




  % omask = logical(interp2(mask,[1:subs:size(mask,2)]',[1:subs:size(mask,1)],'*linear'));
  % osst(omask) = nan;










      warning('off','MISST:NoRawFile');
      warning('on','MISST:NoRawFile');





function stns = get_nmc_winds(stns)
%function stns = get_nmc_winds(stns)
%
% Subset U and V wind components in M/S (and calculate speed in KTS and
% direction degT) from NMC Reanalysis, for all sites in struct(s) STNS.
%
% Last Saved Time-stamp: <Sun 2011-03-06 13:42:02  Lew.Gramer>

  %DEBUG:
  tic,

  datapath = get_ecoforecasts_path('data');
  nmcpath = fullfile(datapath,'nmc');

  vars = {'uwnd','vwnd'};
  flds = {'nmc_wind_u','nmc_wind_v'};

  for ix=1:length(stns)
    for varix = 1:length(vars)
      fld = ['native_' flds{varix}];
      stns(ix).(fld) = struct('date',[],'data',[]);
    end;
  end;

  for yr=2002:2009
    for varix = 1:length(vars)
      var = vars{varix};
      fld = ['native_' flds{varix}];
      %DEBUG:
      disp({yr,var});
      fname = fullfile(nmcpath,[var '.10m.gauss.' num2str(yr) '.nc']);
      nc = mDataset(fname);
      if ( isempty(nc) )
        warning('Skipping bad or missing file "%s"',fname);
        break;
      end;
      try
        lon=nc{'lon'}(:,:);
        lat=nc{'lat'}(:,:);
        hrs = nc{'time'}(:);
        dts = datenum(1,1,1,0,0,0) + (hrs/24);
        ndts = length(dts);
        for ix=1:length(stns)
          lonix=round(interp1(lon,1:length(lon),360+stns(ix).lon));
          latix=round(interp1(lat,1:length(lat),stns(ix).lat));
          mylon = nc{var}(:,latix-2:latix+2,lonix-2:lonix+2);
          mydat = nc{var}(:,latix-2:latix+2,lonix-2:lonix+2);
          stns(ix).(fld).date(end+1:end+ndts) = dts;
          stns(ix).(fld).data(end+1:end+ndts) = interp2(mylat,mylon,;
        end;
      catch
        close(nc); clear nc
        rethrow(lasterror);
      end;
      close(nc); clear nc
    end;
  end;

  for varix = 1:length(vars)
    rawfld = ['native_' flds{varix}];
    fld = flds{varix};
    for ix=1:length(stns)
      stns(ix).(fld).date = ...
          stns(ix).(rawfld).date(1):(1/24):stns(ix).(rawfld).date(end);
      stns(ix).(fld).data = ...
          interp1(stns(ix).(rawfld).date,stns(ix).(rawfld).data,stns(ix).(fld).date,'spline');
    end;
  end;

  for ix=1:length(stns)
    stns(ix).native_nmc_wind_speed.date = stns(ix).native_nmc_wind_u.date;
    stns(ix).native_nmc_wind_speed.data = mps2kts(uv_to_spd(stns(ix).native_nmc_wind_u.data,stns(ix).native_nmc_wind_v.data));
    stns(ix).native_nmc_wind_dir.date = stns(ix).native_nmc_wind_u.date;
    stns(ix).native_nmc_wind_dir.data = uv_to_dir(stns(ix).native_nmc_wind_u.data,stns(ix).native_nmc_wind_v.data);

    stns(ix).nmc_wind_speed.date = stns(ix).nmc_wind_u.date;
    stns(ix).nmc_wind_speed.data = mps2kts(uv_to_spd(stns(ix).nmc_wind_u.data,stns(ix).nmc_wind_v.data));
    stns(ix).nmc_wind_dir.date = stns(ix).nmc_wind_u.date;
    stns(ix).nmc_wind_dir.data = uv_to_dir(stns(ix).nmc_wind_u.data,stns(ix).nmc_wind_v.data);
  end;


  %DEBUG:
  toc,

return;











for ix = [1:4 5 9 6:8]

% [ig,sortix] = sort({sites.name});
[ig,sortix] = sort([sites.lat]);





          ...
          'FFS-32'    'FFS-R30'    'PHR-31'    'PHR-32'    'LIS-R7'    'KUR-17' ...
          'KUR-R36'    'KUR-14'    'LIS-10' ...
          ...
          -166.2306   -166.2059  -175.9735   -175.9396   -173.9706   -178.3662 ...
          -178.3714   -178.3283   -173.9223 ...
          ...
          23.8061    23.8587    27.7758    27.7728    25.9538    28.4321 ...
          28.4203    28.4537    25.9410 ...
          ...
          9    1    9    8    11    3    2    2    16 ...






  % Use default color order if user did not give a plot spec (see PLOT)
  for tix = 1:nts
    if ( isempty(pargs{tix}) )
      pargs{tix} = {pargs{tix}{:} 'Color',co(tix,:)};











    ylbl{end+1} = [sites(ix).name ' (#' num2str(ix) ' Z=' sites(ix).max_depth 'm)'];




sites(1).latix = sites(1).latix-1;
sites(2).lonix = sites(2).lonix:sites(2).lonix+1; sites(2).latix = sites(2).latix:sites(2).latix+1;
sites(5).latix = sites(5).latix-1:sites(5).latix;
sites(6).latix = sites(6).latix:sites(6).latix+1;
sites(7).latix = sites(7).latix-1:sites(7).latix;
sites(8).lonix = sites(8).lonix-1:sites(8).lonix; sites(8).latix = sites(8).latix:sites(8).latix+1;
sites(9).latix = sites(9).latix-1:sites(9).latix;





sites(3).latix = sites(3).latix-1:sites(3).latix;

sites(4).lonix = sites(4).lonix-1:sites(4).lonix;






  if ( ~exist('region','var') || isempty(region) )
    region = 'florida';
  end;






  fname = fullfile(misstpath,sprintf('mw_ir.papa.fusion.%04d.%03d.v03',yr,jd);



          % Algorithm sometimes returns complex numbers??
          result(ix,1:length(res)) = real( res );







function [new_indices, new_data] = gap_expand(indices, data)
%function [new_indices, new_data] = gap_expand(indices, data)
%
% Expand a time series 'data' (e.g., sea temperatures), presumed to be indexed
% by another time series 'indices' (e.g., measurement dates), so that it has
% no indicial gaps: any gaps from the original data are filled with NaNs. The
% 'new_indices' contains a continuous series of equally-spaced datenums from
% the beginning to the end of the vector 'indices' - including datenums for
% all times where a NaN was added to 'new_data'.
%
% Last SavedXXXX: <XXXXSat 2009-02-28 21:32:23 Eastern Standard Time gramerXXXX>

  % Data may contain nans, but dates ('indices') may not!
  if ( any(~isfinite(indices)) )
    error('Indices must not contain NaN or Inf! Please review data...');
  end;
  if ( min(diff(indices)) < 0 )
      error('Indices must increase monotonically! Please review data...');
  end;

  % NOTE: 'indices' may come in as a float vector beginning with any value
  idx = round((indices - indices(1)) ./ min(diff(indices))) + 1;

  % Create a new TS big enough to hold all the filled gaps
  new_n = floor(idx(end) - idx(1) + 1);
  new_indices = linspace(indices(1), indices(end), new_n);

  new_data = repmat(nan, [size(data, 1) new_n]);
  new_data(idx) = data;

return;









  switch ( lower(method) ),
   case {'ssr','rss'},
    result = sum( (resid.^2) );
   case {'l2','l-2','rms','rmsd','rmse'},
    result = sqrt( sum(resid.^2) ) / sqrt(n);
   case {'nrmsd','nrmse'},
    result = sqrt( sum(resid.^2) ) / (max(resid) - min(resid));
   case {'cvrmsd','cvrmse'},
    result = sqrt( sum(resid.^2) ) /  mean(resid);
   case {'l1','l-1','lad','rmsd','rmse'},
    result = sum( abs(resid) );
   otherwise,
    error('Unknown error evaluation method "%s"',method);
  end;





   case {'rsd'},
    result.data = sqrt( ((ts1.data(ix1) - ts2.data(ix2)).^2) );




      % Use default color order if user did not give a plot spec (see PLOT)
      if ( length(pargs{nts}) == 1 )
        if ( ~ischar(arg) || length(arg) > 4 || ...
             isempty(regexp(arg,'^[.ox+*sdv^<>ph:-]*[bgrcmykw][.ox+*sdv^<>ph:-]*$')) )
          pargs{nts} = {pargs{nts}{:} 'Color',co(mod(nts-1,ncos)+1,:)};
        end;
      end;





  npp = get(ax,'NextPlot');
  % Restore previous HOLD state
  set(ax,'NextPlot',npp);







        if ( ~ischar(arg) || length(arg) > 4 || ~ismember(arg(1),'bgrcmykw') )




      havecolor = false;





    if ( ~isempty(pargs{tix}) )
      lh(tix) = plot(ax,ts{tix}.date,ts{tix}.data,pargs{tix}{:});
    else
      lh(tix) = plot(ax,ts{tix}.date,ts{tix}.data);
    end;






  for tix = 1:nts
    % if ( isempty(pargs{tix}) )
      pargs{tix} = {pargs{tix}{:} 'Color',co(tix,:)};
    % end;
  end;







    argix = argix + 1;

    if ( argix <=nargin && ischar(varargin{argix}) )
      pspec(end+1) = varargin{argix};
    end;





    dat = []; clear dat;
    dat = repmat(nan,[length(flds) 1]);




  % Convert station air pressure to sea-level
  p = p ./ exp( -pz ./ ( (a + 273.15) .* 29.263 ) );

  % Convert station air pressure to sea-level
  p = barom_to_surf(p,a,pz);






  sp = p ./ exp( -pz ./ ( (a + 273.15) .* 29.263 ) );





      [wu,wv] = spddir_to_uv(stn.(wfld).data(wsix),stn.(wdfld).data(wdix));
      wu = interp1_stn_fld(stn, wpfld, dts);
      ou = interp1_stn_fld(stn, whfld, dts);








  % Perform Fairall et al. (TOGA-COARE 2.0) calculation
  if ( isempty(dsfld) )
    disp('HFBULKTC 2.0 (no cool-skin)');
    disp(length(w));
    result = hfbulktc(w,wz,a,az,q,az,p,t);

    tau = result(:,4);
    shf = result(:,1);
    %lhf = result(:,2);
    % Be sure to include the 'Webb' correction!
    % (Comment from HFBULKTC: eqn. 22, Fairall et al. (1996), JGR, 101, p3751.)
    lhf = result(:,2) + result(:,3);

  % Perform Fairall et al. (TOGA-COARE 2.6) calculation - NO RAIN
  elseif ( isempty(prfld) )
    disp('HFBULKTC 2.6 (cool-skin, no rain)');
    disp(length(w));
    result = hfbulktc(w,wz,a,az,q,az,p,t,36,dlrf,dsrf,srf);

    tau = result(:,4);
    shf = result(:,1);
    %lhf = result(:,2);
    % Be sure to include the 'Webb' correction!
    % (Comment from HFBULKTC: eqn. 22, Fairall et al. (1996), JGR, 101, p3751.)
    lhf = result(:,2) + result(:,3);

  else

    % Make sure Fairall et al. code is in our path
    % (And remove later, if it wasn't originally... MATLAB namespace pollution)
    if ( isempty(strfind(path, 'fairall')) )
      added_fairall_path = true;
      % FAIRALLHOME = 'c:/Documents and Settings/gramer/My Documents/MATLAB/fairall';
      FAIRALLHOME = get_ecoforecasts_path('../fairall');
      addpath(FAIRALLHOME);
      rehash
    end;

    % Air specific humidity [kg/kg]
    sa = relhumid_to_spechumid(a,q);
    sa = sa .* 1e3; % Fairall expects [g/kg]
    % Saturated ("sea-surface") specific humidity [kg/kg]
    ss = 0.98 .* relhumid_to_spechumid(t,100);
    ss = ss .* 1e3; % Fairall expects [g/kg]


    % (Atmospheric) Planetary Boundary Layer height == inversion height [m]
    pblz = 600;
    %%%% Is 600m a good mean value for the FLORIDA KEYS??
    % Hsu (1979) shows 600m for the TX Gulf coast, while Kara et al (1998)
    % show an *upper bound* 500m for nocturnal conditions over Tallahassee:
    % they show profiles indicating 125m may be a better mean for there.

    % Change this function soon, to include an optional "PBLZFLD" arg, to
    % allow this value to be input, e.g., from ERA Iterim reanalysis...


    % Project ocean currents (if we have them) onto wind direction (ditto)
    if ( isempty(wdfld) || strcmpi(wdfld,'default') || ...
         isempty(oufld) || strcmpi(oufld,'default') || ...
         isempty(ovfld) || strcmpi(ovfld,'default') )
      % Assume surface current projects onto 2% of wind velocity [e.g., Ardhuin et al. 2009]
      warning('Estimating projected ocean currents from wind speed');
      ou = 0.020 .* w;
    else
      [wsix,wdix,ouix,ovix] = ...
          intersect_all_dates([],stn.(wfld).date,stn.(wdfld).date,stn.(oufld).date,stn.(ovfld).date);
      wd = stn.(wdfld).data(wdix);
      ouraw = cosd(wd)
      [wu,wv] = spddir_to_uv(stn.(wfld).data(wsix),stn.(wdfld).data(wdix));
      wu = interp1_stn_fld(stn, wpfld, dts);
      ou = interp1_stn_fld(stn, whfld, dts);
    end;

    if ( isempty(wpfld) )

      % Perform Fairall et al. (TOGA-COARE 2.6) calculation
      disp('COR26 (cool-skin, rain)');
      disp(length(w));
      nby10 = floor(length(w)/10);
      for ix = 1:length(w)
        res = cor30a([w(ix),ou(ix),t(ix),a(ix),ss(ix),sa(ix),dsrf(ix),dlrf(ix),prcp(ix),pblz,p(ix),wz,az,az,stn.lat,1,0,0,0]);
        result(ix,1:length(res)) = res;
        if (mod(ix,nby10)==0); disp(ix); end;
      end;

    else

      if ( strcmpi(wpfld,'default') )
        % Calculate wind-wave period and height from wind speed
        warning('Estimating wave period and height from wind speed');
        [wvper,wvhgt] = wind_to_wave(w);
      else
        wvper = interp1_stn_fld(stn, wpfld, dts);
        wvhgt = interp1_stn_fld(stn, whfld, dts);
      end;

      if ( ~doWarm )

        % Perform Fairall et al. (TOGA-COARE 3.0a) calculation
        disp('COR30A (no warm-layer)');
        disp(length(w));
        nby10 = floor(length(w)/10);
        for ix = 1:length(w)
          res = cor30a([w(ix),ou(ix),t(ix),a(ix),ss(ix),sa(ix),dsrf(ix),dlrf(ix),prcp(ix),pblz,p(ix),wz,az,az,stn.lat,1,1,wvper(ix),wvhgt(ix)]);
          result(ix,1:length(res)) = res;
          if (mod(ix,nby10)==1); disp(ix); end;
        end;

      else

        % Perform Fairall et al. (TOGA-COARE 3.0a + WARM CORRECTION) calculation
        disp('COR30A (warm-layer)');
        result = cor30a_warm(dts,w,ou,t,a,ss,sa,dsrf,dlrf,prcp,pblz,p,stn.lon,stn.lat,wz,az,az,stz,1,1,1,wvper,wvhgt);

      end; %if ~doWarm else

    end; %if isempty(wpfld) else

    tau = result(:,3);
    shf = -result(:,1);
    lhf = -result(:,2);
    rhf = -result(:,14);

    % Diagnostic and error propagation values
    diags = [PFX 'cordiags'];
    disp(['Updating ' diags]);
    stn.(diags).date = dts;
    stn.(diags).ustar = result(:,8);
    stn.(diags).tstar = result(:,9);
    stn.(diags).qstar = result(:,10);
    stn.(diags).Cd = result(:,16);
    stn.(diags).Ch = result(:,17);
    stn.(diags).Ce = result(:,18);
    stn.(diags).Ug = result(:,22);
    if ( size(result,2) >= 24 )
      stn.(diags).dtwarm = result(:,23);
      stn.(diags).dxwarm = result(:,24);
    end;

  end; %if isempty(dsfld) else











        % ??? HACKENSTEIN! Need to estimate Stokes drift here from wave data...
        STOKES_DRIFT = 0;
        % ??? HACKOLA! Assume surface current ~ 2% of wind speed [Ardhuin et al. 2009]
        ou = (0.020 .* w) + STOKES_DRIFT;



        % error('WARM LAYER CALCULATION NOT RECODED JUST YET...');


        % ??? HACKENSTEIN! Need to estimate Stokes drift here from wave data...
        STOKES_DRIFT = 0;
        % HACKOLA! Assume surface current ~ 2% of wind speed [Ardhuin et al. 2009]
        ou = (0.020 .* w) + STOKES_DRIFT;











      if ( ~doWarm )

        % Perform Fairall et al. (TOGA-COARE 3.0a) calculation
        disp('COR30A (no warm-layer)');
        disp(length(w));
        nby10 = floor(length(w)/10);
        % ??? HACKENSTEIN! Need to estimate Stokes drift here from wave data...
        STOKES_DRIFT = 0;
        % ??? HACKOLA! Assume surface current ~ 2% of wind speed [Ardhuin et al. 2009]
        ou = (0.020 .* w) + STOKES_DRIFT;
        for ix = 1:length(w)
          res = cor30a([w(ix),ou(ix),t(ix),a(ix),ss(ix),sa(ix),dsrf(ix),dlrf(ix),prcp(ix),pblz,p(ix),wz,az,az,stn.lat,1,1,wvper(ix),wvhgt(ix)]);
          result(ix,1:length(res)) = res;
          if (mod(ix,nby10)==1); disp(ix); end;
        end;

      else

        % error('WARM LAYER CALCULATION NOT RECODED JUST YET...');

        % Perform Fairall et al. (TOGA-COARE 3.0a + WARM CORRECTION) calculation
        disp('COR30A (warm-layer)');
        % ??? HACKENSTEIN! Need to estimate Stokes drift here from wave data...
        STOKES_DRIFT = 0;
        % HACKOLA! Assume surface current ~ 2% of wind speed [Ardhuin et al. 2009]
        ou = (0.020 .* w) + STOKES_DRIFT;
        result = cor30a_warm(dts,w,ou,t,a,ss,sa,dsrf,dlrf,prcp,pblz,p,stn.lon,stn.lat,wz,az,az,stz,1,1,1,wvper,wvhgt);

      end; %if ~doWarm else







      % MHOME = get_ecoforecasts_path('..');
      MHOME = 'c:/Documents and Settings/gramer/My Documents/MATLAB';
      addpath([MHOME '/fairall']);
      rehash



  if ( exist('added_fairall_path','var') )
    warning('off','MATLAB:rmpath:DirNotFound');
    rmpath([MHOME '/fairall']);
    warning('on','MATLAB:rmpath:DirNotFound');
    rehash
  end;







  % Convert winds from kts. to m/s
  w = w .* 0.5144444444;

  % Convert station air pressure to sea-level
  p = p ./ exp( -pz ./ ( (a + 273.15) .* 29.263 ) );







    % Convert winds from [m/s] to [kts] - be consistent with raw data
    if ( isfield(station,'ndbc_wind1_speed') )
      station.ndbc_wind1_speed.data = station.ndbc_wind1_speed.data ./ 0.5144444444;
    end;
    if ( isfield(station,'ndbc_wind1_gust') )
      station.ndbc_wind1_gust.data = station.ndbc_wind1_gust.data ./ 0.5144444444;
    end;






    %DEBUG:    [yix,xix,yerr,xerr,y,x]
    %DEBUG:    [y1,y,y2,x1,x,x2],
    %DEBUG:    if (y1<=1||y2>=nlats||x1<=1||x2>=nlons); disp('I AM ON THE EDGE!'); end;





  vals(:,1) = squeeze(vals(:,1));






  y = yix + yerr;
  y1 = floor(y-eps); y2 = ceil(y+eps);

  x = xix + xerr;
  x1 = floor(x-eps); x2 = ceil(x+eps);





function vals = interp_field(lats,lons,fld,lat,lon,dlat,dlon)
%function vals = interp_field(lats,lons,fld,lat,lon,[dlat],[dlon])
%
% Return bilinear interpolation on a *time series field*: LATS and LONS are
% monotonic vectors of length M and N, resp.; FLD is a matrix of size DxNxM;
% LAT,LON are the coordinates of the desired location. *Optional* args DLAT,
% DLON are minimum grid-spacing in latitude, longitude, and are calculated if
% needed. VALS is a Dx1 vector (time series) of values interpolated on FLD.
%
% This function is offered as a kindness to fellow dyslexics: the orgy of row
% flipping, dimension permuting, and plaid'ing required by INTERP3, and the
% glacial slowness of interpreted looping on INTERP2, all make it seemingly
% impossible in MATLAB to do cleanly and quickly what this function does.
%
% Last Saved Time-stamp: <Thu 2011-01-13 12:21:18  lew.gramer>

  ulats = unique(lats(:));
  ulons = unique(lons(:));
  nlats = length(ulats);
  nlons = length(ulons);

  if ( ~exist('dlat','var') || isempty(dlat) )
    dlat = mean(diff(ulats));
  end;
  if ( ~exist('dlon','var') || isempty(dlon) )
    dlon = mean(diff(ulons));
  end;

  if ( lats(1) > lats(2) )
    lats = lats(end:-1:1);
    fld = fld(:,end:-1:1,:);
  end;
  if ( lons(1) > lons(2) )
    lons = lons(end:-1:1);
    fld = fld(:,:,end:-1:1);
  end;

  [ig,yix] = min(abs(lat-lats));
  [ig,xix] = min(abs(lon-lons));
  yerr = (lat-lats(yix))/dlat;
  xerr = (lon-lons(xix))/dlon;

  y = yix + yerr;
  y1 = floor(y-eps); y2 = ceil(y+eps);

  x = xix + xerr;
  x1 = floor(x-eps); x2 = ceil(x+eps);

  %DEBUG:  [yix,xix,yerr,xerr,y,x]
  %DEBUG:
  [y1,y,y2,x1,x,x2],

  if ( y1<0 || y2>(nlats+1) || x1<0 || x2>(nlons+1) )
    % This could be considered a caller ERROR, but be kind
    warning('interp_field:OutOfBounds',...
            'Site (%g;%g) outside bounding box (%g,%g;%g%g)',...
            lat,lon,min(ulats),max(ulats),min(ulons),max(ulons));
    vals = repmat(nan,[size(fld,1) 1]);

  else
    %DEBUG:    if (y1<=1||y2>=nlats||x1<=1||x2>=nlons); disp('I AM ON THE EDGE!'); end;

    y1ix=y1; y2ix=y2;
    if (y1ix<1);	y1ix=1;		end;
    if (y2ix>nlats);	y2ix=nlats;	end;
    x1ix=x1; x2ix=x2;
    if (x1ix<1);	x1ix=1;		end;
    if (x2ix>nlons);	x2ix=nlons;	end;


    f11 = fld(:,y1ix,x1ix);
    f12 = fld(:,y2ix,x1ix);
    f21 = fld(:,y1ix,x2ix);
    f22 = fld(:,y2ix,x2ix);
    %DEBUG:    [f11(1),f12(1),f21(1),f22(1),],
    %DEBUG:
    [(f11.*(x2-x).*(y2-y)) , (f21.*(x-x1).*(y2-y)) , ...
           (f12.*(x2-x).*(y-y1)) , (f22.*(x-x1).*(y-y1))],

    vals = (f11.*(x2-x).*(y2-y)) + (f21.*(x-x1).*(y2-y)) + ...
           (f12.*(x2-x).*(y-y1)) + (f22.*(x-x1).*(y-y1));
  end;

return;







    %DEBUG:    [y1,y,y2,x1,x,x2],



    if ( y1 < 1 )
      f11 = fld(:,y1,x1);
      f21 = fld(:,y1,x2);
    else
    end;
    if ( y2 > nlats )
    end;
    if ( x1 < 1 )
    end;
    if ( x2 > nlons )
    end;







    %DEBUG:    [y1,y,y2,x1,x,x2],






  else
    %DEBUG:    if (y1<=1||y2>=nlats||x1<=1||x2>=nlons); disp('I AM ON THE EDGE!'); end;

    if (y1<1);		y1=1;		end;
    if (y2>nlats);	y2=nlats;	end;
    if (x1<1);		x1=1;		end;
    if (x2>nlons);	x2=nlons;	end;

    %DEBUG:    [y1,y,y2,x1,x,x2],

    f11 = fld(:,y1,x1);
    f12 = fld(:,y2,x1);
    f21 = fld(:,y1,x2);
    f22 = fld(:,y2,x2);
    %DEBUG:
    [f11(1),f12(1),f21(1),f22(1),],

    vals = (f11.*(x2-x).*(y2-y)) + (f21.*(x-x1).*(y2-y)) + ...
           (f12.*(x2-x).*(y-y1)) + (f22.*(x-x1).*(y-y1));
  end;








  y = yix + yerr;
  if ( yerr >= 0 )
    y1=yix;
    y2=yix+1;
  else
    y1=yix-1;
    y2=yix;
  end;

  x = xix + xerr;
  if ( xerr >= 0 )
    x1=xix;
    x2=xix+1;
  else
    x1=xix-1;
    x2=xix;
  end;






function vals = interp_field(lats,lons,fld,latix,lonix,laterr,lonerr,dlat,dlon)
%function vals = interp_field(lats,lons,fld,latix,lonix,laterr,lonerr,dlat,dlon)
%
% Return bilinear interpolation on a *time series field*: LATS and LONS are
% monotonic vectors of length M and N, resp.; FLD is a matrix of size DxNxM;
% LATIX,LONIX indices in LONS,LATS closest to desired site; LATERR,LONERR
% both in the interval [0,0.5) are grid distances to the nearest gridpoint.
% *Optional* args DLAT,DLON are minimum grid-spacing in latitude, longitude,
% and are calculated if needed. VALS is a Dx1 vector (time series) of values
% interpolated onto FLD.
%
% Last Saved Time-stamp: <Mon 2011-01-10 17:42:13 Eastern Standard Time gramer>

  if ( ~exist('dlat','var') )
    dlat = mean(diff(unique(lats(:))));
  end;
  if ( ~exist('dlon','var') )
    dlon = mean(diff(unique(lons(:))));
  end;

  y = latix + laterr; x = lonix + lonerr;
  if (laterr>0); y1=latix; y2=latix+1; else; y1=latix-1; y2=latix; end;
  if (lonerr>0); x1=lonix; x2=lonix+1; else; x1=lonix-1; x2=lonix; end;
  f11 = dat(:,y1,x1); f12 = dat(:,y2,x1); f21 = dat(:,y1,x2); f22 = dat(:,y2,x2);

  vals = (f11.*(x2-x).*(y2-y)) + (f21.*(x-x1).*(y2-y)) + ...
         (f12.*(x2-x).*(y-y1)) + (f22.*(x-x1).*(y-y1));

return;










    disp(smry');




function r = chkrnd(a,b,dx)
  r = ( abs(a-b) <= dx  );
  % ar=round(a/dx),
  % br=round(b/dx),
  % r = (abs(ar-br) <= 1);
return;





  %DEBUG:  disp({mfilename,'tol/2 unique',numel(unique(ix1)),numel(unique(ix2))});



  %DEBUG:
  disp({mfilename,'tol unique',numel(unique(ix1)),numel(unique(ix2))});




%DEBUG:
tol=(30.0+(0.009/60.0))/(24.0*60.0);



%DEBUG:
tol=(30.0+(1/60.0))/(24.0*60.0);








  if ( ~exist('tol','var') || isempty(tol) )
    % Default tolerance for timestamp matching is 29 minutes
    tol = 29.0/(24.0*60.0);
  end;







  %%%% ??? HACK! Arithmetic above gives a 991x1321 grid, but SST is 990x1320
%   lons(end) = [];
%   lats(end) = [];
  lons(1) = [];
  lats(1) = [];







function h = subplot_tight(m,n,p)
%function h = subplot_tight(m,n,p)

  r = ceil(p/n);
  c = p - ((r-1)*n);

  % Distance between subplots in normalized units
  gutter = 0.07;

  w = (1/n) - (gutter*2);
  h = (1/m) - (gutter*2);

  l = ((c-1)*w) + ((c-1)*(gutter*2));
  b = ((m-r+1)*h) + ((m-r)*(gutter*2));

  h = subplot('position',[l b w h]);

return;







%
% Load station data if it does not already exist in memory
%

if ( ~exist('station', 'var') || isempty(station) )

  matfname = fullfile(datapath, [stnam '-ndbc.mat']);

  if ( exist(matfname, 'file') )
    disp(['Reloading station data from ' matfname]);
    load(matfname, 'station');
    if ( ~isfield(station, 'station_name') || ...
         ~strcmpi(station.station_name, stnam) )
      fprintf(2, '\n\n');
      warning('Station.station_name DID NOT MATCH "%s"', stnam);
      fprintf(2, '\n\n');
      station.station_name = lower(stnam);
    end;

  else
    disp('Reloading station data from raw NDBC files');
    station = load_all_ndbc_data(station, stnam, 1981:2009);
    save(matfname, 'station');

  end;

end;









  % For most coral reef areas, "bathymetry" above MHHW is just a nuisance
  fld = stn.ngdc_92m_bathy.field;
tic,
  fld(fld > -.2) = nan;
toc,
tic,
  contourf(stn.ngdc_92m_bathy.lon,stn.ngdc_92m_bathy.lat,fld,my_contours);
toc,
  maxigraph;
  colorbar;




    % my_contours = -[ 0:4:40 80 160 ];


  % Remove gaps in original from result
  % HACK ALERT! Should use DTS/DAT(GOODIX)! We don't because CONTOURF hates gaps
  x.ts.date = dts(:);
  x.ts.data = dat(:);
  x = filter_gaps(x,'ts','fts',maxgap,(hlp/24.0)/2);
  newdts = x.fts.date(:);
  newdat = x.fts.data(:);





  badix = find(~isfinite(dat));
  dts(badix) = [];
  dat(badix) = [];






  if ( ischar(vals) )
    if ( iscell(vals) )
      kern = find(~strcmp(vals(1:end-1), vals(2:end)));
    else
      % Stupid MATLAB char arrays
      error('Do not handle char arrays yet!');
    end;
  else
    if ( exist('dim','var') || ~isempty(dim) )
      kern = find(diff(vals,[],dim) ~= 0);
    else
      kern = find(diff(vals) ~= 0);
    end;
  end;








        %DEBUG:        set_more; return;



    for yr = 2000:get_year(now)
      for jd = 1:366
    for yr = 2006:2006
      for jd = 1:2





        % sst = read_misst(fname, xlen, ylen);




    cfgline = fgetl(fid);



  switch ( region ),
   case 'world',	lonoff = 0;	latoff = 0;	lonlen = Inf;	latlen = Inf;
   case 'asam',		lonoff = 0;	latoff = 0;	lonlen = Inf;	latlen = Inf;
   case 'ecarib',	lonoff = 0;	latoff = 0;	lonlen = Inf;	latlen = Inf;
   case 'freef',	lonoff = 0;	latoff = 0;	lonlen = Inf;	latlen = Inf;
   case 'gbr',		lonoff = 0;	latoff = 0;	lonlen = Inf;	latlen = Inf;
   otherwise,		error('Unrecognized region "%s"!',region);
  end;




    % DEBUG:    disp(fname);



        fname = fnames{ix};
        %DEBUG:
        disp(fname);
        % sst = read_misst(fname, xlen, ylen);
        [lon,lat,sst,dx] = read_misst_region(rgn, fname)
        for stnix = 1:length(stns)
          stns.sst.date(end+1) = datenum(yr,1,1) + jd - 1;
          if ( strcmp(rgn,'world') )
            stns.sst.data(end+1) = sst(stns.misst_latix,stns.misst_lonix);
          else
            stns.sst.data(end+1) = sst(stns.(latfld),stns.(lonfld));
          end;
        end;
        %DEBUG:        set_more; return;




    for yr = 2006
      for jd = 1:2





    dataset = 'fusion';
    dataset = [region '.fusion'];





    for ix = 1:length(fnames)
      fname = fnames{ix};
      %DEBUG:
      disp(fname);
      sst = read_misst(fname, xlen, ylen);
      for stnix = 1:length(stns)
        stns.sst.date(end+1) = datenum(yr,1,1) + jd - 1;
        stns.sst.data(end+1) = sst(stns.(latfld),stns.(lonfld));
      end;
    end;



        stns(nstns).misst_lonix = str2num(flds{4}{:}):str2num(flds{5}{:}) - xoff;





  fid = fopen(cfgfname,'r');
  flds = textscan(fid, '%[^,],%[^,],%[^,],%[^,],%[^,\n]\n');
  fclose(fid);

  if ( numel(flds) ~= 3 )
    error('Bad format in config file "%s"!', cfgfname);
  end;

  % Eliminate all blank lines and comments
  commentlns = regexp(flds{1},'^ *[#]');
  goodix = find(cellfun(@isempty,commentlns));
  flds{1} = flds{1}(goodix);
  flds{2} = flds{2}(goodix);
  flds{3} = flds{3}(goodix);

  % Sanity check - may go away after this function is stable
  for ix = 1:length(flds{1})
    stnm = flds{1}{ix};
    stns(ix).misst_lonix = eval(flds{2}{cfgix});
    stns(ix).misst_latix = eval(flds{3}{cfgix});

    if ( isempty(stns(ix).misst_lonix) || ...
         any(size(stns(ix).misst_lonix) ~= size(stns(ix).misst_latix)) )
      warning('MISST:BadConfigIndices', ...
              'Bad indices "%s,%s" loaded for station "%s"!', ...
              flds{2}{cfgix}, flds{3}{cfgix}, stns(ix).name);
    end;
  end;









  dts = stn.(fld).date;
  dat = stn.(fld).data;
  [yrs,mos,dys] = datevec(dts);
  if ( ~isempty(domos) )
    domoix = find(ismember(mos,domos));
    yrs = yrs(domoix);
    dat = dat(domoix);
  end;





  % contourf(stn.ngdc_92m_bathy.lon,stn.ngdc_92m_bathy.lat,stn.ngdc_92m_bathy.field);
  % my_contours = -[ 0:2:16 ];
  % my_contours = -[ 0:1:10 12:20 30 50 80 100 150 200 ];
  % my_contours = -[ 0:1:10 12:20 30 50 80 -nanmin(stn.ngdc_92m_bathy.field(:)) ];





    % This mesh would be HUGE: We don't really seem to need it??
    % [freef_ngdc_topo.lons,freef_ngdc_topo.lats] = meshgrid(freef_ngdc_topo.tlon,freef_ngdc_topo.tlat);




    cs = round(size(stn.ngdc_92m_bathy.field)./2);
    cy = cs(1); cx = cs(2);
%     if ( nanmin(nanmin(stn.ngdc_92m_bathy.field(cx-5:cx+5,cy-5:cy+5))) > -20 )
      my_contours = -[ 0:1:16 -nanmin(stn.ngdc_92m_bathy.field(:)) ];
%     end;




    result.field = griddata(stn.(fldnm).lon,stn.(fldnm).lat,stn.(fldnm).field, ...
                            result.lon,result.lat,'cubic');




  if ( isfield(stn,'station_name') )
    stnm = lower(stn.station_name);
  elseif ( ischar(stn) )
    stnm = lower(stn);
    clear stn;
    stn.station_name = upper(stnm);
  else
    error('First arg should be either station name or struct with .station_name field!');
  end;





    lon = rawxyz.lon;
    lat = rawxyz.lat;
    depth = rawxyz.depth;





[0:1:15]);%


xyz_bathy_subset_station




% Find all points "inside" our bounding ellipse: ugh, geography sucks
smaecc(1) = xrad;
smaecc(2) = axes2ecc(xrad,yrad);
% WGS 84 reference ellipsoid for the Earth, in [km]
wgs84_ell = [6356.752 (1/298.25722356)];
[elllat,elllon] = ellipse1(stn.lat,stn.lon,smaecc,90,[],wgs84_ell,[],100);


    % Find all points "inside" our bounding ellipse: ugh, geography sucks
    smaecc(1) = xrad;
    smaecc(2) = axes2ecc(xrad,yrad);
    % WGS 84 reference ellipsoid for the Earth, in [km]
    wgs84_ell = [6356.752 (1/298.25722356)];
    [elllat,elllon] = ellipse1(stn.lat,stn.lon,smaecc,90,[],wgs84_ell,[],100);
    inix = find(inside(lon,lat, elllon,elllat) > 0);







    ddeg = sqrt(((stn.lon - lon).^2) + ((stn.lat - lat).^2));
    inix = find(ddeg < (rad
    dkm = sw_dist([
    inix = find(





  % nlon=length(unique(LON));
  % nlat=length(unique(LAT));
  % stn.ngdc_92m_bathy.field = reshape(dat',[nlon nlat]);
  stn.ngdc_92m_bathy.field = griddata(lon,lat,dat,LON',LAT');






LOAD_NCDC_WINDS.m:


  emptyts.date = []; emptyts.data = [];

  [stns.ncdc_wind_speed] = deal(emptyts);
  [stns.ncdc_wind_dir] = deal(emptyts);
  [stns.ncdc_wind_u] = deal(emptyts);
  [stns.ncdc_wind_v] = deal(emptyts);

  stns(1).ncdc_wind_speed.date = dts;
  stns(1).ncdc_wind_speed.data = rwnd;

  stns(1).ncdc_wind_dir.date = dts;
  stns(1).ncdc_wind_dir.data = rdir;

  stns(1).ncdc_wind_u.date = dts;
  stns(1).ncdc_wind_u.data = rwnd .* (-sind(rdir));

  stns(1).ncdc_wind_v.date = dts;
  stns(1).ncdc_wind_v.data = rwnd .* (-cosd(rdir));

  save('data/asam_stns_0km.mat', '-append', 'stns');




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DEBUG





  x.t = t;
  x.f = f;
  x = blank_derived_var( x,'t','f',[],(ceil(n/2)/24) );
  f = x.f;






  morestate = get(0, 'More');
  more('off');

  more(morestate);






%   figure;
%   plot(1:12, t);
%   maxigraph;

%   figure;
%   plot(1:12, [n ; rf ; tf]);
%   maxigraph;
%   legend('Q_N_E_T', 'Q_R_A_D', 'Q_T_U_R');







%   set(lh1,'box','off');
%   set(lh2,'box','off');
%   pos1 = get(lh1,'position');
%   pos2 = get(lh2,'position');
%   pos2
%   set(lh2,'position',pos2);






      %DEBUG:      datestr(stn.quasi_eulerian_u.date(uix)),




  % Create colorbar
  % axes('Units','norm', 'Position',[0.05 0.05 0.9 0.9], 'Visible','off');
  % colorbar('FontSize',18);
  axes('Units','norm', 'Position',[0.10 0.10 0.80 0.80], 'Visible','off');
  colorbar;
  caxis([cmin cmax]);

  % cbf=figure;
  % set(gca,'vis','off');
  % cbh=colorbar;
  % caxis([cmin cmax]);
  % set(cbh,'units','norm','outerpos',[0 0 1 1]);
  % set(cbf,'units','norm','outerpos',[.5 .1 .1 .8]);






      if ( ~isempty(strfind(tide_fld,'tide')) )
        % Make sure we have a value in [m] (not FEET!) to compare to...
        tide_fld_m = [tide_fld '_m'];
        stn.(tide_fld_m) = stn.(tide_fld);
        stn.(tide_fld_m).data = stn.(tide_fld_m).data * (12*0.0254);
      end;








function stn = station_tmd_tide(stn,dts)
%function stn = station_tmd_tide(stn,dts)
%
% Use TMD_TOOLBOX (v.) to predict tides for station STN, using both the
% global TPXO7.2 solution, and best available regional solution. Assumes
% STN.station_name contains 5-char code for site; if STN.lon and STN.lat are
% not found containing site coordinates, calls GET_STATION_COORDS for them.
%
% Last Saved Time-stamp: <Mon 2010-07-12 16:23:01 Eastern Daylight Time gramer>

  datapath = get_ecoforecasts_path('data');


  tide_mpath = which('tide_pred');
  if ( isempty(tide_mpath) )
    error('No path found to TIDE_PRED function!');
  end;
  [tide_mpath,ig] = fileparts(tide_mpath);
  tidepath = fullfile(tide_mpath,'DATA');
  if ( ~exist(tidepath,'dir') )
    error('No ./DATA sub-directory in TMD path "%s"!',tide_mpath);
  end;
  %DEBUG:
  dir(tidepath),

  if ( ~isfield(stn,'station_name') )
    error('STN struct must have "station_name" field!');
  end;

  matfname = fullfile(datapath,[stn.station_name '_tmd_tide.mat']);

  if ( exist(matfname,'file') )
    disp(['Loading pre-saved file ' matfname]);
    load(matfname,'station');

    if ( isfield(station,'tmd_tide') )
      stn.tmd_tide = station.tmd_tide;
      % [ix1,ix2] = intersect_dates(dts,station.tmd_tide.date);
      % stn.tmd_tide.date = station.tmd_tide.date(ix2);
      % stn.tmd_tide.data = station.tmd_tide.data(ix2);
    end;
    if ( isfield(station,'tpxo_tide') )
      stn.tpxo_tide = station.tpxo_tide;
      % [ix1,ix2] = intersect_dates(dts,station.tpxo_tide.date);
      % stn.tpxo_tide.date = station.tpxo_tide.date(ix2);
      % stn.tpxo_tide.data = station.tpxo_tide.data(ix2);
    end;
    station = []; clear station;

  else

    if ( ~exist('dts','var') || isempty(dts) )
      tide_fld = 'ctd_deep_i_depth';
      if (~isfield(stn,tide_fld)); tide_fld = 'ctd_shallow_i_depth';	end;
      if (~isfield(stn,tide_fld)); tide_fld = 'tide';			end;
      if (~isfield(stn,tide_fld)); tide_fld = 'ndbc_tide';		end;
      if (~isfield(stn,tide_fld)); tide_fld = 'wave_tide_height';	end;
      if (~isfield(stn,tide_fld)); tide_fld = 'sea_t';		end;
      if (~isfield(stn,tide_fld)); tide_fld = 'ct_shallow_seatemp';	end;
      if (~isfield(stn,tide_fld))
        error('No DTS specified and no appropriate comparison field found!');
      else
        dts = stn.(tide_fld).date(1):(1/24):stn.(tide_fld).date(end);
      end;
    end;

    if ( ~isfield(stn,'lon') )
      [stn.lon,stn.lat,stn.depth] = get_station_coords(stn.station_name);
    end;
    if ( ~isfield(stn,'lon') )
      error('Cannot find coordinates for station "%s"!', stn.station_name);
    end;

    stn.tmd_tide.date = dts;
    stn.tpxo_tide.date = dts;

%     mydir = pwd;
%     cd(tide_mpath);
    stn.tmd_tide.data = tide_pred(fullfile(tidepath,'Model_Mex'),dts,stn.lat,stn.lon,'z');
    stn.tpxo_tide.data = tide_pred(fullfile(tidepath,'Model_tpxo7.2'),dts,stn.lat,stn.lon,'z');
%     cd(mydir);

    disp(['Saving result to ' matfname]);
    station.tmd_tide = stn.tmd_tide; 
    station.tpxo_tide = stn.tpxo_tide;
    save(matfname,'station');
    station = []; clear station;

  end;

return;








  Q0_factor = Q0factor(t,[],stz);

  heat_flux_term = [PFX 'heat_flux_term'];
  stn.(heat_flux_term).date = stn.(net_heat_flux).date;
  stn.(heat_flux_term).data = (stn.(net_heat_flux).data ./ Q0_factor) .* (60*60);

  latent_flux_term = [PFX 'latent_flux_term'];
  stn.(latent_flux_term).date = stn.(latent_heat_flux).date;
  stn.(latent_flux_term).data = (stn.(latent_heat_flux).data ./ Q0_factor) .* (60*60);

  sensible_flux_term = [PFX 'sensible_flux_term'];
  stn.(sensible_flux_term).date = stn.(sensible_heat_flux).date;
  stn.(sensible_flux_term).data = (stn.(sensible_heat_flux).data ./ Q0_factor) .* (60*60);

  if ( ~isempty(sfld) )
    shortwave_flux_term = [PFX 'shortwave_flux_term'];
    stn.(shortwave_flux_term).date = stn.(sfld).date;
    stn.(shortwave_flux_term).data = (stn.(sfld).data ./ Q0_factor) .* (60*60);
  end;
  if ( ~isempty(lfld) )
    longwave_flux_term = [PFX 'longwave_flux_term'];
    stn.(longwave_flux_term).date = stn.(lfld).date;
    stn.(longwave_flux_term).data = (stn.(lfld).data ./ Q0_factor) .* (60*60);
  end;











function wdir = uv_to_dir(u, v)
%function wdir = uv_to_dir(u, v)
% Convert u and v vector components into a direction in degrees True
% NOTE: This version of 'uv-to-dir' is coded for WINDS: "dir" here means
% "source direction", NOT "target direction" as it would for ocean currents.
% Last Saved Time-stamp: <Sun 2010-05-23 11:42:19 Eastern Daylight Time gramer>

    wdir = repmat(nan, size(u));

    spd = uv_to_spd(u, v);

    b0idx = find(u==0 & v==0);
    bpidx = find(u> 0 & v> 0);
    vpidx = find(u<=0 & v> 0);
    upidx = find(u> 0 & v<=0);
    npidx = find(u<=0 & v<=0);

    wdir(b0idx) = 0;
    wdir(bpidx) = ( 180 + asind(u(bpidx) ./ spd(bpidx)) );
    wdir(vpidx) = ( 180 + asind(u(vpidx) ./ spd(vpidx)) );
    wdir(upidx) = ( 360 - asind(u(upidx) ./ spd(upidx)) );
    wdir(npidx) = ( - asind(u(npidx) ./ spd(npidx)));

return;







  if ( isempty(tol) )
    tol = 29/(24*60);
  end;
  if ( ~isnumeric(tol) || ~isscalar(tol) )
    error('First arg should be a numeric TOLERANCE');
  end;

  dts = varargin;
  ixs = {};

  origdts1 = dts{1};
  for ix = 2:length(dts)
    [ixs{1},ig] = intersect_dates(dts{1},dts{ix});
    dts{1} = dts{1}(ixs{1});
  end;
  ixs{1} = find(ismember(origdts1,dts{1}(ixs{1})))';
  for ix = 2:length(dts)
    [ig,ixs{ix}] = intersect_dates(dts{1},dts{ix});
  end;








function stn = station_advect_field(stn,UdotdelQ,centerx,centery,euleru,eulerv,field)
%function stn = station_advect_field(stn,UdotdelQ,centerx,centery,euleru,eulerv,field)
%
% Advect the time series of 2-D fields STN.(FIELD).field across the center
% indices [centerx centery], using the time series of Eulerian current vector
% components STN.(EULERU) and STN.(EULERV). Save the resultant time series as
% STN.(UDOTDELQ). All STN time series must have .date (datenum) subfields;
% STN.(FIELD) should also have .lon and .lat vector subfields.
%
% Last Saved Time-stamp: <Sun 2010-05-16 11:46:05 Eastern Daylight Time gramer>

  lon = stn.(field).lon(centerx-2:centerx+2);
  lat = stn.(field).lat;
  dts = stn.(field).date;
  fld = stn.(field).field(:,:,centerx-2:centerx+2);

  [ig,ix1,ix2] = intersect(floor(dts),floor(stn.(euleru).date));
  dts = dts(ix1);
  fld = fld(ix1,:,:);
  % Goofy INTERP3 (below) assumes 3-d matrix is in LAT,LON,TIME order!
  fld = permute(fld,[2 3 1]);


  % For simplicity, assume for now that FLD is on a 1x1 km grid and that
  % STN.(EULERU).data and STN.(EULERV).data are [m/s]. Convert to [km/h].

  eudts = stn.(euleru).date(ix2);
  euu = stn.(euleru).data(ix2) .* (3600/1e3);
  euv = stn.(eulerv).data(ix2) .* (3600/1e3);


  % U and V are in target direction - subtract
  farx = centerx - euu;
  fary = centery - euv;

  stn.(UdotdelQ).date = eudts;
  stn.(UdotdelQ).data = interp3(lon,lat,dts,fld,farx,fary,eudts,'linear',nan);

return;











  stn.(UdotdelT).data = interp3(lon,lat,dts,fld,farx,fary,eudts,'spline',nan);








  [ig,ig,dyix] = unique(floor(stn.(euleru).date));






  stn.(UdotdelT).data = ...
      interp3(dts(dyix),lon,lat,fld(dyix,:,:),eudts,farx,fary,'spline');







  % val = interp2(lon,lat,squeeze(fld(dix,:,:)),linx(dix,:),liny(dix,:),'spline');




  for dix = 1:length(dts)
    s = [0:23] / 23;
    linx = (centerx.*s) + (farx.*(1-s));
    liny = (centery.*s) + (fary.*(1-s));

    vals = interp2(lon,lat,squeeze(fld(dix,:,:)),linx(dix,:),liny(dix,:),'spline');

  end;








  dts = stn.(field).date;


  [ig,ix1,ix2] = intersect(floor(dts),floor(stn.(euleru).date));
  dts = dts(ix1(1)):(1/24):dts(ix1(end));








    if ( idx > 1 && begidx > 2 && strcmpi(varname(findidx(idx-1):begidx-1),'ndbc_') )
      instfld = ['ndbc_' instfld];
      newfld = ['ndbc_' newfld];
    end;







  sd = (5e-4.*(1.25 - (0.25.*((0.5./fc).^(1/3)))).*w.*basew) + (0.025.*(hs-0.4));








  if ( isnumeric(hsfld) || isnumeric(tpfld) || isnumeric(tdfld) )
    if ( ischar(hsfld) )
      [ix1,ix2] = ...
          intersect_dates(stn.(wfld).date,stn.(hsfld).date);
      w = stn.(wfld).data(ix1);
      hs = stn.(hsfld).data(ix2);
      tp = tpfld;
    elseif ( ischar(tpfld) )
      [ix1,ix2] = ...
          intersect_dates(stn.(wfld).date,stn.(tpfld).date);
      w = stn.(wfld).data(ix1);
      hs = hsfld;
      tp = stn.(tpfld).data(ix2);
    else
      ix1 = 1:length(stn.(wfld).date);
      w = stn.(wfld).data;
      hs = hsfld;
      tp = tpfld;
    end;
  else
    [ix1,ix2,ix3] = ...
        intersect_all_dates(stn.(wfld).date,stn.(hsfld).date,stn.(tpfld).date);
    w = stn.(wfld).data(ix1);
    hs = stn.(hsfld).data(ix2);
    tp = stn.(tpfld).data(ix3);
  end;

















      if ( strcmp(wpfld,'default') )
        % [wvper,wvhgt] = wind_to_wave(w);
        Ca = 0.018;
        Cb = 0.729;
        wvper = Cb .* w;
        wvhgt = Ca .* (w.^2) .* (1 + (0.015 .* w));
      else
        wvper = interp1_stn_fld(stn, wpfld, dts);
        wvhgt = interp1_stn_fld(stn, whfld, dts);
      end;





    rmpath([MHOME '/fairall']);
    rmpath([MHOME '/fairall']);





function stn = station_par_to_insol(stn,parfld,wndfld,cldfld,insfld,uswfld,goodix)
%function stn = station_par_to_insol(stn,parfld,wndfld,cldfld,insfld,uswfld,goodix)
%
% Calculate new fields insolation STN.(INSFLD), upward (reflected) short
% wave radiation STN.(USWFLD), and net shortwave radiation SWFLD, using in
% situ (or other) PAR light field STN.(PARFLD), wind speed STN.(WNDFLD) and cloud cover STN.(CLDFLD). If
% STN.(CLDFLD) does not exist, it is estimated from STN.(PARFLD) using hourly
% solar altitude (see SORADNA1). If GOODIDX is specified, only use the given
% indices from the time series STN.(PARFLD) when estimating insolation.
%
% Last Saved Time-stamp: <Fri 2010-04-09 15:55:18 Eastern Daylight Time gramer>

  if ( ~exist('goodix','var') || isempty(goodix) )
    goodix = 1:length(stn.(parfld).date);
  end;

  stn.(insfld).date = stn.(parfld).date(goodix);
  stn.(insfld).data = par_to_insol(stn.(parfld).data(goodix));

  % For a while, used empirical relationship derived from regressing NCEP
  % USRF against NCEP DSRF and NCEP wind speed in knots. Before that, just
  % assumed a constant Albedo, chosen in Godfrey et al (1991) as 0.06. Then
  % tried 0.04 based on regression using TOMS GSIP satellite estimates. Now
  % we have finally settled on using SWHF from the AIR_SEA toolbox instead!

  % [ix1,ix2] = intersect_dates(stn.(insfld).date,stn.(wndfld).date);

  % alpha = 0.04;
  % stn.(uswfld).date = stn.(insfld).date;
  % stn.(uswfld).data = stn.(insfld).data .* alpha;

  stn.(uswfld).date = stn.(insfld).date;
  stn.(uswfld).data = stn.(insfld).data .* alpha;


return;








    % Air specific humidity [kg/kg]
    sa = relhumid_to_spechumid(a,q);
    sa = sa .* 1e3;
    % Saturated ("sea-surface") specific humidity [kg/kg]
    % ss = relhumid_to_spechumid(t,100);
    % ss = ss .* 1e3;
    ss = qsee([t p]);






        if (ix==(nby10*2)); break; end;





    % Air specific humidity [kg/kg]
    sa = relhumid_to_spechumid(a,q);
    %sa = sa .* 1e3;
    % Saturated ("sea-surface") specific humidity [kg/kg]
    % ss = relhumid_to_spechumid(t,100);
    % ss = ss .* 1e3;
    ss = qsee([t p]);
    sa = sa ./ 1e3;










  hf = [PFX 'latent_heat_flux'];
  ft = [PFX 'latent_flux_term'];
  stn.(ft).date = stn.(hf).date;
  stn.(ft).data = (stn.(hf).data ./ Q0f) .* (60*60);

  hf = [PFX 'sensible_heat_flux'];
  ft = [PFX 'sensible_flux_term'];
  stn.(ft).date = stn.(hf).date;
  stn.(ft).data = (stn.(hf).data ./ Q0f) .* (60*60);








  % Perform Fairall et al. (TOGA-COARE 2.0) calculation
  if ( ~isempty(dsfld) )
    disp('HFBULKTC');
%     result = hfbulktc(w,wz,a,az,q,az,p,t);
srf = interp1_stn_fld(stn,sfld,dts);
result = hfbulktc(w,wz,a,az,q,az,p,t,36,dlrf,dsrf,srf);

    tau = result(:,4);
    shf = result(:,1);
    lhf = result(:,2);

  else

    % Air specific humidity [kg/kg]
    sa = relhumid_to_spechumid(a,q);
    % Saturated ("sea-surface") specific humidity [kg/kg]
    ss = relhumid_to_spechumid(t,100);

    % (Atmospheric) Planetary Boundary Layer height == inversion height [m]
    pblz = 600; %%%% IS THIS A GOOD MEAN VALUE FOR FLORIDA KEYS??

    if ( isempty(wpfld) )

      % Perform Fairall et al. (TOGA-COARE 2.6) calculation
      disp('COR26');
      % ??? HACKOLA! Assume surface current projects onto 1.5% of wind velocity [Shay et al. 1998]
      ou = 0.015 .* w;
      length(w),
      for ix = 1:length(w)
        res = cor30a([w(ix),ou(ix),t(ix),a(ix),ss(ix),sa(ix),dsrf(ix),dlrf(ix),prcp(ix),pblz,p(ix),wz,az,az,stn.lat,1,0,0,0]);
        result(ix,1:length(res)) = res;
        if (mod(ix,10000)==1); ix, end;
      end;

    else

      if ( ~doWarm )

        % Perform Fairall et al. (TOGA-COARE 3.0a) calculation
        disp('COR30A - NO WARM LAYER');
        % ??? HACKENSTEIN! Need to estimate Stokes drift here from wave data...
        STOKES_DRIFT = 0;
        % ??? HACKOLA! Assume surface current ~ 1.5% of wind speed [Shay et al. 1998]
        ou = (0.015 .* w) + STOKES_DRIFT;
        result = cor30a([w,ou,t,a,ss,sa,dsrf,dlrf,prcp,pblz,p,wz,az,az,stn.lat,1,1,wvper,wvhgt]);

      else

        error('WARM LAYER CALCULATION NOT RECODED JUST YET...');
        % % Perform Fairall et al. (TOGA-COARE 3.0a + WARM CORRECTION) calculation
        % disp('COR30A with Warm Layer');
        % % HACKOLA! Assume surface current ~ 1.5% of wind speed [Shay et al. 1998]
        % ou = (0.015 .* w) + STOKES_DRIFT;
        % result = cor30a_warm([w,ou,t,a,ss,sa,dsrf,dlrf,prcp,pblz,p,wz,az,az,stn.lat,1,1,wvper,wvhgt,???]);

      end; %if doWarm else

    end; %if isempty(wpfld) else

    tau = result(:,1);
    shf = result(:,2);
    lhf = result(:,4);
    rhf = result(:,14);

  end; %if isempty(dsfld) else









    tau = result(:,4);
    shf = result(:,1);
    lhf = result(:,2);





  stn.(wind_stress).data = result(:,4);
  % stn.(wind_stress).data = result(:,1);

  stn.(sensible_heat_flux).data = result(:,1);
  % stn.(sensible_heat_flux).data = result(:,2);

  stn.(latent_heat_flux).data = result(:,2);
  % stn.(latent_heat_flux).data = result(:,4);






  % Conversion factor from Insolation to PAR
  % See Papaioannou, Papanikolaou and Retalis, 1993
  PAR_PER_INSOL = 0.473;

  % Conversion factor from W/m^2 to micromole quanta/m^2.s
  % See Morel and Smith, 1974
  UMOL_PER_WATT = 4.1513469579;
  % NOTE: Dye (JGR-A, 2004) recommends 4.56 mumol/J instead!

  stn.(insfld).date = stn.(parfld).date(goodix);
  stn.(insfld).data = (stn.(parfld).data(goodix) ./ UMOL_PER_WATT) ./ PAR_PER_INSOL;





MULTIPLOT_STATION:
  if (nargin < 5 || isempty(ylbl)); ylbl = strrep(lower(fldnms),'_','\_'); end;
  if (nargin < 6); xlm = []; end;
  if (nargin < 7); doverify = true; end;
  if (nargin < 8); linspec = []; end;

  X = {};
  Y = {};
  max_ndates = 0;
  for ix = 1:length(fldnms)
    fldnm = fldnms{ix};
    if ( doverify )
      stn = verify_variable(stn, fldnm);
    end;
    if ( ~isfield(stn, fldnm) )
      warning('No field "%s" in station struct!', fldnm);
      ylbl(ix) = [];
    else
      X{end+1} = stn.(fldnm).date;
      Y{end+1} = stn.(fldnm).data;
      max_ndates = max(max_ndates,length(stn.(fldnm).date));
    end;
  end;






  % Calculation total cloud fraction (if it was not supplied in DSRFLD)
  if ( ~strcmp(dsrfld,cfld) )
    if ( isfield(stn,cfld) ); stn = rmfield(stn,cfld); end;

    [yr,mo,dy,hr] = datevec(stn.(dsrfld).date);
    yd = stn.(dsrfld).date - datenum(yr,1,1);
    [theta,cscdsr] = soradna1(yd,yr,-stn.lon,stn.lat);
    % This method produces unreliable results near sunrise/set
    eveningix = find(theta < 10);

    peakdsr = nanmax(stn.(dsrfld).data);
    Cdts = stn.(dsrfld).date;
    Craw = 100 - (100 .* stn.(dsrfld).data ./ (sind(theta).*peakdsr));

    Cdts(eveningix) = [];
    Craw(eveningix) = [];
    Craw(Craw < 0) = 0;
    Craw(isnan(Craw)) = 0;

    % stn.(cfld).date = stn.(dsrfld).date;
    % % This method produces NO RESULT (thus interpolation) during night hours
    % stn.(cfld).data = interp1(Cdts, Craw, stn.(cfld).date, 'linear');
    % DEBUG:
    stn.(cfld).date = Cdts;
    stn.(cfld).data = Craw;
  end;








  legend( [strrep(xname,'_','\_') ' vs. ' strrep(yname,'_','\_')], ...
           sprintf('%.5g + %.5g*X, R^2~%0.2g, \n RMSE=%g, p=%0.5g, N=%d', ...
                   B(1), B(2), (Stats.coeffcorr(2,1)^2), ...
                   Stats.s, roundn(Stats.p(2),-5), length(x)), ...
           'Location','NorthWest');





  dtP = stn.(qfld).date(qidx);
  [ig,ix] = intersect_dates(dtTa,dtQ);

  dtC = stn.(qfld).date(qidx);
  [ig,ix] = intersect_dates(dtTa,dtQ);

  dtS = stn.(qfld).date(qidx);
  [ig,ix] = intersect_dates(dtTa,dtQ);







STATION_PAR_TO_INSOL:
  % This needs some big improvements probably!
  % Matches assumption by Godfrey et al (1991)
  alpha = 0.06;

  stn.(uswfld).date = stn.(insfld).date;
  stn.(uswfld).data = stn.(insfld).data .* alpha;





FROM RSTOOL_STATION:
       'yhat', ...
       'leverage', ...
       'hatmat', ...
       's2_i', ...
       'beta_i', ...
       'dfbetas', ...





FROM RSTOOL_STATION:
  if ( ~exist('ylbl','var') || isempty(ylbl) )
    ylbl = strrep(fldy, '_', '\_');
  end;
  if ( ~exist('xlbls','var') || isempty(xlbls) )
    xlbls = strrep(flds, '_', '\_');
  end;

  if ( ~exist('fixy','var') || isempty(fixy) )
    fixy = 1:numel(stn.(fldy).data);
  elseif ( isa(fixy,'function_handle') )
    ylbl = [ylbl ' (' upper(strrep(func2str(fixy),'_','\_')) ')'];
    fixy = fixy(fldy);
  end;
  if ( ~exist('fixes','var') || isempty(fixes) )
    for ix = 1:nx
      fixes{ix} = 1:numel(stn.(flds{ix}).data);
    end;
  else
    for ix = 1:nx
      if ( isa(fixes{ix},'function_handle') )
        xlbls{ix} = [xlbls{ix} ' (' upper(strrep(func2str(fixes{ix}),'_','\_')) ')'];
        fixes{ix} = fixes{ix}(flds{ix});
      end;
    end;
  end;








      else
        warning('Result was NOT a struct: Did NOT save "%s"!', matfname);




    % Remove fields not needed for EFs - SAVES MEMORY
    stn = remove_diag_fields(stn);
    stn = remove_intrahour_winds(stn);





  % Retrieve data from a path relative to this M-file's local directory
  [pathroot, ig, ig, ig] = fileparts(mfilename('fullpath'));
  if ( ~exist('datapath', 'var') || isempty(datapath) )
    datapath = fullfile(pathroot, 'data', '');
    datapath = get_ecoforecasts_path('data');
  end;
  if ( ~exist(datapath, 'dir') )
    error('The data source path "%s" does not exist!', datapath);
  end;








    elseif ( iscellstr(arg) )
      for ix = 1:length(arg)
        [X,Y,ylbl] = multiplot_ts_get_arg(curstn,curstnm,arg{ix});
        argn = arg{ix};
        X = { X{:} curstn.(argn).date };
        Y = { Y{:} curstn.(argn).data };
        ylbl = { ylbl{:} sprintf('%s.%s', curstn.station_name, argn) };
      end;








FROM SCATTER_FIT:
  % Doing regress here just to estimate an R^2
  Y = y;
  X = [x , repmat(1,size(x))];
  [r_B,r_BINT,r_R,r_RINT,r_STATS] = regress(Y,X);

  plot(x,y,'b.');
  plot(x,(B(1)+(B(2).*x)),'r-');
  legend( [strrep(xname,'_','\_') ' vs. ' strrep(yname,'_','\_')], ...
           sprintf('%.5g + %.5g*X, R^2~%0.2g, \n RMSE=%g, p=%0.5g, N=%d', ...
                   B(1), B(2), r_STATS(1), ...
                   Stats.s, roundn(Stats.p(2),-5), length(x)), ...
           'Location','NorthWest');
  xlabel(xname); ylabel(yname);




% EXTREME TWEAKING...
  sri = 3;
  switch ( fuzzy ),
   case {'drastic-low', 'drastic-high', 'very-conducive'},
    sri = sri * 3;
   case {'very-low', 'very-high', 'conducive'},
    sri = sri * 2.5;
   case {'low', 'high', 'somewhat-conducive'},
    sri = sri * 2;
  end;





function sri = fuzzy_sri(fuzzy)

% % PER ICON/G2
%   sri = 3;
%   switch ( fuzzy ),
%    case {'drastic-low', 'drastic-high', 'conducive'},
%     sri = sri * 2.5;
%    case {'very-low', 'very-high', 'somewhat-conducive'},
%     sri = sri * 2;
%   end;

% % PER Jim's (and Jank's) original CLIPS code
%   sri = 3;
%   switch ( fuzzy ),
%    case {'drastic-low', 'drastic-high', 'conducive'},
%     sri = sri * 2;
%   end;

% TWEAKING...
  sri = 3;
  switch ( fuzzy ),
   case {'drastic-low', 'drastic-high', 'very-conducive'},
    sri = sri * 2.5;
   case {'very-low', 'very-high', 'conducive'},
    sri = sri * 2;
   case {'somewhat-low', 'somewhat-high', 'marginal-conducive'},
    sri = sri * 0.5;
  end;

return;







          %DEBUG:
          if ( ishandle(fhs(fi)) );          fhname = get(fhs(fi),'Name');
          else          fhname = sprintf('Figure %g', fhs(fi));        end;
          disp(sprintf('UIWAIT: "%s"',fhname));








        if ( ishandle(fhs(fi)) );          fhname = get(fhs(fi),'Name');
        else          fhname = sprintf('Figure %g', fhs(fi));        end;
        warndlg(sprintf('Checking "%s"',fhname),'Check','modal');






      txt = sprintf('%s (%g)', rng(1), rng{2}(1));



    % PRCTILES = [0.1 1 3 7 15 30 70 85 93 97 99 99.9];
    % PRCTILES = [00.13 00.62 02.27 06.68 15.87 30.85 69.15 84.13 93.32 97.72 99.38 99.87];







function [cutoffs,ranges,FUZZIES] = create_prctile_fuzzy(stn, varname, PRCTILES, FUZZIES)
%function [cutoffs,ranges,FUZZIES] = create_prctile_fuzzy(stn, varname, PRCTILES, FUZZIES)
%
% Use percentiles (v. PRCTILE, MATLAB Statistics Toolbox) to generate simple
% fact fuzzy ranges for variable VARNAME in station STN. PRCTILES: optional
% 2nd arg to PRCTILE, DEFAULT: [0.01 1 3 7 15 30 70 85 93 99 99.99].
%
% Last Saved Time-stamp: <Tue 2010-02-16 11:02:08 Eastern Standard Time gramer>

% SYMBOLIC VALUE                  LOWER BOUND PERCENTILE
% 'unbelievably-low',             0
% 'drastic-low',                  1
% 'very-low',                     3
% 'low',                          7
% 'somewhat-low',                 15
% 'average',                      30
% 'somewhat-high',                70
% 'high',                         85
% 'very-high',                    93
% 'drastic-high',                 97
% 'unbelievably-high',            99


  if ( ~exist('PRCTILES','var') || isempty(PRCTILES) )
    % PRCTILES = [0.1 1 3 7 15 30 70 85 93 97 99 99.9];
    % PRCTILES = [00.13 00.62 02.27 06.68 15.87 30.85 69.15 84.13 93.32 97.72 99.38 99.87];
    PRCTILES = [00.62 02.27 06.68 15.87 30.85 69.15 84.13 93.32 97.72 99.38];
  end;

  if ( ~exist('FUZZIES','var') || isempty(FUZZIES) )
    FUZZIES = { ...
        'unbelievably-low', ...
        'drastic-low', ...
        'very-low', ...
        'low', ...
        'somewhat-low', ...
        'average', ...
        'somewhat-high', ...
        'high', ...
        'very-high', ...
        'drastic-high', ...
        'unbelievably-high', ...
              };
  end;

  cutoffs = prctile(stn.(varname).data, PRCTILES);
  for cix = 2:(length(cutoffs)-1)
    ranges{cix-1} = { FUZZIES{cix} , cutoffs([cix-1 cix]) };
  end;

return;












      % begdix = find(dts >= floor(dts(1) + (perdays*npers)), 'first');





  if ( ~exist('climvar','var') || isempty(climvar) )
    if ( perdays == 30 )
      climvar = [basevar '_monthly_clim'];
    elseif ( perdays == 1 )
      climvar = [basevar '_anom'];
    else
      climvar = [basevar '_weekly_clim'];
    end;
  end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  climvar = [basevar '_monthly_clim'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%







      % % EXPERIMENTING: Various ways of characterizing "trivial anomaly"
      % dat(abs(dat) < 0.5) = nan;
      dat(dat < 0.5) = nan;
      % dat(abs(dat) < 1.0) = nan;
      % dat(dat < 1.0) = nan;









function stns = calc_cum_anom(stns, basevar, anomvar, cumvar, npers, perdays)
%function stns = calc_cum_anom(stns, basevar, anomvar, cumvar, npers, perdays)
%
% Use climatological anomaly of value BASEVAR (a string) for each struct in
% vector STNS, to calculate a cumulative anomaly of that variable. A field
% ANOMVAR must already exist in STNS and a new field CUMVAR is added to each
% struct. Commonly used to calculate, e.g., Degree Heating Weeks on a sea
% temperature: for this, DEFAULTS of ANOMVAR=[BASEVAR '_weekly_anom'],
% CUMVAR=[BASEVAR '_dhw'], NPERS=12 (weeks), and PERDAYS=7 are handy.
%
% Last Saved Time-stamp: <Tue 2010-02-09 21:15:57 Eastern Standard Time gramer>

  if ( ~exist('npers','var') || isempty(npers) )
    npers = 12;
  end;
  if ( ~exist('perdays','var') || isempty(perdays) )
    perdays = 7;
  end;

  if ( ~exist('anomvar','var') || isempty(anomvar) )
    if ( perdays == 30 )
      anomvar = [basevar '_monthly_anom'];
    elseif ( perdays == 1 )
      anomvar = [basevar '_anom'];
    else
      anomvar = [basevar '_weekly_anom'];
    end;
  end;
  if ( ~exist('cumvar','var') || isempty(cumvar) )
    if ( perdays == 30 )
      cumvar = [basevar '_dhm'];
    elseif ( perdays == 1 )
      cumvar = [basevar '_dhd'];
    else
      cumvar = [basevar '_dhw'];
    end;
  end;


  for ix = 1:length(stns)
    stns(ix).(cumvar) = [];

    if ( ~isempty(stns(ix).(anomvar)) && ~isempty(stns(ix).(anomvar).date) )
      dts = stns(ix).(anomvar).date;
      rawdat = stns(ix).(anomvar).data;

      % Round to nearest 0.1oC and blank out trivial anomalies
      dat = roundn(rawdat, -1);

      % % EXPERIMENTING: Various ways of characterizing "trivial anomaly"
      % dat(abs(dat) < 0.5) = nan;
      % dat(dat < 0.5) = nan;
      % dat(abs(dat) < 1.0) = nan;
      dat(dat < 1.0) = nan;

      % begdix = find(dts >= floor(dts(1) + (perdays*npers)), 'first');
      begdix = (perdays*npers);

      stns(ix).(cumvar).date = dts(begdix:end);
      for dix = 1:length(stns(ix).(cumvar).date)
        stns(ix).(cumvar).data(dix,1) = nansum(dat(dix:(dix+begdix-1)));
      end;
    end;
  end;

return;













      [yrs,mos,dys] = datevec(dts);
      jds = datenum(yrs,mos,dys) - datenum(yrs,1,1) + 1;
      wks = floor((jds-1)./perdays) + 1;
      % A few ugly special cases
      if ( perdays == 7 )
        wks(wks > 52) = 52;
      elseif ( perdays == 30 )
        wks(wks > 12) = 12;
      end;





      % dts(abs(dat) < 0.5) = [];
      % dat(abs(dat) < 0.5) = [];








      cumvar = [basevar '_dhw'];



      datmtx = reshape(dat, [length(stns(ix).(cumvar).date) begdix])';
      stns(ix).(cumvar).data = sum(datmtx);





      stns(ix).(cumvar).date = dts(begdix:end);
      for dix = begdix:length(dts)
        stns(ix).(cumvar).data(dix,1) = sum(dat((dix-begdix+1):dix));
      end;






      nvals = length(dts(begdix:end));
      stns(ix).(cumvar).date(1:nvals,1) = dts(begdix:end);






      badix = unique( [ (badix(:)-cutoffidx) , badix(:) , (badix(:)+cutoffidx) ] );
      newts(badix) = nan;









  for ix = 2:length(vars)

    % This is in here to keep from overheating ESRL's servers!
    pause(0.2);

    var = vars{ix};
    varstub = varstubs{ix};
    fld = flds{ix};
    if ( ~isfield(stn, fld) || isempty(stn.(fld)) )
      stn.(fld).date = [];
      stn.(fld).data = [];
    end;

    url = sprintf('%s/%s.gauss.%04d.nc', BASEURL, var, yr);
    nc = mDataset(url);
    if ( isempty(nc) )
      error('Error opening URL "%s"!', url);
    end;
    if ( ~isempty(jds) )
      dat = nc{varstub}(tidx,latix,lonix);
    else
      dat = nc{varstub}(:,latix,lonix);
    end;
    close(nc); clear nc;

    stn.(fld).date(end+1:end+length(dat),1) = dts(:);
    stn.(fld).data(end+1:end+length(dat),1) = dat(:);

  end;











function stn = get_ncep_reanalysis(stn, yr, jds, vars, flds, BASEURL)
%function stn = get_ncep_reanalysis(stn, yr, jds, vars, flds, BASEURL)
%
% Load NCEP Reanalysis of atmosphere-ocean variables VARS (cellstr), directly
% from NOAA/ESRL NCEP Reanalysis netCDF file for YR, range of Julian days JD.
% Arg STN must be a struct with (at least) scalar fields STN.lon and STN.lat:
% fields are added to STN for each value in cellstr FLDS (DEFAULT: variable
% names returned by netCDF THREDDS/dodsC ncep.reanalysis query on VARS.)
%
% For a list of NCEP Reanalysis variables accessible by this routine, see:
%  http://www.esrl.noaa.gov/psd/thredds/catalog/Datasets/ncep.reanalysis/surface_gauss/catalog.html
%
% If optional BASEURL is a valid THREDDS URL, use that dataset instead of
% 'http://www.esrl.noaa.gov/psd/thredds/catalog/Datasets/ncep.reanalysis/surface_gauss'.
%
% EXAMPLE:
%  >> % Get Specific Humidity at z=2m, Molasses Reef, all days in 1987 and 1988
%  >> s.lat =  25.01;
%  >> s.lon = -80.38;
%  >> s = get_ncep_reanalysis(s, 1987, [], 'shum.2m');
%  >> s = get_ncep_reanalysis(s, 1988, [], 'shum.2m');
%  >> mean(s.ncep_reanalysis_shum.data),
%  ans =
%      0.0154
%  >> plot(s.ncep_reanalysis_shum.date,s.ncep_reanalysis_shum.data);
%  >> datetick;
%
% Last Saved Time-stamp: <Sat 2010-02-06 15:39:13 Eastern Standard Time gramer>


  if ( ~iscell(vars) )
    vars = { vars };
  end;
  % CHEAT! To get actual netCDF var names, strip everything after initial Dot
  varstubs = regexprep(vars,'([^.]*)[.].*','$1');

  if ( ~exist('flds','var') || isempty(flds) )
    flds = strcat('ncep_reanalysis_', varstubs);
  end;

  if ( ~exist('BASEURL','var') || ~ischar(BASEURL) )
    BASEURL = 'http://www.esrl.noaa.gov/psd/thredds/dodsC/Datasets/ncep.reanalysis/surface_gauss';
  end;


  lon = stn.lon;
  lon(lon < 0) = 360 + lon(lon < 0);
  lat = stn.lat;

  if ( ~isempty(jds) )
    % NCEP is every six hours, daily
    tidx = (( (jds(1)-1) * 4 ) + 1):(jds(end)*4);
  end;


  ix = 1;
  var = vars{ix};
  varstub = varstubs{ix};
  fld = flds{ix};
  if ( ~isfield(stn, fld) || isempty(stn.(fld)) )
    stn.(fld).date = [];
    stn.(fld).data = [];
  end;

  url = sprintf('%s/%s.gauss.%04d.nc', BASEURL, var, yr);
  nc = mDataset(url);
  if ( isempty(nc) )
    error('Error opening URL "%s"!', url);
  end;
  [ig,latix] = min(abs(nc{'lat'}(:) - lat));
  [ig,lonix] = min(abs(nc{'lon'}(:) - lon));
  if ( ~isempty(jds) )
    t = nc{'time'}(tidx);
    dat = nc{varstub}(tidx,latix,lonix);
  else
    t = nc{'time'}(:);
    dat = nc{varstub}(:,latix,lonix);
  end;
  close(nc); clear nc;

  dts = datenum(1,1,1) + (t./24);

  stn.(fld).date(end+1:end+length(dat),1) = dts(:);
  stn.(fld).data(end+1:end+length(dat),1) = dat(:);


  for ix = 2:length(vars)

    % This is in here to keep from overheating ESRL's servers!
    pause(0.2);

    var = vars{ix};
    varstub = varstubs{ix};
    fld = flds{ix};
    if ( ~isfield(stn, fld) || isempty(stn.(fld)) )
      stn.(fld).date = [];
      stn.(fld).data = [];
    end;

    url = sprintf('%s/%s.gauss.%04d.nc', BASEURL, var, yr);
    nc = mDataset(url);
    if ( isempty(nc) )
      error('Error opening URL "%s"!', url);
    end;
    if ( ~isempty(jds) )
      dat = nc{varstub}(tidx,latix,lonix);
    else
      dat = nc{varstub}(:,latix,lonix);
    end;
    close(nc); clear nc;

    stn.(fld).date(end+1:end+length(dat),1) = dts(:);
    stn.(fld).data(end+1:end+length(dat),1) = dat(:);

  end;

  % This is in here to keep from overheating ESRL's servers!
  pause(0.2);

return;













function [lon,lat,dep] = get_station_coords(stnm)
%function [lon,lat,dep] = get_station_coords(stnm)
%
% Return station coordinates and depth for all known reef monitoring stations
%
% Last Saved Time-stamp: <Wed 2010-01-20 08:19:59 Eastern Standard Time Lew.Gramer>

  switch ( upper(stnm) ),
   case 'LKWF1',    lon = -80.033;    lat = 26.612;    dep = 1.0;
   case 'FWYF1',    lon = -80.100;    lat = 25.590;    dep = 3.0;
   case 'CRYF1',    lon = -80.212;    lat = 25.222;    dep = 11.0;
   case 'MLRF1',    lon = -80.380;    lat = 25.010;    dep = 2.0;
   case 'LONF1',    lon = -80.860;    lat = 24.840;    dep = 1.0;
   case 'TNRF1',    lon = -80.7833;   lat = 24.750;    dep = 6.0;
   case 'SMKF1',    lon = -81.110;    lat = 24.630;    dep = 2.0;
   case 'AMSF1',    lon = -81.520;    lat = 24.525;    dep = 2.0;
   case 'SANF1',    lon = -81.880;    lat = 24.460;    dep = 1.0;
   case 'PLSF1',    lon = -82.770;    lat = 24.690;    dep = 1.0;
   case 'DRYF1',    lon = -82.862;    lat = 24.638;    dep = 1.0;
   case '42003',    lon = -85.594;    lat = 25.966;    dep = 3282.7;

   case 'CMRC3',    lon = -76.139;    lat = 23.791;    dep = 6.0;
   case 'SRVI2',    lon = -64.761;    lat = 17.784;    dep = 6.0;
   case 'LPPR1',    lon = -67.052;    lat = 17.939;    dep = 6.0;
   case 'DBJM1',    lon = -77.42;     lat = 18.47;     dep = 6.0;
   case 'LCIY2',    lon = -80.06;     lat = 19.699;    dep = 5.5;
   % NDBC buoy, Christiansted, St. Croix
   case '41140',    lon = -64.723;    lat = 17.769;    dep = 244.0;

   % Later: add AIMS, AIMS-GBROOS, NMFS-CRED, CRW, etc.
   % OR download/read 'stations.txt' file from ICON data-integration server

   otherwise, error('Do not know coordinates for station "%s"!', stnm);
  end;

return;








% error('This function HAS BEEN replaced by GET_STATION_COORDS.m!');



warning('This function SHOULD BE replaced by GET_STATION_COORDS.m!');





  % Remove any data in VAR that matches gaps in BASEVAR longer than MAXGAP days
  dtdiff = diff(station.(basevar).date);
  gapix = find(dtdiff >= maxgap);
  for ixix = length(gapix):-1:1
    ix = gapix(ixix);
    [ig, preix] = min(abs(station.(var).date - station.(basevar).date(ix) + tol));
    [ig, postix] = min(abs(station.(var).date - (station.(basevar).date(ix+1) + ringlen - tol)));
    % Ensure we only remove elements that are INSIDE the time gap (plus TOL and RINGLEN)...
    if ( preix < (postix-1) )
      station.(var).date((preix+1):(postix-1)) = [];
      station.(var).data((preix+1):(postix-1)) = [];
    end;
  end;





    if ( station.(var).date(preix) < (station.(basevar).date(ix) - tol) )
      preix = preix + 1;
    end;






  % Cp and rho from [Fofonoff and Millard, 1983]
  representative_mean_salinity = 38;
  rho = sw_dens0( repmat(representative_mean_salinity,size(t)), t );
  Cp = sw_cp( representative_mean_salinity, mean(t(:)), stz );
  h = stz;
  stn.heat_flux_term.date = dts;
  stn.heat_flux_term.data = stn.net_heat_flux.data ./ (rho*Cp*h);





function station = filter_gaps(station, basevar, var, tol)
% function station = filter_gaps(station, basevar, var, tol)
%
% Restrict the time series 'station.(var)' to only those dates that lie
% within 'tol' (default 29 min) of some datapoint in 'station.(basevar)'.
%
...
  [ix1,ix2] = intersect_dates(station.(basevar).date(:), station.(var).date(:));
  station.(var).date = station.(var).date(ix2(:));
  station.(var).data = station.(var).data(ix2(:));










  disp(['Start ' mfilename]);
  tic;



  for ixix = 1:length(gapix)
    ix = gapix(ixix);
    [ig, preix] = min(abs(station.(var).date - station.(basevar).date(ix)));
    [ig, postix] = min(abs(station.(var).date - (station.(basevar).date(ix+1) + echolen)));
    datestr(station.(basevar).date([ix ix+1])),
    datestr(station.(var).date([preix postix])),
    disp('--');
  end;

  disp(['End ' mfilename]);
  toc;










  % fid = fopen(fname,'r');



function [clim,climstd,climn] = read_9km_clim(fname)
%function [clim,climstd,climn] = read_9km_clim(fname)

  clim = [];
  climstd = [];
  climn = [];

  if ( ~exist('fname','var') || isempty(fname) )
    fname = '//cygnus/gramer/home/coral/misst/climatology/misst_climatology.dat';
  end;

  type='real  ';
  bchar='float32';
  ibit=4;
  dchar=deblank(bchar);

  arg1=4096;
  arg2=2048;
  arg3=2;

  idim = [arg1 arg2 arg3];
  ilen = prod(idim);
  args = cat(3, arg1, arg2, arg3);

  fid = fopen(fname,'r');
  % fid = fopen(fname,'r','ieee-be');
  if ( fid < 0 )
    error('UNABLE TO OPEN CLIMATOLOGY "%s"!', fname);
  end;

  for ix=1:nargout
    a = fread(fid, ilen, dchar);
    varargout{ix}=reshape(a, idim);
  end;

  fclose(fid);

return;











%  % DEBUG: Try filtering out responses that are TOO extreme?!
%  andts(andat > 0.8) = [];
%  andat(andat > 0.8) = [];



 % DEBUG: Try filtering out responses that are TOO extreme?!
 andts(andat > 1.0) = [];
 andat(andat > 1.0) = [];





 % Find all 'events' - periods when response variable is above cutoff value
 % cutoff = prctile(andat, cutoff_prctile);
 % DEBUG: First, though, try filtering out responses that are TOO extreme?!
 cutoff = prctile(andat(andat <= 1.0), cutoff_prctile);





  disp(dervar);
  disp('intersect 1:'); tic;
  [six,dix] = intersect_dates(sdts, ddts, tol);
  rmix = setdiff(1:length(ddts), dix);
  toc;


  rmdts = ddts(rmix);
  if ( ~isempty(rmix) && ~all(dys == 0) )
    disp('intersect 2:'); tic;
    [brmix, ig] = intersect_dates(ddts, rmdts, dys(1));
    toc;
    brmix = union(brmix, rmix);
  end;










  ix1 = sort(ix1);
  ix2 = sort(ix2);



method = 3
if ( method == 1 )

  while ( ~isempty(preix1) && ~isempty(preix2) )
    if ( (dts1(preix1(1)) + tol) < dts2(preix2(1)) )
      badix = find( (abs(dts1(preix1)-dts2(preix2(1))) <= tol), 1 );
      if ( isempty(badix) )
        break;
      end;
      preix1(1:(badix-1)) = [];
    elseif ( (dts2(preix2(1)) + tol) < dts1(preix1(1)) )
      badix = find( (abs(dts1(preix1(1))-dts2(preix2)) <= tol), 1 );
      if ( isempty(badix) )
        break;
      end;
      preix2(1:(badix-1)) = [];
    else
      % If we hit a single, go for a grand slam!
      maxix = min(length(preix1), length(preix2));
      x1 = preix1(1:maxix);
      x2 = preix2(1:maxix);
      goodix = find( abs(dts1(x1) - dts2(x2)) < tol );
      ix1(end+1:end+length(goodix)) = x1(goodix);
      preix1(goodix) = [];
      ix2(end+1:end+length(goodix)) = x2(goodix);
      preix2(goodix) = [];
    end;
  end;

elseif ( method == 2 )

  ini1 = 1;
  ini2 = 1;
  while ( ini1 <= length(preix1) && ini2 <= length(preix2) )
    if ( abs(dts1(preix1(ini1))-dts2(preix2(ini2))) <= tol )
      % If we hit a single, go for a grand slam!
      maxix = min(length(preix1)-ini1, length(preix2)-ini2);
      x1 = preix1(ini1:ini1+maxix);
      x2 = preix2(ini2:ini2+maxix);
      goodix = find(abs(dts1(x1) - dts2(x2)) < tol);
      ix1(end+1:end+length(goodix)) = x1(goodix);
      ini1 = ini1 + length(goodix);
      ix2(end+1:end+length(goodix)) = x2(goodix);
      ini2 = ini2 + length(goodix);
    else
      if ( dts1(preix1(ini1)) < dts2(preix2(ini2)) )
        ini1 = ini1 + 1;
      else
        ini2 = ini2 + 1;
      end;
    end;
  end;

elseif ( method == 3 )

  % Ugly but somewhat fast code to find all matching dates
  ini1 = 1;
  ini2 = 1;
  while ( ini1 <= length(preix1) && ini2 <= length(preix2) )
    if ( abs(dts1(preix1(ini1))-dts2(preix2(ini2))) <= tol )
      ix1(end+1) = preix1(ini1);
      ix2(end+1) = preix2(ini2);
      ini1 = ini1 + 1;
      ini2 = ini2 + 1;
      % If we hit a double, try to make it - not a triple - a grand slam!
      if ( ini1 <= length(preix1) && ini2 <= length(preix2) && ...
           abs( dts1(preix1(ini1)) - dts2(preix2(ini2)) ) <= tol )
        maxix = min(length(preix1)-ini1, length(preix2)-ini2);
        x1 = preix1(ini1:ini1+maxix);
        x2 = preix2(ini2:ini2+maxix);
        goodix = find(abs(dts1(x1) - dts2(x2)) < tol);
        ix1(end+1:end+length(goodix)) = x1(goodix);
        ini1 = ini1 + length(goodix);
        ix2(end+1:end+length(goodix)) = x2(goodix);
        ini2 = ini2 + length(goodix);
      end;
    else
      if ( dts1(preix1(ini1)) < dts2(preix2(ini2)) )
        ini1 = ini1 + 1;
      else
        ini2 = ini2 + 1;
      end;
    end;
  end;

elseif ( method == 4 )

  % SLOOOWOWOWOW! Need a clever MATLAB-ish way to do this...
  for ixi = 1:length(preix1)
    ix = preix1(ixi);
    goodix = find(abs(dts1(ix)-dts2) <= tol);
    if ( ~isempty(goodix) )
      ix1(end+1) = ix;
      ix2(end+1) = goodix(1);
    end;
  end;

end;











    badix = find( (dts1(preix1) + tol) < dts2(preix2(1)) );
    if ( (dts1(preix1(1)) + tol) < dts2(preix2(1)) )
      preix1(1) = [];
    elseif ( (dts2(preix2(1)) + tol) < dts1(preix1(1)) )
      preix2(1) = [];
    else
      % If we hit a single, go for a grand slam!
      maxix = min(length(preix1), length(preix2));
      x1 = preix1(1:maxix);
      x2 = preix2(1:maxix);
      goodix = find( abs(dts1(x1) - dts2(x2)) < tol );
      ix1(end+1:end+length(goodix)) = x1(goodix);
      preix1(goodix) = [];
      ix2(end+1:end+length(goodix)) = x2(goodix);
      preix2(goodix) = [];
    end;
  end;







    if ( abs( dts1(preix1(1)) - dts2(preix2(1)) ) <= tol )





  % First get all the low-hanging fruit - exact timestamp matches
fprintf(1, 'Exact\n');
tic;
  [exactix1,exactix2] = find(ismember(dts1(preix1), dts2(preix2)));
  ix1 = preix1(exactix1); ix2 = preix2(exactix2);
  preix1(exactix1) = []; preix2(exactix2) = [];
toc;
fprintf(1,'Found %d\n', length(ix1));

fprintf(1, 'Searching...\n');
tic;
method = 2;
if ( method == 1 )

  % SLOOOWOWOWOW! Need a clever MATLAB-ish way to do this...
  for ixi = 1:length(preix1)
    ix = preix1(ixi);
    goodix = find(abs(dts1(ix)-dts2) <= tol);
    if ( ~isempty(goodix) )
      ix1(end+1) = ix;
      ix2(end+1) = goodix(1);
    end;
  end;

else


end;
toc;
fprintf(1,'Total found %d\n', length(ix1));


















 disp('blank an:'); tic;
 toc; disp('blank fc1:'); tic;




    disp('union:'); tic;
    toc;



 for ix = 1:nvars
   station = verify_variable(station, fcfld{ix});
   % station = blank_derived_var(station, 'sea_t', anfld);
   % HACK!!!!!!!
   station.(fcfld{ix}).data(station.(fcfld{ix}).data < 0.1) = nan;
 end;





 anmn = nanmean(station.(anfld).data(ix1));
 ansd = nanstd(station.(anfld).data(ix1));



 for ix = 1:nvars
   fcdts{ix} = station.(fcfld{ix}).date(ix2);
   fcdat{ix} = station.(fcfld{ix}).data(ix2);

   % Normalize our independent variables
   fcmn(ix) = nanmean(station.(fcfld{ix}).data(ix2));
   fcsd(ix) = nanstd(station.(fcfld{ix}).data(ix2));
   fcdat{ix} = (fcdat{ix} - fcmn(ix)) ./ fcsd(ix);
 end;



 % Re-de-normalize our independent variables
 for ix = 1:nvars
   fcdat{ix} = (fcdat{ix} .* fcsd(ix)) + fcmn(ix);
 end;




 interpmthd = 'pchip';
 sominit = 'lininit';







method = 2;
tic;
if ( method == 1 )

  % SLOOOWOWOWOW! Need a clever MATLAB-ish way to do this...
  for ixi = 1:length(preix1)
    ix = preix1(ixi);
    goodix = find(abs(dts1(ix)-dts2) <= tol);
    if ( ~isempty(goodix) )
      ix1(end+1) = ix;
      ix2(end+1) = goodix(1);
    end;
  end;

else

end;
toc;






    % DEBUG:
    if ( preix1(ini1) >= 116980 )
      keyboard;
    end;






function [ix1,ix2] = intersect_dates(dts1, dts2, tol)
%function [ix1,ix2] = intersect_dates(dts1, dts2, tol)
%
% Return indices of all the elements in two series of timestamps (datenums)
% 'dts1' and 'dts2', that match: a "matching" timestamp is defined as one
% where the difference between the elements of the two series is less than
% some tolerance level 'tol'. DEFAULT 'tol' is 29 minutes, i.e., 29/(24*60).
%
% Last Saved Time-stamp: <Sun 2009-08-09 00:19:57 Eastern Daylight Time gramer>

  ix1 = [];
  ix2 = [];

  if ( any(diff(dts1) <= 0) || any(diff(dts2) <= 0) )
    error('Both date series must be monotonically increasing!');
  end;
  if ( ~exist('tol','var') || isempty(tol) )
    % Default tolerance for timestamp matching is 29 minutes
    tol = 29.0/(24.0*60.0);
  end;
  if ( tol <= 0 )
    error('Tolerance argument should be greater than zero!');
  end;

  preix1 = find( (dts2(1) <= dts1) & (dts1 <= dts2(end)) );
  preix2 = find( (dts1(1) <= dts2) & (dts2 <= dts1(end)) );

  dts1 = dts1(preix1);
  dts2 = dts2(preix2);

  % SLOOOWOWOWOW! Need a clever MATLAB-ish way to do this...
tic;
%   ix1 = repmat(nan, size(preix1));
%   ix2 = repmat(nan, size(preix2));
%   curix = 1;
%   for ixi = 1:length(preix1)
%     ix = preix1(ixi);
%     goodix = [];
%     %goodix = find(abs(dts1(ix)-dts2) <= tol);
%     if ( ~isempty(goodix) )
%       ix1(curix) = ix;
%       ix2(curix) = goodix(1);
%       curix = curix + 1;
%     end;
%   end;
%   ix1(isnan(ix1)) = [];
%   ix2(isnan(ix2)) = [];
  for ix = 1:length(dts1)
    goodix = find(abs(dts1(ix)-dts2) <= tol);
    if ( ~isempty(goodix) )
      ix1(end+1) = ix;
      ix2(end+1) = goodix(1);
%       dts2(1:goodix(1)) = [];
    end;
  end;
%   for ixi = 1:length(preix1)
%     ix = preix1(ixi);
%     goodix = find(abs(dts1(ix)-dts2(preix2)) <= tol);
%     if ( ~isempty(goodix) )
%       ix1(end+1) = ix;
%       ix2(end+1) = preix2(goodix(1));
%     end;
%   end;
toc;

  ix1 = sort(ix1);
  ix2 = sort(ix2);

return;









        flds = fieldnames(stn);
        for ix = 1:length(flds)
          fld = flds{ix};
          if ( isfield(station,fld) )
            warning('Field %s of station struct was overwritten!', 
          else
            station.(fld) = stn.(fld);
          end;
        end;



  % Retrieve data from a path relative to this M-file's local directory
  [pathroot, ig, ig, ig] = fileparts(mfilename('fullpath'));
  if ( ~exist('datapath', 'var') || isempty(datapath) )
    datapath = fullfile(pathroot, 'data', '');
  end;



  STATIONS.codes = ...
      { ...
          'CMRC3', ...
          'SRVI2', ...
          'LPPR1', ...
          'DBJM1', ...
          'LCCI1', ...
          ...
          'LKWF1', ...
          'FWYF1', ...
          'MLRF1', ...
          'LONF1', ...
          'SMKF1', ...
          'SANF1', ...
          '42003', ...
      };
  STATIONS.lons = ...
      [ ...
          -76.139, ...
          -64.761, ...
          -67.052, ...
          -77.420, ...
          -80.060, ...
          ...
          -80.033, ...
          -80.100, ...
          -80.380, ...
          -80.860, ...
          -81.110, ...
          -81.880, ...
          -85.594, ...
      ];
  STATIONS.lats = ...
      [ ...
          23.791, ...
          17.784, ...
          17.939, ...
          18.470, ...
          19.699, ...
          ...
          26.612, ...
          25.590, ...
          25.010, ...
          24.840, ...
          24.630, ...
          24.460, ...
          25.966, ...
      ];







  % lkwf1,26.612,-80.033
  % fwyf1,25.590,-80.100
  % cryf1,25.222,-80.212
  % mlrf1,25.010,-80.380
  % lonf1,24.840,-80.860
  % tnrf1,24.750,-80.78333
  % smkf1,24.630,-81.110
  % amsf1,24.525,-81.520
  % sanf1,24.460,-81.880
  % plsf1,24.690,-82.770
  % 42003,25.966,-85.594







  station = stn2; save(['data/' stnm2 '-ndbc.mat'], 'station'); station = []; clear station;




function [dts1,dat1,dts2,dat2] = extract_matching_dates(dts1,dat1,dts2,dat2)
  [ix1, ix2] = intersect_dates(dts1, dts2);
  dts1 = dts1(ix1); dat1 = dat1(ix1);
  dts2 = dts2(ix2); dat2 = dat2(ix2);
return;







function stn = station_heat_flux(stn,wfld,afld,qfld,pfld,tfld)
%function stn = station_heat_flux(stn,wfld,afld,qfld,pfld,tfld)
%
% Use Fairall et al. (1996) method to calculate wind stress, and sensible and
% latent surface heat flux, from station data struct 'stn', using in situ
% sensor data named in field names 'wfld' (wind speed), 'afld' (air temp),
% 'qfld' (RELATIVE humidity), 'pfld' (air pressure), and 'tfld' (sea temp).
% Arg 'qfld' may optionally be a scalar value, or the name of a dew-point
% temperature field (as in many SEAKEYS stations) from which relative humidity
% can be calculated. Fields 'wind_stress', 'sensible_flux', 'latent_flux' are
% added to 'stn', each with 'date' subfields equal to the intersection of all
% the input fields' timestamps (time matching tolerance defaults to 29 min).
%
% Last Saved Time-stamp: <Mon 2009-07-13 08:10:33 Eastern Daylight Time Lew.Gramer>

  % Get height/depth [m] of meteorology and 'shallow' sea temperature sensors
  [wz,az,pz,tz] = station_instrument_heights(stn.station_name);

  wd = stn.(wfld).date; w = stn.(wfld).data;
  ad = stn.(afld).date; a = stn.(afld).data;
  pd = stn.(pfld).date; p = stn.(pfld).data;
  td = stn.(tfld).date; t = stn.(tfld).data;

  [wd,w,ad,a] = extract_matching_dates(wd,w,ad,a);
  [wd,w,pd,p] = extract_matching_dates(wd,w,pd,p);
  [wd,w,td,t] = extract_matching_dates(wd,w,td,t);

  % Convert station air pressure to sea-level
  [ad,a,pd,p] = extract_matching_dates(ad,a,pd,p);
  p = p ./ exp(-pz./(a + 273.15));

  % Figure out something for relative humidity
  if ( isnumeric(qfld) )
    q = repmat(qfld, size(w));

  else

    % Is the named field dewpoint temperature or relative humidity?
    if ( strfind(lower(qfld), 'dew') )
      [td,t,qd,d] = extract_matching_dates(td,t,stn.(qfld).date,stn.(qfld).data);
      q = dewp_to_relhumid(t,d);

      stn.relative_humidity_calc.date = qd;
      stn.relative_humidity_calc.data = q;
    else
      qd = stn.(qfld).date;
      q = stn.(qfld).data;
      [wd,w,qd,q] = extract_matching_dates(wd,w,qd,q);
    end;

  end;

  result = hfbulktc(w,wz,a,az,q,az,p,t);

  stn.wind_stress_calc.date = wd;
  stn.wind_stress_calc.data = result(1);
  stn.sensible_flux_calc.date = wd;
  stn.sensible_flux_calc.data = result(2);
  stn.latent_flux_calc.date = wd;
  stn.latent_flux_calc.data = result(4);

return;


function [dts1,dat1,dts2,dat2] = extract_matching_dates(dts1,dat1,dts2,dat2)
  [ix1, ix2] = intersect_dates(dts1, dts2);
  dts1 = dts1(ix1); dat1 = dat1(ix1);
  dts2 = dts2(ix2); dat2 = dat2(ix2);
return;

% August-Roche-Magnus approximation
function d = relhumid_to_dewp(t,q)
  a = 17.271; b = 237.7;
  gamma = ((a.*t)/(b+t)) + log(q./100);
  d = b.*gamma/(a - gamma);
return;

% August-Roche-Magnus approximation - inverted
function q = dewp_to_relhumid(t,d)
  a = 17.271; b = 237.7;
  gamma = (a.*d) / (d + b);
  q = 100.*exp(gamma - ((a.*t)/(b + t)));
return;











  multiplot_datetick({tdts,bdts},{tdat,bdat},'Good data');
  keyboard;






function [hlines, haxes, hfig] = multiplot_datetick(X, Y, ttl, xlbl, ylbl, xlm, datetick_fmt)
  for ix = 1:length(ylbl)
    ylbl{ix} = strrep(lower(ylbl{ix}),'_','\_');
  end;






  ddts = tdts;
  ddat = repmat(nan, size(ddts));
  kdts = tdts;
  kdat = repmat(nan, size(kdts));
  for dix = 1:length(ddts)

    % Only pick days with reasonable amount of data
    if ( length(tix) >= 16 )
      mbdat(dix) = mean(bdat(bix));
      ddat(dix) = mtdat(dix) - mbdat(dix);
      if ( mtdat(dix) > 0 && mbdat(dix) > 0 )
        kdat(dix) = -log(mbdat(dix) / mtdat(dix));
      end;
    end;
  end;









function stn = calc_kd_daily(stn, topfld, btmfld, zfld)
%function stn = calc_kd_daily(stn, topfld, btmfld, zfld)
%
% From the ICON/G2 definitions/documentation for the 'Kd' functions:
% // Used by Kd attenuation coefficient: assumes "surface" and "shallow" light sensors have
% //  a CONSTANT 1m of water between them. Iz = I0 * exp(-Kd*z), Kd = -(1/z)ln(Iz/I0)
% negln1m(val) = ( if ((val)<=0) then -9.0 else (- ln(val)) )
% // Used by Kd attenuation coefficient: assumes "shallow" and "deep" light sensors 2m apart
% //negln2m(val) = ( if (val<=0) then -9.0 else - ln((val)/2) ) - According to E. Stabenau - ?
% negln2m(val) = ( if ((val)<=0) then -9.0 else (- (ln(val)/2.0)) )

  stn = verify_variable(stn, topfld);
  stn = verify_variable(stn, btmfld);
  if ( ischar(zfld) )
    stn = verify_variable(stn, zfld);
  end;

  tdts = stn.(topfld).date(stn.(topfld).data >= 0);
  tdat = stn.(topfld).data(stn.(topfld).data >= 0);
  bdts = stn.(btmfld).date(stn.(btmfld).data >= 0);
  bdat = stn.(btmfld).data(stn.(btmfld).data >= 0);

  [tix,bix] = intersect_dates(tdts, bdts);
  tdts = tdts(tix); tdat = tdat(tix);
  bdts = bdts(bix); bdat = bdat(bix);

  % Simplify life: just pick the 1700UT (i.e., near-peak) sample for each day!
  hrs = round((tdts - floor(tdts)) * 24);
  peakix = find(hrs == 17);

  tdts = tdts(peakix);  tdat = tdat(peakix);
  bdts = bdts(peakix);  bdat = bdat(peakix);

  newfld = ['kd_1_day_' topfld '_' btmfld];
  % Make sure we start with a clean slate...
  if ( isfield(stn, newfld) )
    stn = rmfield(stn, newfld);
  end;
  stn.(newfld).date = tdts;
  stn.(newfld).data = repmat(-9.0, size(stn.(newfld).date));

  ix = find(tdat > 0);
  val = stn.(newfld).data;
  val(ix) = bdat(ix) ./ tdat(ix);

  if ( isnumeric(zfld) )
    z = zfld;
    ix = find(val > 0);
    stn.(newfld).data(ix) = -log(val(ix)) ./ z;
  else
    [nix,zix] = intersect_dates(stn.(newfld).date, stn.(zfld).date);
    z = stn.(zfld).data;
    ix = find(val(nix) >= 0);
    stn.(newfld).data(nix(ix)) = -log(val(nix(ix))) ./ z(zix(ix));
  end;

return;














  for ix = 2:length(haxes)
    set(haxes(ix), 'XTickLabel', []);
  end;




  % Simplify life: just pick the 1700UT (i.e., near-peak) sample for each day!
  tdvec = datevec(tdts);
  peakix = find(tdvec(:,4) == 17);







  [tix,bix] = intersect_dates(stn.(topfld).date, stn.(btmfld).date);





function stn = calc_kd_daily(stn, topfld, btmfld, zfld)
%function stn = calc_kd_daily(stn, topfld, btmfld, zfld)
%
% From the ICON/G2 definitions/documentation for the 'Kd' functions:
% // Used by Kd attenuation coefficient: assumes "surface" and "shallow" light sensors have
% //  a CONSTANT 1m of water between them. Iz = I0 * exp(-Kd*z), Kd = -(1/z)ln(Iz/I0)
% negln1m(val) = ( if ((val)<=0) then -9.0 else (- ln(val)) )
% // Used by Kd attenuation coefficient: assumes "shallow" and "deep" light sensors 2m apart
% //negln2m(val) = ( if (val<=0) then -9.0 else - ln((val)/2) ) - According to E. Stabenau - ?
% negln2m(val) = ( if ((val)<=0) then -9.0 else (- (ln(val)/2.0)) )

  stn = verify_variable(stn, topfld);
  stn = verify_variable(stn, btmfld);
  stn = verify_variable(stn, zfld);

  tdts = stn.(topfld).date(stn.(topfld).data >= 0);
  tdat = stn.(topfld).data(stn.(topfld).data >= 0);
  bdts = stn.(btmfld).date(stn.(btmfld).data >= 0);
  bdat = stn.(btmfld).data(stn.(btmfld).data >= 0);

  [tix,bix] = intersect_dates(stn.(topfld).date, stn.(btmfld).date);
  tdts = tdts(tix); tdat = tdat(tix);
  bdts = bdts(bix); bdat = bdat(bix);

  firstix = find((tdat < 50 & bdat < 50), 1);
  tdts = tdts(firstix:end);
  tdat = tdat(firstix:end);
  bdts = bdts(firstix:end);
  bdat = bdat(firstix:end);

  badtix = find(diff(tdts) <= ((1/24.0) - eps));
  tdts(badtix) = [];
  tdat(badtix) = [];
  bdts(badtix) = [];
  bdat(badtix) = [];

  [tdts,tdat] = gap_expand(tdts, tdat);
  [bdts,bdat] = gap_expand(bdts, bdat);

  ndays = floor(length(tdat) / 24);
  endix = ndays * 24;

  tdts = tdts(1:endix);  tdat = tdat(1:endix);
  bdts = bdts(1:endix);  bdat = bdat(1:endix);

  tdat = reshape(tdat, [24 ndays]);
  bdat = reshape(bdat, [24 ndays]);

  [tdat, tix] = nanmax(tdat); tdts = tdts(tix);
  [bdat, bix] = nanmax(bdat); bdts = bdts(bix);


  newfld = ['kd_1_day_' topfld '_' btmfld];
  % Make sure we start with a clean slate...
  if ( isfield(stn, newfld) )
    stn = rmfield(stn, newfld);
  end;
  stn.(newfld).date = tdts;
  stn.(newfld).data = repmat(-9.0, size(stn.(newfld).date));

  ix = find(tdat ~= 0);
  val = stn.(newfld).data;
  val(ix) = bdat(ix) ./ tdat(ix);

  [nix,zix] = intersect_dates(stn.(newfld).date, stn.(zfld).date);
  z = stn.(zfld).data';

  ix = find(val(nix) >= 0);
  stn.(newfld).data(nix(ix)) = -log(val(nix(ix))) ./ z(zix(ix));

return;












    disp(['Reloading station data from MAT file ' fname '...']);
    x = load(fname, 'station');
    stn = x.station;
    clear x;
%%%%
    stn.station_name = stanm;
    station = stn;
    save(fname, 'station');
    station = [];
    clear station;
%%%%
    clear fname;






ax = suplabel(sprintf('PCA Modes (%s %s), 1991-2009 %s, %d-day frames (%3.1f%% of %s)', ...
                      upper(stnam), strrep(upper(fldnm),'_','\_'), period, ...
                      ndays, totpct, '\sigma^2'), ...
              't', [0.075 0.075 0.850 0.870]);




    station = [];
    yrs = 1981:2009;
    for ix = 1:length(yrs)
      yr = yrs(ix);
      fname = sprintf('%s/%sh%d.txt', datapath, lower(stnam), yr);
      if ( exist(fname, 'file') )
        station = load_ndbc_data(station, fname);
      end;
    end;





  preix1 = find( (abs(dts2(1)-dts1) <= tol) & (abs(dts1-dts2(end)) <= tol) );
  preix2 = find( (abs(dts1(1)-dts2) <= tol) & (abs(dts2-dts1(end)) <= tol) );





  if ( ~exist('datpath', 'var') )
    datpath = './data/';
  end;
  if ( ~exist(datpath, 'dir') )
    error('The data source path "%s" does not exist!', datpath);
  end;







  [sea_yrs,mos,dys] = datevec(station.sea_t.date);
  figure;
  set(gcf, 'units','normalized', 'outerposition',[0 0 1 1]);
  boxplot(station.sea_t.data, sea_yrs, 'notch','on', 'whisker',1.5);
  xlabel('Year'); title([upper(stnm) ' T_s_e_a']);
  ylim([5 35]);
  print('-dpng', [stnm '-mpo624-year-sea_t-boxplot.png']);






%%% Do some monthly boxplots - Tsea, Tair
[yrs,air_mos,dys] = datevec(station.air_t.date);
figure;
set(gcf, 'units','normalized', 'outerposition',[0 0 1 1]);
boxplot(station.air_t.data, air_mos, 'notch','on', 'whisker',1.5);
xlabel('Year Month'); title([upper(stnm) ' T_a_i_r']);
ylim([5 35]);
print('-dpng', [stnm '-mpo624-month-air_t-boxplot.png']);

[yrs,sea_mos,dys] = datevec(station.sea_t.date);
figure;
set(gcf, 'units','normalized', 'outerposition',[0 0 1 1]);
boxplot(station.sea_t.data, sea_mos, 'notch','on', 'whisker',1.5);
xlabel('Year Month'); title([upper(stnm) ' T_s_e_a']);
ylim([5 35]);
print('-dpng', [stnm '-mpo624-month-sea_t-boxplot.png']);









mo_air = {};
mo_sea = {};
for mo = 1:12
  [y,m,d] = datevec(station.air_t.date);
  fidx = find(m == mo);
  yr_air = unique([yr_air y]);
  mo_air{mo} = station.air_t.data(fidx);
  [y,m,d] = datevec(station.sea_t.date);
  fidx = find(m == mo);
  mo_sea{mo} = station.sea_t.data(fidx);
end;
clear y; clear m; clear d;
boxplot(mo_air, 'notch','on', 'whisker',6.0, 





station = verify_variable(station, 'sea_t_qc_1_day_deviation_3_day_average');
station = verify_variable(station, 'sea_t_qc_3_day_average_0_day_asof_diff_sea_t_qc_1_day_minimum');

%%% Do some monthly boxplots - Tsea, Tair


fidx = 0;
clear multidat;
clear fldabbr;

fidx = fidx + 1;
multidts{fidx} = station.sea_t_qc.date;
multidat{fidx} = station.sea_t_qc.data;
fldabbr{fidx} = 'T_s_e_a';

fidx = fidx + 1;
multidts{fidx} = station.sea_t_qc_1_day_deviation_3_day_average.date;
multidat{fidx} = station.sea_t_qc_1_day_deviation_3_day_average.data;
fldabbr{fidx} = '\mu_3_d\sigma_1_d(T_s_e_a)';

fidx = fidx + 1;
multidts{fidx} = station.sea_t_qc_3_day_average.date;
multidat{fidx} = station.sea_t_qc_3_day_average.data;
fldabbr{fidx} = '\mu_3_d(T_s_e_a)';

fidx = fidx + 1;
multidts{fidx} = station.sea_t_qc_3_day_average_0_day_asof_diff_sea_t_qc_1_day_minimum.date;
multidat{fidx} = station.sea_t_qc_3_day_average_0_day_asof_diff_sea_t_qc_1_day_minimum.data;
fldabbr{fidx} = 'anom^m^i^n_3_d(T_s_e_a)';








% station = verify_variable(station, 'sea_t_qc_3_day_maximum_0_day_current_diff_sea_t_qc_3_day_average');




  diary([nm '-stats.txt']);
  diary off;




%   % Datetick labels spill over edges on all plots above lowest one
%   for ix = 1:(length(haxes)-1)
%     set(haxes(ix),'XTickLabel','')
%   end;





    % [ig,ig,YYYY,MM,DD,hh,mm] = parseusfurl(pname);
    [yr, mo, dy, hr, mn, sc] = datevec(dts(ix))


    title(sprintf( 'BMU = %d : %04d%02d%02d %02d:%02d', ...
                   sm.bmus(ix), yr, mo, dy, hr, mn ));




nhrs = 168;
mxel = 10000;
nmel = mxel - rem(mxel,nhrs);
x = rand([1 nmel]);
x = reshape(x, [nhrs (numel(x)/nhrs)])';


ix = reshape([1:numel(x)]', [nhrs (numel(x)/nhrs)])';

[ig,mx] = max(x,[],2);

modeix = ceil(2*nhrs/3) + 1;

% First and last frames are special cases - we don't want to zero-pad these
% two frames, nor to have to invent data! So limit how they can shift...
[ig, mx(1)] = max(x(1,modeix:end),[],2);
mx(1) = mx(1) + modeix - 1;
[ig, mx(end)] = max(x(end,1:(modeix-1)),[],2);

shiftix = repmat((modeix - mx), [1 size(x,2)]);

resultix = ix - shiftix;






XXXXXXXXXXXXXXXXXXXX

[ig,ix] = max(fda,[],2);
ix2 = repmat((ceil(nhrs/2) - ix;
newix = repmat([1:(nhrs*2)], [size(fda,1) 1]);

XXXXXXXXXXXXXXXXXXXX

station = verify_variable(station, [varname '_7_day_lp']);








%%%% ???
% %     set(gca, 'clim', [1.1 3.5]);
%     set(gca, 'clim', [-2.1 4.7]);
%%%% ???




%%%% ???
%     set(gca, 'clim', [-30 +70]);
%%%% ???






   [yy mm dd] = datevec(fdt(:, 1));
   jdays = datenum(yy,mm,dd) - datenum(yy,1,1) + 1;






%%%% ???
    set(gca, 'clim', [1.2 3.2]);
%%%% ???




  fprintf('Missing %d weekly means (%d total)\n', length(find(rmrows)), size(sst,1));




% ???
%      framedts(:,2:(end-1)) = [];



  print('-dtiff', '-r300', ...
        fullfile(figspath, ...
                 sprintf('%s-PCA-NDBC-%s-%dx%d-%s.tiff', stanm, ...
                         period, mapdims(1), mapdims(2), framestr)));



  % % Linear init. - faster SOM runs, but MORE MEMORY (and must coddle NaNs)
  % sm = som_make(goodsst, 'msize',mapdims, 'lininit', 'shape','sheet', ...
  %               'mask',double(~landmask), 'training','long', 'neigh','ep');





 saveas(gcf, [pfname '.fig']);
   saveas(gcf, [pfname '.fig']);




%   % Add annotations
%   axes(axs(end));
%   hold on;
%   for ix = 1:length(frameix)
%     annotline(framedts(frameix(ix)));
%   end;
%   hold off;







  multidts{end+1} = framedts(frameix);
  multidat{end+1} = repmat(1, size(framedts(frameix)));








format short; format compact;

more_status = get(0, 'More');
more('off');

% Done
more(more_status);








 % covm(~isfinite(covm)) = 0;




%%%% ???
   cmin = []; cmax = [];





 % For each variable, remove TIME AVERAGE over the WHOLE RECORD
 for fidx = 1:length(fldnms)
   begix = ((fidx-1) * nhrs) + 1;
   endix = begix + nhrs - 1;
   varmean{fidx} = mean(mean(fdat(:,begix:endix)));
   fanom(:,begix:endix) = fdat(:,begix:endix) - varmean{fidx};
 end;




% fldnms{1} = 'wind1_speed_40_hour_lp';
% fldnms{2} = 'wind1_u_40_hour_lp_1_day_deviation_sum_wind1_v_40_hour_lp';
% fldnms{3} = 'air_t_qc_40_hour_lp_1_day_deviation';
% fldnms{4} = 'sea_t_qc_40_hour_lp_1_day_deviation';






% fldnms{1} = 'sea_t_qc_1_day_deviation_3_day_average';
% fldnms{2} = 'air_t_qc_1_day_deviation_3_day_average';
% fldnms{3} = 'wind1_u_3_day_average';
% fldnms{4} = 'wind1_v_3_day_average';

% fldnms{1} = 'sea_t_qc';
% fldnms{2} = 'air_t_qc';
% fldnms{3} = 'wind1_u';
% fldnms{4} = 'wind1_v';

% fldnms{1} = 'sea_t_qc_40_hour_lp';
% fldnms{2} = 'air_t_qc_40_hour_lp';
% fldnms{3} = 'wind1_u_40_hour_lp';
% fldnms{4} = 'wind1_v_40_hour_lp';

% fldnms{1} = 'sea_t_qc_1_day_deviation_40_hour_lp';
% fldnms{2} = 'air_t_qc_1_day_deviation_40_hour_lp';
% fldnms{3} = 'wind1_u_40_hour_lp';
% fldnms{4} = 'wind1_v_40_hour_lp';

% fldnms{1} = 'wind1_v_40_hour_lp';
% fldnms{2} = 'wind1_u_40_hour_lp';
% fldnms{3} = 'air_t_qc_40_hour_lp_1_day_deviation';
% fldnms{4} = 'sea_t_qc_40_hour_lp_1_day_deviation';







    % ???
    % Add in the frame just BEFORE each date in 'keepdts'
    keepix = unique([keepix, (keepix - 1)]);
    keepix(keepix == 0) = [];






  % % Remove the daily temperature cycle from each day of time series
  % clim = [ ...
  %     0.030 0.036 0.058  0.106 0.169 0.246 ...
  %     0.314 0.367 0.410  0.426 0.427 0.398 ...
  %     0.355 0.295 0.263  0.234 0.207 0.177 ...
  %     0.142 0.112 0.096  0.070 0.051 0.036 ...
  %        ];
  % clim = repmat(clim, [??? ??? nrows ndays]);
  % dat = dat - clim;

  % % Remove diurnal cycle from interpolated time series with a bandstop filter
  % if ( strcmpi(fldnm, 'sea_t') || strcmpi(fldnm, 'tide') )
  %   disp('Removing Fourier components of 10-14 and 18-30 hour periods');
  %   [B,A] = butter(filtOrder, [(2/14) (2/10)], 'stop');
  %   dat = filtfilt(B, A, dat);
  %   [B,A] = butter(filtOrder, [(2/30) (2/18)], 'stop');
  %   dat = filtfilt(B, A, dat);
  % end;










   fdat(:,begix:endix) - mean(mean(fdat(:,begix:endix)));


 [nrows, ncols] = size(fdat);
 for fidx = 1:length(fldnms)
   for rix = 1:nrows
     fanom(rix,1:ncols) = fdat(rix,:) - nanmean(fdat(rix,:));
   end;
 end;





 %%%% ???
 cmin = repmat(0, size(cmn)); cmax = cmx;




      yr = 2009;
      fname = sprintf('%s/%sh%d.txt', datapath, lower(stnam), yr);
      if ( exist(fname, 'file') )
        station = load_ndbc_data(station, fname);
      end;
    station.station_name = lower(stnam);
    save(matfname, 'station');









  if ( iscell(period) )

    varnam = period{1};
    varrng = period{2};

  else

  end;






  clear datapath yrs ix yr fname matfname;
clear more_status;
clear graphvis_status;




&& ...
         numel(cmin) == numel(hl) && numel(cmax) == numel(hl) )



 print('-dtiffnocompression', '-adobecset', pfname);




  yposn = repmat(strvcat('left','right'), [6 1]);
    set(ha(ix), 'YAxisLocation', yposn(ix,:));



    %ylim(ax(1), [-0.4 0.4]);




  x = 1:nhrs;
  ax = plotyy(x, cdbk(:,[1]), x, cdbk(:,[2 3 4]));
  xlim(ax(1), [1 nhrs]);
  xlim(ax(2), [1 nhrs]);
%   ylim(ax(1), [-0.02 0.02]);


     x = 1:nhrs;
     pcadats = reshape(pcadat(ix,:), [nhrs length(fldnms)]);
     ax = plotyy(x, pcadats(:,[1]), x, pcadats(:,[2 3 4]));
     xlim(ax(1), [1 nhrs]);
     xlim(ax(2), [1 nhrs]);





   case 'winter',
    %wks = [1:12 50:52];
    jdays = [1:84 (344-ndays):366];
   case 'spring',
    %wks = [11:25];
    jdays = [(71-ndays):175];
   case 'summer',
    %wks = [24:38];
    jdays = [(162-ndays):266];
   case 'autumn',
    %wks = [37:51];
    jdays = [(253-ndays):357];



  % ylim(ax(2), [-2 2]);




  for fidx = 1:length(fldnms)
    begcol = ((fidx-1)*nhrs) + 1;
    endcol = begcol + (nhrs-1);
    cdbk{fidx} = sm.codebook(ix,begcol:endcol);
%     cdbk{fidx} = cdbk{fidx} / std(cdbk{fidx});
%     cdbk{fidx} = cdbk{fidx} - min(cdbk{fidx});
    plot(1:nhrs, cdbk{fidx}, color_order(fidx));
  end;

  xlim([1 nhrs]);



  for fidx = 1:length(fldnms)
    begcol = ((fidx-1)*nhrs) + 1;
    endcol = begcol + (nhrs-1);
    pcadats{fidx} = pcadat(ix,begcol:endcol);
%     pcadats{fidx} = pcadats{fidx} / std(pcadats{fidx});
%     pcadats{fidx} = pcadats{fidx} - min(pcadats{fidx});
    plot(1:nhrs, pcadats{fidx}, color_order(fidx));
  end;

  xlim([1 nhrs]);










% fldnms{1} = 'sea_t_qc_1_day_deviation_3_day_average';
% % fldnms{2} = 'wind1_speed_3_day_average';
% % fldnms{2} = 'wind1_u_3_day_average';
% % fldnms{3} = 'wind1_v_3_day_average';
% fldnms{2} = 'air_t_qc_1_day_deviation_3_day_average';
% fldnms{3} = 'wind1_u_7_day_deviation_sum_wind1_v';




% ax = suplabel(sprintf('SOM "Modes" (%s %s), 1991-2009 %s, %d-day frames (%d frames)', ...
%                       upper(stnam), strrep(upper(fldnmstr),'_','\_'), period, ...
%                       ndays, size(fdat,1)), ...
%               't', [0.075 0.075 0.850 0.840]);




% ax = suplabel(sprintf('PCA Modes (%s %s), 1991-2009 %s, %d-day frames (%3.1f%% of %s)', ...
%                       upper(stnam), strrep(upper(fldnmstr),'_','\_'), period, ...
%                       ndays, totpct, '\sigma^2'), ...
%               't', [0.075 0.075 0.850 0.840]);







  %%%% 
  %%%% ??? NEED TO CLEAR ANY VARIABLES DERIVED FROM THIS QC VARIABLE ALSO!
  %%%% 
  flds = fieldnames(stn);
  findix = strmatch(qcname, flds);
  for fidx = 1:length(findix)
    if ( ~strcmpi(flds{fidx}, qcname) )
      fprintf('Should recalculate %s here also!\n', flds{fidx});
    end;
  end;






              't', [0.075 0.075 0.850 0.870]);







% (For now, use same scale as SOM outputs) NOT
minc = min(pcadat(:));
maxc = max(pcadat(:));




  legend(fldnm1, fldnm2, 'Location', 'Best');





  dist = nonnan;
  for idx = 1:sds
    dav = mean(dist);
    sd = std(dist);
    dist = dist(dav - dsd > dist | dist > dav + dsd);
  end;
  clear dist;







  % ---
  ef.rules.hwh.desc = 'Winds bringing atypical water onto the reef (high sea temp. variability + weak air temp. var. + strong wind var.)';
  ef.rules.hwh.facts = {'seandbc_variability','airt_variability','windsd_7day'};
  ef.rules.hwh.fuzzies = {{'high','very-high','drastic-high'}, ...
                      {'drastic-low','very-low','low'}, ...
                      {'high','very-high','drastic-high'}};

  % ---
  ef.rules.hha.desc = 'Rapid cooling possibly coincident with heat flux on the reef (high sea temp. variability + high air temp. var. + any wind var.)';
  ef.rules.hha.facts = {'seandbc_variability','airt_variability','windsd_7day'};
  ef.rules.hha.fuzzies = {{'high','very-high','drastic-high'}, ...
                      {'high','very-high','drastic-high'}, ...
                      {'drastic-low','very-low','low','somewhat-low','average','somewhat-high','high','very-high','drastic-high'}};

  % ---
  ef.rules.hna.desc = 'Heat flux and possibly some air cooling on the reef (high sea temp. variability + normal air temp. var. + any wind var.)';
  ef.rules.hna.facts = {'seandbc_variability','airt_variability','windsd_7day'};
  ef.rules.hna.fuzzies = {{'high','very-high','drastic-high'}, ...
                      {'somewhat-low','average','somewhat-high','high','very-high','drastic-high'}, ...
                      {'drastic-low','very-low','low','somewhat-low','average','somewhat-high','high','very-high', 'drastic-high'}};













cutoff = ceil(hlp/2);




% Remove seasonal and longer-period cycles with a high-pass filter
[B,A] = butter((filtOrder*2), (1/(40*24)), 'high');
datHP = filtfilt(B, A, rawdat);









clear dts dat fdts;
clear covm;
clear evecs;





% Remove daily cycle from the interpolated time series
ncols = (12);
nrows = floor(length(dat) / ncols);
nelts = (nrows * ncols);
dts = reshape(dts(1:nelts)', [ncols nrows])';
dat = reshape(dat(1:nelts)', [ncols nrows])';
% Remove the semidiurnal "trend" (i.e., diurnal cycle)
dat = detrend(dat', 'linear')';







% Remove daily cycle from time series
dvec = datevec(dts);
detrend_bps = find(mod(dvec(:,4), 24) == 0);
dat = detrend(dat', 'linear')';
keyboard;
dat = detrend(dat', 'linear', detrend_bps)';








%              [.075 .075 .85 .85]);






  title('Sea temperature time series (2008)');






% Normalize by variance
ntdat = tdat / (std(tdat));
nwdat = wdat / (std(wdat));

% Remove any mean and trend from normalized sea temperature
ntdat = detrend(ntdat, 'linear');


disp('Plotting wavelet spectra');

nyears = 1;
yearhours = (nyears*365*24) - 1;
wtw = [wdts(end-yearhours:end) ; nwdat(end-yearhours:end)];
wtt = [tdts(end-yearhours:end) ; ntdat(end-yearhours:end)];







% NOW "Remove any mean and trend from normalized sea temperature"
ntdat = detrend(ntdat, 'linear');
wtw = [wdts(end-yearhours:end) ; nwdat(end-yearhours:end)];
wtt = [tdts(end-yearhours:end) ; ntdat(end-yearhours:end)];

figure;
set(gcf, 'units', 'normalized', 'outerposition', [0 0 1 1]);
xwt(wtw, wtt);
datetick('x', 'mmmyy', 'keeplimits','keepticks');
title('Wavelet cross-spectrum: wind vs. DETRENDED temperature (2008)');








figure;
xwt(wtw, wtt);
datetick('x', 'mmmyy', 'keeplimits','keepticks');
title('Wavelet cross-spectrum: wind vs. temperature (2008)');


rwtw = [wdts(end-yearhours:end) ; (1./nwdat(end-yearhours:end))];
figure;
xwt(rwtw, wtt);
datetick('x', 'mmmyy', 'keeplimits','keepticks');
title('Wavelet cross-spectrum: RECIPRICOL-wind vs. temperature (2008)');


figure;
xwt(real(nwdat(end-yearhours:end)), ntdat(end-yearhours:end));
datetick('x', 'mmmyy', 'keeplimits','keepticks');
title('Wavelet cross-spectrum: wind "U" vs. temperature (2008)');


figure;
xwt(imag(nwdat(end-yearhours:end)), ntdat(end-yearhours:end));
datetick('x', 'mmmyy', 'keeplimits','keepticks');
title('Wavelet cross-spectrum: wind "V" vs. temperature (2008)');








wtw = [wdts(end-yearhours:end) ; wdat(end-yearhours:end)];
wtt = [tdts(end-yearhours:end) ; tdat(end-yearhours:end)];
figure;
xwt(wtw, wtt);
datetick('x', 'mmmyy', 'keeplimits','keepticks');
title('Wavelet cross-spectrum: wind vs. temperature (NON-NORMALIZED, 2008)');







% clim = [ ...
%     0.030 0.036 0.058  0.106 0.169 0.246 ...
%     0.314 0.367 0.410  0.426 0.427 0.398 ...
%     0.355 0.295 0.263  0.234 0.207 0.177 ...
%     0.142 0.112 0.096  0.070 0.051 0.036 ...
%        ];





% Now remove the DAILY CYCLE:
% First, build a daily "climatology"
for ix = 1:24
  clim(ix) = 0;
  totix = 0;
  for wkix = ix:7:size(fdat,2)
    gdix = find(~isnan(fdat(:, wkix)));
    totix = totix + length(gdix);
    clim(ix) = clim(ix) + sum(fdat(gdix, wkix));
  end;
  if ( totix > 0 )
    clim(ix) = clim(ix) / totix;
  end;
end;







% Now remove the DAILY CYCLE:
% First, build a daily "climatology"
for ix = 1:24
  clim(ix) = 0;
  totix = 0;
  for wkix = ix:7:size(fdat,2)
    gdix = find(~isnan(fdat(:, wkix)));
    totix = totix + length(gdix);
    clim(ix) = clim(ix) + sum(fdat(gdix, wkix));
  end;
  if ( totix > 0 )
    clim(ix) = clim(ix) / totix;
  end;
end;

% Then subtract it from all elements
clim = repmat(clim, [7 nrows]);
fdat = fdat - clim;








% Remove the mean for each week from that week
fdat = mean(fdat, 2);



mufdat = mean(fdat, 2);
mufdat = repmat(mufdat, [1 size(fdat,2)]);
fdat = fdat - mufdat;
clear mufdat;


% Replace any NaNs in the original, with the value of the new grand mean
fdat(badix) = mean(fdat(~badix));











begwk = 900;

figure;
set(gcf, 'units', 'normalized', 'outerposition', [0 0 1 1]);
for ix = 1:nmaps
  subplot(mapdims(2), mapdims(1), ix);
  bix = (begwk*ncols) + ((ix-1)*ncols);
  plot(dts((bix+1):(bix+ncols)), dat((bix+1):(bix+ncols)));
  datetick('keepticks', 'keeplimits'); set_datetick_cursor;
%   xlim([1 ncols]);
%   ylim([cmin cmax]);
end;


fdts = reshape(dts(1:nelts), [ncols nrows])';
fdat = reshape(dat(1:nelts), [ncols nrows])';

badix = (~isfinite(fdat));

figure;
set(gcf, 'units', 'normalized', 'outerposition', [0 0 1 1]);
for ix = 1:nmaps
  subplot(mapdims(2), mapdims(1), ix);
  plot(fdts(begwk+ix,:), fdat(begwk+ix,:));
  datetick('keepticks', 'keeplimits'); set_datetick_cursor;
%   xlim([1 ncols]);
%   ylim([cmin cmax]);
end;
keyboard;







begix = find(~isnan(station.sea_t.data), 1);
begdt = station.sea_t.date(begix);
begix = begix + 18000;  % Skip ahead a long way
figure;
set(gcf, 'units', 'normalized', 'outerposition', [0 0 1 1]);
for ix = 1:nmaps
  subplot(mapdims(2), mapdims(1), ix);
  bix = begix + ((ix-1)*72);
  plot(station.sea_t.date((bix+1):(bix+72)), station.sea_t.data((bix+1):(bix+72)));
  datetick('keepticks', 'keeplimits'); set_datetick_cursor;
%   xlim([1 72]);
%   ylim([cmin cmax]);
end;



[dts, dat] = gap_expand(station.sea_t.date, station.sea_t.data);


begix = find((dts >= begdt), 1);
begix = begix + 18000;  % Skip ahead a long way
figure;
set(gcf, 'units', 'normalized', 'outerposition', [0 0 1 1]);
for ix = 1:nmaps
  subplot(mapdims(2), mapdims(1), ix);
  bix = begix + ((ix-1)*72);
  plot(dts((bix+1):(bix+72)), dat((bix+1):(bix+72)));
  datetick('keepticks', 'keeplimits'); set_datetick_cursor;
%   xlim([1 72]);
%   ylim([cmin cmax]);
end;
keyboard;









% % If we wanted to try to make this frame-independent...
% fdts = [];
% fdat = [];
% for ix = 1:72
%   dt = dts(ix:end);
%   dd = dat(ix:end);
%   nrows = floor(length(dd) / 72);
%   fdts(end+1:end+nrows, 72) = reshape(dt, [nrows 72]);
%   fdat(end+1:end+nrows, 72) = reshape(dd, [nrows 72]);
% end;
% [ign,dtsi] = sort(fdts(:,1));
% ...












  % Idiotic IMPORTDATA. Idiotic NDBC. (If a field appears in the header, but
  % only begins having non-blank data in it midway through the file, then
  % importdata will parse that extra column as a line all by itself, with all
  % other columns NaN-filled. This is meant to try to fix that, if possible.
  bogusix = find(isnan(dat(:,YYix)) | isnan(dat(:,MMidx)) | isnan(dat(:,DDix)) | isnan(dat(:,hhix)));






function stn = load_ndbc_data(stn, yr, datapath)
%function stn = load_ndbc_data(stn, yr, datapath)
%
% Last Saved Time-stamp: <Sat 2009-02-28 12:56:07 Eastern Standard Time gramer>


  if ( ~exist('datapath', 'var') || isempty(datapath) )
    datapath = './data';
  end;

  fname = sprintf('%s/%sh%d.txt', datapath, lower(stnam), yr);

  if ( ~exist(fname, 'file') )
    warning('Skipping non-existent year file "%s"...', fname);
    return;
  end;




function stn = load_ndbc_data(stn, stnam)


  yrs = 1984:2008;

  for ix = 1:length(yrs)

    yr = yrs(ix);
    fname = sprintf('data/%sh%d.txt', lower(stnam), yr);

    if ( ~exist(fname, 'file') )
      warning('Skipping non-existent year file "%s"...', fname);
      continue;
    end;

    x = importdata(fname);
    if ( ~isfield(x, 'textdata') )
      error('Unable to find header in "%s"!', fname);
    end;
    if ( ~isfield(x, 'data') )
      error('Unable to find data in "%s"!', fname);
    end;

    % NDBC headers have been remarkably consistent for >20 years! From...
    %YY MM DD hh WD   WSPD GST  WVHT  DPD   APD  MWD  BAR    ATMP  WTMP  DEWP  VIS
    %  to ...
    %#YY  MM DD hh mm WDIR WSPD GST  WVHT   DPD   APD MWD   PRES  ATMP  WTMP  DEWP  VIS  TIDE
    %#yr  mo dy hr mn degT m/s  m/s     m   sec   sec degT   hPa  degC  degC  degC   mi    ft

    % But sometimes, flaky IMPORTDATA() fails to parse header lines
    if ( size(x.textdata) == [1 1] )
      hdr = textscan(x.textdata{:}, '%s');
      hdr = hdr{:};
    elseif ( size(x.textdata,1) > 1 )
      hdr = textscan(x.textdata{1,1}, '%s');
      hdr = hdr{:};
    else
      hdr = x.textdata;
    end;

    dat = x.data;


    % Parse timestamp columns

    YYix = find((strcmpi(hdr, 'YEAR') | strcmpi(hdr, 'YR') | strncmpi(hdr, 'YY', 2)), 1);
    MMix = find((strcmp(hdr, 'MM')), 1);
    DDix = find((strcmpi(hdr, 'DD')), 1);
    hhix = find((strcmpi(hdr, 'HH')), 1);
    if ( isempty(MMix) || isempty(DDix) || isempty(hhix) )
      error('Needed timestamp header fields not found in "%s"!', fname);
    end;

    dts = datenum(yr, dat(:,MMix), dat(:,DDix), dat(:,hhix), 0, 0);


    % Parse data columns (of interest)

    dirix = find((strcmpi(hdr, 'DIR') | strncmpi(hdr, 'WD', 2)), 1);
    spdix = find(strcmpi(hdr, 'WSPD'), 1);
    gstix = find((strcmpi(hdr, 'GST') | strncmpi(hdr, 'WG', 2)), 1);
    baromix = find((strcmpi(hdr, 'PRES') | strncmpi(hdr, 'BAR', 3)), 1);
    airtix = find((strcmpi(hdr, 'ATMP') | strncmpi(hdr, 'AIR', 3)), 1);
    seatix = find((strcmpi(hdr, 'WTMP') | strncmpi(hdr, 'SEA', 3)), 1);

    stn = update_ndbc_field(stn, 'wind1_dir', dts, dat, dirix);
    stn = update_ndbc_field(stn, 'wind1_speed', dts, dat, spdix);
    stn = update_ndbc_field(stn, 'wind1_gust', dts, dat, gstix);
    stn = update_ndbc_field(stn, 'barom', dts, dat, baromix);
    stn = update_ndbc_field(stn, 'air_t', dts, dat, airtix);
    stn = update_ndbc_field(stn, 'sea_t', dts, dat, seatix);

  end;

return;


function stn = update_ndbc_field(stn, fld, dts, dat, fldix)

  if ( ~isempty(fldix) )
    npts = numel(dts);

    if ( ~isfield(stn, fld) )
      stn.(fld).date = [];
      stn.(fld).data = [];
    end;

    stn.(fld).date(end+1:end+npts) = dts;
    qc = dat(:,fldix);
    qc(qc == 99.0) = nan;
    qc(qc == 999.0) = nan;
    stn.(fld).data(end+1:end+npts) = qc;
  end

return;







stn.dts = [];
stn.wind1_dir = [];
stn.wind1_speed = [];
stn.air_t = [];
stn.barom = [];
stn.sea_t = [];




    if ( ~isempty(spdix) )
      stn.wind1_speed.date(end+1:end+npts) = dts;
      stn.wind1_speed.data(end+1:end+npts) = dat(:,spdix);
    end
    if ( ~isempty(gstix) )
      stn.wind1_gust.date(end+1:end+npts) = dts;
      stn.wind1_gust.data(end+1:end+npts) = dat(:,gstix);
    end
    if ( ~isempty(baromix) )
      stn.barom.date(end+1:end+npts) = dts;
      stn.barom.data(end+1:end+npts) = dat(:,baromix);
    end
    if ( ~isempty(airtix) )
      stn.air_t.date(end+1:end+npts) = dts;
      stn.air_t.data(end+1:end+npts) = dat(:,airtix);
    end
    if ( ~isempty(seatix) )
      stn.sea_t.date(end+1:end+npts) = dts;
      stn.sea_t.data(end+1:end+npts) = dat(:,seatix);
    end




  fid = fopen(fname, 'r');
  if ( fid < 1 )
    error('Failed to open "%s"!', fname);
  end;

  hdr = fgetl(fid);
  hix = strfind(hdr, ' hh');
  wix = strfind(hdr, ' WTMP');
  wix = wix - hix - 4;

  fmt = ['%d %d %d %d %*' num2str(wix) 'c %f %*c\n'];
  vals = textscan(fid, fmt);
  dts = datenum(vals{1},vals{2},vals{3},vals{4},0,0);
  wtmp = vals{5};
  fclose(fname);

  stn.dts(end+1:end+npts) = dts;
  stn.wtmp(end+1:end+npts) = wtmp;









f.dts = [];
f.dir = [];
f.spd = [];
f.atmp = [];
f.barom = [];
f.wtmp = [];

for ix = 1:length(yrs)

  yr = yrs(ix);
  fname = sprintf('data/fwyf1h%d.txt', yr);

  fid = fopen(fname, 'r');
  if ( fid < 1 )
    error('Failed to open "%s"!', fname);
  end;

  hdr = fgetl(fid);
  hix = strfind(hdr, ' hh');
  wix = strfind(hdr, ' WTMP');
  wix = wix - hix - 4;

  fmt = ['%d %d %d %d %*' num2str(wix) 'c %f %*c\n'];
  vals = textscan(fid, fmt);
  dts = datenum(vals{1},vals{2},vals{3},vals{4},0,0);
  wtmp = vals{5};
  fclose(fname);

  npts = numel(dts);
  f.dts(end+1:end+npts) = dts;
  f.wtmp(end+1:end+npts) = wtmp;

end;






  if ( isfield(x, 'colheaders') )
    wtmpix = find(strcmpi(x.colheaders, 'WTMP'));
    fprintf(1, '%d - %d\n', yr, wtmpix);
  end;




  str='name'; strl=length(str); nameix = find(strncmp(lower(hdrs), str, strl), 1);
  % Avoid problems with stations being "renamed" over the years...
  if ( ~isempty(nameix) )
    dat.textdata{2:end, nameix} = translate_station_name(dat.textdata{2,nameix});
  end;



    % Avoid problems with stations being "renamed" over the years...
    if ( nameidx ~= 0 )
      new_station_name = translate_station_name(vals{nameidx}{1});
      vals(nameidx) = { new_station_name };
    end;




    % Now read all data lines at one blow
    fmt = '';
    for idx = 1:length(hdr)
        if ( idx == stmpidx )
          fmt = [ fmt '%20s' ];
        elseif ( idx ~= stmpidx+1 )
          fmt = [ fmt '%10s' ];
        end;
    end





    if ( jday == 180 )
      disp(jday);
    end;



        if ( strfind(lower(hdr{idx}), 'timestamp') )
        end;



  % Create the button group.
  bgh = uibuttongroup('Visible', 'off', 'Units', 'normalized', ...
                      'Position', [0.00 0.25 0.10 0.50]);




  x = [p1(1) p1(1)+offset(1) p1(1)+offset(1) p1(1) p1(1)];
  y = [p1(2) p1(2) p1(2)+offset(2) p1(2)+offset(2) p1(2)];






  finalRect = rbbox;

  x=finalRect(1); y=finalRect(2); w=finalRect(3); h=finalRect(4);







%   disp(src);
%   disp([eventdata.EventName,'  ',... 
%         get(eventdata.OldValue,'String'),'  ', ...
%         get(eventdata.NewValue,'String')]);
%   disp(get(get(src, 'SelectedObject'), 'String'));
  disp(['New val is ' newval]);
    disp('GOOD BtnDown set');
    disp('BAD BtnDown set');






  qaflgs = stn.(flgname).data;




function stn = qa_ts_edit_graph(stn, varname)
%function stn = qa_ts_edit_graph(stn, varname)
% 
% Allow user to edit a graph by modifying 'stn.([varname '_qc_flags'])', a
% time series of 4-byte bit-vectors, QC flags for 'stn.(varname)'. Flags are:
% 
%     Bit 16 = marked as BAD by human reviewer
%     Bit 17 = marked as GOOD (despite above flags) by human reviewer
% 
% Last Saved Time-stamp: <Fri 2009-02-06 16:17:32 Eastern Standard Time gramer>
% 

  qcname = [varname '_qc'];
  flgname = [varname '_qc_flags'];

  figure; hold on;
  ax = gca;
  set(ax, 'OuterPosition', [0.03 0.00 0.94 1.00]);
  rawph = plot(stn.(varname).date, stn.(varname).data, 'r.');
  set(rawph, 'MarkerSize', 5);
  qcph = plot(stn.(qcname).date, stn.(qcname).data, 'b.');
  set(qcph, 'MarkerSize', 5);
  set(qcph, 'HitTest', 'off');
  datetick;
  hold off;

  qaflgs = stn.(flgname).data;
  h = datacursormode(gcf);
  set(h, 'UpdateFcn', @qa_ts_datacursor);

  rawuic = uicontrol('Style', 'checkbox', 'String', 'Raw', ...
                     'Value', 1, ...
                     'Position', [10 10 60 20]);
  set(rawuic, 'Callback', {@qa_ts_uicb, ax, rawph});
  qcuic = uicontrol('Style', 'checkbox', 'String', 'QC', ...
                    'Value', 1, ...
                    'Position', [10 30 60 20]);
  set(qcuic, 'Callback', {@qa_ts_uicb, ax, qcph});
  set(gcf,'Toolbar','figure');

  % Create the button group.
  % bgh = uibuttongroup('Visible', 'off', 'Position', [20 80 80 120]);
  bgh = uibuttongroup('Visible', 'off', 'Position', [0 0 0.1 0.5]);
  % Create three radio buttons in the button group.
  exO = uicontrol('Style', 'Radio', 'String', 'Explore', ...
                  'Position',[5 300 60 30], 'Parent', bgh);
  mBd = uicontrol('Style', 'Radio', 'String', 'GOOD', ...
                  'Position',[5 200 60 30], 'Parent', bgh);
  mGd = uicontrol('Style', 'Radio', 'String', 'BAD', ...
                  'Position',[5 100 60 30], 'Parent', bgh);
  % Initialize some button group properties. 
  set(bgh, 'SelectionChangeFcn', @qa_ts_goodbadcbk);
  set(bgh, 'SelectedObject', exO);
  set(bgh, 'Visible', 'on');

  uiwait;
  % ...

return;








                  'pos',[5 150 100 30],'parent',bgh,'HandleVisibility','off');






            facts = rule.facts{fi};
            % Allow for conditions that may match one of several 'similar'
            % facts, e.g., a sea temperature condition that can match fact
            % sensors 'seandbc' and/or 'sea1m'.
            if ( ~iscell(facts) )
              facts = { facts };
            end;
            foundFactory = false;
            for fii = 1:length(facts)
              fact = facts{fii};
            if ( isfield(stn.factories, fact) )
              foundFactory = true;






    rules = cellstr(sortrows(char(rulearr)', -1)');





    endcol = find(isnumeric(rawcells{hdrrow,:}) & isnan(rawcells{hdrrow,:}), 1);
    if ( isempty(endcol) )
      error('Cannot find last data column in Excel file "%s"?', fname);
    end;
    endcol = endcol - 1;





  % Do checks in order, from least- to most-stringent






    dstr = [ rawcells{begrow:endrow, datecol} ' ' ...
             num2str(rawcells{begrow:endrow, datecol}) ':0:0' ];
    dnum = datenum(dstr);










    for col = begcol:endcol
      if ( isnumeric(
      result.(hdr{col}).data = rawcells{begrow:endrow, col};
    end;










        %set(ah(1), 'XTick', []);




                if ( edays(end) ~= length(hits(fi).datenum) )
                    edays(end+1) = ;
                end;




    for vi = 1:length(vars)

        figure;

        % Plot S/RI for all events above the actual environmental data
        [ah, p1, p2] = ...
            plotyy(station.(vars{vi}).date, station.(vars{vi}).data, ...
                   uniqmos, sri);
        set(p2, 'LineStyle', 'none');
        set(p2, 'Marker', 'v');
        set(p2, 'MarkerEdgeColor', 'r');

        axes(ah(1));
        datetick;
        axes(ah(2));
        datetick;
        legend('EF S/RI');
        syncxaxes(ah);

        th = title(sprintf('%s vs. time: Station %s', vars{vi}, stanm));
        set(th, 'Interpreter', 'none');
        print('-dpng', sprintf('%s-%s-%s.png', stanm, ef.name, vars{vi}));

    end;







    % Plot all antecedent time series, annotated by monthly cumulative S/RIs
    vars = associate_variables(station, events);
    for vi = 1:length(vars)
        rng = [ min(station.(vars{vi}).data), (max(station.(vars{vi}).data)*1.1) ];
        sriadj = interp1(linspace(0, max(sri)), linspace(rng(1), rng(2)), sri);

        figure;
        hold on;
        plot(station.(vars{vi}).date, station.(vars{vi}).data);
        % Plot S/RI for all events above the actual environmental data
        plot(uniqmos, sriadj, 'rv');
        legend('Data', 'EF S/RI');
        datetick;
        ylim([rng(1) (rng(2)*1.1)]);
        th = title(sprintf('%s vs. time: Station %s', vars{vi}, stanm));
        set(th, 'Interpreter', 'none');
        print('-dpng', sprintf('%s-%s-%s.png', stanm, ef.name, vars{vi}));
    end;





function [ uniqmos, sri ] = accumulate_monthly_sri(stn, events)
%function [ uniqmos, sri ] = accumulate_monthly_sri(stn, events)
%
% Return a monthly cumulative S/RI for *all events* in the struct 'events'
%
% Last Saved Time-stamp: <Tue 2008-08-19 15:30:06 Eastern Daylight Time gramer>
%

    [Y M D h m s] = datevec(events.jday);
    uniqmos = datenum(Y, M, 1, 0, 0, 0);
    for idx = 1:length(uniqmos)
        uniqmo = uniqmos(idx);
        sri(idx) = length(find(uniqmos == uniqmo));
    end;

return;








            fprintf(1, '\t EF.%s: %d matching days\n', ...
                    upper(rules{ri}), length(events.rules.(rules{ri}).jday));

            events.jday = union(events.jday, events.rules.(rules{ri}).jday);







                % And then remove that matching fact from the agenda
                %rmidx = find(ismember(agenda.(rule.facts{fi}).datenum, ...
                %                      hits(fi).datenum));








            events.rules.(rules{ri}).jday = unique(floor([hits(:).datenum]));








            for fi = 1:length(events.rules.(rules{ri}).facts)
                dn = find(floor(hits(:).datenum) == events.rules.(rules{ri}).jday(ji));
                events.rules.(rules{ri}).sri = sum([hits(:).sri(dn)]);
            end;
            events.rules.(rules{ri}).sri = sum(hits(:).sri(hi));









                    epers = [1 (find(diff(dn) < (4/24)) + 1)];
                    perlen = 0;
                    % If no extended periods, just pick the first match
                    maxidx = edays(di);
                    for pi = 1:length(epers)-1
                        if ( perlen < length(dn(epers(pi))) )
                            perlen = length(dn(epers(pi)));
                            maxidx = edays(di) + epers(pi) - 1;
                        end;
                    end;
                    
                    epers = [ 1 (find(diff(floor(hits(fi).datenum()) == 0) + 1)];






    x = { 'vl', 'lo', 'lo', 'vl',  'vl', 'sl', 'av', 'sl' };
    rl = diff([ 0 find(~strcmp(x(1:end-1), x(2:end))) length(x) ]);
    ev = x([ find(~strcmp(x(1:end-1), x(2:end))) length(x) ]);







                    x = hits(fi).datenum(edays(di):edays(di+1)-1);
                    % Run-length encoding algorithm courtesy of "MATLAB array
                    % manipulation tips and tricks", Peter J. Acklam, Norway
                    rl = diff([ 0 find(x(1:end-1) ~= x(2:end)) length(x) ]);
                    ev = x([ find(x(1:end-1) ~= x(2:end)) length(x) ]);







                edays = [1 (find(diff(floor(hits(fi).datenum)) == 0) + 1)];
                for di = 1:length(edays)-1
                    drng = edays(di):edays(di+1)-1;
                    contig = [1 (find(diff(hits(fi).datenum(drng)) < (4/24)) + 1)];
                    hitidx = [hitidx contigidx];
                    hits(fi).datenum(drng);
                    hits(fi).fuzzies(drng);
                end;





        % Were there "events" - i.e., any days when ALL conditions matched?
        if ( ~isempty(rulehits(ri).jday) )

            for fi = 1:length(rule.facts)
                % Keep only hits from jdays when ALL conditions matched
                hitidx = find(ismember(floor(hits(fi).datenum), rulehits(ri).jday));
                hits(fi).datenum = hits(fi).datenum(hitidx);
                hits(fi).fuzzies = hits(fi).fuzzies(hitidx);

% Ignore same-fuzzy criterion for now!
%                 for zi = 1:length(rule.fuzzies)
%                     fz = [1 (find(strcmp(hits(fi).fuzzies, rule.fuzzies{zi})) + 1)];
%                     contig = [1 (find(diff(hits(fi).datenum(fz)) < (4/24)) + 1)];
%                     sameday = [1 (find(diff(floor(hits(fi).datenum(fz))) == 0) + 1)];
%                 end;
                contig = [1 (find(diff(hits(fi).datenum) < (4/24)) + 1)];
                edays = [1 (find(diff(floor(hits(fi).datenum)) == 0) + 1)];
                for di = 1:length(edays)-1
                    drng = edays(di):edays(di+1)-1;
                    hits(fi).datenum(drng);
                end;
                synerg = intersect(contig, sameday);

%                 % Just pick the FIRST matching fact for each jday
%                 hitidx = [ 1 ( find(diff(floor(hits(fi).datenum))~=0) + 1 ) ];
%                 hits(fi).datenum = hits(fi).datenum(hitidx);
%                 hits(fi).fuzzies = hits(fi).fuzzies(hitidx);

                % Store each matching fact on the event...
                events.rules.(rules{ri}).(rule.facts{fi}) = hits(fi);

                % And remove that matching fact from the agenda
                rmidx = find(ismember(agenda.(rule.facts{fi}).datenum, hits(fi).datenum));
                agenda.(rule.facts{fi}).datenum(rmidx) = [];
                agenda.(rule.facts{fi}).fuzzies(rmidx) = [];
                agenda.(rule.facts{fi}).averages(rmidx) = [];
            end;









%     mdl = ( a1 * sin((b1 * (t/365.0) * (8.0*pi)) + c1) ) ...
%         + ( a2 * sin((b2 * (t/365.0) * (8.0*pi)) + c2) );








1;
% SCRIPT par_model
%
% Parameters for composited "two sine-curve" model of daily total insolation,
% at various Caribbean and Gulf of Mexico sites [van Woesik et al., 2006]
%
% Last Saved Time-stamp: <Wed 2008-08-13 14:45:29 Eastern Daylight Time gramer>
%

mlrf1 = [ 6.443     0.1815  0.4659  0.221   1.519   1.737   ];
fgbl1 = [ 6.611     0.2084  0.2675  -0.2992 0.8668  0.2826  ];
cmrc3 = [ -2.209    0.356   2.741   7.802   0.01513 6.847   ];
srvi2 = [ -0.3057   1.372  -0.2698  -6.495  0.1299  4.0     ];
lppr1 = [ 6.456     0.1298  0.8587  0.3093  1.378	-3.453  ];
dbjm1 = [ 6.476     0.1292  0.8537  0.2422  1.386   -3.475  ];
lcci1 = [ 0.3586    1.374   -3.546  6.404   0.1371  0.8115  ];
% Belize numbers are actually from Veracruz, Mexico!
grbz1 = [ 0.3362    0.8195  1.404   6.181   0.178   0.4923  ];

% t = 1:365;
t = 1:12;

muMol_PER_W = 4.1513469579;

W_PER_kWh_PER_D = 41.666667;

muMol_PER_kWh_PER_D = (muMol_PER_W * W_PER_kWh_PER_D);

mI = insolation_model(t, mlrf1, muMol_PER_kWh_PER_D);
fI = insolation_model(t, fgbl1, muMol_PER_kWh_PER_D);
bI = insolation_model(t, grbz1, muMol_PER_kWh_PER_D);

% figure;
% plot(t, [ mI ; fI ; bI ]);
% xlim([t(1) t(end)]);
% legend('MLRF1','FGBL1','GRBZ1', 'Location', 'East');
% xlabel('Month');
% %ylabel('kW m^-^2 day^-^1');
% ylabel('Avg. \mu-mole quanta m^-^2 s^-^1');
% title('van Woesik et al. insolation model');


mFit = fit(t', mI', 'pchipinterp');
mDer1 = differentiate(mFit, t);
fFit = fit(t', fI', 'pchipinterp');
fDer1 = differentiate(fFit, t);
bFit = fit(t', bI', 'pchipinterp');
bDer1 = differentiate(bFit, t);

keyboard;
figure;
plot(t, [ mDer1 ; fDer1 ; bDer1 ]);
xlim([t(1) t(end)]);
legend('dMLRF1','dFGBL1','dGRBZ1', 'Location', 'East');
xlabel('Month');
%ylabel('kW m^-^2 day^-^1');
ylabel('\Delta \mu-mole quanta m^-^2 s^-^1 month^-^1');
title('van Woesik et al. DERIVATIVES');


return;







function mdl = insolation_model(t, site, unit_conv)
%function mdl = insolation_model(t, site, unit_conv)
%
% Implement van Woesik et al (2006) two-sine curve insolation model
%
% Last Saved Time-stamp: <Wed 2008-08-13 08:21:31 Eastern Daylight Time gramer>
%

    a1 = site(1); b1 = site(2); c1 = site(3);
    a2 = site(4); b2 = site(5); c2 = site(6);

%     mdl = ( a1 * sin((b1 * (t/365.0) * (8.0*pi)) + c1) ) ...
%         + ( a2 * sin((b2 * (t/365.0) * (8.0*pi)) + c2) );
    mdl = ( a1 * sin((b1 * t) + c1) ) ...
        + ( a2 * sin((b2 * t) + c2) );

    if ( exist('unit_conv', 'var') && ~isempty(unit_conv) && isnumeric(unit_conv) )
      mdl = mdl * unit_conv;
    end;

return;
















    disp('Loading more station data from BB files! Please wait...');
    stafiles = dir([stanm '-*.bb']);
    for fidx = 1:length(stafiles)
      disp(sprintf('Loading data file %s', stafiles(fidx).name));
      station = load_bb_data(stafiles(fidx).name, station);
    end;
    clear stafiles;






fprintf(1, 'Agenda.(%s).datenums = %d\n', rule.facts{fi}, length(agenda.(rule.facts{fi}).datenum));





  % There was a check right here to see if the variable already existed: this
  % has been REMOVED - this m-func is FAST, and you never know when things
  % may have changed somewhere down the line of dependencies for a given
  % variable: so why not build it over and over again when needed? (If the
  % variable building starts to slow things down later, a "smart check" for
  % existing data may have to come back in here somewhere...)





    % There was a check right here to see if the variable already existed: this
    % has been REMOVED - this m-func is FAST, and you never know when things
    % may have changed somewhere down the line of dependencies for a given
    % variable: so why not build it over and over again when needed? (If the
    % variable building starts to slow things down later, a "smart check" for
    % existing data may have to come back in here somewhere...)







  % Don't rebuild if we already have everything in place!
  if ( isfield(stn, varname) )
    return;
  end;







    for idx = 1:length(flds{9})
        result.([ flds{7}{idx} '_' flds{8}{idx} ]).date(idx) = ...
            datenum(flds{4}(idx), flds{2}(idx), ...
                    flds{3}(idx), flds{5}(idx), flds{6}(idx), 0);

        result.([ flds{7}{idx} '_' flds{8}{idx} ]).data(idx) = ...
            flds{9}(idx);
    end;










    % Construct result structure from loaded (station-matching) data
    for idx = 1:length(flds{9})
        result.([ flds{7}{idx} '_' flds{8}{idx} ]).date(idx) = ...
            datenum(flds{4}(idx), flds{2}(idx), ...
                    flds{3}(idx), flds{5}(idx), flds{6}(idx), 0);

        result.([ flds{7}{idx} '_' flds{8}{idx} ]).data(idx) = ...
            flds{9}(idx);
    end;










station = rmfield(station, 'factories');
flds = fieldnames(station);
for idx = 1:length(flds)
  % Try to 'translate' field into the equivalent ICON/G2 variable name
  newfld = translate_var_name(flds{idx});
  if ( ~strcmp(flds{idx}, newfld) )
    station.(newfld) = station.(flds{idx});
    station = rmfield(station, flds{idx});
  end;
end;
disp('Saving station data to MAT for future runs...');
save([stanm '.mat'], 'station');









% if ( ~isfield(station, 'factories') )
% end;







1;

% stanm = 'fwyf1';
stanm = 'mlrf1';
% stanm = 'smkf1';
% stanm = 'sanf1';
% stanm = 'dryf1';

% Make sure our data is already loaded

changed = false;
if ( ~exist('station', 'var') )

  % If possible, load from a MAT file, or...
  if ( exist([stanm '.mat'], 'file') )

    fname = [stanm '.mat'];
    disp(['Reloading station data from MAT file ' fname '...']);
    load(fname, 'station');
    changed = false;


  % If we absolutely have to, load everything from scratch!
  else

    disp('Loading station data from raw files! Please wait...');
    station = [];

    stafiles = dir([stanm '-*.txt']);
    for fidx = 1:length(stafiles)
      disp(sprintf('Loading data file %s', stafiles(fidx).name));
      station = load_10col_data(stafiles(fidx).name, station);
    end;
    clear stafiles;

    changed = true;

  end;

end;

if ( ~isfield(station, 'factories') )
  disp('Loading fact factories from CSV file...');
  station = load_factories([stanm '-factories.csv'], station);
  changed = true;
end;

if ( ~isfield(station, 'Name') || ~isfield(station.Name, 'data') || ...
     ~iscell(station.Name.data) || any(~strcmpi(strtrim(station.Name.data), stanm)) )
    error('Did we load the wrong MAT file, or is station name wrong??');
end;

if ( changed )
  disp('Saving station data to MAT for future runs...');
  save([stanm '.mat'], 'station');
end;


% Build our ecological forecast and evaluate it with this station's data


clear ef;
efnam = 'ibf';
ef = feval(['load_ef_' efnam]);


disp('Evaluating ecoforecast...');
% events = evaluate_ecoforecast(ef, station, datenum(2005,1,1), now);
events = evaluate_ecoforecast(ef, station);




% Show the caller our results onscreen


if ( isempty(events) )
    warning('No events produced!');

else

    fprintf(1, '%d event days, %d periods...\n', ...
            length(events.jday), length(events.period));

    [ uniqmos, sri ] = accumulate_monthly_sri(station, events);


    % Plot all antecedent time series, annotated by monthly cumulative S/RIs
    vars = associate_variables(station, events);
    for vi = 1:length(vars)
        rng = [ min(station.(vars{vi}).data), (max(station.(vars{vi}).data)*1.1) ];
        sriadj = interp1(linspace(0, max(sri)), linspace(rng(1), rng(2)), sri);

        figure;
        hold on;
        % Plot S/RI for all events above the actual environmental data
        plot(uniqmos, sriadj, 'rv');
        legend('S/RI');

        plot(station.(vars{vi}).date, station.(vars{vi}).data);
        datetick;
        ylim(rng);
        th = title(sprintf('%s vs. time: Station %s', vars{vi}, stanm));
        set(th, 'Interpreter', 'none');
        print('-dpng', sprintf('%s-%s-%s.png', stanm, ef.name, vars{vi}));
    end;


    [Y M D] = datevec(events.jday);
    yrs = unique(Y);
    nyrs = length(yrs);


    % Make a plot of seasonality in events per year
    clear hits;
    mos = 1:12;
    for mo = mos
        hits(mo) = (length(find(M == mo)) / nyrs);
    end;
    figure;
    plot(mos, hits);
    ylim([-1 max(hits)]);
    xlabel('Month #');
    ylabel('Events / year');
    th = title(sprintf('Event seasonality for "%s" at "%s"', ef.name, stanm));
    set(th, 'Interpreter', 'none');
    print('-dpng', sprintf('%s-%s-climatology.png', stanm, ef.name));


    % Make a plot of events by year
    clear hits;
    for iyr = 1:nyrs
        yr = yrs(iyr);
        hits(iyr) = length(find(Y == yr));
    end;
    figure;
    plot(yrs, hits);
    ylim([-1 max(hits)]);
    xlabel('Year');
    ylabel('Events');
    th = title(sprintf('Events per year for "%s" at "%s"', ef.name, stanm));
    set(th, 'Interpreter', 'none');
    print('-dpng', sprintf('%s-%s-per-year.png', stanm, ef.name));


    % Plot average length of periods, and length of longest period, by year
    clear hits maxes;
    [pY pM pD] = datevec(events.period);
    for iyr = 1:nyrs
        yr = yrs(iyr);
        findx = find(Y == yr);
        hits(iyr) = length(events.jday(findx)) / length(find(pY == yr));
        % Matlab... it's magic!
        maxes(iyr) = max(diff(find(diff(events.jday(findx)) > 2)));
    end;
    figure;
    plot(yrs, hits, yrs, maxes);
    ylim([-1 max(maxes)]);
    xlabel('Year');
    ylabel('Event avg/max duration (days)');
    th = title(sprintf('Event duration per year for "%s" at "%s"', ef.name, stanm));
    set(th, 'Interpreter', 'none');
    print('-dpng', sprintf('%s-%s-event-duration.png', stanm, ef.name));


    % Plot a histogram of period lengths over whole sample record
    figure;
    hist(diff(find(diff(events.jday) > 2)), max(maxes));
    th = title(sprintf('Event duration histogram for "%s" at "%s"', ef.name, stanm));
    set(th, 'Interpreter', 'none');
    print('-dpng', sprintf('%s-%s-duration-histo.png', stanm, ef.name));

end;
















%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
station = load_factories([stanm '-factories.csv'], station);
disp('Saving station data to MAT for future runs...');
save([stanm '.mat'], 'station');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%






    
disp('Parse header'); tic;
toc; disp('Textscan file'); tic;
toc; disp('Parse fields'); tic;
toc; disp('Build datenum'); tic;
toc; disp('Assign timestamp to all'); tic;
toc; disp('Remove invalid values'); tic;
toc;






toc; disp('Parse fields'); tic;

    % Ensure 'datenum' is always the first field in result
    result.datenum = [];

    nGoodCols = 0;
    for idx = 1:length(hdr)
        % Only process non-blank columns
        if ( length(hdr{idx}) > 0 )
            nGoodCols = nGoodCols + 1;
            % Return a vector of numbers where ever possible
            result.(hdr{idx}) = str2double(vals{idx});
            total = length( result.(hdr{idx}) );
            nans = sum( isnan(result.(hdr{idx})) );
            emps = sum( isempty(vals{idx}) );
            if ( (nans > (emps * 1.1)) && (nans > (total * 0.5)) )
                % If most of the values in the column are non-numeric but
                % non-empty, then return a cell array of strings instead.
                result.(hdr{idx}) = vals{idx};
            end;
        end;
    end;


toc; disp('Build datenum'); tic;











        hits(iyr) = length(find(diff(events.jday(findx)) <= 2) + 1) ...
            / length(find(pY == yr));






    % Plot average length of periods by year
    clear hits;
    for iyr = 1:nyrs
        yr = yrs(iyr);
        findx = find(Y == yr);
        hits(iyr) = length(find(diff(events.jday(findx)) <= 2) + 1);
    end;
    keyboard;
    figure;
    plot(yrs, hits);
    ylim([-1 max(hits)]);
    xlabel('Year');
    ylabel('Event duration (days)');
    title(sprintf('Average event duration per year for "%s" at "%s"', ef.name, stanm));
    print('-dpng', sprintf('%s-%s-event-duration.png', stanm, ef.name));








  figure;
  ah1 = gca;
  hold on;
  plot(station.sea_t_1_day_deviation_3_day_average.date, ...
       station.sea_t_1_day_deviation_3_day_average.data, 'r');
  ylabel('\DeltaT_1_d_(_3_d_)');
  % This should be a "reasonable" range for T variability!
  ylim([-0.1 2]);
  set(gca, 'YColor', 'r');


  ah2 = newaxis;
  hold on;


  plot(station.wind1_speed_3_day_average.date, ...
       station.wind1_speed_3_day_average.data, 'b');
  ylabel('|U^w^i^n^d_3_d|');
  set(gca, 'YColor', 'b');

  %% Annotate starting date of each "extended event"
  %axis(axis);
  %maxdat = max(station.wind1_speed_3_day_average.data);
  %arrow( [events.period ; repmat((maxdat*1.1), [1 length(events.period)])], ...
  %       [events.period ; repmat(maxdat, [1 length(events.period)])], ...
  %       'FaceColor', 'y' );


  % Make sure x-axes line up!
  mindt = min(min(station.sea_t_1_day_deviation_3_day_average.date), ...
              min(station.wind1_speed_3_day_average.date));
  maxdt = max(max(station.sea_t_1_day_deviation_3_day_average.date), ...
              max(station.wind1_speed_3_day_average.date));
  axes(ah1);
  xlim([mindt maxdt]);
  % Stupid matlab can't do the same thing the same way twice...
  % So just TURN OFF display of X-ticks for one of our two axes.
  set(ah1, 'XTick', []);
  set(ah1, 'XTickLabel', []);
  axes(ah2);
  xlim([mindt maxdt]);
  datetick('keeplimits', 'keepticks');
  hold off;

  title(sprintf('Events from model "%s" at site "%s"', ef.name, stanm));
  print('-dpng', sprintf('%s-%s-events.png', stanm, ef.name));




  % Make a plot of any seasonality in the dates of our events
  [Y M D h m s] = datevec(events.jday);
  clear hits;
  yrs = unique(Y);
  nyrs = length(yrs);
  mos = 1:12;
  for mo = mos
      hits(mo) = (length(find(M == mo)) / nyrs);
  end;
  figure;
  plot(mos, hits);
  ylim([-1 max(hits)]);
  xlabel('Month #');
  ylabel('Events / year');
  title(sprintf('Event seasonality for "%s" at "%s"', ef.name, stanm));
  print('-dpng', sprintf('%s-%s-climatology.png', stanm, ef.name));


  % Make a plot of events by year
  clear hits;
  for iyr = 1:nyrs
      yr = yrs(iyr);
      hits(iyr) = length(find(Y == yr));
  end;
  figure;
  plot(yrs, hits);
  ylim([-1 max(hits)]);
  xlabel('Year');
  ylabel('Events');
  title(sprintf('Events per year for "%s" at "%s"', ef.name, stanm));
  print('-dpng', sprintf('%s-%s-per-year.png', stanm, ef.name));









      for fi = 1:length(facts)
        events.rules.(rules{ri}).(facts{fi});
      end;





            events.(rules{ri}).variables. = agenda.variables.(fact) = stn.factories.(fact).variable;





    station = load_factories([stanm '-factories.csv'], station);
    disp('Saving station data to MAT for future runs...');
    save([stanm '.mat'], 'station');





AND THEN UNDONE!

In evaluate_factory.m, THIS:
    if ( ~isfield(stn, varname) )
        error('Variable name %s is not a field in "stn"!', varname);
    end;

CHANGED TO THIS:
    % Make sure all required 'derived' variables have been calculated
    stn = verify_variable(stn, varname);








    % Calculate "derived-variable" time series from raw loaded data
    disp('Calculating derived time series from raw data...');

    goodidx = find(0 < station.sea_t.data & station.sea_t.data < 38);
    [station.sea_t_1_day_deviation.date, ...
     station.sea_t_1_day_deviation.data] = ...
        window_func(station.sea_t.date(goodidx), ...
                    station.sea_t.data(goodidx), 'std', 24, 1);

    [station.sea_t_1_day_deviation_3_day_average.date, ...
     station.sea_t_1_day_deviation_3_day_average.data] = ...
        window_func(station.sea_t_1_day_deviation.date, ...
                    station.sea_t_1_day_deviation.data, 'avg', 72, 1);

    if ( strcmp(stanm, 'mlrf1') )
      [station.wind2_speed_3_day_average.date, ...
       station.wind2_speed_3_day_average.data] = ...
          window_func(station.wind2_speed.date, station.wind2_speed.data, 'avg', 72, 1);

      station.wind2_u.date = station.wind2_speed.date;
      station.wind2_u.data = station.wind2_speed.data .* (-sind(station.wind2_dir.data));
      station.wind2_v.date = station.wind2_speed.date;
      station.wind2_v.data = station.wind2_speed.data .* (-cosd(station.wind2_dir.data));

      [dts, u] = window_func(station.wind2_u.date, station.wind2_u.data, ...
                             'avg', (7*24), 24);
      [dts, v] = window_func(station.wind2_v.date, station.wind2_v.data, ...
                             'avg', (7*24), 24);
      station.wind2_u_7_day_average_dir_wind2_v.date = dts;
      station.wind2_u_7_day_average_dir_wind2_v.data = uv_to_dir(u, v);

    else
      [station.wind1_speed_3_day_average.date, ...
       station.wind1_speed_3_day_average.data] = ...
          window_func(station.wind1_speed.date, station.wind1_speed.data, 'avg', 72, 1);

      station.wind1_u.date = station.wind1_speed.date;
      station.wind1_u.data = station.wind1_speed.data .* (-sind(station.wind1_dir.data));
      station.wind1_v.date = station.wind1_speed.date;
      station.wind1_v.data = station.wind1_speed.data .* (-cosd(station.wind1_dir.data));

      [dts, u] = window_func(station.wind1_u.date, station.wind1_u.data, ...
                             'avg', (7*24), 24);
      [dts, v] = window_func(station.wind1_v.date, station.wind1_v.data, ...
                             'avg', (7*24), 24);
      station.wind1_u_7_day_average_dir_wind1_v.date = dts;
      station.wind1_u_7_day_average_dir_wind1_v.data = uv_to_dir(u, v);
    end;








      fprintf(1, 'Building %s\n', var1);
      fprintf(1, 'Building %s\n', var2);

fprintf(1, 'per %s, derop %s\n', per, derop);
return;



    fprintf(1, 'Var %s, Inst %s, Comp %s\n', ...
            varname(begidx:endidx+1), varname(begidx:endidx-1), varname(endidx+1));




%     switch( fun )
%      case {'avg', 'average'},
%      otherwise,
%     end;






    flds = fieldnames(station);
    for idx = 1:length(flds)
        % Try to 'translate' field into the equivalent ICON/G2 variable name
        newfld = translate_var_name(flds{idx});
        if ( ~strcmp(flds{idx}, newfld) )
            station.(newfld) = station.(flds{idx});
            station = rmfield(station, flds{idx});
        end;
    end;


    save([stanm '.mat'], 'station');










            fprintf(1, '\t %s.%s: %d matches\n', ...
                    rules{ri}, rule.facts{fi}, length(hits(fi).datenum));

            % If ANY combination of conditions is not met, rule FAILS
            if ( isempty(rulehits(ri).jday) )
                break;
            end;







  set(ah1, 'Color', [.5 .5 .5])




  mindat = min( min(station.sea_t_1_day_deviation_3_day_average.data), ...
                min(station.wind2_speed_3_day_average.data) );
  maxdat = max( max(station.sea_t_1_day_deviation_3_day_average.data), ...
                max(station.wind2_speed_3_day_average.data) );








    % If input was empty, just return loaded data
    if ( isempty(stn) )
        stn = result;

    % Otherwise merge new data intelligently into 'stn' using timestamps...
    % NOTE: Where no timestamp exists, WIPE OUT the old field in 'stn'!
    else

        % If we get this far - we MUST enforce uniqueness on all timestamps
        stn = unique_10col_data(stn);
        result = unique_10col_data(result);

        resflds = fieldnames(result);
        resflds = resflds(~strcmp(resflds, 'datenum'));
        for fi = 1:length(resflds)

            fld = resflds{fi};

            % If this is a new field, our job is easy
            if ( ~isfield(stn, fld) )
                stn.(fld) = result.(fld);

            % NOTE: If existing field has no timestamps, we *TRASH* it here!
            % (Also, this only trashes a field if its name happens to match a
            % loaded COLUMN HEADER. User-defined fields should be preserved.)
            elseif ( ~isfield(stn.(fld), 'date') || isempty(stn.(fld).date) || ...
                     ~isfield(stn.(fld), 'data') || isempty(stn.(fld).data) )
                warning('Field %s cannot be merged: REPLACING old value!', fld);
                stn = rmfield(stn, fld);
                stn.(fld) = result.(fld);

            % Otherwise, merge all values - overwriting old data with new
            else
                % First do the delicate merge of datenums
                newdatenum = union(stn.(fld).date, result.(fld).date);
                si = find(ismember(newdatenum, stn.(fld).date));
                ri = find(ismember(newdatenum, result.(fld).date));
                stn.(fld).date = newdatenum;
                stn.(fld).data(si) = stn.(fld).data;
                % If the input field was a cell array, keep it that way
                if ( iscell(stn.(fld).data) == iscell(result.(fld).data) )
                    stn.(fld).data(ri) = result.(fld).data;
                elseif ( iscell(stn.(fld).data) )
                    stn.(fld).data(ri) = cellstr(num2str(result.(fld).data(:)));
                else
                    stn.(fld).data(ri) = str2double(result.(fld).data);
                end;
                if ( length(stn.(fld).date) ~= length(stn.(fld).data) )
                    warning('Field %s date and data do not match!', fld);
                end;

            end;

        end;

    end;
















    disp('Loading station data from raw files! Please wait...');
    station = [];

    station = load_10col_data('mlrf1-2001-01.txt', mlrf1);
    mlrf1 = load_10col_data('mlrf1-2001-02.txt', mlrf1);
    mlrf1 = load_10col_data('mlrf1-2001-03.txt', mlrf1);

    mlrf1 = load_10col_data('mlrf1-2002-01.txt', mlrf1);
    mlrf1 = load_10col_data('mlrf1-2002-02.txt', mlrf1);

    mlrf1 = load_10col_data('mlrf1-2003-01.txt', mlrf1);
    mlrf1 = load_10col_data('mlrf1-2003-02.txt', mlrf1);

    mlrf1 = load_10col_data('mlrf1-2004.txt', mlrf1);

    mlrf1 = load_10col_data('mlrf1-2005.txt', mlrf1);

    mlrf1 = load_10col_data('mlrf1-2006.txt', mlrf1);

    mlrf1 = load_10col_data('mlrf1-2007-01.txt', mlrf1);
    mlrf1 = load_10col_data('mlrf1-2007-02.txt', mlrf1);

    mlrf1 = load_10col_data('mlrf1-2008.txt', mlrf1);








  figure;

  subplot(2, 1, 1);
  hold on;
  plot(mlrf1.mlrf1_sea_t_1_day_deviation_3_day_average.date, ...
       mlrf1.mlrf1_sea_t_1_day_deviation_3_day_average.data, 'r');
  ylabel('\DeltaT_1_d_(_3_d_)');
  % Annotate starting date of each "extended event"
  mindat = min(mlrf1.mlrf1_sea_t_1_day_deviation_3_day_average.data);
  maxdat = max(mlrf1.mlrf1_sea_t_1_day_deviation_3_day_average.data);
  axis(axis);
  arrow( [events.period ; repmat((maxdat*1.1), [1 length(events.period)])], ...
         [events.period ; repmat(maxdat, [1 length(events.period)])], ...
         'FaceColor', 'y' );
  % Make sure all x-axes line up!
  xlim([mindt maxdt]);
  datetick('keeplimits', 'keepticks');
  hold off;


  subplot(2, 1, 2);
  hold on;
  plot(mlrf1.mlrf1_wind2_speed_3_day_average.date, ...
       mlrf1.mlrf1_wind2_speed_3_day_average.data, 'b');
  ylabel('|U^w^i^n^d_3_d|');
  % Annotate starting date of each "extended event"
  mindat = min(mlrf1.mlrf1_wind2_speed_3_day_average.data);
  maxdat = max(mlrf1.mlrf1_wind2_speed_3_day_average.data);
  axis(axis);
  arrow( [events.period ; repmat((maxdat*1.1), [1 length(events.period)])], ...
         [events.period ; repmat(maxdat, [1 length(events.period)])], ...
         'FaceColor', 'y' );
  % Make sure all x-axes line up!
  xlim([mindt maxdt]);
  datetick('keeplimits', 'keepticks');
  hold off;







  mindat = min( min(mlrf1.mlrf1_sea_t_1_day_deviation_3_day_average.data), ...
                min(mlrf1.mlrf1_wind2_speed_3_day_average.data) );
  maxdat = max( max(mlrf1.mlrf1_sea_t_1_day_deviation_3_day_average.data), ...
                max(mlrf1.mlrf1_wind2_speed_3_day_average.data) );

  figure;
  ah1 = gca;
  hold on;
  plot(mlrf1.mlrf1_sea_t_1_day_deviation_3_day_average.date, ...
       mlrf1.mlrf1_sea_t_1_day_deviation_3_day_average.data, 'r');
  ylabel('\DeltaT_1_d_(_3_d_)');
  set(gca, 'YColor', 'r');

  ah2 = newaxis;
  hold on;
  plot(mlrf1.mlrf1_wind2_speed_3_day_average.date, ...
       mlrf1.mlrf1_wind2_speed_3_day_average.data, 'b');
  ylabel('|U^w^i^n^d_3_d|');
  set(gca, 'YColor', 'b');

  % Annotate starting date of each "extended event"
  axis(axis);
  arrow( [events.period ; repmat((maxdat*1.1), [1 length(events.period)])], ...
         [events.period ; repmat(maxdat, [1 length(events.period)])], ...
         'FaceColor', 'y' );

  % Make sure x-axes line up!
  axes(ah1);
  axis(axis);
  xlim([mindt maxdt]);
  datetick;
  axes(ah2);
  axis(axis);
  xlim([mindt maxdt]);
  datetick;
  hold off;


















fprintf(1, 'BEFORE hits(%d).datenum %d\n', fi, length(hits(fi).datenum));
fprintf(1, 'hitidx(%d) %d\n', fi, length(hitidx));
fprintf(1, 'AFTER  hits(%d).datenum %d\n', fi, length(hits(fi).datenum));

fprintf(1, 'events.%s.jday %d\n', rules{ri}, length(events.(rules{ri}).jday));


        fprintf(1, 'Rule %d (%s): %d hit dates, %d total hits\n', ...
                ri, rules{ri}, ...
                length(rulehits(ri).jday), length(events.(rules{ri}).jday));

        datestr(events.(rules{ri}).jday(~ismember(events.(rules{ri}).jday, rulehits(ri).jday)))













    % line(repmat(events.period, [2 1]), ...
    %      repmat([mindat; maxdat], [1 length(events.period)]));







%     x = 0:101;
%     enddidx = lcm(24,max(x));
%     d = ((enddidx-length(x)+1):enddidx)/24.0; d(end),





    if ( gcd(samplehrs, nhrs) == samplehrs )

        newdts = repmat(nan, [1 floor(length(tdts)/samplehrs)]);
        newts = repmat(nan, [1 floor(length(tdts)/samplehrs)]);

        % "Reshape" a copy of time series, once for each sample-per-period.
        % The reason we build a vector, then "reshape" it into a rectangular
        % matrix is so that the highly optimized row-wise versions of Matlab
        % funs 'mean', 'std', 'min, 'max', 'trapz', 'diff', ... can be used.
        ncopies = lcm(samplehrs, nhrs) / samplehrs;
        for idx = 1:ncopies

            % End each copy based on a "samplehrs" reading frame: so the
            % first copy will end at the true end of the time series, the
            % second will end 'samplehrs' earlier, and so forth.
            endidx = length(tsrc) - ((idx-1)*samplehrs);
            % How many windows of size 'nhrs' can *this copy* accommodate?
            nwnds = floor(endidx / nhrs);
            % Start at the beginning of the earliest possible reading frame
            begidx = endidx - (nwnds * nhrs) + 1;

            begidx = endidx - (nwnds * nhrs) + 1;
% assert(begidx == (nhrs-idx+1), 'beg');
assert(idx ~= 1 || max(tdts(begidx:nhrs:endidx)) == tdts(end), 'end==end');

            newendidx = length(newts) - (idx-1);
            newnwnds = floor(newendidx / ncopies);
            newbegidx = newendidx - (newnwnds * ncopies) + 1;

            % Assign new time-stamps based on the END of each "nhrs" period
            newdts(newbegidx:ncopies:newendidx) = tdts(begidx:nhrs:endidx);

            % Rearrange vector into nhrs-sized columns, starting from END
            newtsrc = reshape(tsrc(begidx:endidx), [nhrs nwnds]);
            % Take a row-wise mean/std/min/max/etc. - one value per column!
            newts(newbegidx:ncopies:newendidx) = feval(func, newtsrc);

        end;
















%     n = w;
%     [dts1, baseline] = run_win_test('baseline', d, x, w, n, false);
%     [dts2, baseline2] = run_win_test('baseline2', d, x, w, n, true);
%     if ( ~isequal(dts1, dts2) || ~isequal(baseline, baseline2) )
%         disp([ baseline(end), baseline2(end) ; ...
%                dts1(end), dts2(end)]);
%         warning('Old/new baselines do not match?? (len %d,%d) ==%d?', ...
%               length(baseline), length(baseline2), isequal(baseline, baseline2));
%     end;

%     n = 1;
%     [dts1, res1] = run_win_test('old form', d, x, w, n, false);
%     [dts2, res2] = run_win_test('new form', d, x, w, n, true);
%     if ( ~isequal(dts1, dts2) || ~isequal(res1, res2) )
%         disp([ res1(end), res2(end) ; ...
%                dts1(end), dts2(end)]);
%         warning('Old/new results for n=%d do not match?? (len %d,%d) ==%d?', ...
%               n, length(res1), length(res2), isequal(res1, res2));
%     end;

%     n = 2;
%     [dts1, res1] = run_win_test('old form', d, x, w, n, false);
%     [dts2, res2] = run_win_test('new form', d, x, w, n, true);
%     if ( ~isequal(dts1, dts2) || ~isequal(res1, res2) )
%         disp([ res1(end), res2(end) ; ...
%                dts1(end), dts2(end)]);
%         warning('Old/new results for n=%d do not match?? (len %d,%d) ==%d?', ...
%               n, length(res1), length(res2), isequal(res1, res2));
%     end;

%     n = 12;
%     [dts1, res1] = run_win_test('old form', d, x, w, n, false);
%     [dts2, res2] = run_win_test('new form', d, x, w, n, true);
%     if ( ~isequal(dts1, dts2) || ~isequal(res1, res2) )
%         disp([ res1(end), res2(end) ; ...
%                dts1(end), dts2(end)]);
%         warning('Old/new results for n=%d do not match?? (len %d,%d) ==%d?', ...
%               n, length(res1), length(res2), isequal(res1, res2));
%     end;


function [dts, result] = run_win_test(testid, d, x, w, n, doNew)

    func = 'avg';

    disp(sprintf('%s: %d by %d (New? %d) time:', testid, w, n, doNew));
    if ( doNew )
        tic;
          [dts, result] = window_func(d, x, func, w, n);
        toc;
    else
        tic;
          [dts, result] = window_func_v0(d, x, func, w, n);
        toc;
    end;

return;


































            % Rearrange vector into nhrs-sized columns, starting from END
            newtsrc = tsrc( ((end-startidx+1) - (nwnds*nhrs) + nhrs):(end-startidx+1) );
            newtsrc = reshape(newtsrc, [nhrs nwnds]);

            newdts(idx:ncopies:end) = ...
                tdts(((end-idx+1) - (nwnds*nhrs) + nhrs):nhrs:(end-idx+1));
            % Take a row-wise mean/std/min/max/etc. - one value per column!
            newts(idx:ncopies:end) = feval(func, newtsrc);










        % How many windows of size 'nhrs' can our time series accommodate?
        nwnds = floor(length(tdts) / nhrs);






                tdts(((end-idx+1) - (nwnds*nhrs)):nhrs:(end-idx+1));
%                tdts(((end-idx+1) - (nwnds*nhrs) + 1):nhrs:(end-idx+1));





    % SPECIAL CASE: Much faster, if we just want one sample per window
    if ( nhrs == samplehrs )
        % How many windows of size 'nhrs' can our time series accommodate?
        nwnds = floor(length(tdts) / nhrs);

        % Truncate time series based on even multiples of window-size
        tsrc = tsrc(1:(nwnds * nhrs));

        % Rearrange our vector into nhrs-sized columns
        tsrc = reshape(tsrc, [nhrs nwnds]);

        newdts = tdts((nhrs-1):nhrs:end);
        newts = feval(func, tsrc);










fwyf1.factories.spawning_seatemp.variable = 'WaterT_1d_avg';
fwyf1.factories.photo_accum.variable = 'PAR0_30d_avg';

if ( isfield(fwyf1, 'WaterT_1d_avg') )
    fwyf1 = rmfield(fwyf1, {'WaterT_1d_avg_dts', 'WaterT_1d_avg'});
end;
if ( isfield(fwyf1, 'PAR0_30d_avg') )
    fwyf1 = rmfield(fwyf1, {'PAR0_30d_avg_dts','PAR0_30d_avg'});
end;

[fwyf1.WaterT_1d_avg_dts, fwyf1.WaterT_1d_avg] = ...
    window_func(fwyf1.datenum, fwyf1.WaterT, 'avg', 24, 1);
[fwyf1.PAR0_30d_avg_dts, fwyf1.PAR0_30d_avg] = ...
    window_func(fwyf1.datenum, fwyf1.PAR0, 'avg', (24*30), 24);




clear pas;
pas.tl.facts = {'spawning_seatemp', 'photo_accum'};
pas.tl.fuzzies = {{'conducive'}, {'high','very-high','drastic-high'}};








        % Then merge every other field in loaded 'result' into 'stn'
        resflds = fieldnames(result);
        resflds = resflds(~strcmp(resflds, 'datenum'));
        for fi = 1:length(resflds)
            fld = resflds{fi};
            if ( isfield(stn, fld) )
                % Note this should work whether field is numeric or cellstr
                stn.(fld)(si) = stn.(fld);
            else
                if ( isnumeric(result.(fld)) )
                    stn.fld = repmat(nan, size(newdatenum));
                else
                    stn.fld = repmat({[]}, size(newdatenum));
                end;
            end;

            % If the input field was a cell array, keep it that way
            if ( iscell(stn.(fld)) == iscell(result.(fld)) )
                stn.(fld)(ri) = result.(fld);
            elseif ( iscell(stn.(fld)) )
                stn.(fld)(ri) = cellstr(num2str(result.(fld)(:)));
            else
                stn.(fld)(ri) = str2double(result.(fld));
            end;
        end;

    else

        error('File "%s" lacked valid date-times: Struct "stn" unchanged!', ...
              fname);








    % CAVEAT CALLER...
    % NOTE: If existing struct 'stn' had no datenums, WE *TRASH* 'stn' HERE!
    elseif ( ~isfield(stn, 'datenum') || isempty(stn.datenum) )
        warning('All existing values of second argument "stn" have been lost!');
        clear stn;
        stn = [];






fwyf1.factories.spawning_seatemp.variable = 'WaterT_1d_avg';
fwyf1.factories.photo_accum.variable = 'PAR0_30d_avg';

rmfield(fwyf1, {'WaterT_1d_avg_dts', 'WaterT_1d_avg', ...
                'PAR0_30d_avg_dts','PAR0_30d_avg'});

[fwyf1.WaterT_1d_avg_dts, fwyf1.WaterT_1d_avg] = ...
    window_func(fwyf1.datenum, fwyf1.WaterT, 'avg', 24, 1);
[fwyf1.PAR0_30d_avg_dts, fwyf1.PAR0_30d_avg] = ...
    window_func(fwyf1.datenum, fwyf1.PAR0, 'avg', (24*30), 24);



pas.t.facts = {'spawning_seatemp', 'photo_accum'};
pas.t.fuzzies = {{'conducive'}, {'high','very-high','drastic-high'}};


events = evaluate_ecoforecast(pas, fwyf1);

if ( isempty(events) )
  error('No events produced!');

else

  figure;
  hold on;
  plot(fwyf1.WaterT_1d_avg_dts, fwyf1.WaterT_1d_avg);
  newaxis;
  plot(fwyf1.PAR0_30d_avg_dts, fwyf1.PAR0_30d_avg);
  line(repmat(events.datenum, 
  hold off;

end;








%         % Parse 'bound' string for each fuzzy value-range
%         fuzzies = strfind(bound, '''');
%         for fuzzidx = fuzzies
%             strtok
%         end;
    end;
    rest = bounds;
    while ( ~isempty([rest{:}]) )
      [bd, rest] = strtok(rest, ',');
      stn.(sensors) = 
    end;





        for vi = 1:length(stn.(ef.(rules{ri}).var))
            var = stn.(ef.(rules{ri}).variable){vi};
            lbd = ef.(rules{ri}).lower{vi};
            ubd = ef.(rules{ri}).upper{vi};
            if ( isempty(hits(ri)) )
              hits(ri) = find(lbd <= var & var <= ubd));
            else
              hits(ri) = intersect(hits(ri), ...
                                   find(lbd <= var & var <= ubd));
            end;
        end;







  msg = ferror(fid);
  if ( ~isempty(msg) )
    fprintf(1, 'ERROR in textscan: "%s"\n', msg);
  end


  % flds = fscanf(fid, '%*s,%*s,%d/%d/%d %d:%d,%f\n');
  while 1
    ln = fgetl(fid);
    if ( ~ischar(ln) ); break; end;
    fprintf(1, '%s\n', ln);
    flds = sscanf(ln, '%*s,%*s,%d/%d/%d %d:%d,%f');
    keyboard;
    break;
  end;
