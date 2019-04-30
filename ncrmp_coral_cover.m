1;
%% SCRIPT ncrmp_coral_cover.m - Extract coral cover estimates from NOAA NCRMP datasets
%
% Expects: CSV files converted from original RDA (e.g., via R package 'rio')
%
% Calls: NANMAX, NANMEAN, NANMIN, NANSUM, HIST, TITLE
%
% See Also:
%   https://github.com/shgroves/NCRMP_benthics
%     Especially ./ncrmp.benthics.analysis/R/NCRMP_calculate_cover.R
%   https://cran.r-project.org/web/packages/rio/index.html
%
% Author: Lew.Gramer@noaa.gov, Lew.Gramer@gmail.com, 2018-12-20
% 
% Last Saved Time-stamp: <Mon 2019-02-11 17:01:36 Eastern Standard Time gramer>


ncrmppath = get_ecoforecasts_path(fullfile('data','NCRMP_benthics-master','ncrmp.benthics.analysis','data'));

regions_years = {...
    'FGBNMS',		2013;...
    'FGBNMS',		2015;...
    'FLK',		2014;...
    'FLK',		2016;...
    'PRICO',		2014;...
    'PRICO',		2016;...
    'SEFCRI',		2014;...
    'SEFCRI',		2016;...
    'TortugasMarq',	2014;...
    'TortugasMarq',	2016;...
    'USVI',		2013;...
    'USVI',		2015;...
    'USVI',		2017;...
};


for ryix=1:size(regions_years,1)
  region = regions_years{ryix,1};
  year = regions_years{ryix,2};
  
  rgnyr = sprintf('%s_%04d',region,year);
  %DEBUG:
  disp(rgnyr);
  
  ncrmp.(rgnyr).basename = rgnyr;
  ncrmp.(rgnyr).two_stage = false;
  if ( year == 2014 && (region == string('SEFCRI') || region == string('FLK')) )
    ncrmp.(rgnyr).basename = [ncrmp.(rgnyr).basename,'_2stage'];
    ncrmp.(rgnyr).two_stage = true;
  end;
  
  ncrmp.(rgnyr).csvfname = fullfile(ncrmppath,[ncrmp.(rgnyr).basename,'_benthic_cover.csv']);
  
  
  x = importdata(ncrmp.(rgnyr).csvfname);
  
  dat = string(x.textdata(2:end,:));
  % 1: REGION
  % 2: PRIMARY_SAMPLE_UNIT
  % 3: STATION_NR
  % 4: YEAR
  % 5: MONTH
  % 6: DAY
  % 7: HABITAT_CD
  % 8: STRAT
  % 9: RUGOSITY_CD
  %10: WTD_RUG
  %11: LAT_DEGREES
  %12: LON_DEGREES
  %13: MAPGRID_NR
  %14: SUB_REGION_NAME
  %15: SUB_REGION_NR
  %16: ZONE_NAME
  %17: ZONE_NR
  %18: MPA_NAME
  %19: MPA_NR
  %20: ADMIN
  %21: PROT
  %22: DEPTH_STRAT
  %23: MIN_DEPTH
  %24: MAX_DEPTH
  %25: METERS_COMPLETED
  %26: COVER_CAT_CD
  %27: COVER_CAT_NAME
  %28: HARDBOTTOM_P
  %29: SOFTBOTTOM_P
  %30: RUBBLE_P
  
  % X.data: 1: HARDBOTTOM_P, 2: SOFTBOTTOM_P, 3: RUBBLE_P
  
  ncrmp.(rgnyr).stn = double(dat(:,2));			% 2: PRIMARY_SAMPLE_UNIT
  ncrmp.(rgnyr).sno = double(dat(:,3));			% 3: STATION_NR [*proxy for Replicate Number*]
  ncrmp.(rgnyr).dts = datenum(double(dat(:,[4,5,6])));	% 4: YEAR, 5: MONTH, 6: DAY
  ncrmp.(rgnyr).lat = double(dat(:,11));		%11: LAT_DEGREES
  ncrmp.(rgnyr).lon = double(dat(:,12));		%12: LON_DEGREES
  ncrmp.(rgnyr).rug = double(dat(:,10));		%10: WTD_RUG [weighted rugosity]
  ncrmp.(rgnyr).mnh = double(dat(:,23));		%23: MIN_DEPTH
  ncrmp.(rgnyr).mxh = double(dat(:,24));		%24: MAX_DEPTH
  ncrmp.(rgnyr).cod = dat(:,26);			%26: COVER_CAT_CD
  ncrmp.(rgnyr).nam = dat(:,27);			%27: COVER_CAT_NAME
  ncrmp.(rgnyr).hdp = x.data(:,1);			%28: HARDBOTTOM_P
  ncrmp.(rgnyr).sfp = x.data(:,2);			%29: SOFTBOTTOM_P
  ncrmp.(rgnyr).rup = x.data(:,3);			%30: RUBBLE_P
  ncrmp.(rgnyr).cov = nansum(x.data,2);			%dplyr::mutate(Percent_Cvr = rowSums(.[28:30])) %>%
  
  clear x
  
  ncrmp.(rgnyr).meh = nanmean([ncrmp.(rgnyr).mnh,ncrmp.(rgnyr).mxh],2);
  
  
  %% ALL VALUES (NAM and COD)
  %     "Acropora cervicornis"                          "ACR CERV"
  %     "Acropora palmata"                              "ACR PALM"
  %     "Agaricia agaricites"                           "AGA AGAR"
  %     "Bare Substrate"                                "BAR SUB."
  %     "Cliona spp"                                    "CLI SPE."
  %     "Colpophyllia natans"                           "COL NATA"
  %     "Cyanophyta spp"                                "CYA SPE."
  %     "Dichocoenia stokesii"                          "DIC SPE."
  %     "Dictyota spp"                                  "DIC STOK"
  %     "Diploria labyrinthiformis"                     "DIP LABY"
  %     "Encrusting gorgonian"                          "ENCR GORG
  %     "Erythropodium caribaeorum"                     "ERY CARY"
  %     "Eusmilia fastigiata"                           "EUS FAST"
  %     "Gorgonians"                                    "GOR GORG"
  %     "Halimeda spp"                                  "HAL SPE."
  %     "Helioceris cucullata"                          "HEL CUCU"
  %     "Lobophora spp"                                 "LOB SPE."
  %     "MacroOtherCalcareous"                          "MAC CALC"
  %     "MacroOtherFleshy"                              "MAC FLES"
  %     "Madracis decactis"                             "MAD DECA"
  %     ""                                              "MAD MIRA"
  %     "Magnoliophyta spp"                             "MAG SPE."
  %     "Meandrina meandrites"                          "MEA MEAN"
  %     "Millepora spp"                                 "MIL SPE."
  %     "Montastraea cavernosa"                         "MON CAVE"
  %     "Mussa angulosa"                                "MUS ANGU"
  %     "Mycetophyllia danaana"                         "MYC DANA"
  %     "Mycetophyllia spp"                             "MYC SPE."
  %     "Oculina diffusa"                               "OCU DIFF"
  %     "Orbicella annularis"                           "ORB ANNU"
  %     "Orbicella faveolata"                           "ORB FAVE"
  %     "Orbicella franksi"                             "ORB FRAN"
  %     "Other species"                                 "OTH SPE."
  %     "Palythoa spp"                                  "PAL SPE."
  %     "Peysonnellia"                                  "PEY SPE."
  %     "Porifera spp"                                  "POF SPE."
  %     "Porites astreoides"                            "POR ASTR"
  %     "Porites divaricata"                            "POR DIVA"
  %     "Porites furcata"                               "POR FURC"
  %     "Porites porites"                               "POR PORI"
  %     "Pseudodiploria strigosa"                       "PSE STRI"
  %     "Rhodophyta cru. spp"                           "RHO CRUS"
  %     "Siderastrea radians"                           "SID RADI"
  %     "Siderastrea siderea"                           "SID SIDE"
  %     "Solenastrea bournoni"                          "SOL BOUR"
  %     "Solenastrea hyades"                            "SOL HYAD"
  %     "Stephanocoenia intersepta"                     "STE INTE"
  %     "Turf Algae Free of Sediment"                   "TUR FREE"
  
  %% Scleractinia
  %     "Acropora cervicornis"                          "ACR CERV"
  %     "Acropora palmata"                              "ACR PALM"
  %     "Agaricia agaricites"                           "AGA AGAR"
  %     "Colpophyllia natans"                           "COL NATA"
  %     "Dichocoenia stokesii"                          "DIC SPE."
  %     "Diploria labyrinthiformis"                     "DIP LABY"
  %     "Eusmilia fastigiata"                           "EUS FAST"
  %     "Helioceris cucullata"                          "HEL CUCU"
  %     "Madracis decactis"                             "MAD DECA"
  %     ""                                              "MAD MIRA"
  %     "Meandrina meandrites"                          "MEA MEAN"
  %     "Montastraea cavernosa"                         "MON CAVE"
  %     "Mussa angulosa"                                "MUS ANGU"
  %     "Mycetophyllia danaana"                         "MYC DANA"
  %     "Mycetophyllia spp"                             "MYC SPE."
  %     "Orbicella annularis"                           "ORB ANNU"
  %     "Orbicella faveolata"                           "ORB FAVE"
  %     "Orbicella franksi"                             "ORB FRAN"
  %     "Porites astreoides"                            "POR ASTR"
  %     "Porites divaricata"                            "POR DIVA"
  %     "Porites furcata"                               "POR FURC"
  %     "Porites porites"                               "POR PORI"
  %     "Pseudodiploria strigosa"                       "PSE STRI"
  %     "Siderastrea radians"                           "SID RADI"
  %     "Siderastrea siderea"                           "SID SIDE"
  %     "Solenastrea bournoni"                          "SOL BOUR"
  %     "Solenastrea hyades"                            "SOL HYAD"
  %     "Stephanocoenia intersepta"                     "STE INTE"
  
  
  %% "Massive" (vs. branching) Scleractinia
  % Per https://reefdivers.io/quick-guide-mound-boulder-caribbean-corals/5277
  %     "Colpophyllia natans"                           "COL NATA"
  %     "Dichocoenia stokesii"                          "DIC SPE."
  %     "Diploria labyrinthiformis"                     "DIP LABY"
  %     "Montastraea cavernosa"                         "MON CAVE"
  %     "Orbicella annularis"                           "ORB ANNU"
  %     "Orbicella faveolata"                           "ORB FAVE"
  %     "Orbicella franksi"                             "ORB FRAN"
  %     "Porites astreoides"                            "POR ASTR"
  %     "Siderastrea radians"                           "SID RADI"
  %     "Siderastrea siderea"                           "SID SIDE"
  %     "Solenastrea bournoni"                          "SOL BOUR"
  %     "Stephanocoenia intersepta"                     "STE INTE"
  
  % Unique stations
  
  [ucor,uix] = unique([ncrmp.(rgnyr).lon,ncrmp.(rgnyr).lat],'rows');
  ncrmp.(rgnyr).ustn = ncrmp.(rgnyr).stn(uix);
  ncrmp.(rgnyr).ulon = ucor(:,1);
  ncrmp.(rgnyr).ulat = ucor(:,2);
  clear ucor uix
  
  %[ULON,ULAT] = meshgrid(ulon,ulat);
  %[d,a]=distance_wgs84(ULAT,ULON,ULAT',ULON');
  
  %% Scleractinian COVER_CAT_CD list
  % Good codes, not sacred cod
  goodcods = string({'ACR CERV','ACR PALM','AGA AGAR','COL NATA','DIC SPE.','DIP LABY','EUS FAST','HEL CUCU','MAD DECA','MAD MIRA','MEA MEAN','MON CAVE','MUS ANGU','MYC DANA','MYC SPE.','ORB ANNU','ORB FAVE','ORB FRAN','POR ASTR','POR DIVA','POR FURC','POR PORI','PSE STRI','SID RADI','SID SIDE','SOL BOUR','SOL HYAD','STE INTE',});
  goodix = find(ismember(ncrmp.(rgnyr).cod,goodcods));
  
  % Massive coral spp.
  masscods = string({'COL NATA','DIC SPE.','DIP LABY','MON CAVE','ORB ANNU','ORB FAVE','ORB FRAN','POR ASTR','SID RADI','SID SIDE','SOL BOUR','STE INTE',});
  massix = find(ismember(ncrmp.(rgnyr).cod,masscods));
  
  for ix=1:numel(ncrmp.(rgnyr).ulon)
    ncrmp.(rgnyr).uix{ix} = find(ncrmp.(rgnyr).lon==ncrmp.(rgnyr).ulon(ix) ...
                                 & ncrmp.(rgnyr).lat==ncrmp.(rgnyr).ulat(ix));
    
    ncrmp.(rgnyr).gix{ix} = intersect(ncrmp.(rgnyr).uix{ix},goodix);
    ncrmp.(rgnyr).mix{ix} = intersect(ncrmp.(rgnyr).uix{ix},massix);
    
    if ( isempty(ncrmp.(rgnyr).gix{ix}) )
      ncrmp.(rgnyr).unh(ix) = nan;
      ncrmp.(rgnyr).ueh(ix) = nan;
      ncrmp.(rgnyr).uxh(ix) = nan;
      
      ncrmp.(rgnyr).uhd(ix) = nan;
      ncrmp.(rgnyr).usf(ix) = nan;
      ncrmp.(rgnyr).uru{ix} = nan;
      ncrmp.(rgnyr).utc(ix) = nan;
      
      ncrmp.(rgnyr).mhd(ix) = nan;
      ncrmp.(rgnyr).msf(ix) = nan;
      ncrmp.(rgnyr).mru{ix} = nan;
      ncrmp.(rgnyr).mtc(ix) = nan;
      
    
    elseif ( ncrmp.(rgnyr).two_stage )
      ncrmp.(rgnyr).unh(ix) = nanmin(ncrmp.(rgnyr).mnh(ix,:));
      ncrmp.(rgnyr).ueh(ix) = nanmean(ncrmp.(rgnyr).meh(ix,:));
      ncrmp.(rgnyr).uxh(ix) = nanmax(ncrmp.(rgnyr).mxh(ix,:));
      
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
      usno = unique(ncrmp.(rgnyr).sno(ncrmp.(rgnyr).gix{ix}));
      for snoix=1:numel(usno)
        snoixen = find(ncrmp.(rgnyr).sno(ncrmp.(rgnyr).gix{ix})==usno(snoix));
        ncrmp.(rgnyr).shd(ix,snoix) = nansum(ncrmp.(rgnyr).hdp(ncrmp.(rgnyr).gix{ix}(snoixen)));
        ncrmp.(rgnyr).ssf(ix,snoix) = nansum(ncrmp.(rgnyr).sfp(ncrmp.(rgnyr).gix{ix}(snoixen)));
        ncrmp.(rgnyr).sru(ix,snoix) = nansum(ncrmp.(rgnyr).rup(ncrmp.(rgnyr).gix{ix}(snoixen)));
        ncrmp.(rgnyr).stc(ix,snoix) = nansum(ncrmp.(rgnyr).cov(ncrmp.(rgnyr).gix{ix}(snoixen)));
        
        snomixen = find(ncrmp.(rgnyr).sno(ncrmp.(rgnyr).mix{ix})==usno(snoix));
        ncrmp.(rgnyr).mshd(ix,snoix) = nansum(ncrmp.(rgnyr).hdp(ncrmp.(rgnyr).mix{ix}(snomixen)));
        ncrmp.(rgnyr).mssf(ix,snoix) = nansum(ncrmp.(rgnyr).sfp(ncrmp.(rgnyr).mix{ix}(snomixen)));
        ncrmp.(rgnyr).msru(ix,snoix) = nansum(ncrmp.(rgnyr).rup(ncrmp.(rgnyr).mix{ix}(snomixen)));
        ncrmp.(rgnyr).mstc(ix,snoix) = nansum(ncrmp.(rgnyr).cov(ncrmp.(rgnyr).mix{ix}(snomixen)));
      end;
      clear usno snoix snoixen
      
      ncrmp.(rgnyr).uhd(ix) = nanmean(ncrmp.(rgnyr).shd(ix,:));
      ncrmp.(rgnyr).usf(ix) = nanmean(ncrmp.(rgnyr).ssf(ix,:));
      ncrmp.(rgnyr).uru{ix} = nanmean(ncrmp.(rgnyr).sru(ix,:));
      ncrmp.(rgnyr).utc(ix) = nanmean(ncrmp.(rgnyr).stc(ix,:));
      
      ncrmp.(rgnyr).mhd(ix) = nanmean(ncrmp.(rgnyr).mshd(ix,:));
      ncrmp.(rgnyr).msf(ix) = nanmean(ncrmp.(rgnyr).mssf(ix,:));
      ncrmp.(rgnyr).mru{ix} = nanmean(ncrmp.(rgnyr).msru(ix,:));
      ncrmp.(rgnyr).mtc(ix) = nanmean(ncrmp.(rgnyr).mstc(ix,:));
    
    else
      ncrmp.(rgnyr).unh(ix) = nanmin(ncrmp.(rgnyr).mnh(ix,:));
      ncrmp.(rgnyr).ueh(ix) = nanmean(ncrmp.(rgnyr).meh(ix,:));
      ncrmp.(rgnyr).uxh(ix) = nanmax(ncrmp.(rgnyr).mxh(ix,:));
      
      % dplyr::group_by(YEAR, REGION, SUB_REGION_NAME, PRIMARY_SAMPLE_UNIT, STATION_NR, LAT_DEGREES, LON_DEGREES,
      %                 MIN_DEPTH, MAX_DEPTH, ANALYSIS_STRATUM, STRAT, HABITAT_CD, PROT, cover_group) %>%
      % dplyr::summarise(Percent_Cvr = sum(Percent_Cvr)) %>%
      ncrmp.(rgnyr).uhd(ix) = nansum(ncrmp.(rgnyr).hdp(ncrmp.(rgnyr).gix{ix}));
      ncrmp.(rgnyr).usf(ix) = nansum(ncrmp.(rgnyr).sfp(ncrmp.(rgnyr).gix{ix}));
      ncrmp.(rgnyr).uru{ix} = nansum(ncrmp.(rgnyr).rup(ncrmp.(rgnyr).gix{ix}));
      ncrmp.(rgnyr).utc(ix) = nansum(ncrmp.(rgnyr).cov(ncrmp.(rgnyr).gix{ix}));
      
      ncrmp.(rgnyr).mhd(ix) = nansum(ncrmp.(rgnyr).hdp(ncrmp.(rgnyr).mix{ix}));
      ncrmp.(rgnyr).msf(ix) = nansum(ncrmp.(rgnyr).sfp(ncrmp.(rgnyr).mix{ix}));
      ncrmp.(rgnyr).mru{ix} = nansum(ncrmp.(rgnyr).rup(ncrmp.(rgnyr).mix{ix}));
      ncrmp.(rgnyr).mtc(ix) = nansum(ncrmp.(rgnyr).cov(ncrmp.(rgnyr).mix{ix}));
    end; %if ( isempty(ncrmp.(rgnyr).gix{ix}) ) elseif ( ncrmp.(rgnyr).two_stage ) else    
    
  end; %for ix=1:numel(ncrmp.(rgnyr).ulon)
  
  
  dat=[]; clear ans dat ix region year
  
  %{
  fmg; hist(ncrmp.(rgnyr).utc,100); title(textize(ncrmp.(rgnyr).basename)); xlim([0,100]);
  %}

end; %for ryix=1:size(regions_years,1)

clear ans rgnyr ryix
