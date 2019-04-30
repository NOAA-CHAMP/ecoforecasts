1;
%% SCRIPT COMPARISON_PLOT_SPATIAL_DT_HC.m:
%
% Adapted from COMPARE_SPATIAL_DT_HC.m, lew.gramer@gmail.com, 2017-Jul-22
%
% Plot comparisons of predicted sea temperature extremes (CALC_SPATIAL_DT_HC)
% to observed seafloor temperatures from FKNMS, SEAKEYS, Lirman data sets,
% based on each SCENARIO: 1=2010 Cold Snap, 2=Normal Winter, 3=Normal Summer,
% 4=1998 Bleaching), using Horizontal Convection SCALING (SS - DEFAULT, US,
% SU, UU, or for the extreme cases, SS_SU, or SS_UU). If DOPRINT (DEFAULT:
% true), print figures in FIGSPATH (DEFAULT: GET_ECOFORECASTS_PATH('figs')).
%
% Last Saved Time-stamp: <Wed 2017-08-16 11:34:43 Eastern Daylight Time gramer>

if ( ~exist('scaling','var') || isempty(scaling) )
  scaling = 'SS';
end;
if ( ~exist('doFigs','var') || isempty(doFigs) )
  doFigs = false;
end;
if ( ~exist('doPrint','var') || isempty(doPrint) )
  doPrint = false;
end;
if ( ~exist('figspath','var') )
  figspath = get_ecoforecasts_path('figs');
end;

figsbasenm = mfilename;

if ( ~exist('lon','var') && ~exist('lat','var') && ~exist('h','var') )
  if ~exist('subrgn','var'); subrgn='FRT'; end;
  if ~exist('use_habitat_map','var'); use_habitat_map=false; end;
  if ~exist('allow_hard_bottom','var'); allow_hard_bottom=true; end;
  calc_spatial_dt_hc;
end;

dlat = min(diff(unique(lat)));
[LON,LAT] = meshgrid(lon,lat);

if ( ~exist('stns','var') )
  anfknms;
end;
if ( ~isfield(stns,'BNPIN') )
  lstns = get_lirman_thermistors;
  stnms = fieldnames(lstns);
  for stix=1:numel(stnms);
    stnm = stnms{stix};
    if isfield(stns,stnm); error('Lirman already in STNS??'); end;
    STNM = upper(stnm);
    stns.(STNM) = lstns.(stnm);
    stns.(STNM).fknms_seatemp = lstns.(stnm).seatemp;

    sites.lons(end+1) = lstns.(stnm).lon;
    sites.lats(end+1) = lstns.(stnm).lat;
    sites.stnms{end+1} = upper(stnm);
    sites.depths(end+1) = lstns.(stnm).depth;

    [err,ix] = min(abs(LON(:)-stns.(STNM).lon) + abs(LAT(:)-stns.(STNM).lat));
    [jix,iix] = ind2sub(size(LON),ix);
    if ( 1<=jix && jix<=size(bet,2) && 1<=iix && iix<size(bet,1) )
      stns.(STNM).ngdc_depth = h(iix,jix);
      stns.(STNM).ngdc_beta = bet(iix,jix);
      stns.(STNM).ngdc_beta_ang = ang(iix,jix);

      sites.betas(end+1) = bet(iix,jix);
      sites.beta_angs(end+1) = ang(iix,jix);
    else
      stns.(STNM).ngdc_depth = nan;
      stns.(STNM).ngdc_beta = nan;
      stns.(STNM).ngdc_beta_ang = nan;

      sites.betas(end+1) = nan;
      sites.beta_angs(end+1) = nan;
    end;
    sites.iso_angs(end+1) = nan;
  end;
  lstns=[]; clear lstns stnm STNM stnms stix 
end;


% Specify warming rate based on HC scaling
switch ( scaling ),
 case 'SS',
  dTC = squeeze(dTdt_SS(1,:,:));
  dTW = squeeze(dTdt_SS(2,:,:));
  dTS = squeeze(dTdt_SS(3,:,:));
  dTB = squeeze(dTdt_SS(4,:,:));
  dTdthc = dTdthc_SS;
  dTdt = dTdt_SS;

 case 'US',
  dTC = squeeze(dTdt_US(1,:,:));
  dTW = squeeze(dTdt_US(2,:,:));
  dTS = squeeze(dTdt_US(3,:,:));
  dTB = squeeze(dTdt_US(4,:,:));
  dTdthc = dTdthc_US;
  dTdt = dTdt_US;

 case 'SU',
  dTC = squeeze(dTdt_SU(1,:,:));
  dTW = squeeze(dTdt_SU(2,:,:));
  dTS = squeeze(dTdt_SU(3,:,:));
  dTB = squeeze(dTdt_SU(4,:,:));
  dTdthc = dTdthc_SU;
  dTdt = dTdt_SU;

 case 'UU',
  dTC = squeeze(dTdt_UU(1,:,:));
  dTW = squeeze(dTdt_UU(2,:,:));
  dTS = squeeze(dTdt_UU(3,:,:));
  dTB = squeeze(dTdt_UU(4,:,:));
  dTdthc = dTdthc_UU;
  dTdt = dTdt_UU;

 case 'SS_SU',
  dTC = squeeze(dTdt_SU(1,:,:));
  dTW = squeeze(dTdt_SS(2,:,:));
  dTS = squeeze(dTdt_SS(3,:,:));
  dTB = squeeze(dTdt_SU(4,:,:));
  dTdthc = dTdthc_SS;
  dTdt = dTdt_SS;

 case 'SS_UU',
  dTC = squeeze(dTdt_UU(1,:,:));
  dTW = squeeze(dTdt_SS(2,:,:));
  dTS = squeeze(dTdt_SS(3,:,:));
  dTB = squeeze(dTdt_UU(4,:,:));
  dTdthc = dTdthc_SS;
  dTdt = dTdt_SS;

 otherwise,
  error('Scaling was %s, should be a recognized combination of SS, US, SU, UU',scaling);
end;

dTC(maxdepth<h | h<mindepth) = nan;
dTW(maxdepth<h | h<mindepth) = nan;
dTS(maxdepth<h | h<mindepth) = nan;
dTB(maxdepth<h | h<mindepth) = nan;


winter_min = 12;
winter_mean = 26;
summer_mean = 26;
summer_max = 34;

% % % %TC = winter_mean + dTW*7 + dTC*3;
% % % TC = winter_mean + dTW*23 + dTC*7;
% % % % TC = winter_mean + (4*dTW*5) + (2*dTC*1);
% % % % %TC = winter_mean + dTC*4;
% % TC = winter_mean + dTW*14 + dTC*7;
% TC = winter_mean + dTC*14 + dTW*16;
TC = winter_mean + dTC*4 + dTW*3;
%TC = winter_mean + dTC*4 + dTW*10;
TC = max(TC,(Tas(1)),'includenan');

% % %TW = winter_mean + dTW*10;
% % TW = winter_mean + dTW*21;
% TW = winter_mean + dTW*30;
TW = winter_mean + dTW*14;
TW = max(TW,(Tas(2)),'includenan');

% %TS = summer_mean + dTS*10;
% TS = summer_mean + dTS*21;
%TS = summer_mean + dTS*30;
TS = summer_mean + dTS*14;
TS = min(TS,(Tas(3)+2),'includenan');


% %TB = summer_mean + dTS*7 + dTB*3;
% TB = summer_mean + dTS*23 + dTB*7;
% %%TB = summer_mean + dTB*6;
%TB = summer_mean + dTS*23 + dTB*7;
TB = summer_mean + dTS*7 + dTB*7;
TB = min(TB,(Tas(4)+4),'includenan');


%cldcmp = @min;  	wrmcmp = @max;
cldcmp = @prctile3;	wrmcmp = @prctile97;
%cldcmp = @prctile25;	wrmcmp = @prctile75;
%cldcmp = @median;	wrmcmp = @median;

T{1} = TC;	dts{1} = datenum(2010,[1,2],1);	cmp{1}=cldcmp;	scenario_desc{1} = 'Mass Cold Snap';
%T{2} = TW;	dts{2} = datenum(2009,[1,2],1);	cmp{2}=cldcmp;	scenario_desc{2} = 'Normal Winter';
T{2} = TW;	dts{2} = datenum(1991,[1,2],1);	cmp{2}=cldcmp;	scenario_desc{2} = 'Normal Winter';
T{3} = TS;	dts{3} = datenum(2009,[8,9],1);	cmp{3}=wrmcmp;	scenario_desc{3} = 'Normal Summer';
T{4} = TB;	dts{4} = datenum(1998,[8,9],1);	cmp{4}=wrmcmp;	scenario_desc{4} = 'Mass Bleaching';

for scenario=1:4;
  rmse{scenario} = 0;
  N{scenario} = 0;
end;

ts = []; clear ts
validix = [];
stnms = fieldnames(stns);
for stix=1:numel(stnms)
  stnm = stnms{stix};
  % %if ( isfield(stns.(stnm),'fknms_seatemp') )
  % if ( ~isfield(stns.(stnm),'fknms_seatemp') )
  if ( ~isfield(stns.(stnm),'seatemp') )
    keyboard;
  end;
  % HACK HACK HACK
  % % if ( isfield(stns.(stnm),'fknms_seatemp') && stns.(stnm).lat < 26.0 )
  % if ( isfield(stns.(stnm),'seatemp') && stns.(stnm).lat < 26.0 )
  if ( isfield(stns.(stnm),'seatemp') )
    ts(stix).station_name = upper(stnm);
    ts(stix).lon = stns.(stnm).lon;
    ts(stix).lat = stns.(stnm).lat;
    ts(stix).depth = stns.(stnm).depth;
    if ( isfield(stns.(stnm),'ngdc_depth') )
      ts(stix).ngdc_depth = stns.(stnm).ngdc_depth;
    end;
    ts(stix).ngdc_beta = stns.(stnm).ngdc_beta;
    ts(stix).ngdc_beta_ang = stns.(stnm).ngdc_beta_ang;
    ts(stix).date = stns.(stnm).seatemp.date;
    ts(stix).data = stns.(stnm).seatemp.data;

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
        extC = feval(cmp{scenario},ts(stix).data(ts_date_range(ts(stix),dts{scenario})));
        if ( isempty(extC) )
          extC = nan;
        end;
        aExtC{scenario}(stix) = extC;
        aT{scenario}(stix) = nanmean(nanmean(T{scenario}(iix,jix)));
        aQ0{scenario}(stix) = nanmean(dTdtq0_simple(scenario,iix,jix)*24);
        ahc{scenario}(stix) = nanmean(dTdthc(scenario,iix,jix));
        adT{scenario}(stix) = nanmean(dTdt(scenario,iix,jix));

        tdif = aT{scenario}(stix)-aExtC{scenario}(stix);
        if ( ~isnan(tdif) )
          rmse{scenario} = rmse{scenario} + (tdif^2);
          N{scenario} = N{scenario} + 1;
        end; %if ( ~isnan(tdif) )
      end; %for scenario=1:4;
    end;
  end;
end;

%save(fullfile(get_ecoforecasts_path('data'),'T_tses_1.mat'),'ts');

for scenario=1:4;
  if ( N{scenario} > 0 )
    rmse{scenario} = sqrt(rmse{scenario}/N{scenario});
  else
    rmse{scenario} = nan;
  end;
end;

% % Find the "warmest January" - 1991!
% for yr=1990:2010; 
%   dts=datenum(yr,[1,2],1); 
%   n=0; mn=0; 
%   for stix=1:numel(ts); 
%     ix=ts_date_range(ts(stix),dts); 
%     if numel(ix)>(24*15); 
%       n=n+1; mn=mn+nanmean(ts(stix).data(ix)); 
%     end; 
%   end; 
%   disp({yr,n,roundn(mn/n,-1)}); 
% end;


if ( doFigs )

  % goodix = find( ~isnan(aExtC{1}) & aExtC{1}>0 & ~isnan(aExtC{2}) & aExtC{2}>0 & ...
  %                ~isnan(aExtC{3}) & aExtC{3}>0 & ~isnan(aExtC{4}) & aExtC{4}>0 & ...
  %                ~isnan(ah) & ah>0 & ah<20 & ~isnan(abet) & abet>0 );
  goodix = find( ~isnan(ah) & ah>0 & ah<20 & ~isnan(abet) & abet>0 );
  
  [sorth,sorthix] = sort(ah(goodix)); 
  sorthix = goodix(sorthix);
  scenix{1} = sorthix( ~isnan(aExtC{1}(sorthix)) & aExtC{1}(sorthix)>0 );
  scenix{2} = sorthix( ~isnan(aExtC{2}(sorthix)) & aExtC{2}(sorthix)>0 );
  scenix{3} = sorthix( ~isnan(aExtC{3}(sorthix)) & aExtC{3}(sorthix)>0 );
  scenix{4} = sorthix( ~isnan(aExtC{4}(sorthix)) & aExtC{4}(sorthix)>0 );
  
  fmg;
  plot(ah(scenix{1}),aExtC{1}(scenix{1}),'mo-',ah(scenix{1}),aT{1}(scenix{1}),'ms:');
  plot(ah(scenix{2}),aExtC{2}(scenix{2}),'bo-',ah(scenix{2}),aT{2}(scenix{2}),'bs:');
  plot(ah(scenix{3}),aExtC{3}(scenix{3}),'ro-',ah(scenix{3}),aT{3}(scenix{3}),'rs:');
  plot(ah(scenix{4}),aExtC{4}(scenix{4}),'ko-',ah(scenix{4}),aT{4}(scenix{4}),'ks:');
  titlename('Seafloor depth (h, [m]) vs. sea temperature extremes');
  if doPrint; print('-dpng',fullfile(figspath,[figsbasenm,'-h-vs-extrema.png'])); end;
  
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
  titlename('Seafloor slope (\beta) vs. sea temperature extremes');
  if doPrint; print('-dpng',fullfile(figspath,[figsbasenm,'-beta-vs-extrema.png'])); end;
  set(gca,'xscale','log'); xlim([5e-4,5e-2]);
  if doPrint; print('-dpng',fullfile(figspath,[figsbasenm,'-beta-vs-extrema-log.png'])); end;

end; %if ( doFig )
