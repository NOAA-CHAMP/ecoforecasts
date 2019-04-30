1;
%% SCRIPT COMPARE_SPATIAL_DT_HC.m:
%
% Calculate table comparing predicted sea temperature extremes calculated by 
% CALC_SPATIAL_DT_HC to observed sea temperatures from FKNMS, SEAKEYS, Lirman
% data sets, based on SCENARIO # (1=2010 Cold Snap, 2=Normal Winter, 3=Normal
% Summer, 4=1998 Bleaching), using Horizontal Convection SCALING (SS, US, SU,
% UU, or for the extreme cases, SS_SU, or SS_UU). If DOREPORT (DEFAULT: true)
% DISP table of warming rates and ground-truth predicted temperature extreme.
%
% Last Saved Time-stamp: <Thu 2017-11-16 21:39:57 Eastern Standard Time gramer>

if ( ~exist('scenario','var') || isempty(scenario) )
  scenario = 1;
end;
if ( ~exist('doReport','var') || isempty(doReport) )
  doReport = true;
end;
if ( ~exist('scaling','var') || isempty(scaling) )
  scaling = 'SS';
end;

if ( ~exist('lon','var') && ~exist('lat','var') && ~exist('h','var') )
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
      sites.betas(end+1) = bet(iix,jix);
      sites.beta_angs(end+1) = ang(iix,jix);
    else
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
 case 'US',
  dTC = squeeze(dTdt_US(1,:,:));
  dTW = squeeze(dTdt_US(2,:,:));
  dTS = squeeze(dTdt_US(3,:,:));
  dTB = squeeze(dTdt_US(4,:,:));
 case 'SU',
  dTC = squeeze(dTdt_SU(1,:,:));
  dTW = squeeze(dTdt_SU(2,:,:));
  dTS = squeeze(dTdt_SU(3,:,:));
  dTB = squeeze(dTdt_SU(4,:,:));
 case 'UU',
  dTC = squeeze(dTdt_UU(1,:,:));
  dTW = squeeze(dTdt_UU(2,:,:));
  dTS = squeeze(dTdt_UU(3,:,:));
  dTB = squeeze(dTdt_UU(4,:,:));
 case 'SS_SU',
  dTC = squeeze(dTdt_SU(1,:,:));
  dTW = squeeze(dTdt_SS(2,:,:));
  dTS = squeeze(dTdt_SS(3,:,:));
  dTB = squeeze(dTdt_SU(4,:,:));
 case 'SS_UU',
  dTC = squeeze(dTdt_UU(1,:,:));
  dTW = squeeze(dTdt_SS(2,:,:));
  dTS = squeeze(dTdt_SS(3,:,:));
  dTB = squeeze(dTdt_UU(4,:,:));
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


% cldcmp = @min;  	wrmcmp = @max;
cldcmp = @prctile3;	wrmcmp = @prctile97;

switch (scenario)
 case 1,	T = TC;	dts = datenum(2010,[1,2],1);	cmp=cldcmp;	scenario_desc = 'Mass Cold Snap';
 case 2,	T = TW;	dts = datenum(2009,[1,2],1);	cmp=cldcmp;	scenario_desc = 'Normal Winter';
 case 3,	T = TS;	dts = datenum(2009,[8,9],1);	cmp=wrmcmp;	scenario_desc = 'Normal Summer';
 case 4,	T = TB;	dts = datenum(1998,[8,9],1);	cmp=wrmcmp;	scenario_desc = 'Mass Bleaching';
end;

ts = []; clear ts
validix = [];
stnms = fieldnames(stns);
for stix=1:numel(stnms)
  stnm = stnms{stix};
  %if ( isfield(stns.(stnm),'fknms_seatemp') )
  % HACK HACK HACK
  if ( isfield(stns.(stnm),'fknms_seatemp') && stns.(stnm).lat < 26.0 )
    ts(stix) = stns.(stnm).fknms_seatemp;
    extC = feval(cmp,ts(stix).data(ts_date_range(ts(stix),dts)));
    if ( ~isempty(extC) && ~isnan(extC) )
      aExtC(stix) = extC;
      [err,ix] = min(abs(LON(:)-stns.(stnm).lon)+abs(LAT(:)-stns.(stnm).lat));
      [jix,iix] = ind2sub(size(LON),ix);
      if ( 1<=jix && jix<=size(h,2) && 1<=iix && iix<size(h,1) )
        if ( isnan(h(iix,jix)) )
          % iix = iix-1:iix+1;
          % jix = jix-1:jix+1;
        end;
        ah(stix) = nanmean(nanmean(h(iix,jix)));
        abet(stix) = nanmean(nanmean(bet(iix,jix)));
        ahcrng(stix) = nanmean(nanmean(hcrng(iix,jix)));
        ahch(stix) = nanmean(nanmean(hch(iix,jix)));
        aT(stix) = nanmean(nanmean(T(iix,jix)));
        aQ0(stix) = nanmean(dTdtq0_simple(scenario,iix,jix)*24);
        ahc(stix) = nanmean(dTdthc_SS(scenario,iix,jix));
        adT(stix) = nanmean(dTdt_SS(scenario,iix,jix));
        if ( ~isnan(ah(stix)) )
          validix(end+1) = stix;
        end;
      end;
    end;
  end;
end;

if ( doReport )
  [ig,sortedix] = sort(aExtC(validix)); sortedix = validix(sortedix);
  disp([scenario_desc,' (',datestr(dts(1),'mmm yyyy'),')']);
  disp(sprintf('%-15s: % 6s  % 5s  % 6s  % 6s  % 7s  % 7s  % 5s  % 5s  % 5s','Station',...
               'Depth','Beta','HCdep','HCrng','Q0','HC','PredT','ObsT','TDif'));

  rmse=0; N=0; gts=0; lts=0;
  for stix=sortedix(:)'
    stnm = stnms{stix};
    tdif = aT(stix)-aExtC(stix);
    if ( tdif > 1 ); ltGT=sprintf('%+5.1f',tdif); gts=gts+1;
    elseif ( tdif < -1 ); ltGT=sprintf('%+5.1f',tdif); lts=lts+1;
    elseif ( ~isnan(tdif) ); ltGT=sprintf('% 5s','---');
    else; ltGT=sprintf('% 5s','NaN'); end;
    if ( ~isnan(tdif) )
      rmse = rmse + (tdif^2);
      N = N + 1;
    end;
    disp(sprintf('%-15s: % 6.1f  % 5.3f % 6.1f  % 6.0f  % 7.3f  % 7.3f  % 5.1f  % 5.1f  %s',...
                 stnm,ah(stix),abet(stix),ahch(stix),ahcrng(stix),aQ0(stix),ahc(stix),aT(stix),aExtC(stix),ltGT));
  end;
  rmse = sqrt(rmse/numel(sortedix));
  disp(sprintf('%g >, %g <, %g ~ : RMSE=%g',gts,lts,numel(validix)-gts-lts,rmse));
end;
