1;
% SCRIPT PLOT_RANGE_VS_CONTROL.m
%
% Compare scalar control variable CTLVAR (e.g., seafloor slope, 'ngdc_beta')
% vs. time series response variable RESVAR (e.g., 'fknms_seatemp') for all
% sites in STRUCT dataset STNS, whose site names match PATT (v. GREPSTRUCT).
%
% DEFAULTS: Compare stns.FKNMS_.*.beta vs. whole-distribution Jun-Sep 93rd
% minus 7th PRCTILE ranges in time series stns.FKNMS_*.fknms_seatemp. If
% STRUCT STNS does not exist, call ANFKNMS (v.)
%
% Last Saved Time-stamp: <Tue 2017-04-25 14:09:35 Eastern Daylight Time gramer>


set_more off;

if ( ~exist('doPrint','var') || isempty(doPrint) )
  doPrint = false;
end;
figspath = get_ecoforecasts_path('figs');
if ( doPrint )
  disp(['Will print plots to ',figspath]);
end;

if ( ~exist('doWaves','var') || isempty(doWaves) )
  % NOTE: If DOWAVES is set False, it may be SWITCHED BACK ON below, if
  % caller specified a CTLVAR or a RESVAR containing the substring 'wave'
  doWaves = false;
end;

%% Load set of moored time series (STNS) to process
if ( ~exist('stns','var') )
  anfknms;
end;
stnms = fieldnames(stns);
for cstnm=stnms(:)';
  if ( isfield(stns.(cstnm{:}),'fknms_seatemp') && ~isfield(stns.(cstnm{:}),'seatemp') )
    stns.(cstnm{:}).seatemp = stns.(cstnm{:}).fknms_seatemp;
  end;
end;
clear stnms cstnm

% PATTern to match against STNS fieldnames
if ( ~exist('patt','var') )
  %patt = 'FKNMS_.*';
  % Also include the C-MAN (SEAKEYS) sites in our analysis
  patt = '.*';
end;


%% Control and Response variables

% ConTroL VARiable
if ( ~exist('ctlvar','var') || isempty(ctlvar) )
  % %ctlvar = 'ngdc_beta';
  % %ctlvar = {'ngdc_beta',@nanmean,3,3,8};
  % %ctlvar = {'ngdc_beta',@nanmean,3,3,5};
  ctlvar = {'ngdc_beta',@nanmean,4,4,5};
  %ctlvar = {'ngdc_beta',@nanmean,5,5,5};

  %ctlvar = 'ngdc_depth';
  %ctlvar = {'ngdc_depth',@nanmean,4,4,5};

  %ctlvar = 'raw_nwps_sigwavehgt';
end;
% RESponse VARiable
if ( ~exist('resvar','var') || isempty(resvar) )
  resvar = 'seatemp';
end;

if ( ~doWaves )
  if ( (ischar(ctlvar) && ~isempty(strfind(ctlvar,'wave'))) || ...
       (iscell(ctlvar) && ~isempty(strfind(ctlvar{1},'wave'))) || ...
       (ischar(resvar) && ~isempty(strfind(resvar,'wave'))) || ...   
       (iscell(resvar) && ~isempty(strfind(resvar{1},'wave'))) )
    doWaves = true;
  end;
end;

if ( ~doWaves )
  disp(['NOT using NWPS wave data']);
elseif ( ~isfield(stns.(sites.stnms{1}),'raw_nwps_sigwavehgt') )
  try,
    w = warning('OFF','Ecoforecasts:NWPS:NoFile');
    %stns = get_nwps_stations(stns,[],[],false);
    stns = get_nwps_stations(stns,[],[],true);
    warning(w);
    disp(['Using NWPS wave data']);
  catch,
    disp(['No NWPS wave data AVAILABLE!']);
  end;
end;


% Significant Wave Height cutoff range
if ( ~exist('HsRange','var') )
  %HsRange = [0,0.35];
  %HsRange = [0.35,+Inf];
  HsRange = [];
end;

if ( isempty(HsRange) )
  wstr = '';
  wtxt = '';
elseif ( isnumeric(HsRange) && numel(HsRange) == 2 && HsRange(1) < HsRange(2) )
  wstr = ['_Hs_',num2str(HsRange(1)),'_',num2str(HsRange(2))];
  wtxt = [', Hs\in[',num2str(HsRange(1)),',',num2str(HsRange(2)),']'];
else
  error('HsRange if given must be a numeric range (2-vector)');
end;


%% Processing parameters

% PERiod ACCumulator function
if ( ~exist('peracc','var') )
  peracc = @get_jday;
end;
% PERiod SUBsetting function
if ( ~exist('persub','var') )
  persub = @ts_isfinite;
  %persub = @ts_jfm;
  %persub = @ts_amj;
  %persub = @ts_jas;
  %persub = @ts_ond;
end;

perstr = lower(char(persub));
switch ( perstr ),
 case 'ts_isfinite',   perstr = 'all';
end;

speracc = upper(regexprep(char(peracc),'^(get|ts)_',''));
%spersub = upper(regexprep(char(persub),'^(get|ts)_',''));
spersub = upper(regexprep(perstr,'^(get|ts)_',''));


if ( ~exist('doAccDist','var') )
  % Do Distribution over entire Period?
  doAccDist = false;
  % % Or just do Distribution within JDAY (or other PERACC)
  % doAccDist = true;
end;

if ( ~exist('fillBathGaps','var') || isempty(fillBathGaps) )
  % % For smoothed bathymetry, fill gaps using higher-resolution products
  % fillBathGaps = true;
  % Don't bother
  fillBathGaps = false;
end;


if ( ~exist('bat','var') || ~isfield(bat,'field') )
  % Default F010 bathymetry resolution for south Florida is 30 m: 3 pts. ~ 90 m
  [x,bath] = find_ngdc_slope_sites(stns,[],3);
  x=[]; clear x
  bat = bath.ngdc_hires_bathy;
  bath=[]; clear bath
end;

if 0;
  fmg;
  plot_hires_coastline(bat);
  daspect([1,cosd(25),1]);
  lh=plot(sites.lons,sites.lats,'r.');
  print('-dpng',fullfile(figspath,'plot_range_vs_control_FRT_sites.png'));
end;

%% Control variable, which may be a result of smoothing or other statistical manipulation

% E.g., CTLVAR = {'ngdc_beta',@nanmean,3,3,3}
if ( iscell(ctlvar) )
  smoothvar = ctlvar{1};
  smoothparms = ctlvar(2:end);
  if ( numel(smoothparms) >= 1 )
    smoothvar = [smoothvar,'_',matlab.lang.makeValidName(char(smoothparms{1}))];
    hires_smoothparms{1} = smoothparms{1};
  end;
  if ( numel(smoothparms) >= 2 )
    smoothvar = [smoothvar,'_',num2str(smoothparms{2})];
    hires_smoothparms{2} = smoothparms{2}*3;
  end;
  if ( numel(smoothparms) >= 3 )
    smoothvar = [smoothvar,'_',num2str(smoothparms{3})];
    hires_smoothparms{3} = smoothparms{3}*3;
  end;
  if ( numel(smoothparms) >= 4 )
    smoothvar = [smoothvar,'_',num2str(smoothparms{4})];
    hires_smoothparms{4} = smoothparms{4}*3;
  end;
  if ( numel(smoothparms) >= 5 )
    error('CTLVAR if a cell must have one to four elements (see INTERP_FIELD)');
  end;
  switch ( ctlvar{1} ),
   case 'ngdc_beta',
    fldvar = 'beta';
   case 'ngdc_depth',
    fldvar = 'field';
   otherwise,
    error('CTLVAR if a cell must have first elt. "ngdc_depth" or "ngdc_beta"');
  end;

  % Default F010 bathymetry resolution for south Florida is 30 m: 3 pts. ~ 90 m
  if ( fillBathGaps )
    % The 30 m F010 bathymetry has a rectangular "hole" around Key West
    % Fill that hole in F010 with higher-resolution bathymetry
    kwbath = read_hires_bathymetry('sanf1',[20e3,30e3]);
    % Bathymetry resolution near Key West is 10 m: 9 pts. ~ 90 m
    [x,kwbath] = find_ngdc_slope_sites(stns,kwbath,9);
    x=[]; clear x
    kwbat = kwbath.ngdc_hires_bathy;
    kwbath=[]; clear kwbath
  end;
  for ix = 1:numel(sites.stnms)
    fld = sites.stnms{ix};
    [ig,lonix]=min(abs(bat.lon-sites.lons(ix)));
    lonixen = lonix-35:lonix+35; lonixen(1 > lonixen | lonixen > numel(bat.lon)) = [];
    [ig,latix]=min(abs(bat.lat-sites.lats(ix)));
    latixen = latix-34:latix+34; latixen(1 > latixen | latixen > numel(bat.lat)) = [];
    sbat.lon = bat.lon(lonixen);
    sbat.lat = bat.lat(latixen);
    [sLAT,sLON] = meshgrid(sbat.lat,sbat.lon);
    sbat.field = bat.field(latixen,lonixen);
    sbat.beta = bat.beta(latixen,lonixen);
    sbat.(smoothvar) = interp_field(sbat.lat,sbat.lon,sbat.(fldvar),sLAT,sLON,smoothparms);
    sbat.(smoothvar) = reshape(sbat.(smoothvar),size(sbat.beta'))';
    sbat.(smoothvar)(sbat.field >= -1.0) = nan;
    stns.(fld).(smoothvar) = interp2(sbat.lon,sbat.lat,sbat.(smoothvar),stns.(fld).lon,stns.(fld).lat,'nearest');
    sbat=[]; sLAT=[]; sLON=[]; clear sbat sLAT sLON lonix lonixen latix latixen 

    if ( fillBathGaps && isnan(stns.(fld).(smoothvar)) )
      [ig,lonix]=min(abs(kwbat.lon-sites.lons(ix)));
      lonixen = lonix-35:lonix+35; lonixen(1 > lonixen | lonixen > numel(kwbat.lon)) = [];
      [ig,latix]=min(abs(kwbat.lat-sites.lats(ix)));
      latixen = latix-34:latix+34; latixen(1 > latixen | latixen > numel(kwbat.lat)) = [];
      sbat.lon = kwbat.lon(lonixen);
      sbat.lat = kwbat.lat(latixen);
      [sLAT,sLON] = meshgrid(sbat.lat,sbat.lon);
      sbat.field = kwbat.field(latixen,lonixen);
      sbat.beta = kwbat.beta(latixen,lonixen);
      sbat.(smoothvar) = interp_field(sbat.lat,sbat.lon,sbat.(fldvar),sLAT,sLON,hires_smoothparms);
      sbat.(smoothvar) = reshape(sbat.(smoothvar),size(sbat.beta'))';
      sbat.(smoothvar)(sbat.field >= -1.0) = nan;
      stns.(fld).(smoothvar) = interp2(sbat.lon,sbat.lat,sbat.(smoothvar),stns.(fld).lon,stns.(fld).lat,'nearest');
      sbat=[]; sLAT=[]; sLON=[]; clear sbat sLAT sLON lonix lonixen latix latixen 
    end;
  end; %for ix = 1:numel(sites.stnms)

  ctlvar = smoothvar;
  %disp(['CTLVAR is ',ctlvar]);
end;

% if ( doSave )
%   save('FRT_depth_and_beta.mat','-struct','bat');
% end;


if ( ~doWaves || 1 )
  disp(['NOT using WW3 wave data']);
else
  ww3matfname = [mfilename,'_ww3.mat'];
  if ( exist(ww3matfname,'file') )
    disp(['Loading WW3 wave data from ',ww3matfname]);
    load(ww3matfname);

  else
    disp(['Extracting WW3 wave data']);
    sites.ww3_hs_clim = repmat(nan,[365,numel(sites.lons)]);
    sites.ww3_wu_clim = repmat(nan,[365,numel(sites.lons)]);
    sites.ww3_wv_clim = repmat(nan,[365,numel(sites.lons)]);
    sites.ww3_tp_clim = repmat(nan,[365,numel(sites.lons)]);
    % Process three WW3 regions to get wave climatology at as many sites as we can
    for ww3rgn = {'fks','bsc','sef'};
      % Climatologize NOAA Wave Watch III output surrounding a region RGN
      ww3 = seasonalize_ww3_region(ww3rgn{:},false,false,false);

      % ... Linearly attenuate (scale) field for depths in range 0 to -20 m ...
      % ... And finally, land-mask gridpoints in result as NaN
      [hs,fldlat,fldlon] = oversample_attenuate_field({ww3.lon,ww3.lat,ww3.hs.date,ww3.hs.field},3,[],bat,'depth',true);
      sitehs = interp_field(fldlat,fldlon,hs,sites.lats,sites.lons); hs=[]; clear hs
      wvixen = find(any(~isnan(sitehs)));
      for wvix=wvixen(:)'
    sites.ww3_dates = ;
    stns.(fld).ww3_sigwavehgt.data = sites.ww3_hs(:,ix);
    stns.(fld).ww3_peakwaveu_clim = nanmean(sites.ww3_wu_clim(:,ix));
    stns.(fld).ww3_peakwaveu.date = sites.ww3_dates;
    stns.(fld).ww3_peakwaveu.data = sites.ww3_wu(:,ix);
    stns.(fld).ww3_peakwavev_clim = nanmean(sites.ww3_wv_clim(:,ix));
    stns.(fld).ww3_peakwavev.date = sites.ww3_dates;
    stns.(fld).ww3_peakwavev.data = sites.ww3_wv(:,ix);
    stns.(fld).ww3_peakwaveper_clim = nanmean(sites.ww3_tp_clim(:,ix));
    stns.(fld).ww3_peakwaveper.date = sites.ww3_dates;
    stns.(fld).ww3_peakwaveper.data = sites.ww3_tp(:,ix);


        [cum,tid] = grp_ts(sitehs(:,wvix),ww3.hs.date,[],[],7);
        sites.ww3_hs_clim(tid,wvix) = cum;
      end;
      sitehs=[]; clear sitehs

      wu = oversample_attenuate_field({ww3.lon,ww3.lat,ww3.hs.date,ww3.u},3,[],bat,'depth',true);
      sitewu = interp_field(fldlat,fldlon,wu,sites.lats,sites.lons); wu=[]; clear wu
      wvixen = find(any(~isnan(sitewu)));
      for wvix=wvixen(:)'
        [cum,tid] = grp_ts(sitewu(:,wvix),ww3.hs.date,[],[],7);
        sites.ww3_wu_clim(tid,wvix) = cum;
      end;
      sitewu=[]; clear sitewu

      wv = oversample_attenuate_field({ww3.lon,ww3.lat,ww3.hs.date,ww3.v},3,[],bat,'depth',true);
      sitewv = interp_field(fldlat,fldlon,wv,sites.lats,sites.lons); wv=[]; clear wv
      wvixen = find(any(~isnan(sitewv)));
      for wvix=wvixen(:)'
        [cum,tid] = grp_ts(sitewv(:,wvix),ww3.hs.date,[],[],7);
        sites.ww3_wv_clim(tid,wvix) = cum;
      end;
      sitewv=[]; clear sitewv

      tp = oversample_attenuate_field({ww3.lon,ww3.lat,ww3.tp.date,ww3.tp.field},3,[],bat,'depth',true);
      sitetp = interp_field(fldlat,fldlon,tp,sites.lats,sites.lons); tp=[]; clear tp
      wvixen = find(any(~isnan(sitetp)));
      for wvix=wvixen(:)'
        [cum,tid] = grp_ts(sitetp(:,wvix),ww3.tp.date,[],[],7);
        sites.ww3_tp_clim(tid,wvix) = cum;
      end;
      sitetp=[]; clear sitetp

      ww3=[]; fldlat=[]; fldlon=[]; clear ww3 fldlat fldlon
    end;

    disp(['Saving ',ww3matfname]);
    save(ww3matfname,'sites');
  end;

  for ix = 1:numel(sites.stnms)
    fld = sites.stnms{ix};
    stns.(fld).ww3_sigwavehgt_clim = nanmean(sites.ww3_hs_clim(:,ix));
    stns.(fld).ww3_sigwavehgt.date = sites.ww3_dates;
    stns.(fld).ww3_sigwavehgt.data = sites.ww3_hs(:,ix);
    stns.(fld).ww3_peakwaveu_clim = nanmean(sites.ww3_wu_clim(:,ix));
    stns.(fld).ww3_peakwaveu.date = sites.ww3_dates;
    stns.(fld).ww3_peakwaveu.data = sites.ww3_wu(:,ix);
    stns.(fld).ww3_peakwavev_clim = nanmean(sites.ww3_wv_clim(:,ix));
    stns.(fld).ww3_peakwavev.date = sites.ww3_dates;
    stns.(fld).ww3_peakwavev.data = sites.ww3_wv(:,ix);
    stns.(fld).ww3_peakwaveper_clim = nanmean(sites.ww3_tp_clim(:,ix));
    stns.(fld).ww3_peakwaveper.date = sites.ww3_dates;
    stns.(fld).ww3_peakwaveper.data = sites.ww3_tp(:,ix);
  end;

end;

% Descriptive strings (raw, and textized for MATLAB TeX)
sctlvar = upper(ctlvar);
tctlvar = textize(sctlvar);
if ( doAccDist )
  sresvar = upper([spersub,' MEDIAN(',resvar,' ',speracc,' ',')']);
  tresvar = upper([textize(spersub),' MEDIAN(',textize(resvar),'_{',textize(speracc),'}',')']);
  dstr = 'diurnal_';
else
  sresvar = [upper(spersub),' ',upper(resvar),wtxt];
  tresvar = textize(sresvar);
  dstr = '';
end;



%% Process control and response variables for each site (field) in STNS (struct)

disp([sctlvar,' VS. ',sresvar]);

lons = [];
lats = []; 
stnms = {};
ctlvars = [];
resvars = [];
alt_ctlvars = [];
alt2_ctlvars = [];

cstnms = grepstruct(stns,patt)';
Hses = repmat(nan,[1,numel(cstnms)]);
for stnmix = 1:numel(cstnms)
  stnm = cstnms{stnmix};

  %DEBUG:  if ( ~strcmp(stnm,'FKNMS_BROAD_CRK') ); disp(['DEBUG: SKIPPING ',stnm]); continue; end;

  %DEBUG:
  if ( strcmp(stnm,'FKNMS_WELLWOOD') ); disp(['DEBUG: SKIPPING ',stnm]); continue; end;
  %DEBUG:
  if ( strcmp(stnm,'LONF1') ); disp(['DEBUG: SKIPPING ',stnm]); continue; end;

  if ( ~isfield(stns.(stnm),resvar) || ~is_valid_ts(stns.(stnm).(resvar)) )
    disp(['SKIPPING ',stnm,': No valid ',resvar]);
  else
    if ( isfield(sites,'hs') && any(~isnan(sites.ww3_hs_clim(:,stnmix))) )
      %DEBUG:      disp([stnm,' =? ',sites.stnms{stnmix}]);
      Hses(stnmix) = nanmean(sites.ww3_hs_clim(:,stnmix));
    elseif ( isfield(stns.(stnm),'raw_nwps_sigwavehgt') )
      Hses(stnmix) = nanmean(stns.(stnm).raw_nwps_sigwavehgt.data);
    end;
    if ( ~isnan(Hses(stnmix)) && ~isempty(HsRange) && isnumeric(HsRange) )
      if ( HsRange(1) > Hses(stnmix) || Hses(stnmix) > HsRange(2) )
        disp(['EXCLUDING ',stnm,': Hs ',num2str(Hses(stnmix)),' outside [',num2str(HsRange(1)),',',num2str(HsRange(2)),']']);
        %%%% EARLY LOOP TERMINATION
        continue;
      elseif ( ~isfield(stns.(stnm),'seatemp') )
        disp(['EXCLUDING ',stnm,': No SEATEMP field']);
        %%%% EARLY LOOP TERMINATION
        continue;
      elseif ( ~isfield(stns.(stnm),'ngdc_beta') || isnan(stns.(stnm).ngdc_beta) )
        disp(['EXCLUDING ',stnm,': No BETA field']);
        %%%% EARLY LOOP TERMINATION
        continue;
      end; %if HsRange elseif elseif
    end; %if ~isnan

    %disp(stnm);
    lons(end+1) = stns.(stnm).lon;
    lats(end+1) = stns.(stnm).lat;
    stnms{end+1} = textize(stnm);

    try,
      stns.(stnm) = station_ngdc_offshore_slope(stns.(stnm));
      alt_ctlvars(end+1) = stns.(stnm).ngdc_offshore_slope;
    catch,
      %disp(['No precalculated slope for ',stnm]);
      alt_ctlvars(end+1) = nan;
    end;
    alt2_ctlvars(end+1) = stns.(stnm).ngdc_beta;

    ts = subset_ts(stns.(stnm).(resvar),persub);
    if ( strcmpi(ctlvar,'depth') )
      ctlvars(end+1) = -stns.(stnm).(ctlvar);
    elseif ( is_ts(stns.(stnm).(ctlvar)) )
      ctlvars(end+1) = nanmean(stns.(stnm).(ctlvar).data);
    elseif ( isscalar(stns.(stnm).(ctlvar)) )
      ctlvars(end+1) = stns.(stnm).(ctlvar);
    else
      error('CTLVAR must be a scalar or time-series field of station %s',stnm);
    end;

    if ( ~doAccDist )
      p07d = prctile(ts.data, 7) - mean(ts.data);
      p93d = prctile(ts.data,93) - mean(ts.data);
      resvars(end+1) = p93d - p07d;
    else
      %DEBUG:      if ( strcmp(stnm,'FKNMS_BROAD_CRK') ); keyboard; end;
      p07d = []; p93d = [];
      pers = unique(peracc(ts.date));
      for ix=1:numel(pers);
        per = pers(ix);
        perix = find(peracc(ts.date)==per);
        p07d(ix) = prctile(ts.data(perix), 7) - mean(ts.data(perix));
        p93d(ix) = prctile(ts.data(perix),93) - mean(ts.data(perix));
      end; %for ix=1:numel(pers);
      clear ix
      %resvars(end+1) = nanmax(p93d - p07d);
      %resvars(end+1) = prctile(p93d - p07d,93);
      resvars(end+1) = nanmedian(p93d - p07d);
      pers = []; clear per perix pers
    end; %if ( ~doAccDist ) else
    p07d=[]; p93d=[]; ts=[]; clear p07d p93d ts

  end; %if ~isfield else

end; %for stnmix = 1:numel(cstnms)


%% Plot statistical relationships

%plot_diurnal_range_vs_control(tses,ctlvars,peracc,persub,tsnms,'beta')

printname =  [upper(ctlvar),'_vs_',perstr,'_',resvar,'_',dstr,'range',wstr];


trim_stnms = capitalize(lower(strrep(stnms,'FKNMS\_',' ')));
%trim_stnms = strvcat(trim_stnms,' (',num2str(alt_ctlvars),')');
for ix = 1:numel(trim_stnms)
  if ( ~isnan(alt_ctlvars(ix)) )
    trim_stnms{ix} = [trim_stnms{ix},' (',num2str(alt_ctlvars(ix)),')'];
  end;
end;

minctl = nanmin(ctlvars(~isnan(resvars))) * 0.9;
maxctl = nanmax(ctlvars(~isnan(resvars))) * 1.1;
minvar = nanmin(resvars(~isnan(ctlvars))) * 0.9;
maxvar = nanmax(resvars(~isnan(ctlvars))) * 1.1;

scatter_fit(ctlvars,resvars,tctlvar,tresvar);
text(ctlvars,resvars,trim_stnms);
axis([minctl,maxctl,minvar,maxvar]);
legend('Location','Best');
if ( doPrint )
  print('-dpng',fullfile(figspath,[printname,'.png']));
end;


%% Relationship need not be LINEAR

% scatter_curve_fit(ctlvars',resvars','a/(x-b)',tctlvar,tresvar);
% scatter_curve_fit(ctlvars',resvars','a*(x-b)^n',tctlvar,tresvar); % ERROR: NEEDS BOUNDS

% % Net sea-surface buoyancy flux B_0
% res.B0 = (g.*alph.*abs(q0))./(rhoCp);
% % "Characteristic" convective velocity (Sturmann)
% res.uf = (res.B0.*h).^(1/3);

% % Volumetric flow [m^2/s]: Panel letter refers to Fig. 10, Monismith et al (2006)
% % Advective inertial (steady) momentum balance, steady thermal balance: Panel (c)
% res.Qv_SS = res.uf.*h./(bet.^(1/3));
try,
  %scatter_curve_fit(ctlvars',resvars','a*(x-b)^n',tctlvar,tresvar,[],'StartPoint',[1,0,-1/3]); %SS
  scatter_curve_fit(ctlvars',resvars','a*(x-b)^n',tctlvar,tresvar,[],'ConfInt',0.67,'StartPoint',[1,0,-1],'Robust','LAR','MaxFunEvals',5e4,'MaxIter',5e3,'TolFun',1e-7); %SS
  th = text(ctlvars,resvars,trim_stnms);
  axis([minctl,maxctl,minvar,maxvar]);
  legend('Location','Best');
  if ( doPrint )
    print('-dpng',fullfile(figspath,[printname,'_hc_SS.png']));
  end;
catch ME,
  catchwarn(ME,'Exponential -1/3 fit failure');
end;

% % % Viscous (unsteady) inertial balance, balanced thermal forcing: Panel (a)
% % res.Qv_US = sqrt((res.uf.^3) .* (24.*dt) .* h);

% % Advective inertial balance, unbalanced thermal forcing: Panel (f)
% res.Qv_SU = (bet.^(2/3)).*res.uf.*((res.uf.*(24.*dt)./h).^(3/2));
%{
try,
  %scatter_curve_fit(ctlvars',resvars','a*(x-b)^n',tctlvar,tresvar,[],'StartPoint',[1,0,+2/3]); %SU
  scatter_curve_fit(ctlvars',resvars','a*(x-b)^n',tctlvar,tresvar,[],'ConfInt',0.67,'StartPoint',[1,0,+2/3],'Robust','LAR','MaxFunEvals',1e5); %SU
  th = text(ctlvars,resvars,trim_stnms);
  axis([minctl,maxctl,minvar,maxvar]);
  if ( doPrint )
    print('-dpng',fullfile(figspath,[printname,'_hc_SU.png']));
  end;
catch ME,
  catchwarn(ME,'Exponential +2/3 fit failure');
end;
%}

% % % Viscous (unsteady) inertia, unbalanced thermal forcing: Panel (d)
% % res.Qv_UU = bet.*(res.uf.^3).*((24.*dt)^2)./h;

% % And NOTE: Beta and depth are independent of one another!
% scatter_fit(sites.depths,sites.betas)
% scatter_fit(sites.depths(sites.depths>-10),sites.betas(sites.depths>-10))

clear ans cstnms stnm stnmix

set_more;
