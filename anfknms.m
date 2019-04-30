1;
% SCRIPT ANFKNMS.m
%
% Load all available sea temperature time series for Florida Keys National
% Marine Sanctuary network of bottom thermistors (H. Hudson et al.) Download
% bathymetry for all sites possible, and calculate seafloor slope ('beta'),
% slope angle, and local isobath angle, i.e., "alongshore direction". (Useful
% in calls to PLOT_RANGE_VS_CONTROL, to load and prepare data for analysis.)
%
% Last Saved Time-stamp: <Fri 2017-08-11 17:52:04 Eastern Daylight Time gramer>

matfname = fullfile(get_ecoforecasts_path('data'),'anfknms.mat');

if ( exist(matfname,'file') )
  disp(['Loading ',matfname]);
  load(matfname);

else
  if ( ~exist('stns','var') )
    stns = get_fknms_thermistors;
    stns = get_langdon_thermistors(stns);
    stns = rmfield(stns,{'lons','lats'});
    stnms = fieldnames(stns);
    for cstnm=stnms(:)';
      if ( isfield(stns.(cstnm{:}),'fknms_seatemp') )
        stns.(cstnm{:}).seatemp = stns.(cstnm{:}).fknms_seatemp;
      end;
    end;
    clear stnms cstnm

    stn = read_fdep_stevens_data('k');
    stnm = 'FDEPK';
    stns.(stnm).station_name = stn.station_name;
    stns.(stnm).lon = stn.lon;
    stns.(stnm).lat = stn.lat;
    stns.(stnm).depth = 2; % Just a guess
    stns.(stnm).seatemp = stn.fdep_seatemp_shallow;
    stn = []; clear stn

    for cst = {'LKWF1','PVGF1','FWYF1','MLRF1','LONF1','SMKF1','SANF1','PLSF1','DRYF1'};
      stnm = cst{:};
      if ( isfield(stns,stnm) )
        warning('Station %s already in STNS??',stnm);
      else
        %disp(stnm);
        stn = get_station_from_station_name(stnm);
        stn = load_all_ndbc_data(stn);
        stns.(stnm).station_name = stn.station_name;
        stns.(stnm).lon = stn.lon;
        stns.(stnm).lat = stn.lat;
        stns.(stnm).depth = stn.depth;
        stns.(stnm).seatemp = stn.ndbc_sea_t;
        stn = []; clear stn
      end;
    end;
    clear cst stnm

  end;

  if ( ~exist('sites','var') )
    % Get depth and slope according to bathymetric data
    % Default F010 bathymetry resolution for south Florida is 30 m: 3 pts. ~ 90 m
    sites = find_ngdc_slope_sites(stns,[],7); 

    % The 30 m F010 bathymetry has a rectangular "hole" around Key West
    nanix = find(isnan(sites.depths) | isnan(sites.betas));
    if ( ~isempty(nanix) )
      % Fill that hole in F010 with higher-resolution bathymetry
      kwbath = read_hires_bathymetry('sanf1',[20e3,30e3]);
      kwdepths = interp2(kwbath.ngdc_hires_bathy.lon,kwbath.ngdc_hires_bathy.lat,...
                         kwbath.ngdc_hires_bathy.field,sites.lons,sites.lats,'*linear',nan);
      % Bathymetry resolution near Key West is 10 m: 9 pts. ~ 90 m
      [kwbetas,kwbeta_angs,kwiso_angs,kwbath.ngdc_hires_bathy] = ...
          find_ngdc_slope(kwbath.ngdc_hires_bathy,sites.lons,sites.lats,9);
      
      sites.depths(nanix) = kwdepths(nanix);
      sites.betas(nanix) = kwbetas(nanix);
      sites.beta_angs(nanix) = kwbeta_angs(nanix);
      sites.iso_angs(nanix) = kwiso_angs(nanix);
      
      kwbath=[]; clear kwbath
      clear ans kwbetas kwbeta_angs kwdepths kwiso_angs nanix 
    end; %if ( ~isempty(nanix) )
    
    % The 30 m F010 bathymetry also has a restrictive upper latitude bound
    nanix = find(isnan(sites.depths) | isnan(sites.betas));
    if ( ~isempty(nanix) )
      % Fill in beyond that upper bound with higher-resolution bathymetry
      mibath = read_hires_bathymetry('fwyf1',[15e3,20e3]);
      midepths = interp2(mibath.ngdc_hires_bathy.lon,mibath.ngdc_hires_bathy.lat,...
                         mibath.ngdc_hires_bathy.field,sites.lons,sites.lats,'*linear',nan);
      % Bathymetry resolution near Key West is 10 m: 9 pts. ~ 90 m
      [mibetas,mibeta_angs,miiso_angs,mibath.ngdc_mires_bathy] = ...
          find_ngdc_slope(mibath.ngdc_hires_bathy,sites.lons,sites.lats,9);
      
      sites.depths(nanix) = midepths(nanix);
      sites.betas(nanix) = mibetas(nanix);
      sites.beta_angs(nanix) = mibeta_angs(nanix);
      sites.iso_angs(nanix) = miiso_angs(nanix);
      
      mibath=[]; clear mibath
      clear ans mibetas mibeta_angs midepths miiso_angs nanix 
    end; %if ( ~isempty(nanix) )
    
    % If we are STILL missing betas, try the old tried and true...
    nanix = find(isnan(sites.depths) | isnan(sites.betas));
    for ix = nanix(:)'
      lobath.lon = sites.lons(ix);
      lobath.lat = sites.lats(ix);
      lobath = read_hires_bathymetry(lobath,[2e3,2e3],[],false);
      sites.depths(ix) = interp2(lobath.ngdc_hires_bathy.lon,lobath.ngdc_hires_bathy.lat,...
                                 lobath.ngdc_hires_bathy.field,sites.lons(ix),sites.lats(ix),'*linear',nan);
      % NGDC 3" Bathymetry resolution is 92 m: 2 pts. ~ 180 m
      [sites.betas(ix),sites.beta_angs(ix),sites.iso_angs(ix),lobath.ngdc_lores_bathy] = ...
          find_ngdc_slope(lobath.ngdc_hires_bathy,sites.lons(ix),sites.lats(ix),2);
      
      lobath=[]; clear lobath
    end;
    clear ix nanix
    
  end;
  
  for ix = 1:numel(sites.stnms)
    fld = sites.stnms{ix};
    stns.(fld).ngdc_depth = sites.depths(ix);
    stns.(fld).ngdc_beta = sites.betas(ix);
    stns.(fld).ngdc_beta_ang = sites.beta_angs(ix);
    stns.(fld).ngdc_iso_ang = sites.iso_angs(ix);
  end;
  
  % sites=[]; clear sites;
  clear ix fld

  disp(['Saving ',matfname]);
  save(matfname,'-v7.3');
end;

clear ans matfname
