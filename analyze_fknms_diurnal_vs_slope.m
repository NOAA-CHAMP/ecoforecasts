1;

if ( ~exist('stns','var') || ~isfield(stns,'FKNMS_MOLASSES') )
  stns = []; clear stns;
  stns = get_fknms_thermistors;
end;

cumfun = @get_season;

% Local inertial period in hours
[IP_day,IP_hr] = inertial_period(stns.FKNMS_MOLASSES.lat);

% High-pass limit - the local inertial period
hplim = ceil(IP_hr);
hpfld = sprintf('fknms_seatemp_%d_h_hp',hplim);
disp(['High-pass filter: ',hpfld]);

nstds=0;
csts=fieldnames(stns);
for stix=1:numel(csts); 
  stnm = csts{stix}; 
  stn = stns.(stnm); 
  if (isfield(stn,'fknms_seatemp')); 
    stn = verify_variable(stn,hpfld); 
    dt = median(diff(stn.fknms_seatemp.date));
    minN = 0.90 * min_n_ts_filter(cumfun,dt);
    [st,ti] = grp_ts(stn.(hpfld).data,stn.(hpfld).date,cumfun,@nanstd,minN); 
    if ( ~isempty(ti) )
      nstds = nstds + 1;

      stds.nm{nstds} = stnm; 
      stds.st(nstds,1:4) = NaN; 
      stds.st(nstds,ti) = st; 
      stns.(stnm).(hpfld) = stn.(hpfld);
      stns.(stnm).fknms_seatemp_seasonal_diurnal = stds.st(nstds,:);

      bathfld = 'ngdc_hires_bathy';
      if ( ~isfield(stn,bathfld) )
        stn = read_hires_bathymetry(stn,[1e3,1e3]);
      end;
      slope_pts = 7; slope_method = [];
      [stn.slope,stn.slope_orientation,stn.isobath_orientation,stn.(bathfld)] = ...
          find_ngdc_slope(stn.(bathfld),stn.lon,stn.lat,slope_pts,slope_method);
      stds.sl(nstds) = stn.slope;
      stds.io(nstds) = stn.isobath_orientation;
      stns.(stnm).bottom_slope = stn.slope;
      stns.(stnm).isobath_orientation = stn.isobath_orientation;
    end;
  end; 
  stn=[]; clear stn; 
end; 
clear ans cst csts st stix stnm ti
