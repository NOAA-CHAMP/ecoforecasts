function stn = station_ekman_flux(stn_or_stnm,orifld,varargin)
%function stn = station_ekman_flux(stn_or_stnm,orifld,[ISRC,QSRC|afld,qfld,wfld,dfld,pfld,tfld,ffld,sfld])
%
% Calculate cross-shore component of Ekman flux from surface wind stress.
% Useful in predicting (or disproving) wind-driven upwelling along a coast.
% Isobath orientation ORIFLD (degTrue) may be a fieldname or a scalar. ISRC
% is a data source for meteorological variables, possibly excluding Relative
% Humidity; QSRC is a data source for Relative Humidity (or dew-point temp.)
% Alternatively individual field names may be specified for air temp. (AFLD),
% Rel. Hum. (QFLD), wind speed (WFLD) and direction (DFLD), air pressure
% (PFLD), calculated wind-stress (TFLD), calculated Ekman mass flux [kg/ms]
% (FFLD) and Ekman volumetric flux ([FFLD,'_volume']) [m^2/s]. If field SFLD
% is given or if the default sea temperature field for ISRC exists, use it to
% calculate volumetric flux from mass flux; otherwise, assume T=15oC. 
%
% If QFLD is a CELLSTR (or QSRC is 'ndbc'), the 1st element is a fieldname
% for dew-point temp., the 2nd a fieldname for calculated Relative Humidity.
%
% Also calculates orthogonal components for Ekman mass and volume flux, i.e.,
% those corresponding to "cross-shore" wind stress as well as "alongshore".
% Orthogonal components are labeled, e.g., STN.([FFLD,'_orthog']).
%
% Last Saved Time-stamp: <Tue 2018-03-13 16:44:36 EDT lew.gramer>

  stn = get_station_from_station_name(stn_or_stnm);
  if ( ~isfield(stn,'lat') )
    error('First arg should be a valid station code or STRUCT with .lat field');
  end;
  clear stn_or_stnm;

  if ( ~exist('orifld','var') || isempty(orifld) )
    orifld = 'isobath_orientation';
  end;

  args = varargin;
  nargs = numel(varargin);

  ISRC = '';
  QSRC = '';

  switch (nargs)
   case 0,
   case 1,
    ISRC = args{1};
    QSRC = ISRC;
   case 2,
    ISRC = args{1};
    QSRC = args{2};
   case {3,4,5,6,7}
    afld = args{1};
    qfld = args{2};
    wfld = args{3};
    if ( nargs >= 4); dfld = args{4}; end;
    if ( nargs >= 5); pfld = args{5}; end;
    if ( nargs >= 6); tfld = args{6}; end;
    if ( nargs >= 7); ffld = args{7}; end;
   otherwise,
    error('Either args 3-4 must be two data source names, or args 3-7 must be STRUCT field names');
  end;

  % Data source for most meteorology (possibly excluding Relative Humidity)
  switch (ISRC),
   case '',
    % This is the default case
   case 'ndbc',
    % This is the default case
   case 'erai',
    afld = 'erai_air_t';
    wfld = 'erai_wind_speed';
    dfld = 'erai_wind_dir';
    pfld = 'erai_barom';
    tfld = 'erai_bulk_windstress';
    ffld = 'erai_ekman_flux_x';
    fofld = 'erai_ekman_flux_l';
    sfld = 'erai_sea_t';
   case 'fdep',
    afld = 'fdep_airtemp';
    wfld = 'fdep_wind_speed';
    dfld = 'fdep_wind_dir';
    pfld = 'fdep_airpres';
    tfld = 'fdep_bulk_windstress';
    ffld = 'fdep_ekman_flux_x';
    fofld = 'fdep_ekman_flux_l';
    sfld = 'fdep_seatemp_shallow';
   case 'icon',
    afld = 'airtemp';
    wfld = 'wind1_speed';
    dfld = 'wind1_dir';
    pfld = 'barom';
    tfld = 'icon_bulk_windstress';
    ffld = 'icon_ekman_flux_x';
    fofld = 'icon_ekman_flux_l';
    sfld = 'ctd_shallow_seatemp';
   case 'rsmas',
    afld = 'rsmas_air_t';
    wfld = 'rsmas_wind1_speed';
    dfld = 'rsmas_wind1_dir';
    pfld = 'rsmas_barom';
    tfld = 'rsmas_bulk_windstress';
    ffld = 'rsmas_ekman_flux_x';
    fofld = 'rsmas_ekman_flux_l';
    sfld = 'ndbc_sea_t';
   otherwise,
    error('Do not recognize meteorology data source (ISRC) %s',ISRC);
  end;
  % Data source for Relative Humidity (may differ from meteorology data source)
  switch (QSRC),
   case '',
   case 'ndbc',
    qfld = {'ndbc_dew_t','ndbc_relhumid'};
   case 'erai',
    qfld = 'erai_relhumid';
   case 'fdep',
    qfld = 'fdep_relhumid';
   case 'icon',
    qfld = 'relhumid';
   case 'rsmas',
    qfld = 'rsmas_relhumid';
   otherwise,
    error('Do not recognize humidity data source (QSRC) %s',QSRC);
  end;

  % Assign defaults as needed
  if ( ~exist('afld','var') )	afld = 'ndbc_air_t';			end;
  if ( ~exist('qfld','var') )	qfld = {'ndbc_dew_t','ndbc_relhumid'};	end;
  if ( ~exist('wfld','var') )	wfld = 'ndbc_wind1_speed';		end;
  if ( ~exist('dfld','var') )	dfld = 'ndbc_wind1_dir';		end;
  if ( ~exist('pfld','var') )	pfld = 'ndbc_barom';			end;
  if ( ~exist('tfld','var') )   tfld = 'ndbc_bulk_windstress';          end;
  if ( ~exist('ffld','var') )   ffld = 'ndbc_ekman_flux_x';             end;
  if ( ~exist('fofld','var') )  fofld = 'ndbc_ekman_flux_l';            end;

  if ( ~exist('sfld','var') )   sfld = 'ndbc_sea_t';                    end;
  if ( ~exist('vfld','var') )   vfld = [ffld,'_volume'];                end;
  if ( ~exist('vofld','var') )  vofld = [fofld,'_volume'];              end;

  txfld = [tfld,'_x'];
  tyfld = [tfld,'_y'];
  txsfld = [tfld,'_xs'];
  tlsfld = [tfld,'_ls'];

  dewfld = '';
  if ( iscell(qfld) )
    dewfld = qfld{1};
    qfld = qfld{2};
  end;


  nonqflds = {afld,wfld,dfld,pfld};

  ndbcflds = nonqflds(~cellfun(@isempty,strfind(nonqflds,'ndbc')));
  if ( ~isempty(ndbcflds) && any(~isfield(stn,ndbcflds)) )
    stn = load_all_ndbc_data(stn);
  end;


  allflds = {afld,qfld,dewfld,wfld,dfld,pfld};

  eraiflds = allflds(~cellfun(@isempty,strfind(allflds,'erai')));
  if ( ~isempty(eraiflds) && any(~isfield(stn,eraiflds)) )
    stn = get_erai_station(stn);
  end;

  fdepflds = allflds(~cellfun(@isempty,strfind(allflds,'fdep')));
  if ( ~isempty(fdepflds) && any(~isfield(stn,fdepflds)) )
    stn = read_fdep_stevens_data(stn);
  end;


  if ( ~isempty(dewfld) && ~isfield(stn,qfld) )
    stn = station_dewp_to_relhumid(stn,afld,dewfld,qfld);
  end;

  % Calculate sea-surface wind stress
  if ( ~isfield(stn,tfld) )
    stn = station_bulk_windstress(stn,tfld,wfld,[],dfld,afld,qfld,pfld);
  end;

  % Find prevailing shoreline (or isobath) orientation to true North
  if ( ~isempty(orifld) && ~isfield(stn,orifld) )
    if ( isnumeric(orifld) && isscalar(orifld) )
      x = orifld;
      orifld = 'ekman_flux_isobath_orientation';
      stn.(orifld) = x;
    else
      try,
        stn = station_optimal_isobath_orientation(stn);
      catch,
      end;
    end;
  end;

  % Find cross- and alongshore components of wind stress
  if ( ~isfield(stn,tlsfld) )
    if ( isfield(stn,orifld) )
      stn = station_reorient_vectors(stn,orifld,txfld,tyfld,txsfld,tlsfld);
    else
      if ( isfield(stn,'station_name') )
        disp(stn.station_name);
      end;
      warning('Defaulting alongshore to y direction');
      tlsfld = tyfld;
    end;
  end;

  % Estimate Ekman mass flux from alongshore wind stress and Coriolis
  stn.(ffld).date = stn.(tlsfld).date;
  stn.(ffld).data = stn.(tlsfld).data ./ sw_f(stn.lat);

  stn.(fofld).date = stn.(txsfld).date;
  stn.(fofld).data = stn.(txsfld).data ./ sw_f(stn.lat);

  % Sanity checks
  badix = find(stn.(ffld).data ~= real(stn.(ffld).data) | stn.(fofld).data ~= real(stn.(fofld).data));
  stn.(ffld).date(badix) = [];
  stn.(ffld).data(badix) = [];
  stn.(fofld).date(badix) = [];
  stn.(fofld).data(badix) = [];

  % Estimate Ekman volumetric flux per meter of coastline [m^2/s]
  if ( isfield(stn,sfld) )
    % For what dates do we have both Ekman flux and sea temperature?
    [f,t] = intersect_tses(stn.(ffld),stn.(sfld));
    stn.(vfld).date = f.date;
    stn.(vfld).data = f.data ./ sw_dens0(repmat(35,size(t.data)),t.data);

    [fo,t] = intersect_tses(stn.(fofld),stn.(sfld));
    stn.(vofld).date = fo.date;
    stn.(vofld).data = fo.data ./ sw_dens0(repmat(35,size(t.data)),t.data);
  else
    stn.(vfld).date = stn.(ffld).date;
    stn.(vfld).data = stn.(ffld).data ./ sw_dens0(35,15);

    stn.(vofld).date = stn.(fofld).date;
    stn.(vofld).data = stn.(fofld).data ./ sw_dens0(35,15);
  end;

return;
