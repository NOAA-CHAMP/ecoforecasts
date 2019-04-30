function [stn,Q0_factor] = station_heat_flux(stn,wfld,afld,qfld,pfld,tfld,sfld,lfld,PFX,dsfld,dlfld,prfld,wdfld,oufld,ovfld,wpfld,whfld,pblzfld,doWarm,max_wl,sal,doPlot)
%function [stn,Q0_factor] = station_heat_flux(stn,wfld,afld,qfld,pfld,tfld,sfld,lfld,PFX,dsfld,dlfld,prfld,wdfld,oufld,ovfld,wpfld,whfld,pblzfld,doWarm,max_wl,sal,doPlot)
%
% Use Fairall et al. (1996, 2003) methods to calculate wind stress, sensible
% and latent surface heat flux, from station data struct STN, using in situ
% data from fields WFLD (wind speed, KTS), AFLD (air temp), QFLD (*RELATIVE*
% Humidity), PFLD (air pressure), and TFLD (sea temp). Optionally, arg QFLD
% may be a scalar value, or the name of a dew-point temperature field (as in
% existing SEAKEYS stations) to approximate relative humidity. New fields are
% added to STN, with optional fieldname prefix PFX if specified, e.g., [PFX
% '_wind_stress'], and so on for 'sensible_flux', 'latent_flux', 'rain_flux',
% each with 'date' subfields interpolated to hourly time series covering the
% intersection of all the input fields' timestamps. Similarly, if shortwave
% (SFLD) and longwave (LFLD) fields with valid data are specified, a field
% STN.net_flux (or STN.[PFX '_net_flux']) is calculated with these terms
% added; and field 'flux_term' is added containing net heat flux scaled
% [3600*Q0/rho.Cp.h, K/hr] as a heat budget term (Wilson-Diaz et al., 2009).
%
% By default, the TOGA-COARE 2.0 bulk algorithm is used to calculate latent
% and sensible fluxes, by calling AIR_SEA toolbox function HFBULKTC (qv). If
% field names for downward short- (insolation, DSFLD) and downward longwave
% (DLFLD) flux are specified, use TOGA-COARE 2.6 algorithm with a cool-skin
% correction (via appropriate call to HFBULKTC). If precipation rate [mm/hr]
% (PRFLD) is also specified and a non-empty string, TOGA-COARE 2.6 algorithm
% is called via COR30A function (Fairall et al.) instead. In this case, new
% 'rain_flux' field is also added, and summed into 'net_flux' (and '_term')
% field. If fields WDFLD, OUFLD and OVFLD are given, ocean U and V currents
% are projected onto wind direction WD in call to COR30A; otherwise, surface
% current in wind direction is just estimated at 2.0% of wind speed. Also, if
% dominant wave period (WPFLD) and height [m] (WHFLD) data are given, use the
% full TOGA-COARE 3.0a algorithm for fluxes (also via COR30A), with roughness
% scale z0 derived from wave data as per Oost et al. (2001). If no planetary
% boundary layer height (PBLZFLD) is given, DEFAULT is 600m. Finally, if arg
% DOWARM is TRUE, apply an iterative warm-layer correction to fluxes (calls
% COR30A_WARM, code derived from Fairall script COR3_0AH). If arg MAX_WL is
% non-empty, use it as the max warm-layer depth in COR30a.
%
% If optional arg SAL (DEFAULT: 35) is a scalar or STN time series fieldname,
% use that value(s) as salinity (PSU) in cold-skin calculation (v. HFBULKTC).
%
% Last Saved Time-stamp: <Sat 2016-04-09 16:45:56 Eastern Daylight Time gramer>

  if ( ~exist('sfld','var') || isempty(sfld) )
    sfld = '';
  end;
  if ( ~exist('lfld','var') || isempty(lfld) )
    lfld = '';
  end;
  if ( ~exist('PFX','var') || isempty(PFX) )
    PFX = '';
  elseif ( PFX(end) ~= '_' )
    PFX = [ PFX '_' ];
  end;

  if ( ~exist('dsfld','var') || isempty(dsfld) )
    dsfld = '';
  end;
  if ( ~exist('dlfld','var') || isempty(dlfld) )
    dlfld = '';
  end;
  if ( ~exist('prfld','var') || isempty(prfld) )
    prfld = '';
  end;
  if ( ~exist('wdfld','var') || isempty(wdfld) )
    wdfld = '';
  end;
  if ( ~exist('oufld','var') || isempty(oufld) )
    oufld = '';
  end;
  if ( ~exist('ovfld','var') || isempty(ovfld) )
    ovfld = '';
  end;
  if ( ~exist('wpfld','var') || isempty(wpfld) )
    wpfld = '';
  end;
  if ( ~exist('whfld','var') || isempty(whfld) )
    whfld = '';
  end;
  if ( ~exist('pblzfld','var') || isempty(pblzfld) )
    % DEFAULT applied further down
    pblzfld = '';
  end;
  if ( ~exist('doWarm','var') || isempty(doWarm) )
    doWarm = false;
  end;
  if ( ~exist('max_wl','var') || isempty(max_wl) )
    max_wl = [];
  end;
  if ( ~exist('sal','var') || isempty(sal) )
    % ASSUMES CONSTANT SALINITY 35psu
    sal = 35;
  end;
  if ( ~exist('doPlot','var') || isempty(doPlot) )
    doPlot = false;
  end;

  onehour = (1.0/24.0);

  % Get height [m] of each meteorology sensor at the station, and depth of
  % 'shallow' (or default) sea temperature sensor STZ
  try
    [wz,az,pz,stz] = station_instrument_heights(stn.station_name);
  catch
    warning('No in situ instrument heights/depths for %s',stn.station_name);
  end;

  % REANALYSIS and SATELLITE dataset fields are at ~standard heights/depth
  if ( ~isempty(strfind(wfld,'ncep_')) || ~isempty(strfind(wfld,'erai_')) )
    wz = 10;
  end;
  if ( ~isempty(strfind(afld,'ncep_')) || ~isempty(strfind(afld,'erai_')) )
    az = 2;
  end;
  if ( ~isempty(strfind(pfld,'ncep_')) || ~isempty(strfind(pfld,'erai_')) )
    pz = 2;
  end;
  if ( ~isempty(strfind(tfld,'misst_')) )
    stz = 5;
  elseif ( ~isempty(strfind(tfld,'seatemp')) )
    stz = stn.depth;
  elseif ( ~isempty(strfind(tfld,'sst_')) || ~isempty(strfind(tfld,'modis_')) || ~isempty(strfind(tfld,'avhrr_')) )
    stz = 0.5;
  end;



  if (isnumeric(qfld)) minqdt = 0; maxqdt = inf; else
    minqdt = stn.(qfld).date(1); maxqdt = stn.(qfld).date(end); end;

  if (isnumeric(sal)) minsaldt = 0; maxsaldt = inf; else
    minsaldt = stn.(sal).date(1); maxsaldt = stn.(sal).date(end); end;

  mindt = max([ ...
      stn.(wfld).date(1) ...
      stn.(afld).date(1) ...
      minqdt ...
      minsaldt ...
      stn.(pfld).date(1) ...
      stn.(tfld).date(1) ...
              ]);
  maxdt = min([ ...
      stn.(wfld).date(end) ...
      stn.(afld).date(end) ...
      maxqdt ...
      maxsaldt ...
      stn.(pfld).date(end) ...
      stn.(tfld).date(end) ...
              ]);

  % Construct perfect once-hourly series of timestamps
  dts = [mindt:onehour:maxdt]';

  % Interpolate all fields to hourly time series
  w = interp1_stn_fld(stn, wfld, dts);
  a = interp1_stn_fld(stn, afld, dts);
  p = interp1_stn_fld(stn, pfld, dts);
  t = interp1_stn_fld(stn, tfld, dts);

  if (~isnumeric(sal));  sal = interp1_stn_fld(stn,   sal, dts); end;

  if (~isempty(sfld));   srf = interp1_stn_fld(stn,  sfld, dts); end;

  if (~isempty(dsfld)); dsrf = interp1_stn_fld(stn, dsfld, dts); end;
  if (~isempty(dlfld)); dlrf = interp1_stn_fld(stn, dlfld, dts); end;
  if (~isempty(prfld))
    prcp = interp1_stn_fld(stn, prfld, dts);
  else
    prcp = repmat(0, size(dts));
  end;

  % Convert winds from kts. to m/s
  w = kts2mps(w);

  % Convert station air pressure and temperature to sea-level
  p = barom_to_surf(p,a,pz);

  % Use no unphysical air pressures!
  badix = (960>p | p>1040);
  p(badix) = interp1(dts(~badix),p(~badix),dts(badix),'spline');

  % Figure out something for relative humidity
  if ( isnumeric(qfld) )
    warning('Using CONSTANT Relative Humidity %g', qfld);
    q = repmat(qfld, size(w));

  else

    % Is the named field dewpoint temperature or relative humidity?
    if ( strfind(lower(qfld), 'dew') )
      warning('Calculating Relative Humidity from dew-temperature field %s', qfld);
      d = interp1_stn_fld(stn, qfld, dts);
      q = dewp_to_relhumid(a,d);

      relhumid = [PFX 'relhumid'];
      stn.(relhumid).date = dts;
      stn.(relhumid).data = q;

      spechumid = [PFX 'spechumid'];
      if ( ~isfield(stn,spechumid) )
        stn.(spechumid).date = dts;
        stn.(spechumid).data = relhumid_to_spechumid(a,q);
      end;

    elseif ( strfind(lower(qfld), 'hum') )
      q = interp1_stn_fld(stn, qfld, dts);

      spechumid = strrep(qfld,'rel','spec');
      if ( strcmp(spechumid,qfld) )
        warning('Could not calculate specific humidity: Field "%s" is not Relative Humidity??',qfld);
      elseif ( ~isfield(stn,spechumid) )
        stn.(spechumid).date = stn.(qfld).date;
        stn.(spechumid).data = relhumid_to_spechumid(a,stn.(qfld).data);
      end;

    else
      error('Does field "%s" contain dew-point or Relative Humidity??',qfld);

    end;

  end;

  tau = [];
  shf = [];
  lhf = [];
  rhf = [];


  result = [];

  if ( isempty(dsfld) || isempty(prfld) )
    % Perform Fairall et al. (TOGA-COARE 2.0) calculation
    if ( isempty(dsfld) )
      disp('HFBULKTC 2.0 (no cool-skin)');
      disp(length(w));
      % Algorithm sometimes returns complex numbers??
      result = real( hfbulktc(w,wz,a,az,q,az,p,t) );

    % Perform Fairall et al. (TOGA-COARE 2.6) calculation - NO RAIN
    elseif ( isempty(prfld) )
      disp('HFBULKTC 2.6 (cool-skin, no rain)');
      disp(length(w));
      % Algorithm sometimes returns complex numbers??
      result = real( hfbulktc(w,wz,a,az,q,az,p,t,sal,dlrf,dsrf,srf) );
    end;

    tau = result(:,4);
    shf = result(:,1);
    %lhf = result(:,2);
    % Be sure to include the 'Webb' correction!
    % (Comment from HFBULKTC: eqn. 22, Fairall et al. (1996), JGR, 101, p3751.)
    lhf = result(:,2) + result(:,3);

    % A=HFBULKTC(ur,zr,Ta,zt,rh,zq,Pa,Ts,sal,dlw,dsw,nsw) computes the following:
    % 01: Hs      = sensible heat flux INTO ocean [W/m^2]
    % 02: Hl      = latent heat flux INTO ocean [W/m^2]
    % 03: Hl_webb = Webb correction to latent heat flux INTO ocean [W/m^2]
    % 04: stress  = wind stress [N/m^2]
    % 05: U_star  = velocity friction scale [m/s]
    % 06: T_star  = temperature scale [deg C]
    % 07: Q_star  = humidity scale [kg/kg]
    % 08: L       = Monin-Obukhov length [m]
    % 09: zetu    = zr/L
    % 10: CD      = drag coefficient
    % 11: CT      = temperature transfer coefficient (Stanton number)
    % 12: CQ      = moisture transfer coefficient (Dalton number)
    % 13: RI      = bulk Richardson number
    % 14: Dter    = cool-skin temperature difference (optional output) [C]; 
    %               positive if surface is cooler than bulk (presently no 
    %               warm skin permitted by model)

    % Diagnostic and error propagation values
    diags = [PFX 'cordiags'];
    disp(['Updating ' diags]);
    stn.(diags).date = dts;

    stn.(diags).webb = result(:,3);             % (For 3.0a we get "Webb mean w (m/s)")
    stn.(diags).ustar = result(:,5);            % turbulent friction velocity (m/s)
    stn.(diags).tstar = result(:,6);            % temperature scaling parameter (K)
    stn.(diags).qstar = result(:,7);            % humidity scaling parameter (g/g)
    stn.(diags).monin = result(:,8);            % Monin_Obukhov stability length
    stn.(diags).zetu = result(:,9);             % (For 3.0a we get "velocity roughness length (m)")
    stn.(diags).Cd = result(:,10);              % velocity drag coefficient at zu, referenced to u
    stn.(diags).Ch = result(:,11);              % heat transfer coefficient at zt (Stanton number)
    stn.(diags).Ce = result(:,12);              % moisture transfer coefficient at zq (Dalton number)
    stn.(diags).Ri = result(:,13);              % bulk Richardson number
    if ( ~isempty(dsfld) )
      stn.(diags).dtcool = result(:,14);        % cool skin temperature depression (K)                   
    end;

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
    if ( isempty(pblzfld) || strcmpi(pblzfld,'default') )
      warning('Using constant Planetary Boundary Layer height 600m');
      pblz = repmat(600, size(w));
      %%%% Is 600m a good mean value for the FLORIDA KEYS??
      % Hsu (1979) shows 600m for the TX Gulf coast, while Kara et al (1998)
      % show an *upper bound* 500m for nocturnal conditions over Tallahassee:
      % they show profiles indicating 125m may be a better mean for there...
    elseif ( isnumeric(pblzfld) )
      warning('Using constant Planetary Boundary Layer height %g', pblzfld);
      pblz = repmat(pblzfld, size(w));
    else
      pblz = interp1_stn_fld(stn, pblzfld, dts);
    end;


    % Project ocean currents (if we have them) onto wind direction (ditto)
    if ( isempty(wdfld) || strcmpi(wdfld,'default') || ...
         isempty(oufld) || strcmpi(oufld,'default') || ...
         isempty(ovfld) || strcmpi(ovfld,'default') )
      % Assume surface current projects onto 2% of wind velocity [e.g., Ardhuin et al. 2009]
      warning('Estimating projected ocean currents from wind speed');
      ou = 0.020 .* w;
    else
      ou = interp1_stn_fld(stn, oufld, dts);
      ov = interp1_stn_fld(stn, ovfld, dts);
      [wix,oix] = intersect_dates(stn.(wdfld).date,dts);
      wd = stn.(wdfld).data(wix);
      ouraw = -sind(wd).*ou(oix) - cosd(wd).*ov(oix);
      ou = interp1(dts(oix(isfinite(ouraw))),ouraw(isfinite(ouraw)),dts);
    end;

    if ( isempty(wpfld) )

      % Perform Fairall et al. (TOGA-COARE 2.6) calculation
      disp('COR26 (cool-skin, rain)');
      disp(length(w));
      nby10 = floor(length(w)/10);
      for ix = 1:length(w)
        res = cor30a([w(ix),ou(ix),t(ix),a(ix),ss(ix),sa(ix),dsrf(ix),dlrf(ix),prcp(ix),pblz(ix),p(ix),wz,az,az,stn.lat,1,0,0,0]);
        % Preallocate for speed
        if ( isempty(result) )
          result = repmat(nan,[length(w) length(res)]);
        end;
        % Algorithm sometimes returns complex numbers??
        result(ix,:) = real( res );
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

      % Use no unphysical wave periods!
      badix = (1>wvper | wvper>20 | ~isfinite(wvper));
      %DEBUG:      disp(length(find(badix)));
      wvper(badix) = interp1(dts(~badix),wvper(~badix),dts(badix),'spline');


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
        %result = real( cor30a_warm(dts,w,ou,t,a,ss,sa,dsrf,dlrf,prcp,pblz,p,stn.lon,stn.lat,wz,az,az,stz,1,1,1,wvper,wvhgt) );
        % Lew.Gramer@noaa.gov, 2013 Feb 06:
        % JWave==2: Use Taylor & Yelland z0
        result = real( cor30a_warm(dts,w,ou,t,a,ss,sa,dsrf,dlrf,prcp,pblz,p,stn.lon,stn.lat,wz,az,az,stz,1,1,2,wvper,wvhgt,max_wl) );
      end; %if ~doWarm else

    end; %if isempty(wpfld) else

    tau = result(:,3);                          % stress (nt/m^2)
    shf = -result(:,1);                         % sensible heat flux (w/m^2)
    lhf = -result(:,2);                         % latent heat flux (w/m^2)
    rhf = -result(:,14);                        % rain heat flux(w/m^2)

    % Diagnostic and error propagation values
    diags = [PFX 'cordiags'];
    disp(['Updating ' diags]);
    stn.(diags).date = dts;
    stn.(diags).z0u = result(:,4);              % velocity roughness length (m)                          
    stn.(diags).z0t = result(:,5);              % temperature roughness length (m)                       
    stn.(diags).z0q = result(:,6);              % moisture roughness length (m)                          
    stn.(diags).monin = result(:,7);            % Monin_Obukhov stability length                         
    stn.(diags).ustar = result(:,8);            % turbulent friction velocity (m/s), including gustiness 
    stn.(diags).tstar = result(:,9);            % temperature scaling parameter (K)                      
    stn.(diags).qstar = result(:,10);           % humidity scaling parameter (g/g)                       
    stn.(diags).dtcool = result(:,11);          % cool skin temperature depression (K)                   
    stn.(diags).dqcool = result(:,12);          % cool skin humidity depression (g/g)                    
    stn.(diags).dxcool = result(:,13);          % cool skin thickness (m)                                

    stn.(diags).wbar = result(:,15);            % Webb mean w (m/s)
    stn.(diags).Cd = result(:,16);              % velocity drag coefficient at zu, referenced to u       
    stn.(diags).Ch = result(:,17);              % heat transfer coefficient at zt                        
    stn.(diags).Ce = result(:,18);              % moisture transfer coefficient at zq                    
    stn.(diags).Cdn = result(:,19);             % neutral 10-m velocity drag coefficient, including gustiness    
    stn.(diags).Chn = result(:,20);             % neutral 10-m heat transfer coefficient, including gustiness    
    stn.(diags).Cen = result(:,21);             % neutral 10-m humidity transfer coefficient, including gustiness
    stn.(diags).Ug = result(:,22);		% wind "gustiness" parameter
    if ( size(result,2) >= 25 )
      stn.(diags).dtwarm = result(:,23);        % warm layer deltaT
      stn.(diags).dxwarm = result(:,24);        % warm layer thickness (Price-Weller-Pinkel)
      stn.(diags).twarm = result(:,25);         % warm layer corrected Tsea
    end;

  end; %if isempty(dsfld) || isempty(prfld) else


  % Update our station data structure

  wind_stress = [PFX 'wind_stress'];
  if ( isfield(stn,wind_stress) ); stn = rmfield(stn,wind_stress); end;
  stn.(wind_stress).date = dts;
  stn.(wind_stress).data = tau;
  if ( ~isempty(wdfld) && ~strcmpi(wdfld,'default') && isfield(stn,wdfld) )
    [wix,oix] = intersect_dates(stn.(wdfld).date,dts);
    wd = stn.(wdfld).data(wix);
    stn.([wind_stress '_u']).date = dts(oix);
    stn.([wind_stress '_v']).date = dts(oix);
    [stn.([wind_stress '_u']).data,stn.([wind_stress '_v']).data] = spddir_to_uv(tau(oix),wd);
  end;

  sensible_flux = [PFX 'sensible_flux'];
  if ( isfield(stn,sensible_flux) ); stn = rmfield(stn,sensible_flux); end;
  stn.(sensible_flux).date = dts;
  stn.(sensible_flux).data = shf;

  latent_flux = [PFX 'latent_flux'];
  if ( isfield(stn,latent_flux) ); stn = rmfield(stn,latent_flux); end;
  stn.(latent_flux).date = dts;
  stn.(latent_flux).data = lhf;

  net_flux = [PFX 'net_flux'];
  if ( isfield(stn,net_flux) ); stn = rmfield(stn,net_flux); end;
  stn.(net_flux).date = stn.(sensible_flux).date;
  stn.(net_flux).data = stn.(sensible_flux).data + stn.(latent_flux).data;


  % Add in rain heat flux, if we were given rain data
  if (~isempty(rhf))
    rain_flux = [PFX 'rain_flux'];
    if ( isfield(stn,rain_flux) ); stn = rmfield(stn,rain_flux); end;
    stn.(rain_flux).date = dts;
    stn.(rain_flux).data = rhf;

    stn.(net_flux).data = stn.(net_flux).data + stn.(rain_flux).data;
  end;


  %%%%
  %% REMOVE any timestamps corresponding to big gaps in original data

  %wfld,afld,qfld,pfld,tfld
  stn = filter_gaps(stn,wfld,wind_stress);
  stn = filter_gaps(stn,afld,wind_stress);
  if ~isnumeric(qfld); stn = filter_gaps(stn,qfld,wind_stress); end;
  stn = filter_gaps(stn,pfld,wind_stress);
  stn = filter_gaps(stn,tfld,wind_stress);
  %wfld,afld,qfld,pfld,tfld
  stn = filter_gaps(stn,wfld,sensible_flux);
  stn = filter_gaps(stn,afld,sensible_flux);
  if ~isnumeric(qfld); stn = filter_gaps(stn,qfld,sensible_flux); end;
  stn = filter_gaps(stn,pfld,sensible_flux);
  stn = filter_gaps(stn,tfld,sensible_flux);
  %wfld,afld,qfld,pfld,tfld
  stn = filter_gaps(stn,wfld,latent_flux);
  stn = filter_gaps(stn,afld,latent_flux);
  if ~isnumeric(qfld); stn = filter_gaps(stn,qfld,latent_flux); end;
  stn = filter_gaps(stn,pfld,latent_flux);
  stn = filter_gaps(stn,tfld,latent_flux);
  if (~isempty(rhf))
    %wfld,afld,qfld,pfld,tfld
    stn = filter_gaps(stn,wfld,rain_flux);
    stn = filter_gaps(stn,afld,rain_flux);
    if ~isnumeric(qfld); stn = filter_gaps(stn,qfld,rain_flux); end;
    stn = filter_gaps(stn,pfld,rain_flux);
    stn = filter_gaps(stn,tfld,rain_flux);
  end;
  %wfld,afld,qfld,pfld,tfld
  stn = filter_gaps(stn,wfld,net_flux);
  stn = filter_gaps(stn,afld,net_flux);
  if ~isnumeric(qfld); stn = filter_gaps(stn,qfld,net_flux); end;
  stn = filter_gaps(stn,pfld,net_flux);
  stn = filter_gaps(stn,tfld,net_flux);


  % If specified, add in NET insolation and longwave radiation terms also...
  if ( isempty(sfld) )
    warning('Ecoforecasts:Heat:NoRadiative',...
            '"Net" heat flux DOES NOT INCLUDE RADIATIVE FLUXES!');
  else
    [ix1,ix2] = intersect_dates(stn.(net_flux).date(:), stn.(sfld).date(:));
    if ( isempty(ix1) )
      warning('No matching dates between STN.%s and STN.%s!', net_flux, sfld);
      stn.(net_flux).date = [];
      stn.(net_flux).data = [];
    else
      stn.(net_flux).data(ix1(:)) = stn.(net_flux).data(ix1(:)) + stn.(sfld).data(ix2(:));
      stn = filter_gaps(stn, sfld, net_flux);

      if ( isempty(lfld) )
        warning('"Net" heat flux DOES NOT INCLUDE NET LONG-WAVE FLUX!');
      else
        [ix1,ix2] = intersect_dates(stn.(net_flux).date(:), stn.(lfld).date(:));
        if ( isempty(ix1) )
          warning('No matching dates between STN.%s and STN.%s!', net_flux, lfld);
          stn.(net_flux).date = [];
          stn.(net_flux).data = [];
        else
          stn.(net_flux).data(ix1(:)) = stn.(net_flux).data(ix1(:)) + stn.(lfld).data(ix2(:));
          stn = filter_gaps(stn, lfld, net_flux);
        end;
      end;
    end;

  end;


  % Scale using form from heat budget equation HF./(rho*Cp*h), and
  % change from MKS units of [K/s], to more user-friendly [K/hr]

  Q0_factor = Q0factor(t,[],stz);

  flux_term = [PFX 'flux_term'];
  stn.(flux_term).date = stn.(net_flux).date;
  stn.(flux_term).data = (stn.(net_flux).data ./ Q0_factor) .* (60*60);

  latent_flux_term = [PFX 'latent_flux_term'];
  stn.(latent_flux_term).date = stn.(latent_flux).date;
  stn.(latent_flux_term).data = (stn.(latent_flux).data ./ Q0_factor) .* (60*60);

  sensible_flux_term = [PFX 'sensible_flux_term'];
  stn.(sensible_flux_term).date = stn.(sensible_flux).date;
  stn.(sensible_flux_term).data = (stn.(sensible_flux).data ./ Q0_factor) .* (60*60);

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


  if ( doPlot )
    X = { dts, dts, dts, dts, dts, dts, dts, dts };
    Y = { w, a, q, p, t, stn.(wind_stress).data, ...
          stn.(sensible_flux).data, stn.(latent_flux).data };
    multiplot(X, Y, 'Title',['Heat Flux: ' stn.station_name], ...
              'YLabel',{'Wind','Air_T','RH','Baro','Sea_T','\tau','H_s','H_l'});
    datetick3;
    maxigraph;
  end;

  if ( exist('added_fairall_path','var') )
    warning('off','MATLAB:rmpath:DirNotFound');
    rmpath(FAIRALLHOME);
    warning('on','MATLAB:rmpath:DirNotFound');
    rehash
  end;

return;


%%%%%%%%%%
%%%%%%%%%% PRIVATE FUNCTIONS
%%%%%%%%%%

function d = interp1_stn_fld(stn, fld, dts)
  d = interp1(stn.(fld).date(isfinite(stn.(fld).data)), ...
              stn.(fld).data(isfinite(stn.(fld).data)), ...
              dts);
return;
