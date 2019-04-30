function stn = station_ndbc_hfbulk(stn,srf,lrf,rh,wz,az)
%function stn = station_ndbc_hfbulk(stn,srf,lrf,[rh[,wz[,az]]])
%
% Calculate turbulent heat fluxes from station NDBC (Quality Controlled)
% data, using simple bulk formulae with bulk coefficients per Smith (1988).
% RH is relative humidity fieldname; if no in situ RH is available, specify
% another fieldname (e.g., 'erai_relhumid'), or a constant scalar. Optional
% WZ, AZ are sensor heights [m] for air temp. and wind speed (DEFAULT: 30m).
%
% Last Saved Time-stamp: <Wed 2012-03-28 13:47:38  Lew.Gramer>


  if ( ~exist('rh','var') || isempty(rh) )
    rh = 'ndbc_relhumid';
  end;

  % "Representative" values from Smith (1988)
  Cdr = 1.1e-3;
  Cth = 1.0e-3;
  Cqs = 1.3e-3;

  Q0f = Q0factor(stn.ndbc_sea_t.data,[],2);

  if ( ~exist('wz','var') || isempty(wz) )
    try
      [wz,ig,ig,ig] = station_instrument_heights(stn.station_name);
    catch
      warning('No in situ wind sensor height for STN: assuming 30m');
      wz = 30;
    end;
  end;
  if ( ~exist('az','var') || isempty(az) )
    try
      [ig,az,ig,ig] = station_instrument_heights(stn.station_name);
    catch
      warning('No in situ air temp. sensor height for STN: assuming 30m');
      az = 30;
    end;
  end;


  % Theta_a_AZ - Theta_t
  [ix1,ix2] = intersect_dates(stn.ndbc_air_t.date, stn.ndbc_sea_t.date);
  stn.ndbc_air_sea_t.date = stn.ndbc_air_t.date(ix1);
  stn.ndbc_air_sea_t.data = stn.ndbc_air_t.data(ix1) - ...
      stn.ndbc_sea_t.data(ix2);

  % Needed for U^2 + V^2
  if ( ~isfield(stn,'ndbc_wind1_u') )
    stn = vector_func(stn, 'ndbc_wind1', 'u');
    % Convert winds from kts. to m/s
    stn.ndbc_wind1_u.data = kts2mps(stn.ndbc_wind1_u.data);
  end;
  if ( ~isfield(stn,'ndbc_wind1_v') )
    stn = vector_func(stn, 'ndbc_wind1', 'v');
    % Convert winds from kts. to m/s
    stn.ndbc_wind1_v.data = kts2mps(stn.ndbc_wind1_v.data);
  end;

  % Convert data to standard 10m height assuming a log-layer
  % Uses simple algorithms from COR30a by Fairall et al. (2003)
  if ( ~isfield(stn,'ndbc_air_t10m') || ~isfield(stn,'ndbc_air_sea_t10m') )
    % Theta_a_10m - Theta_t
    stn.ndbc_air_sea_t10m.date = stn.ndbc_air_sea_t.date;
    stn.ndbc_air_sea_t10m.data = stn.ndbc_air_sea_t.data - (0.0098*az);

    [dix,six] = intersect_dates(stn.ndbc_air_sea_t10m.date, stn.ndbc_sea_t.date);
    stn.ndbc_air_t10m.date = stn.ndbc_air_sea_t10m.date(dix);
    stn.ndbc_air_t10m.data = stn.ndbc_air_sea_t10m.data(dix) + ...
        stn.ndbc_sea_t.data(six);
  end;
  if ( ~isfield(stn,'ndbc_wind1_u10m') || ~isfield(stn,'ndbc_wind1_v10m') )
    stn.ndbc_wind1_u10m.date = stn.ndbc_wind1_u.date;
    stn.ndbc_wind1_u10m.data = stn.ndbc_wind1_u.data .* log(10/1e-4)/log(wz/1e-4);
    stn.ndbc_wind1_v10m.date = stn.ndbc_wind1_v.date;
    stn.ndbc_wind1_v10m.data = stn.ndbc_wind1_v.data .* log(10/1e-4)/log(wz/1e-4);
  end;


  % Representative constants from AIR_SEA toolbox and SMKF1.ndbc_sea_t
  % Density of moist air
  rhoa = 1.2;
  % Specific heat capacity of moist air at constant pressure [J/kg/K]
  Cpa = 1.01e3;
  % Latent heat of evaporation [J/kg]
  %Le = 2.4e6;
  Le = vapor(nanmean(stn.ndbc_air_t10m.data));


  % SENSIBLE HEAT FLUX:
  % rhoa x Cpa x Cdr^.5 x Cth ^.5 x (U^2 + V^2)^.5 x (th_s - th_a)
  [ix1,ix2] = intersect_dates(stn.ndbc_air_sea_t10m.date, stn.ndbc_wind1_u10m.date);
  stn.ndbc_hfbulk_sensible_heat_flux.date = stn.ndbc_air_sea_t10m.date(ix1);
  stn.ndbc_hfbulk_sensible_heat_flux.data = sqrt(Cdr) .* sqrt(Cth) .* ...
      rhoa .* Cpa .* ...
      sqrt((stn.ndbc_wind1_u10m.data(ix2).^2) + (stn.ndbc_wind1_v10m.data(ix2).^2)) .*...
      stn.ndbc_air_sea_t10m.data(ix1);

  stn.ndbc_hfbulk_sensible_flux_term.date = stn.ndbc_hfbulk_sensible_heat_flux.date;
  stn.ndbc_hfbulk_sensible_flux_term.data = ...
      (stn.ndbc_hfbulk_sensible_heat_flux.data ./ Q0f) .* (60*60);


  % qs_s
  stn.ndbc_sea_spechumid.date = stn.ndbc_sea_t.date;
  % 0.98 factor from Stommel - accounts for salinity
  stn.ndbc_sea_spechumid.data = 0.98 .* relhumid_to_spechumid(stn.ndbc_sea_t.data,100);

  % qs_a - qs_s
  if ( ~isfield(stn,rh) )
    warning('Cannot calculate NDBC HFBulk fluxes accurately: no humidity data! Using constant RH=80');
    stn.const_relhumid.date = stn.ndbc_air_t10m.date;
    stn.const_relhumid.data = repmat(80, size(stn.ndbc_air_t10m.data));
    stn = station_relhumid_to_spechumid(stn,'ndbc_air_t10m','const_relhumid','ndbc_spechumid10m');
  else
    stn = station_relhumid_to_spechumid(stn,'ndbc_air_t10m',rh,'ndbc_spechumid10m');
  end;

  [ix1,ix2] = intersect_dates(stn.ndbc_spechumid10m.date, stn.ndbc_sea_spechumid.date);
  stn.ndbc_air_sea_spechumid10m.date = stn.ndbc_spechumid10m.date(ix1);
  stn.ndbc_air_sea_spechumid10m.data = (stn.ndbc_spechumid10m.data(ix1) - ...
                                        stn.ndbc_sea_spechumid.data(ix2));

  % LATENT HEAT FLUX:
  % rhoa x Le x Cdr^.5 x Cth ^.5 x (U^2 + V^2)^.5 x (qs_a - qs_s)
  [ix1,ix2] = intersect_dates(stn.ndbc_air_sea_spechumid10m.date, stn.ndbc_wind1_u10m.date);
  stn.ndbc_hfbulk_latent_heat_flux.date = stn.ndbc_air_sea_spechumid10m.date(ix1);
  stn.ndbc_hfbulk_latent_heat_flux.data = sqrt(Cdr) .* sqrt(Cqs) .* ...
      rhoa .* Le .* ...
      sqrt((stn.ndbc_wind1_u10m.data(ix2).^2) + (stn.ndbc_wind1_v10m.data(ix2).^2)) .*...
      stn.ndbc_air_sea_spechumid10m.data(ix1);

  stn.ndbc_hfbulk_latent_flux_term.date = stn.ndbc_hfbulk_latent_heat_flux.date;
  stn.ndbc_hfbulk_latent_flux_term.data = ...
      (stn.ndbc_hfbulk_latent_heat_flux.data ./ Q0f) .* (60*60);


  if ( exist('srf','var') && exist('lrf','var') )
    [ix1,ix2] = intersect_dates(stn.ndbc_hfbulk_latent_heat_flux.date, stn.ndbc_hfbulk_sensible_heat_flux.date);
    dts = stn.ndbc_hfbulk_latent_heat_flux.date(ix1);
    dat = stn.ndbc_hfbulk_latent_heat_flux.data(ix1) + stn.ndbc_hfbulk_sensible_heat_flux.data(ix2);

    [ix1,ix2] = intersect_dates(dts, stn.(srf).date);
    dts = dts(ix1);
    dat = dat(ix1) + stn.(srf).data(ix2);
    [ix1,ix2] = intersect_dates(dts, stn.(lrf).date);
    dts = dts(ix1);
    dat = dat(ix1) + stn.(lrf).data(ix2);
    stn.ndbc_hfbulk_net_heat_flux.date = dts;
    stn.ndbc_hfbulk_net_heat_flux.data = dat;

    stn.ndbc_hfbulk_heat_flux_term.date = stn.ndbc_hfbulk_net_heat_flux.date;
    stn.ndbc_hfbulk_heat_flux_term.data = ...
        (stn.ndbc_hfbulk_net_heat_flux.data ./ Q0f) .* (60*60);
  end;


return;
