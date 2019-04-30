function [stn,Q0_factor] = station_heat_flux_term_inverse(stn,nffld,htfld,tfld,sfld,dfld)
%function [stn,Q0_factor] = station_heat_flux_term_inverse(stn,nffld,htfld,tfld,sfld,dfld)
%
% This function does the *inverse* of function STATION_HEAT_FLUX_TERM (v.).
% Calculate divisor Q0_FACTOR for heat budget equations (rho_w*Cp*h), and
% divide result by 60*60. Multiply heating rate term, STN.(HTFLD) (rate of
% temperature change from surface fluxes in [K/hr]) by this value to estimate
% STN.(NFFLD), a heat flux from the temperature change in [W/m2]. STN.(TFLD)
% is sea temperature time series; args SFLD (salinity [psu]), and DFLD (site
% depth [m]) may be field names, scalar, or [] ("use a reasonable default").
%
% Last Saved Time-stamp: <Sun 2010-10-03 17:06:17 Eastern Daylight Time gramer>

  % DEFAULTS
  representative_mean_salinity = 36;
  try
    [wz,az,pz,stz] = station_instrument_heights(stn.station_name);
    representative_depth_meters = stz;
  catch
    representative_depth_meters = 2;
  end;

  s = [];
  if ( ~exist('sfld','var') || isempty(sfld) )
    s = representative_mean_salinity;
  elseif ( isscalar(sfld) && isnumeric(sfld) )
    s = sfld;
  end;

  d = [];
  if ( ~exist('dfld','var') || isempty(dfld) )
    d = representative_depth_meters;
  elseif ( isscalar(dfld) && isnumeric(dfld) )
    d = dfld;
  end;

  if ( ~isempty(s) && ~isempty(d) )
    [htix,stix] = intersect_dates(stn.(htfld).date,stn.(tfld).date);
    s = repmat(s,size(stn.(htfld).date(htix)));
    d = repmat(d,size(stn.(htfld).date(htix)));
  elseif ( ~isempty(d) )
    [htix,stix,ssix] = intersect_all_dates([],stn.(htfld).date,stn.(tfld).date,stn.(sfld).date);
    s = stn.(sfld).data(ssix);
    d = repmat(d,size(stn.(htfld).date(htix)));
  elseif ( ~isempty(s) )
    [htix,stix,sdix] = intersect_all_dates([],stn.(htfld).date,stn.(tfld).date,stn.(dfld).date);
    s = repmat(s,size(stn.(htfld).date(htix)));
    d = stn.(dfld).data(sdix);
  else
    [htix,stix,ssix,sdix] = intersect_all_dates([],stn.(htfld).date,stn.(tfld).date,stn.(sfld).date,stn.(dfld).date);
    s = stn.(sfld).data(ssix);
    d = stn.(dfld).data(sdix);
  end;

  dts = stn.(htfld).date(htix);
  ht = stn.(htfld).data(htix);
  t = stn.(tfld).data(stix);

  % Cp and rho from [Fofonoff and Millard, 1983]
  rho = sw_dens0( s, t );
  Cp = sw_cp( s, t, d );
  h = d;

  Q0_factor = (rho.*Cp.*h);

  if ( isfield(stn,nffld) )
    stn = rmfield(stn,nffld);
  end;

  stn.(nffld).date = dts;
  stn.(nffld).data = (ht .* Q0_factor) ./ (60*60);

return;
