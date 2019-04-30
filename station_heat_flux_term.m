function [stn,Q0_factor] = station_heat_flux_term(stn,nffld,htfld,tfld,sfld,dfld)
%function [stn,Q0_factor] = station_heat_flux_term(stn,nffld,htfld,tfld,sfld,dfld)
%
% Calculate divisor Q0_FACTOR for STN.(NFFLD)in heat budget equation, and
% multiply result by 60*60 to calculate a rate of temperature change from
% surface fluxes in [K/hr]. Store result in new field STN.(HTFLD). STN.(TFLD)
% is sea temperature time series; args SFLD (salinity [psu]), and DFLD (site
% depth [m]) may be field names, scalar, or [] ("use a reasonable default").
%
% Last Saved Time-stamp: <Fri 2010-06-25 18:26:13 Eastern Daylight Time gramer>

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
    [nfix,stix] = intersect_dates(stn.(nffld).date,stn.(tfld).date);
    s = repmat(s,size(stn.(nffld).date(nfix)));
    d = repmat(d,size(stn.(nffld).date(nfix)));
  elseif ( ~isempty(d) )
    [nfix,stix,ssix] = intersect_all_dates([],stn.(nffld).date,stn.(tfld).date,stn.(sfld).date);
    s = stn.(sfld).data(ssix);
    d = repmat(d,size(stn.(nffld).date(nfix)));
  elseif ( ~isempty(s) )
    [nfix,stix,sdix] = intersect_all_dates([],stn.(nffld).date,stn.(tfld).date,stn.(dfld).date);
    s = repmat(s,size(stn.(nffld).date(nfix)));
    d = stn.(dfld).data(sdix);
  else
    [nfix,stix,ssix,sdix] = intersect_all_dates([],stn.(nffld).date,stn.(tfld).date,stn.(sfld).date,stn.(dfld).date);
    s = stn.(sfld).data(ssix);
    d = stn.(dfld).data(sdix);
  end;

  dts = stn.(nffld).date(nfix);
  nf = stn.(nffld).data(nfix);
  t = stn.(tfld).data(stix);

  % Cp and rho from [Fofonoff and Millard, 1983]
  rho = sw_dens0( s, t );
  Cp = sw_cp( s, t, d );
  h = d;

  Q0_factor = (rho.*Cp.*h);

  if ( isfield(stn,htfld) )
    stn = rmfield(stn,htfld);
  end;

  stn.(htfld).date = dts;
  stn.(htfld).data = (nf ./ Q0_factor) .* (60*60);

return;
