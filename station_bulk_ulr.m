function stn = station_bulk_ulr(stn,sfld,ulrf,dlrf,lrf,epsilon)
%function stn = station_bulk_ulr(stn,sfld,ulrf,dlrf,lrf,epsilon)
%
% Calculate upward longwave radiative flux STN.(ULRF) with bulk formula
% using sea temperature STN.(sfld). Also adjusts for re-radiated longwave
% flux if optional fieldname for downward longwave flux STN.(DLRF) is given,
% and returns net longwave flux STN.(LRF) if that fieldname is given.
%
% NOTE: Unlike other Ecoforecasts air-sea fluxes, ULRF is *POSITIVE UPWARD*
% to match the signs in TOMS GSIP, NCEP, ERA-I, OAFlux, and other products.
%
% ALSO NOTE: STN.(SFLD) should be a *surface* sea temperature: this function
% makes no attempt to adjust STN.(SFLD) from bulk to cool-skin temperature.
%
% Last Saved Time-stamp: <Fri 2012-08-03 12:54:25  lew.gramer>

  if ( isfield(stn,ulrf) ); stn = rmfield(stn,ulrf); end;

  sigma = 5.67e-8;  % Stefan-Boltzman constant

  if ( ~exist('epsilon','var') || isempty(epsilon) )
    % Was 0.96 - why?
    epsilon = 0.97;   % Ocean surface emissivity (see Anderson, 1952)
                      % Xue et al (1998) say: "0.97 after Anderson (1952) and Reed (1976)"
                      % Kraus and Businger suggest 0.98
  end;

  if ( ~exist('dlrf','var') || isempty(dlrf) )
    Ts = stn.(sfld);
    Ts.data = Ts.data + 273.14; % convert to [K]

    % Simple black-body radiation model of ULRF
    stn.(ulrf).date = Ts.date;
    stn.(ulrf).data = (epsilon.*sigma.*(Ts.data.^4));

  else

    if ( ~isfield(stn,dlrf) || ~is_valid_ts(stn.(dlrf)) )
      error('No valid DLR time series STN.%s found',dlrf);
    end;

    [Ts,dlr] = intersect_tses(stn.(sfld),stn.(dlrf));
    Ts.data = Ts.data + 273.14; % convert to [K]

    % Black-body radiation model of ULRF - with tweak
    % Formula below suggested by reading of Weller et al (2008), and Fung et al (1984)
    stn.(ulrf).date = Ts.date;
    stn.(ulrf).data = (epsilon.*sigma.*(Ts.data.^4)) - ((1-epsilon).*dlr.data);

    if ( exist('lrf','var') && ~isempty(lrf) )
      stn.(lrf) = ts_op(stn.(dlrf),stn.(ulrf),'-');
    end;

  end;

return;
