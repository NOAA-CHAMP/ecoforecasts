function stn = station_bulk_albedo(stn,albfld,Wfld,cfld,chl_a,whfld)
%function stn = station_bulk_albedo(stn,albfld,Wfld,cfld,chl_a,whfld)
%
% Calculate bulk shortwave ocean albedo (SOA) time series STN.(ALBFLD) based
% on wind speed STN.(WFLD) and other environmental conditions. Simple model
% 6-15% based on -SQRT(wind speed), or empirical model of Jin et al. (2011)
% based on STN.(WFLD) and cloud fraction STN.(CFLD). Field names CHL_A and
% WHFLD are reserved for later improvements based on water type and waves.
%
% Last Saved Time-stamp: <Mon 2011-12-19 15:44:40  lew.gramer>

  if ( ~exist('cfld','var') || isempty(cfld) )
    cfld = '';
  end;
  if ( ~exist('chl_a','var') || isempty(chl_a) )
    chl_a = '';
  end;
  if ( ~exist('whfld','var') || isempty(whfld) )
    whfld = '';
  end;

  if ( isempty(cfld) )
    haveClouds = false;
  elseif ( ~isfield(stn,cfld) || ~is_valid_ts(stn.(cfld)) )
    haveClouds = false;
    warning('Found no valid cloud fraction STN.(%s)',cfld);
  else
    haveClouds = true;
  end;

  if ( ~haveClouds )
    dts = stn.(Wfld).date;
    W = kts2mps(stn.(Wfld).data);
    C = 0.40; % Mean cloud cover in the Florida Keys ~40%
  else
    [wix,cix] = intersect_dates(stn.(Wfld).date,stn.(cfld).date);
    dts = stn.(Wfld).date(wix);
    W = kts2mps(stn.(Wfld).data(wix));
    C = stn.(cfld).data(cix);
  end;

  %%%%%DEBUG???
  % simpleModel = ~haveClouds;
  simpleModel = false;

  if ( simpleModel )

    minA = 0.06;
    maxA = 0.15;
    facA = (maxA-minA);

    maxW = 30;
    W(W>maxW) = maxW;

    alb = minA + ( facA.*(sqrt(maxW)-sqrt(W))./sqrt(maxW) );

  else

    % Empirical relationship of Jin et al. (2011)

    % Angle of incidence - approximated by solar zenith angle
    if ( ~isfield(stn,'lon') || ~isfield(stn,'lat') )
      error('Solar zenith angle requires longitude STN.lon and latitude STN.lat');
    end;
    % Center hourly averages on the half-hour
    [yds,yrs] = get_yearday(dts-(0.5/24));
    [sun_alt,ig] = soradna1(yds,yrs,-stn.lon,stn.lat);
    theta = 90-sun_alt;

    mu = cosd(theta);

    % Baseline relative refractive index of water and air (see Jin)
    n0 = 1.34;

    % Mean fraction PAR [W] of insolation for overcast skies (Jacovides et
    % al. 2003), *PLUS* mean fraction UV(A+B) [W] of insolation, for annual
    % range of Keys total cloud cover (40%) from ERAI (Leal et al. 2011).
    % NOTE this could be improved somewhat by estimate based on conditions?
    Ppen = 0.549;

    % Refractive indices of air and seawater
    na = 1.005; % Reasonable guess
    nw_PAR = 1.339; % Leyendekkers (1977): T~25, S~35, p~1, lambda~550nm
    nw_NIR = 1.256; % Filipiak (2008): T~25, S~35, lambda~935nm, k~1700/cm
    % Full reference: Filipiak, Mark (2008): Refractive indices (500-3500
    % cm-) and emissivity (600-3350 cm-1) of pure water and seawater.
    % http://datashare.is.ed.ac.uk/handle/10283/17

    nws = [nw_PAR,nw_NIR];
    %DEBUG:    nws = nw_PAR;
    %DEBUG:    nws = nw_NIR;
    %DEBUG:    nws = 1.1315; %Extreme NIR n


    albSdir = repmat(nan,[numel(nws),numel(dts)]);
    albSdif = repmat(nan,[numel(nws),numel(dts)]);

    for nwix = 1:numel(nws)
      nw = nws(nwix);

      % Relative refractive index of water and air
      n = nw ./ na;
      % Fresnel reflectance
      rf = fresnel_reflectance(theta,nw,na);
      rf0 = fresnel_reflectance(theta,n0,1.00);

      % Cox & Munk (1954) surface roughness
      sigma = sqrt(0.003 + (0.00512.*W));

      % Direct reflectance at sea surface
      % (Jin et al. (2011) empirical relationship)
      p0 = 0.0152;  p1 = -1.7873; p2 = 6.8972;   p3 = -8.5778;
      p4 = 4.071;   p5 = -7.6446; p6 = 0.1643;   p7 = -7.8409;
      p8 = -3.5639; p9 = -2.3588; p10 = 10.0538;
      f = (p0 + (p1.*mu) + (p2.*(mu.^2)) + (p3.*(mu.^3)) + (p4.*sigma) + (p5.*mu.*sigma)) ...
          .* exp(p6 + (p7.*mu) + (p8.*(mu.^2)) + (p9.*sigma) + (p10.*mu.*sigma));

      albSdir(nwix,:) = rf - ( (rf./rf0) .* f );
      albSdir(nwix,albSdir(nwix,:)<0) = 0;
      albSdir(nwix,albSdir(nwix,:)>1) = 1;

      % Diffuse reflectance at sea surface
      % (Jin et al. (2011) empirical relationship)
      albSdif(nwix,:) = -0.1482 - (0.012.*sigma) + (0.1608.*n) - (0.0244.*n.*sigma);
      albSdif(nwix,albSdif(nwix,:)<0) = 0;
      albSdif(nwix,albSdif(nwix,:)>1) = 1;

      % Volume scattering due to water column (direct+diffuse)
      %albW = 0.006; % Suggested by Jin et al. (2011) for Case 1 waters
      %albW = 0.008; % Estimated from Jin et al. (2011) Fig. 7 for chl_a ~ 0.4
      albW = 0.012; % From Jin but for an "effective chl_a" (water column+shallow bottom) ~ 1.0
    end;

    % Fraction of direct vs. total (direct+diffuse) insolation
    % Clear sky
    fdir_clear_sky = 0.85; % From Jin et al. (2011) Fig. 8, ignoring zenith angle

    % Varying cloud conditions (e.g., from airport reports or reanalysis)
    fdir = fdir_clear_sky .* (1 - C);

    % Fraction of diffuse vs. total (direct+diffuse) insolation
    fdif = 1 - fdir;

    % Weighted average of direct FDIR vs. diffuse FDIF, and PAR+UV
    % ("penetrative") and NIR insolation fractions, resp.
    albSdirmn = (Ppen.*albSdir(1,:)) + ((1-Ppen).*albSdir(2,:));
    albSdifmn = (Ppen.*albSdif(1,:)) + ((1-Ppen).*albSdif(2,:));
    alb = (fdir.*albSdirmn(:)) + (fdif.*albSdifmn(:)) + albW;

    % Foam correction (op cit, "proposed by Koepke" (1984))
    fwc = 2.95e-6 .* (W.^3.52);
    albwc = 0.55;
    alb = (fwc.*albwc) + ((1-fwc).*alb);

    % Calculations for Florida Keys were also consistent with Chang & Dickey (2004)
    alb(alb<0) = 0;
    alb(alb>1) = 1;

  end;

  stn.(albfld).date = dts(:);
  stn.(albfld).data = alb(:);

return;


%%%%%%%%%%%%%%%%%%%%%
% INTERNAL FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%

function [rf,theta_t] = fresnel_reflectance(theta,nw,na)
%function [rf,theta_t] = fresnel_reflectance(theta,nw,na)

  % Angle of refraction
  theta_t = asind( sind(theta).*na ./ nw );
  % Reflectance
  rf = ( (na.*cosd(theta) - nw.*cosd(theta_t)) ./ (na.*cosd(theta) + nw.*cosd(theta_t)) ).^2;

return;
