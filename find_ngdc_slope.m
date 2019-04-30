function [bet,ang,iso,bat] = find_ngdc_slope(bat,lon,lat,N,method,extrapval)
%function [bet,ang,iso,bat] = find_ngdc_slope(bat,lon,lat,N,method,extrapval)
%
% Find seafloor slope BET from bathymetry struct BAT at location(s) LON and
% LAT. If not present, adds fields .beta_deg and .beta to the returned
% version of BAT (optional fourth output arg). If they are calculated here,
% then the result depends on the number of points N per finite difference
% (DEFAULT: 2): if N<3, calc. finite difference with GRADIENTM (Map Toolkit);
% otherwise, passes N on to GRADIENTN (Ecoforecasts Toolkit). 
%
% Interpolates to LON,LAT by calling INTERP_FIELD (v.) with METHOD (DEFAULT:
% '*nearest') and EXTRAPVAL (DEFAULT: NaN). If LON or LAT is either missing
% or empty, just ensures BAT is populated with the abovementioned fields.
%
% ANG, angle of seafloor slope in degrees clockwise from True North, can also
% be returned; regardless of METHOD, ANG is interpolated using 'nearest'.
% Local isobath angle ISO can also be returned: this is MOD(ANG-90,360).
%
% NOTE: Assumes BAT.lon, BAT.lat in degrees, depths BAT.field in meters.
%
% SEE ALSO: FIND_NGDC_SLOPE_SITES: accepts a variety of forms of input, calls
% FIND_NGDC_SLOPE (this function), and returns a new STRUCT of site STRUCTs.
%
% Last Saved Time-stamp: <Mon 2017-05-01 16:58:39 Eastern Daylight Time gramer>

  bet=[];
  ang=[];
  iso=[];
  if ( iscell(bat) && numel(bat) == 3 && all(cellfun(@isnumeric,bat)) )
    cel = bat; bat=[];
    bat.lon = cel{1};
    bat.lat = cel{2};
    bat.field = cel{3};
    cel=[]; clear cel
  elseif ( ~isfield(bat,'lon') || ~isfield(bat,'lat') || ~isfield(bat,'field') )
    error('BAT must be a bathymetry struct with fields .lon,.lat,.field');
  end;

  if ( ~exist('N','var') || isempty(N) )
    N = 2;
  end;
  if ( ~isscalar(N) || ~isnumeric(N) || (floor(N) ~= N) )
    error('If specified, N must be an integer scalar!');
  end;
  if ( ~exist('method','var') || isempty(method) )
    method = '*nearest';
  end;
  if ( ~exist('extrapval','var') || isempty(extrapval) )
    extrapval = nan;
  end;

  [BLON,BLAT] = meshgrid(bat.lon,bat.lat);

  if ( ~isfield(bat,'beta') )
    if ( N < 3 )
      [aspect_deg,slope_deg,beta_y,beta_x] = gradientm(BLAT,BLON,bat.field,wgs84Ellipsoid);
    else
      hx = 1e3 .* distance_wgs84(bat.lat(1),bat.lon(1),bat.lat(1),bat.lon(2));
      hy = 1e3 .* distance_wgs84(bat.lat(1),bat.lon(1),bat.lat(2),bat.lon(1));
      [beta_x,beta_y] = gradientn(bat.field,N,hx,hy);
      hx=[]; hy=[]; clear hx hy
    end;

    % We want direction of steepest DESCENT - so invert the gradient
    beta_x = -beta_x; beta_y = -beta_y;
    bat.beta = uv_to_spd(beta_x,beta_y);
    bat.beta_deg = uv_to_dir_curr(beta_x,beta_y);
    aspect_deg=[]; slope_deg=[]; beta_y=[]; beta_x=[]; clear aspect_deg slope_deg beta_y beta_x
  end;

  if ( exist('lon','var') && ~isempty(lon) && exist('lat','var') && ~isempty(lat) )
    %bet = interp2(BLON,BLAT,bat.beta,lon,lat,method,NaN);
    bet = interp_field(bat.lat,bat.lon,bat.beta,lat,lon,method,extrapval);
    if ( nargout > 1 )
      %ang = interp2(BLON,BLAT,bat.beta_deg,lon,lat,'nearest',NaN);
      ang = interp_field(bat.lat,bat.lon,bat.beta_deg,lat,lon,'nearest',extrapval);

      if ( nargout > 2 )
        iso = mod(ang - 90, 360);
      end;
    end;

    if ( any(isnan(bet)) )
      warning('Ecoforecasts:find_ngdc_slope',...
              'One or more points not found in bathymetry region!');
    end;
  end;

  BLON=[]; BLAT=[]; clear BLON BLAT

return;
