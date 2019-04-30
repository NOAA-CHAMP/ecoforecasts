function bat = calc_ngdc_slope(bat,N)
%function bat = calc_ngdc_slope(bat,N)
%
% Add seafloor slope and aspect angle fields .beta and .beta_deg, resp., to
% bathymetry struct BAT. Gradient is calculated using an N (DEFAULT: 2) point
% finite difference: if arg N<3, call GRADIENTM with wgs84Ellipsoid (see Map
% Toolkit); otherwise, pass N to GRADIENTN, with HX and HY calculated using
% DISTANCE_WGS84 (see Ecoforecasts Toolkit) on southwesternmost points. 
%
% Aspect angle of seafloor slope is in degrees clockwise from True North.
%
% NOTE: Assumes BAT.lon, BAT.lat in degrees, depths BAT.field in meters.
%
% Last Saved Time-stamp: <Sun 2016-12-11 22:05:42 Eastern Standard Time gramer>

  if ( ~isfield(bat,'lon') || ~isfield(bat,'lat') || ~isfield(bat,'field') )
    error('BAT must be a bathymetry struct with fields .lon,.lat,.field');
  end;
  if ( ~exist('N','var') || isempty(N) )
    N = 2;
  end;
  if ( ~isscalar(N) || ~isnumeric(N) || (floor(N) ~= N) )
    error('If specified, N must be an integer scalar!');
  end;

  [LON,LAT] = meshgrid(bat.lon,bat.lat);

  if ( N < 3 )
    [aspect_deg,slope_deg,beta_y,beta_x] = gradientm(LAT,LON,bat.field,wgs84Ellipsoid);
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

  LON=[]; LAT=[]; clear LON LAT

return;
