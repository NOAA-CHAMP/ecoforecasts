function [rng,az,dlon,dlat] = georuler(varargin)
%function [rng,az,dlon,dlat] = georuler(varargin)
%
% Use GRULER (v.) to draw a line on the current longitude-latitude plot, and
% return the distance (RNG, km) and azimuthal heading (AZ, deg True) between
% the two end-points of the line. Also returns run (DLON) and rise (DLAT).
%
% Last Saved Time-stamp: <Fri 2011-05-06 07:52:36  lew.gramer>

  [ig,ig,x,y,dlon,dlat] = gruler(varargin{:});

  [rng,az] = distance_wgs84(y(1),x(1),y(2),x(2));

  if ( nargout < 1 )
    disp([rng,az]);
  end;

return;
