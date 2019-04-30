function [dst,hdg,lon,lat,dlon,dlat] = gruler_wgs84(varargin)
%function [dst,hdg,lon,lat,dlon,dlat] = gruler_wgs84(varargin)
%
% Use GLINE (v.) to draw a line interactively on current Longitude-Latitude
% plot, and call DISTANCE_WGS84 (v.) to return the length (DST; in KM), and
% heading angle (HDG; degT) between the two end-points of line. Also returns
% coordinates of the two end-points (LON,LAT), and run (DLON) and rise (DLAT)
% of line. The line drawn by GLINE is removed immediately from the figure.
%
% Last Saved Time-stamp: <Tue 2015-01-27 13:08:35 Eastern Standard Time gramer>

  if ( nargin > 0 )
    fh = varargin{1};
  else
    fh = get(0,'CurrentFigure');
  end;
  if ( ~ishandle(fh) )
    error('No valid figure');
  end;


  % GLINE does not work well with any of the Figure Modes
  plotedit(fh,'off');
  pan(fh,'off');
  zoom(fh,'off');
  rotate3d(fh,'off');
  datacursormode(fh,'off');

  h = gline(varargin{:});
  pause;

  lon = get(h,'XData');
  lat = get(h,'YData');
  z = [];
  try,
    z = get(h,'ZData');
  catch,
  end;
  c = [];
  try,
    c = get(h,'CData');
  catch,
  end;

  dlon = diff(lon);
  dlat = diff(lat);

  [dst,hdg] = distance_wgs84(lat(1),lon(1),lat(2),lon(2));

  delete(h);

  if ( nargout < 1 )
    disp([num2str(dst),' km, ',num2str(roundn(hdg,-2)),' deg']);
    if ( ~isempty(z) )
      disp([num2str(z),' m']);
    end;
    if ( ~isempty(c) )
      disp([num2str(c),' ?']);
    end;
    clear dst hdg lon lat dlon dlat
  end;

return;
