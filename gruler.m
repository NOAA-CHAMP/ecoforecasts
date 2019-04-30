function [dst,ang,x,y,dx,dy] = gruler(varargin)
%function [dst,ang,x,y,dx,dy] = gruler(varargin)
%
% Use GLINE (v.) to draw a line on the current X-Y plot, and return the
% length (DST) in plot units, and the angle to the vertical (deg clockwise)
% between the two end-points of the line. Also returns plot coordinates of
% the two end-points (X,Y), and the run (DX) and rise (DY) of the line. The
% line drawn by GLINE is removed immediately from the figure.
%
% Last Saved Time-stamp: <Wed 2012-02-22 15:14:01  lew.gramer>

  h = gline(varargin{:});
  pause;

  x = get(h,'XData');
  y = get(h,'YData');

  dx = diff(x);
  dy = diff(y);

  dst = sqrt((dx^2) + (dy^2));
  ang = rad2deg( angle(dx + (i*dy)) );

  delete(h);

  if ( nargout < 1 )
    disp([dst,ang]);
  end;

return;
