function center_square_axes(ax)
%function center_square_axes(ax)
%
% Set identical limits for X- and Y-axes of axes AX (DEFAULT: GCA) and ensure
% that the origin [0,0] lies at the center of the figure.
%
% Last Saved Time-stamp: <Sat 2012-07-28 21:53:33  lew.gramer>

  if ( ~exist('ax','var') || isempty(ax) )
    ax = gca;
  end;
  xlm = xlim(ax);
  ylm = ylim(ax);
  lm = max(abs([xlm,ylm]));
  xlim(ax,[-lm,lm]);
  ylim(ax,[-lm,lm]);
  axis(ax,'square');

return;
