function square_axes(ax)
%function square_axes(ax)
%
% Set identical limits for X- and Y-axes of axes AX (DEFAULT: GCA)
%
% Last Saved Time-stamp: <Sat 2012-07-28 21:48:42  lew.gramer>

  if ( ~exist('ax','var') || isempty(ax) )
    ax = gca;
  end;
  xlm = xlim(ax);
  ylm = ylim(ax);
  lms=[min([xlm(1),ylm(1)]),max([xlm(2),ylm(2)])];
  xlim(ax,lms);
  ylim(ax,lms);
  axis(ax,'square');

return;
