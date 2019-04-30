function ax = inset_image(fname,pos)
%function ax = inset_image(fname,pos)
%
% Display the image file FNAME in the current FIGURE at position POS
% (DEFAULT: [.5,.5,.5,.5]. Creates a new AXES to contain the inset image, but
% then restores the original current AXES of that figure as drawing target.
%
% CALLS: AXES('Position',POS), IMSHOW(FNAME,...)
%
% Last Saved Time-stamp: <Wed 2018-02-28 16:07:29 EST lew.gramer>

  fh = get(groot,'CurrentFigure');
  if ( isempty(fh) || ~ishandle(fh) )
    error('No CurrentFigure');
  end;
  pax = get(fh,'CurrentAxes');

  if ( ~exist(fname,'file') )
    error('No image file named "%s"!',fname);
  end;
  if ( ~exist('pos','var') || isempty(pos) )
    pos = [0.5,0.5, 0.5,0.5];
  end;

  ax = axes('Position',pos);
  imshow(fname,'Parent',ax,'Border','tight','InitialMagnification','fit');

  %axes(pax); % D'oh!
  % Ensure parent AXES continues to be the drawing target, but DON'T raise it
  set(fh,'CurrentAxes',pax);

return;
