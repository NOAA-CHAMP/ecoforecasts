function xes = brushpoints(varargin)
%error('THIS IS NOT READY FOR PRIME TIME!');

%function xes = brushpoints([xes|fh...])
%
% Use Data Cursor Mode to select the start and end of a range of points in
% figure FH (DEFAULT: GCF), and return a 2x1 matrix of XES (start and end
% ordinates) selected. If XES is passed in as a 2xN matrix, returns XES as a
% 2x(N+1) matrix with selected start and end X-coordinates appended to it.
%
% Last Saved Time-stamp: <Wed 2016-06-08 12:24:34 Eastern Daylight Time gramer>

error('THIS IS NOT READY FOR PRIME TIME!');

  xes = [];
  fh = get(0,'CurrentFigure');

  args = varargin(:);
  nargs = numel(args);
  while ( nargs > 0 )
    if ( ishandle(args{1}) )
      fh = args{1};
    elseif ( ismatrix(args{1}) )
      xes = args{1};
    end;
    args(1) = [];
    nargs = nargs - 1;
  end;

  if ( isempty(fh) || ~ishandle(fh) )
    error('No valid FIGURE found');
  end;
  if ( ~isempty(xes) && size(xes,1) ~= 2 )
    error('XES arg is not a 2xN matrix');
  end;

  h = datacursormode(fh);

  hi = getCursorInfo(h);
  if ( isempty(hi) )
    datacursormode on;
    disp(['Place datacursor and press any key to start deletion region']);
    pause;
    hi = getCursorInfo(h);
  else
    disp(['Using current data cursor as start of deletion']);
  end;

  begx = hi.Position(1);
  disp(['Move datacursor and press any key at end of deletion']);
  pause;
  endx = hi.Position(1);

  xes(:,end+1) = [begx,endx]';

return;
