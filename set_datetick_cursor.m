function set_datetick_cursor(fh, fn)
%function set_datetick_cursor(fh, fn)
%
% For plot calls where 'datetick' (or 'datetick2') has been called, let user
% see date strings in place of raw datenums for X values. If 'fh' is not
% specified, current figure is used. Second arg defaults to @tyz_select_cb.
%
% Last Saved Time-stamp: <Wed 2016-10-26 17:45:32 Eastern Daylight Time lew.gramer>

  if ( ~exist('fh', 'var') || isempty(fh) )
    fh = get(0,'CurrentFigure');
  end;
  if ( ~exist('fn', 'var') || isempty(fn) )
    fn = @tyz_select_cb;
  end;

  if ( exist('isOctave') && isOctave )
    warning('No DATACURSORMODE in Octave');
  else
    if ( ~isempty(fh) )
      h = datacursormode(fh);
      set(h, 'UpdateFcn', fn);
    end;
  end;

return;
