function set_hovmuller_cursor(fh,fn)
%function set_hovmuller_cursor(fh,fn)
%
% For pcolor and similar 3-D plot calls, let user see C-values (fourth
% dimension) when in datacursor mode. If 'fh' is not specified, current
% figure is used, if any. Second arg defaults to @xyzc_select_cb.
%
% Last Saved Time-stamp: <Wed 2016-10-26 16:22:47 Eastern Daylight Time lew.gramer>

  if ( ~exist('fh','var') || isempty(fh) )
    fh = get(0,'CurrentFigure');
  end;
  if ( ~exist('fn', 'var') || isempty(fn) )
    fn = @tyzc_select_cb;
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
