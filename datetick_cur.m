function datetick_cur(fh, fn)
%function datetick_cur(fh, fn)
% Simply calls SET_DATETICK_CURSOR (see).

  % Stupid argument handling in MATLAB...
  if ( ~exist('fh', 'var') || isempty(fh) )
    fh = get(0,'CurrentFigure');
  end;
  if ( ~exist('fn', 'var') || isempty(fn) )
    fn = @tyz_select_cb;
  end;

  set_datetick_cursor(fh, fn);

return;
