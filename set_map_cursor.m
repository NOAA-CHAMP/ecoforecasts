function set_map_cursor(fh,fn)
%function set_map_cursor(fh,fn)
%
% For MAP Toolkit (v.) plots, let user see Longitude and Latitude plus any
% Z-value if present, when in datacursor mode. If 'fh' is not specified,
% current figure is used, if any. Second arg DEFAULT: @XYZC_MAP_SELECT_CB.
% (If indices are a nuisance, specify FN=@XYZC_SELECT_CB instead.)
%
% Last Saved Time-stamp: <Fri 2017-04-07 17:25:39 Eastern Daylight Time gramer>

warning('THIS "TMW recommended" CODE CAUSES SOME WEIRD ERRORS IN MATLAB 2016B!')

  if ( ~exist('fh','var') || isempty(fh) )
    fh = get(0,'CurrentFigure');
  end;
  if ( ~exist('fn', 'var') || isempty(fn) )
    fn = @xyz_map_select_cb;
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
