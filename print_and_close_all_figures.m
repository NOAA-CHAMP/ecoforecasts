function print_and_close_all_figures(filePrefix, fmt)
%FUNCTION print_and_close_all_figures(filePrefix, fmt)
%
% PRINT out to a PS file, and then CLOSE, every figure window currently open
% If 'fmt' is 'png', print to a PNG file instead of PostScript. Ditto 'jpg'.
%
% Note: The 'close()' is NOT optional: this func doesn't work without it...
%
%  RCSID: $Id$
RCSID = '$Id$';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Revision History:
% $Source$
%
% $Log$
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  if ( ~exist('filePrefix', 'var') )
    filePrefix = '';
  end
  if ( ~isempty(filePrefix) )
    filePrefix = [ filePrefix '-' ];
  end

  if ( ~exist('fmt', 'var') || isempty(fmt) )
    fmt = 'ps';
  end

  gidx = get(0, 'currentfigure');
  while ~isempty(gidx)

    printFileName = sprintf('%s%g.%s', filePrefix, gidx, fmt);
    switch ( lower(fmt) )
     case 'ps',
      print(gidx, '-dps2c', printFileName);
     case 'png',
      print(gidx, '-dpng', printFileName);
     case {'jpg', 'jpeg'},
      print(gidx, '-djpeg', printFileName);
     otherwise,
      error('Unknown print format %s!', fmt);
    end;

    disp(['Created ' printFileName]);
    close(gidx);
    gidx = get(0, 'currentfigure');
  end

return;
