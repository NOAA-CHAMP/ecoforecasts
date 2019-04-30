function [lineh,texth] = annotline(X,Y,txt,clr,linstyl,varargin)
%function [lineh,texth] = annotline(X,Y,txt,clr,linstyl[,varargin])
%
% Draw a straight line across the current axes in the given line color (the
% DEFAULT is axes line color), and annotated with the given text. If 'X' is
% empty, the start and end X coords of the line are the 'XLim' of the axes.
% If 'Y' is empty, start and end Y coords are the 'YLim' of axes. Otherwise
% matrices 'X' and 'Y' may be any shape supported by m-function LINE (qv.).
% VARARGIN: Additional name-value pairs for LINE - LineWidth, Marker*, etc.
% (If caller wishes to set TEXT properties, use TEXTH handle return value.)
%
% Last Saved Time-stamp: <Fri 2017-03-17 15:56:59 Eastern Daylight Time gramer>

  lineh = [];
  texth = [];

  ax = gca;

  if ( isempty(X) )
    X = xlim(ax);
  end;
  if ( ~exist('Y', 'var') || isempty(Y) )
    Y = ylim(ax);
  end;
  if ( ~exist('txt', 'var') )
    txt = [];
  end;
  if ( ~exist('clr', 'var') || isempty(clr) )
    % clr = 'default';
    clr = 'r';
  end;
  % Let user specify Matlab's funky '-', etc. line styles
  if ( ~exist('linstyl', 'var') || isempty(linstyl) )
    linstyl = 'default';
  end;

  while ( length(X) < length(Y) )
    X(end+1) = X(end);
  end;
  while ( length(Y) < length(X) )
    Y(end+1) = Y(end);
  end;

  % Preserve axes limits
  xl = xlim(gca); yl = ylim(gca);
  lineh = line(X, Y, 'Color', clr, 'LineStyle', linstyl,varargin{:});
  if ( ~isempty(txt) )
    texth = text(mean(X), mean(Y), txt, 'Color', clr);
  end;
  xlim(gca, xl); ylim(gca, yl);

return;
