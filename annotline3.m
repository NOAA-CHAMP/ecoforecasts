function [lineh,texth] = annotline3(X,Y,Z,txt,clr,linstyl,varargin)
%function [lineh,texth] = annotline3(X,Y,Z,txt,clr,linstyl[,varargin])
%
% Just like "ANNOTLINE" (v.), but accepts a Z coordinate (which may be []).
%
% Draw a straight line across the current axes in the given line color (the
% DEFAULT is axes line color), and annotated with the given text. If 'X' is
% empty, the start and end X coords of the line are the 'XLim' of the axes.
% If 'Y' is empty, start and end Y coords are the 'YLim' of axes. Otherwise
% matrices 'X' and 'Y' may be any shape supported by m-function LINE (qv.).
% VARARGIN: Additional name-value pairs for LINE - LineWidth, Marker*, etc.
% (If caller wishes to set TEXT properties, use TEXTH handle return value.)
%
% Last Saved Time-stamp: <Wed 2018-12-05 12:44:07 Eastern Standard Time gramer>

  lineh = [];
  texth = [];

  ax = gca;

  if ( isempty(X) )
    X = xlim(ax);
  end;
  if ( ~exist('Y', 'var') || isempty(Y) )
    Y = ylim(ax);
  end;
  if ( ~exist('Z', 'var') || isempty(Z) )
    Z = zlim(ax);
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

  maxl = max([length(X),length(Y),length(Z)]);
  while ( length(X) < maxl )
    X(end+1) = X(end);
  end;
  while ( length(Y) < maxl )
    Y(end+1) = Y(end);
  end;
  while ( length(Z) < maxl )
    Z(end+1) = Z(end);
  end;

  % Preserve axes limits
  xl = xlim(gca); yl = ylim(gca); zl = zlim(gca);
  lineh = line(X, Y, Z, 'Color', clr, 'LineStyle', linstyl,varargin{:});
  if ( ~isempty(txt) )
    texth = text(mean(X), mean(Y), mean(Z), txt, 'Color', clr);
  end;
  xlim(gca, xl); ylim(gca, yl); zlim(gca, zl);

return;
