function [hlines,haxes,hfig] = multiplot_datetick(X,Y,ttl,xlbl,ylbl,xlm,ylms,datetick_fmt,linspec)
%function [hlines,haxes,hfig] = multiplot_datetick(X,Y,ttl,xlbl,ylbl,xlm,ylms,datetick_fmt,linspec)
%
% Call Matlab Exchange function MULTIPLOT, then call Matlab Exchange function
% DATETICK2 and local m-function SET_DATETICK_CURSOR on that multiplot. Use
% DATETICK_FMT as date format for multi-axis X. DEFAULT: DATETICK chooses.
% Optional arg YLMS is either a 1xN cell array of 1x2 vectors, or a single
% numeric 1x2 vector, that specifies y-axis limits for each of the N plots.
% Optional cell or char array arg LINSPEC is passed to MULTIPLOT (qv.); the
% DEFAULT is a cell array of strings as follows: {'b.','g.','r.','c.'}. Any
% optional arg may be empty, which implies "use a default for this arg".
%
% Uses EXIST to determine the most recent available incarnation of DATETICK
% (e.g., DATETICK2, DATETICK2_5) and calls that function with DATETICK_FMT.
%
% Last Saved Time-stamp: <Sun 2013-06-09 19:35:10 Eastern Daylight Time gramer>

  if (nargin < 3 || isempty(ttl))
    ttl = 'Time Series Data';
  end;
  if (nargin < 4 || (~ischar(xlbl) && isempty(xlbl))); xlbl = 'Date-Time'; end;
  if (nargin < 5 || isempty(ylbl)); ylbl = {}; end;
  if (nargin < 6); xlm = []; end;
  if (nargin < 7); ylms = []; end;
  if (nargin < 8); datetick_fmt = []; end;
  if (nargin < 9 || isempty(linspec)); linspec = {'b.','g.','r.','c.'}; end;
  if (~iscell(linspec)); linspec={linspec}; end;

  if ( ~isempty(ylms) && ~iscell(ylms) )
    ylms = repmat({ylms},size(X));
  end;

  if ( numel(linspec) > numel(Y) )
    linspec = linspec(1:numel(Y));
  elseif ( numel(linspec) < numel(Y) )
    linspec = repmat(linspec,[1,ceil(numel(Y)/numel(linspec))]);
    linspec = linspec(1:numel(Y));
  end;

  [hlines, haxes] = ...
      multiplot(X, Y, 'LineSpec',linspec, ...
                'Title',ttl, 'XLabel',xlbl, 'YLabel',ylbl, 'XLim',xlm);

  if ( exist('datetick2_5','file') )
    %DEBUG:    disp('Calling DATETICK2_5');
    dtfn = @datetick2_5;
  elseif ( exist('datetick2','file') )
    %DEBUG:    disp('Calling DATETICK2');
    dtfn = @datetick2;
  else
    %DEBUG:    disp('Calling DATETICK');
    dtfn = @datetick;
  end;
  if ( isempty(datetick_fmt) )
    dtfn('keeplimits','keepticks');
  else
    dtfn('x',datetick_fmt,'keeplimits','keepticks');
  end;

  if ( exist('set_datetick_cursor','file') )
    set_datetick_cursor;
  end;

  for ix = 1:length(haxes)
    set(get(haxes(ix),'YLabel'),'FontName','Times');
    % set(get(haxes(ix),'YLabel'),'FontSize',[7]);
    if ( ~isempty(ylms) )
      ylim(haxes(ix), ylms{ix});
    end;
  end;

  hfig = get(haxes(1), 'Parent');
  set(hfig, 'Units','normalized', 'OuterPosition',[0 0 1 1]);

return;
