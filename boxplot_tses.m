function [hs,stats,lh] = boxplot_tses(tses,grpfun,varargin)
%function [hs,stats,lh] = boxplot_tses(tses,grpfun,varargin)
%
% Create a BOXPLOT of each time series struct in the array or cell array of
% structs TSES, grouping elements of each TS by GRPFUN(TS.date). See HELP
% BOXPLOT_TS for a description of GRPFUN and all other optional arguments.
% By default, choose a distinct color (BOXPLOT_TS arg ALLCOLOR) and/or width
% (BOXPLOT arg WIDTHS) for each TS (v. 'widthFirst','widthAndColor' below).
%
% Optional args: 'LEGEND',{'legend_string_1',...} to add a boxplot LEGEND;
% 'WIDTHFIRST' (DEFAULT: false) to vary box widths first; 'WIDTHANDCOLOR'
% (DEFAULT: false) to vary both width and color with every TS (limits number
% of distinguishably drawn TSes). All other args passed to BOXPLOT_TS (v.).
%
% Returns cell arrays of BOXPLOT handle vectors and BOXPLOT_TS STATS arrays.
% If LEGEND is requested, also returns legend handle LH.
%
% NOTE ALSO: Chaos results if all TSES do not share data for the same unique
% values returned by GRPFUN (e.g., if some TSES have no data for January).
%
% Last Saved Time-stamp: <Mon 2013-05-06 19:07:41 Eastern Daylight Time gramer>


  if ( nargin<2 )
    grpfun = [];
  end;

  hs = {};
  stats = {};
  lh = [];

  clrs = {'k','b',[0,.5,0],'r','m','c',[.3,.3,.3],[.8,.8,.8]};
  %wids = [0.5,0.35,0.25,0.1];
  wids = [0.5,0.4,0.3,0.2,0.1];

  [legs,varargin] = getarg(varargin,'legend');
  [widthAndColor,varargin] = getarg(varargin,'widthAndColor','default',false);
  [widthFirst,varargin] = getarg(varargin,'widthFirst','default',false);

  for tix=1:numel(tses)
    if ( widthAndColor )
      clr = clrs{mod(tix-1,numel(clrs))+1};
      wid = wids(mod(tix-1,numel(wids))+1);
    elseif ( widthFirst )
      wid = wids(mod(tix-1,numel(wids))+1);
      clrix = ceil(tix/numel(wids));
      clrix = mod(clrix-1,numel(clrs))+1;
      clr = clrs{clrix};
    else
      clr = clrs{mod(tix-1,numel(clrs))+1};
      widix = ceil(tix/numel(clrs));
      widix = mod(widix-1,numel(wids))+1;
      wid = wids(widix);
    end;
    [h,stat] = boxplot_ts(tses{tix},grpfun,'allcol',clr,'widths',wid,varargin{:});
    hs{end+1} = h;
    stats{end+1} = stat;
  end;

  % Add LEGEND if caller requested it
  if ( ~isempty(legs) )
    for legix=1:numel(hs)
      lh = hs{legix};
      legh(legix) = lh(1,1);
    end;
    lh = legend(legh,legs);
  end;

  if ( nargout < 3 )
    lh = []; clear lh;
    if ( nargout < 2 )
      stats = {}; clear stats;
      if ( nargout < 1 )
        hs = []; clear hs;
      end;
    end;
  end;

return;
