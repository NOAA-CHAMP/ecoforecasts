function [Bs,Stats,ax] = scatter_fit_ts_seasons(fld1,fld2,fix1,fix2,xlbl,ylbl,fh,statToPlot,showPerfectFit,lag,axLms,tol,legendArgs)
%function [Bs,Stats,ax] = scatter_fit_ts_seasons(fld1,fld2,fix1,fix2,xlbl,ylbl,fh,statToPlot,showPerfectFit,lag,axLms,tol,legendArgs)
%
% Create four subplots, showing SCATTER_FIT_TS for each of the four seasons.
% Returns 4x1 cell arrays BS and STATS of each of results of SCATTER_FIT_TS.
% Optional arg AXLMS if non-empty forces all four subplots to have identical
% X-Y axes limits: if AXLMS is a 2-vector, set XLIM and YLIM to that value;
% if a 4-vector, set XLIM to AXLMS(1:2), YLIM to AXLMS(3:4); otherwise, e.g.,
% if AXLMS==True (DEFAULT: False), set limits to include all points on all
% plots. If STRNCMPI(AXLMS,'eq',2), set X and Y limits to same value, i.e.,
% square. Timestamp-matching tolerance TOL, if given, is passed through to
% INTERSECT_DATES (v.)
%
% LEGENDARGS is passed through to SCATTER_FIT (v.) If not given, DEFAULTS are
% season dependent: for Winter, Spring, and Fall {'Location','NorthWest'};
% for summer, LEGENDARGS={'Location','SouthWest'}.
%
% SEE ALSO: SCATTER_FIT_TS, SCATTER_FIT, GET_SEASON.
%
% Last Saved Time-stamp: <Thu 2017-07-06 23:35:39 Eastern Daylight Time gramer>

  if ( ~exist('xlbl','var') || isempty(xlbl) )
    xlbl = inputname(1);
    if ( isempty(xlbl) ); xlbl = 'X'; end;
  end;
  if ( ~exist('ylbl','var') || isempty(ylbl) )
    ylbl = inputname(2);
    if ( isempty(ylbl) ); ylbl = 'Y'; end;
  end;

  if ( ~exist('fix1','var') || isempty(fix1) )
    fix1 = 1:numel(fld1.data);
  elseif ( isa(fix1,'function_handle') )
    xlbl = [xlbl ' (' upper(strrep(func2str(fix1),'_','\_')) ')'];
    fix1 = fix1(fld1);
  elseif ( (~isnumeric(fix1) && ~islogical(fix1)) || ~isvector(fix1) )
    error('Optional third arg FIX1 must be a function handle, numeric- or logical-vector!');
  end;
  if ( ~exist('fix2','var') || isempty(fix2) )
    fix2 = 1:numel(fld2.data);
  elseif ( isa(fix2,'function_handle') )
    ylbl = [ylbl ' (' upper(strrep(func2str(fix2),'_','\_')) ')'];
    fix2 = fix2(fld2);
  elseif ( (~isnumeric(fix2) && ~islogical(fix2)) || ~isvector(fix2) )
    error('Optional fourth arg FIX2 must be a function handle, numeric- or logical-vector!');
  end;

  if ( ~exist('fh','var') || isempty(fh) )
    fh=figure;
    maxigraph;
  end;

  if ( ~exist('statToPlot','var') || isempty(statToPlot) )
    statToPlot = {};
  elseif ( ~iscell(statToPlot) )
    statToPlot = { statToPlot };
  end;

  if ( ~exist('showPerfectFit','var') || isempty(showPerfectFit) )
    showPerfectFit = false;
  end;

  if ( ~exist('lag','var') || isempty(lag) )
    lag = 0;
  end;

  if ( ~exist('tol','var') || isempty(tol) )
    tol = [];
  end;

  if ( ~exist('axLms','var') || isempty(axLms) )
    axLms = false;
  end;

  if ( ~exist('legendArgs','var') || isempty(legendArgs) )
    legendArgs = [];
  end;

  [ix1,ix2] = intersect_dates(fld1.date(fix1),fld2.date(fix2),tol);
  if ( min(diff(ix1)) < 0 || min(diff(ix2)) < 0 )
    % When regressing vs. non-hourly data - override TOLerance if needed
    tol = max(median(diff(fld1.date(fix1))),median(diff(fld2.date(fix2))));
    [ix1,ix2] = intersect_dates(fld1.date(fix1),fld2.date(fix2),tol);
  end;
  fix1 = fix1(ix1);
  fix2 = fix2(ix2);

  mn = [+Inf]; xmn = [+Inf]; ymn = [+Inf];
  mx = [-Inf]; xmx = [-Inf]; ymx = [-Inf];

  ax = repmat(nan,[1,4]);

  % JFM: Jan-Feb-Mar - Winter
  sfix1 = fix1(get_season(fld1.date(fix1))==1);
  if ( isempty(sfix1) )
    warning('Ecoforecasts:Scatter:Empty','No data for JFM');
  else
    if isempty(legendArgs); myLegendArgs = {'Location','NorthWest'};
    else myLegendArgs = legendArgs; end;

    ax(1) = subplot_tight(2,2,1);
    try,
      [Bs{1},Stats{1}] = scatter_fit_ts(fld1,fld2,sfix1,fix2,['JFM: ',xlbl],ylbl,fh,statToPlot,showPerfectFit,lag,tol,myLegendArgs);
      mn  = nanmin([mn(:),xlim(ax(1)),ylim(ax(1))]);
      xmn = nanmin([xmn(:),xlim(ax(1))]);
      ymn = nanmin([ymn(:),ylim(ax(1))]);
      mx  = nanmax([mx(:),xlim(ax(1)),ylim(ax(1))]);
      xmx = nanmax([xmx(:),xlim(ax(1))]);
      ymx = nanmax([ymx(:),ylim(ax(1))]);
    catch,
      catchwarn;
    end;
  end;
  title([]);

  % AMJ: Apr-May-Jun - Spring
  sfix1 = fix1(get_season(fld1.date(fix1))==2);
  if ( isempty(sfix1) )
    warning('Ecoforecasts:Scatter:Empty','No data for AMJ');
  else
    if isempty(legendArgs); myLegendArgs = {'Location','NorthWest'};
    else myLegendArgs = legendArgs; end;

    ax(2) = subplot_tight(2,2,2);
    try,
      [Bs{2},Stats{2}] = scatter_fit_ts(fld1,fld2,sfix1,fix2,['AMJ: ',xlbl],ylbl,fh,statToPlot,showPerfectFit,lag,tol,myLegendArgs);
      mn  = nanmin([mn(:),xlim(ax(2)),ylim(ax(2))]);
      xmn = nanmin([xmn(:),xlim(ax(2))]);
      ymn = nanmin([ymn(:),ylim(ax(2))]);
      mx  = nanmax([mx(:),xlim(ax(2)),ylim(ax(2))]);
      xmx = nanmax([xmx(:),xlim(ax(2))]);
      ymx = nanmax([ymx(:),ylim(ax(2))]);
    catch,
      catchwarn;
    end;
  end;
  title([]);

  % JAS: Jul-Aug-Sep - Summer
  sfix1 = fix1(get_season(fld1.date(fix1))==3);
  if ( isempty(sfix1) )
    warning('Ecoforecasts:Scatter:Empty','No data for JAS');
  else
    if isempty(legendArgs); myLegendArgs = {'Location','SouthWest'};
    else myLegendArgs = legendArgs; end;

    ax(3) = subplot_tight(2,2,3);
    try,
      [Bs{3},Stats{3}] = scatter_fit_ts(fld1,fld2,sfix1,fix2,['JAS: ',xlbl],ylbl,fh,statToPlot,showPerfectFit,lag,tol,myLegendArgs);
      mn  = nanmin([mn(:),xlim(ax(3)),ylim(ax(3))]);
      xmn = nanmin([xmn(:),xlim(ax(3))]);
      ymn = nanmin([ymn(:),ylim(ax(3))]);
      mx  = nanmax([mx(:),xlim(ax(3)),ylim(ax(3))]);
      xmx = nanmax([xmx(:),xlim(ax(3))]);
      ymx = nanmax([ymx(:),ylim(ax(3))]);
    catch,
      catchwarn;
    end;
  end;
  title([]);

  % OND: Oct-Nov-Dec - Fall
  sfix1 = fix1(get_season(fld1.date(fix1))==4);
  if ( isempty(sfix1) )
    warning('Ecoforecasts:Scatter:Empty','No data for OND');
  else
    if isempty(legendArgs); myLegendArgs = {'Location','NorthWest'};
    else myLegendArgs = legendArgs; end;

    ax(4) = subplot_tight(2,2,4);
    try,
      [Bs{4},Stats{4}] = scatter_fit_ts(fld1,fld2,sfix1,fix2,['OND: ',xlbl],ylbl,fh,statToPlot,showPerfectFit,lag,tol,myLegendArgs);
      mn  = nanmin([mn(:),xlim(ax(4)),ylim(ax(4))]);
      xmn = nanmin([xmn(:),xlim(ax(4))]);
      ymn = nanmin([ymn(:),ylim(ax(4))]);
      mx  = nanmax([mx(:),xlim(ax(4)),ylim(ax(4))]);
      xmx = nanmax([xmx(:),xlim(ax(4))]);
      ymx = nanmax([ymx(:),ylim(ax(4))]);
    catch,
      catchwarn;
    end;
  end;
  title([]);

  if ( isfinite(mn) && isfinite(mx) )
    if ( numel(axLms) == 1 )
      if ( axLms )
        subplots_set(fh,'XLim',[xmn,xmx],'YLim',[ymn,ymx]);
      end;
    elseif ( ischar(axLms) && strncmpi(axLms,'eq',2) )
      subplots_set(fh,'XLim',[mn,mx],'YLim',[mn,mx]);
    elseif ( numel(axLms) == 2 )
      subplots_set(fh,'XLim',axLms,'YLim',axLms);
    elseif ( numel(axLms) == 4 )
      subplots_set(fh,'XLim',axLms(1:2),'YLim',axLms(3:4));
    else
      warning('Ignoring arg AXLMS...')
    end;
  else
    warning('Ecoforecasts:Scatter:Empty','No valid data...');
  end;

  % Note: SCATTER_FIT_TS calls TITLENAME (v.) - we must undo that
  ttl = ['Seasons: ',xlbl,' vs. ',ylbl];
  if ( exist('suptitlename') == 2 )
    suptitlename(ttl);
  elseif ( exist('suptitle') == 2 )
    suptitle(ttl);
  else
    nm = strrep(ttl,'\_','_');
    set(fh,'Name',nm);
  end;
  % Return focus to upper-left panel in case caller wishes to assign a TITLE
  if ( ishandle(ax(1)) )
    subplot_tight(ax(1));
  elseif ( ishandle(ax(2)) )
    subplot_tight(ax(2));
  elseif ( ishandle(ax(3)) )
    subplot_tight(ax(3));
  elseif ( ishandle(ax(4)) )
    subplot_tight(ax(4));
  end;

return;
