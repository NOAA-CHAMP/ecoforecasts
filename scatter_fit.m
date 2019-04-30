function [B,Stats,fh,lh] = scatter_fit(x,y,xlbl,ylbl,fh,wfun,tune,showPerfectFit,legendArgs)
%function [B,Stats,fh,lh] = scatter_fit(x,y,xlbl,ylbl,fh,wfun,tune,showPerfectFit,legendArgs)
%
% Create a scatter-plot of X vs. Y (two vectors), do simple linear regression
% using ROBUSTFIT, and superimpose fitted line onto the plot. If XLBL (also
% used in LEGEND) or YLBL are not given, try to figure them out. If valid
% FIGURE handle FH is given, plot in it. If FH is empty, make a new figure.
% If FH is 'none', *do not plot anything*: merely return results. For meaning
% of optional WFUN and TUNE see ROBUSTFIT. Return intercept, slope, a struct
% of regression statistics, FIGURE handle, and a vector of LINESERIES handles
% (see PLOT). If optional arg SHOWPERFECTFIT is present and true, also plot
% "perfect fit" line (slope 1, intercept 0).
%
% LEGENDARGS is passed to LEGEND (v.) DEFAULT: {'Location','NorthWest'}
%
% Last Saved Time-stamp: <Wed 2019-02-20 11:09:03 Eastern Standard Time gramer>

  if ( ~exist('fh','var') || isempty(fh) )
    fh = [];
  end;
  if ( ~exist('showPerfectFit','var') || isempty(showPerfectFit) )
    showPerfectFit = false;
  end;
  if ( ~exist('legendArgs','var') || isempty(legendArgs) )
    legendArgs = {'Location','NorthWest'};
  end;
  % Special option for less gargantuan LEGENDs
  shortLegend = false;
  if ( strncmpi(legendArgs{1},'short',5) )
    shortLegend = true;
    legendArgs(1) = [];
  end;

  B=[];
  Stats=[];
  lh=[];

  if ( isvector(x) && size(x,1) == 1 )
    x = x(:);
  end;
  if ( isvector(y) && size(y,1) == 1 )
    y = y(:);
  end;

  if ( exist('wfun','var') && ~isempty(wfun) )
    if ( exist('tune','var') && ~isempty(tune) )
      [B,Stats] = robustfit(x,y,wfun,tune);
    else
      [B,Stats] = robustfit(x,y,wfun);
    end;
  else
    [B,Stats] = robustfit(x,y);
  end;
  Stats.yhat = B(1) + (B(2).*x);
  Stats.R2_1 = corr(y,Stats.yhat).^2;
  Stats.sse = Stats.dfe .* (Stats.robust_s.^2);
  Stats.ssr = norm(Stats.yhat - nanmean(Stats.yhat)).^2;
  Stats.R2_2 = 1 - (Stats.sse ./ (Stats.sse + Stats.ssr));
  % Rank correlation
  [Stats.R2_3,Stats.p_3] = corr(y,Stats.yhat,'type','Spearman');
  Stats.R2_3 = Stats.R2_3.^2;

  X(1:length(x),1) = 1;
  X(1:length(x),2) = x(:);
  [ig,ig,ig,ig,Stats.regress_stats] = regress(y,X);
  Stats.N = min(length(find(isfinite(x))),length(find(isfinite(y))));

  addTitle = true;

  if ( ischar(fh) && strcmpi(fh,'none') )
    %%%%%%%%%%%%%%%
    % EARLY RETURN
    %%%%%%%%%%%%%%%
    return;
  elseif ( ishandle(fh) )
    % Useful to pass in fh for, e.g., use in subplots
    figure(fh);
    hold on;
    % In that case, do not add our usual title
    addTitle = false;
  elseif ( isempty(fh) )
    % New figure
    %fh = figure;
    %set(fh, 'units','normalized', 'outerposition',[0 0 1 1]);
    %hold on;
    fh = fmg;
  else
    error('Optional arg FH should be empty, FIGURE handle, or ''none''');
  end;

  set(gca,'FontSize',20);  % For publication-ready fonts when printing

  if ( exist('xlbl','var') && ~isempty(xlbl) )
    xname = xlbl;
  else
    xname = strrep(inputname(1),'_','\_');
    if ( isempty(xname) )
      xname = 'X';
    end;
  end;
  if ( exist('ylbl','var') && ~isempty(ylbl) )
    yname = ylbl;
  else
    yname = strrep(inputname(2),'_','\_');
    if ( isempty(yname) )
      yname = 'Y';
    end;
  end;

  % % NOTE TO SELF: Consider using GSCATTER (v.) as optional replacement for
  % % PLOT: user may specify optional grouping variable *or* FUNCTION_HANDLE
  %
  %     GSCATTER(X,Y,G) creates a scatter plot of the vectors X and Y grouped
  %     by G.  Points with the same value of G are shown with the same color
  %     and marker.  G is a grouping variable defined as a categorical
  %     variable, vector, cell array of strings, or string matrix, and it must
  %     have the same number of rows as X and Y.  Alternatively G can be a cell
  %     array of grouping variables (such as {G1 G2 G3}) to group the values in
  %     X by each unique combination of grouping variable values.
  %
  %     H = GSCATTER(...) returns an array of handles to the objects

  yhat = B(1)+(B(2).*x);

  lh(1) = plot(x,y,'b.','MarkerSize',16);
  lh(2) = plot(x,yhat,'r-','LineWidth',2);
  if ( showPerfectFit )
    lh(3) = plot(x,x,'k-');
  end;

  if ( shortLegend )
    legs = sprintf('R^2~%0.2g%s \\pm%0.2g (%d)', ...
                    Stats.regress_stats(1),sig_pstar(Stats.p(2)),Stats.s,Stats.N);
    legend( lh(2), legs, legendArgs{:});
  else
    legs1 = [xname ' vs. ' yname];
    legs2 = sprintf('%.5g + %.5g*X, R^2~%0.2g (%0.2g;%0.2g;%0.2g), \n RMSE=%g, p=%0.5g, N=%d', ...
                    B(1), B(2), (Stats.regress_stats(1)), Stats.R2_1, Stats.R2_2, Stats.R2_2, ...
                    Stats.s, roundn(Stats.p(2),-5), Stats.N);
    if ( ~showPerfectFit )
      legend( legs1, legs2, legendArgs{:});
    else
      legend( legs1, legs2, 'Perfect fit', legendArgs{:});
      minval = nanmin([nanmin(x(:)),nanmin(y(:)),nanmin(yhat(:))]);
      maxval = nanmax([nanmax(x(:)),nanmax(y(:)),nanmax(yhat(:))]);
      axis([minval,maxval,minval,maxval]); axis square
    end;
  end;
  xlabel(xname); ylabel(yname);
  grid on;
  if ( addTitle )
    titlename(legs1);
  end;

return;
