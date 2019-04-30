function [FM,Stats,fh,lh] = scatter_curve_fit(x,y,ftype,xlbl,ylbl,fh,varargin)
%function [FM,Stats,fh,lh] = scatter_curve_fit(x,y,ftype,xlbl,ylbl,fh,varargin)
%
% Create a scatter-plot of X vs. Y (two vectors), do curve fit of type FTYPE
% using FIT (v.) and superimpose fitted curve onto the scatter plot. If XLBL
% (also used in legend) or YLBL are not given, try to figure them out. If
% valid FIGURE handle FH given, plot in it. If FH empty, make new figure. If
% FH 'none', *plot nothing*: merely return results. For a list of optional
% available args see FIT. Return fitted model FM, a struct of goodness-of-fit
% statistics STATS, FIGURE handle, and LINESERIES handle vector (see PLOT).
%
% Last Saved Time-stamp: <Sat 2018-08-25 21:59:06 Eastern Daylight Time gramer>

  if ( ~exist('fh','var') || isempty(fh) )
    fh = [];
  end;

  args = varargin;
  clev = 0.95;
  if ( numel(args) >= 2 && strcmpi(args{1},'confint') )
    clev = args{2};
    args = args(3:end);
  end;

  FM=[];
  Stats=[];
  lh=[];


  goodix = find( isfinite(x) & isfinite(y) );

  %[FM,Stats,Output] = fit(x(goodix),y(goodix),ftype,varargin{:});
  [FM,Stats,Output] = fit(x(goodix),y(goodix),ftype,args{:});

  Stats.N = numel(goodix);

  flds = fieldnames(Output);
  for fldix = 1:numel(flds)
    Stats.Output.(flds{fldix}) = Output.(flds{fldix});
  end;

  % Ensure yhat is the same size as x - NaNs and all
  Stats.yhat = feval(FM,x);
  Stats.R2 = Stats.adjrsquare;
  Stats.R2_1 = corr(y(goodix),Stats.yhat(goodix)).^2;
  % Stats.sse = Stats.dfe .* (Stats.rmse.^2);
  Stats.ssr = norm(Stats.yhat(goodix) - nanmean(Stats.yhat)).^2;
  Stats.R2_2 = 1 - (Stats.sse / (Stats.sse + Stats.ssr));

  % Ensure residuals and Jacobian same length as x - NaNs and all
  resid = repmat(nan,size(x));
  resid(goodix) = Stats.Output.residuals;
  Stats.Output.residuals = resid;
  if ( ~isempty(Stats.Output.Jacobian) )
    jacob = repmat(nan,[numel(x),size(Stats.Output.Jacobian,2)]);
    jacob(goodix,:) = Stats.Output.Jacobian;
    Stats.Output.Jacobian = jacob;
  end;

  if ( ischar(fh) && strcmpi(fh,'none') )
    %%%%%%%%%%%%%%%
    % EARLY RETURN
    %%%%%%%%%%%%%%%
    return;
  elseif ( ishandle(fh) )
    % Useful to pass in fh for, e.g., use in subplots
    figure(fh);
    hold on;
  elseif ( isempty(fh) )
    % fh = figure;
    % set(fh, 'units','normalized', 'outerposition',[0 0 1 1]);
    % hold on;
    fh = fmg;
  else
    error('Optional arg FH should be empty, FIGURE handle, or ''none''');
  end;

  if ( exist('xlbl','var') && ~isempty(xlbl) )
    xname = xlbl;
  else
    xname = inputname(1);
    if ( isempty(xname) )
      xname = 'X';
    end;
  end;
  if ( exist('ylbl','var') && ~isempty(ylbl) )
    yname = ylbl;
  else
    yname = inputname(2);
    if ( isempty(yname) )
      yname = 'Y';
    end;
  end;

  %fitstr = strtrim(evalc('disp(FM)'));
  % NOTE: GENRESULTS is a method of class CFIT (v.)
  resstr = genresults(FM,Stats,Output,[],[],[],clev);
  resstr(cellfun(@isempty,resstr)) = {''};
  resstr = strtrim(resstr);
  if ( Output.exitflag <= 0 )
    fitstr = sprintf('%s\n%s\n%s\n',resstr{1},resstr{6:7});
  else
    fitstr = sprintf('%s\n%s\n',resstr{1:2});
  end;
  modelstr = [strtok(fitstr,sprintf('\n')),' ',xname,' vs. ',yname];

  lh(1) = plot(x,y,'b.','MarkerSize',16);
  [sortx,sortix] = sort(x); 
  sortyhat = Stats.yhat(sortix);
  lh(2) = plot(sortx,sortyhat,'r.-','LineWidth',2);
  % legs1 = sprintf('%s, \n R^2~%0.2g (%0.2g;%0.2g), RMSE=%g, N=%d \n ', ...
  %                 modelstr, Stats.R2, Stats.R2_1, Stats.R2_2, ...
  %                 Stats.rmse, Stats.N);
  legs1 = sprintf('%s, \n R^2~%0.2g, RMSE=%g, N=%d', ...
                  modelstr, Stats.R2, ...
                  Stats.rmse, Stats.N);
  legs2 = sprintf(' \n %s',fitstr);

  legend( legs1, legs2, 'Location','Best' );
  xlabel(xname); ylabel(yname);
  grid on;
  titlename(modelstr);

return;
