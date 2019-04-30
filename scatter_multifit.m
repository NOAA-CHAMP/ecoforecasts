function [B,Stats,fh] = scatter_multifit(x,y,xlbl,ylbl,fh,wfun,tune)
%function [B,Stats,fh] = scatter_multifit(x,y,xlbl,ylbl,fh,wfun,tune)
%
% Create a scatter-plot of 'x' vs. 'y' (matrix vs vector), do a linear regression
% using ROBUSTFIT and superimpose fitted line on the plot. If 'xlbl' (also
% used in legend) or 'ylbl' are not given, try to figure them out. If figure
% handle 'fh' is not specified, create new figure. For meaning of optional
% 'wfun' and 'tune' arguments, see ROBUSTFIT. Return intercept and slope of
% regression, a struct of regression statistics, and figure-handle.
%
% Last Saved Time-stamp: <Fri 2010-06-11 10:43:19 Eastern Daylight Time gramer>

error('INCOMPLETE FUNCTION! Try using RSTOOL_STATION instead.');

% See below for missing lines...

  if ( exist('wfun','var') && ~isempty(wfun) )
    if ( exist('tune','var') && ~isempty(tune) )
      [B,Stats] = robustfit(x,y,wfun,tune);
    else
      [B,Stats] = robustfit(x,y,wfun);
    end;
  else
    [B,Stats] = robustfit(x,y);
  end;

  if ( exist('fh','var') && ~isempty(fh) )
    % Useful to pass in fh for, e.g., use in subplots
    figure(fh);
    hold on;
  else
    fh = figure;
    set(fh, 'units','normalized', 'outerposition',[0 0 1 1]);
    hold on;
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

  maxix = min(2,size(x,1));

  plot(x(1:maxix,:),y,'b.');
  plot(x,(B(1)+(B(2:(maxix+1)).*x(1:maxix,:)),'r-');


%%%%%%% FINISH EDITING ME HERE!
sprintf('%.5g + %.5g*X'
%%%%%%% FINISH EDITING ME HERE!
  for ix = 2:
  legend( [strrep(xname,'_','\_') ' vs. ' strrep(yname,'_','\_')], ...
           sprintf('%.5g + %.5g*X, \n RMSE=%g, p=%.5g, N=%d', ...
                   B(1), B(2), Stats.s, Stats.p(2), length(x)), ...
           'Location','NorthWest');
  xlabel(xname); ylabel(yname);

return;
