function [FM,Stats,fh,fix1,fix2] = scatter_curve_fit_ts(fld1,fld2,ftype,fix1,fix2,xlbl,ylbl,fh,lag,plotResid,varargin)
%function [FM,Stats,fh,fix1,fix2] = scatter_curve_fit_ts(fld1,fld2,ftype,fix1,fix2,xlbl,ylbl,fh,lag,plotResid,varargin)
%
% Fit a curve (with SCATTER_CURVE_FIT, qv.) of FTYPE on matching timestamps
% of two time series structs FLD1 and FLD2 (each with fields .DATE, .DATA).
% If FIX1 or FIX2 given, only compare those indices in FLD1 and/or FLD2. If
% XLBL not given, try INPUTNAME(1) as x-label. Ditto YLBL and INPUTNAME(2).
% If valid FIGURE handle FH is given, plot in it. If FH is empty, make a new
% figure. If FH is 'none', *do not plot anything*: merely return results.
%
% FIX1 and/or FIX2 may also be a FUNCTION_HANDLE (qv.). In this case, the
% function referred to should accept a "time series" (.DATE/.DATA struct) as
% its only required argument and return a vector of indices into it. Function
% return values then index the actual values regressed, and the returned FIX1
% and FIX2 are the date intersection of indices RETURNED by these functions.
%
% Returns fitted model FM, statistics struct STATS (see FIT), figure handle
% FH, and fuzzy time-based intersection of time-series indices in FIX1 and
% FIX2. (See INTERSECT_DATES for timestamp-matching tolerance.)
%
% If LAG given, introduce lag of LAG *days* in data for FLD1. May be < 0.
%
% If PLOTRESID is true, also plot a new figure of residuals time series.
%
% EXAMPLES:
%  >> mlrf1 = load_station_data('mlrf1');
%  >> scatter_curve_fit_ts(mlrf1.wind1_speed,mlrf1.sea_t,@ts_jas,@ts_jas,...
%      'MLRF1 Wind Speed (kts)','MLRF1 Sea Temp (oC)',[],1,true);
%
% The above example regresses and plots an order-3 polynomial fit for Wind
% Speed and the NDBC Sea Temperature one day later (LAG=1), from SEAKEYS
% station "MLRF1" (Molasses Reef), restricting regression to indices returned
% by the function TS_JAS (Jul-Aug-Sep), and does separate plot of residuals.
%
% SEE ALSO: TS_JFM, TS_AMJ, TS_JAS, TS_OND, TS_BOREAL_WARM, TS_BOREAL_COOL,
%  SCATTER_CURVE_FIT, SCATTER_FIT, INTERSECT_DATES, PLOT, FUNCTION_HANDLE.
%
% Last Saved Time-stamp: <Fri 2011-12-09 13:15:59  Lew.Gramer>

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
    fh = [];
  end;

  if ( ~exist('plotResid','var') || isempty(plotResid) )
    plotResid = false;
  end;

  if ( ~exist('lag','var') || isempty(lag) )
    lag = 0;
  end;


  [ix1,ix2] = intersect_dates(fld1.date(fix1)+lag,fld2.date(fix2));

  fix1 = fix1(ix1);
  fix2 = fix2(ix2);

  [FM,Stats,fh] = scatter_curve_fit(fld1.data(fix1),fld2.data(fix2),ftype,...
                                    xlbl,ylbl,fh,varargin{:});

  if ( ischar(fh) && strcmpi(fh,'none') )
    % No figure(s) desired

  elseif ( plotResid )
    % Plot time series of residuals from curve fit
    figure;
    dts = fld2.date(fix2);
    plot(dts,Stats.Output.residuals,'.');
    maxigraph;
    if ( (dts(end)-dts(1)) >= 365 )
      datetick3('x',12,'keeplimits');
    else
      datetick3;
    end;
    grid on;
    fitstr = strtok(evalc('disp(FM)'),sprintf('\n'));
    titlename(sprintf('%s Residual: %s vs. %s', fitstr, xlbl, ylbl));

    % Return focus to scatter-plot
    figure(fh);
  end;

return;
