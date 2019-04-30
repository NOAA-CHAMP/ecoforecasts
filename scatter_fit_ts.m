function [B,Stats,fh,fix1,fix2,lh] = scatter_fit_ts(fld1,fld2,fix1,fix2,xlbl,ylbl,fh,statToPlot,showPerfectFit,lag,tol,legendArgs)
%function [B,Stats,fh,fix1,fix2,lh] = scatter_fit_ts(fld1,fld2,fix1,fix2,xlbl,ylbl,fh,statToPlot,showPerfectFit,lag,tol,legendArgs)
%
% Regress (with SCATTER_FIT, qv.) matching timestamps of two "time series"
% (structs FLD1 and FLD2 each with fields .DATE and .DATA). If FIX1 or FIX2
% given, only compare those indices in FLD1 and/or FLD2. If XLBL not given,
% use INPUTNAME(1) as x-axis label if non-empty. Ditto YLBL and INPUTNAME(2).
% If valid FIGURE handle FH is given, plot in it. If FH is empty, make a new
% figure. If FH is 'none', *do not plot anything*: merely return results.
%
% FIX1 and/or FIX2 may also be a FUNCTION_HANDLE (qv.). In this case, the
% function referred to should accept a "time series" (.DATE/.DATA struct) as
% its only required argument and return a vector of indices into it. Function
% return values then index the actual values regressed, and the returned FIX1
% and FIX2 are the intersection of the indices RETURNED by these functions.
%
% Returns regression coefficients B, statistics struct STATS (see ROBUSTFIT),
% figure handle FH, and fuzzy time-based intersection of time-series indices
% in FIX1 and FIX2. (Timestamp-matching tolerance TOL if given is passed
% through to INTERSECT_DATES, please see.) Optionally returns line handles LH
% returned by SCATTER_FIT (v.) also.
%
% If STATTOPLOT is given and non-empty, also plot that residual or score
% field against the corresponding timestamps given in FLD1.date. FH then is a
% vector containing all figures plotted. STATTOPLOT may be a CELLSTR. 
%
% If SHOWPERFECTFIT is given and true, also plot a 1-to-1 perfect line.
%
% If LAG given, introduce lag of LAG *days* in data for FLD1. May be < 0.
%
% LEGENDARGS (DEFAULT: {}) is passed through to SCATTER_FIT (v.)
%
% EXAMPLES:
%  >> smkf1 = load_station_data('smkf1');
%  >> mlrf1 = load_station_data('mlrf1');
%  >> scatter_fit_ts(smkf1.sea_t,mlrf1.sea_t,@ts_summer,@ts_summer,...
%      'SMKF1 Sea Temp','MLRF1 Sea Temp','resid');
%
% The above example regresses and plots the fit for coincident NDBC Sea
% Temperatures from the two SEAKEYS stations "SMKF1" and "MLRF1" (Sombrero
% and Molasses), restricting regression to only those indices return by the
% (*hypothetical*) function TS_SUMMER, also plotting unnormalized residuals.
%
% SEE ALSO: TS_JFM, TS_AMJ, TS_JAS, TS_OND, TS_BOREAL_WARM, TS_BOREAL_COOL,
%  SCATTER_FIT, INTERSECT_DATES, ROBUSTFIT, DATETICK, PLOT, FUNCTION_HANDLE.
%
% Last Saved Time-stamp: <Thu 2019-02-14 14:27:13 Eastern Standard Time gramer>

  if ( ~exist('xlbl','var') || isempty(xlbl) )
    xlbl = strrep(inputname(1),'_','\_');
    if ( isempty(xlbl) ); xlbl = 'X'; end;
  end;
  if ( ~exist('ylbl','var') || isempty(ylbl) )
    ylbl = strrep(inputname(2),'_','\_');
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

  if ( ~exist('legendArgs','var') || isempty(legendArgs) )
    legendArgs = {};
  end;


  [ix1,ix2] = intersect_dates(fld1.date(fix1)+lag,fld2.date(fix2),tol);

  fix1 = fix1(ix1);
  fix2 = fix2(ix2);

  [B,Stats,fh,lh] = scatter_fit(fld1.data(fix1),fld2.data(fix2),xlbl,ylbl,fh,[],[],showPerfectFit,legendArgs);
  Stats.datenums = fld2.date(fix2);

  if ( ischar(fh) && strcmpi(fh,'none') )
    %%%%%%%%%%%%%%%
    % EARLY RETURN
    %%%%%%%%%%%%%%%
    return;
  end;

  if ( ~isempty(statToPlot) )
      for ix = 1:length(statToPlot)
        st = statToPlot{ix};
        if ( ~isfield(Stats,st) )
          warning('No field "%s" found in ROBUSTFIT Stats!', st);
        else
          fh(end+1) = figure;
          plot(fld2.date(fix2),Stats.(st),'.');
          maxigraph;
          datetick3;
          titlename(sprintf('ROBUSTFIT Stats.%s: %s vs. %s', st, xlbl, ylbl));
        end;
      end;
      figure(fh(1));
  end;

return;
