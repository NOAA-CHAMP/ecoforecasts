function [stn,B,Stats,fh,fix1,fix2]=scatter_fit_station(stn,fld1,fld2,fix1,fix2,xlbl,ylbl,statToPlot,showPerfectFit,lag)
%function [stn,B,Stats,fh,fix1,fix2]=scatter_fit_station(stn,fld1,fld2,fix1,fix2,xlbl,ylbl,statToPlot,showPerfectFit,lag)
%
% Regress (with SCATTER_FIT, qv.) matching timestamps of two "time series"
% (structs STN.FLD1 and STN.FLD2 each with fields .DATE and .DATA). If FIX1
% or FIX2 given, only compare those indices in FLD1/FLD2. If XLBL not given,
% use FLD1 as x-axis label. Ditto YLBL and FLD2.
%
% FIX1 and/or FIX2 may also be a FUNCTION_HANDLE (qv.). In this case, the
% function referred to should accept a "time series" (.DATE/.DATA struct) as
% its only required argument and return a vector of indices into it. Function
% return values then index the actual values regressed, and the returned FIX1
% and FIX2 are the intersection of the indices RETURNED by these functions.
%
% Returns updated STN struct (e.g., in case it was necessary to *create*
% either FLD1 or FLD2 using VERIFY_VARIABLE, qv.), regression coefficients B,
% statistics struct STATS (see ROBUSTFIT), figure handle FH, and fuzzy
% time-based intersection of time-series indices in FIX1 and FIX2. (See
% INTERSECT_DATES for timestamp-matching tolerance.)
%
% If STATTOPLOT is given and non-empty, also plot that residual or score
% field against the corresponding timestamps given in FLD1.date. FH then is a
% vector containing all figures plotted. STATTOPLOT may be a CELLSTR. 
%
% If SHOWPERFECTFIT is given and true, also plot a 1-to-1 perfect line.
%
% If LAG given, introduce lag of LAG *days* in data for FLD1. May be < 0.
%
% EXAMPLES:
%  >> smkf1 = load_station_data('smkf1');
%  >> scatter_fit_ts(smkf1,'sea_t','ct_shallow_seatemp',@ts_summer,...
%      @ts_summer,'SMKF1 NDBC Sea Temp','SMKF1 CT Sea Temp','resid');
%
% The above example regresses and plots the fit for coincident NDBC and CT
% Sea Temperatures from the SEAKEYS station "SMKF1" (Sombrero), restricting
% regression to only those indices return by the function TS_SUMMER; also
% produce a separate figure plotting non-normalized residuals.
%
% SEE ALSO: TS_JFM, TS_AMJ, TS_JAS, TS_OND, TS_BOREAL_WARM, TS_BOREAL_COOL,
%  SCATTER_FIT, INTERSECT_DATES, ROBUSTFIT, DATETICK, PLOT, FUNCTION_HANDLE.
%
% Last Saved Time-stamp: <Fri 2014-01-03 13:58:07 Eastern Standard Time gramer>

  stn = verify_variable(stn,fld1);
  stn = verify_variable(stn,fld2);

  if ( ~exist('xlbl','var') || isempty(xlbl) )
    xlbl = strrep(fld1, '_', '\_');
  end;
  if ( ~exist('ylbl','var') || isempty(ylbl) )
    ylbl = strrep(fld2, '_', '\_');
  end;

  if ( ~exist('fix1','var') || isempty(fix1) )
    fix1 = 1:numel(stn.(fld1).data);
  elseif ( isa(fix1,'function_handle') )
    xlbl = [xlbl ' (' upper(strrep(func2str(fix1),'_','\_')) ')'];
    fix1 = fix1(stn.(fld1));
  end;
  if ( ~exist('fix2','var') || isempty(fix2) )
    fix2 = 1:numel(stn.(fld2).data);
  elseif ( isa(fix2,'function_handle') )
    ylbl = [ylbl ' (' upper(strrep(func2str(fix2),'_','\_')) ')'];
    fix2 = fix2(stn.(fld2));
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


  [ix1,ix2] = intersect_dates(stn.(fld1).date(fix1)+lag,stn.(fld2).date(fix2));

  fix1 = fix1(ix1);
  fix2 = fix2(ix2);

  [B,Stats,fh] = scatter_fit(stn.(fld1).data(fix1),stn.(fld2).data(fix2),xlbl,ylbl,[],[],[],showPerfectFit);
  Stats.datenums = stn.(fld2).date(fix2);

  if ( ~isempty(statToPlot) )
    if ( ~isfield(Stats,statToPlot) )
      warning('No field "%s" found in ROBUSTFIT Stats!', statToPlot);
    else
      for ix = 1:length(statToPlot)
        st = statToPlot{ix};
        fh(end+1) = figure;
        plot(stn.(fld2).date(fix2),Stats.(st),'.');
        maxigraph;
        datetick3;
        titlename(sprintf('ROBUSTFIT Stats.%s: %s vs. %s', st, xlbl, ylbl));
      end;
      figure(fh(1));
    end;
  end;

  if ( nargout < 1 )
    stn = [];
    clear stn;
  end;

return;
