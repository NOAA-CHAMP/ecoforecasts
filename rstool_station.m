function [B,Stats,stn,fixy,fixs]=rstool_station(stn,fldy,flds,fixy,fixs,ylbl,xlbls,statToPlot)
%function [B,Stats,stn,fixy,fixs]=rstool_station(stn,fldy,flds,fixy,fixs,ylbl,xlbls,statToPlot)
%
% Regress (with RSTOOL, qv.) matching timestamps of multiple "time series"
% (structs STN.(FLDY) and STN.(FLDS{:}), each with fields .DATE and .DATA).
% If FIXY or FIXS{:} are given, only compare those indices in FLDY/FLDS{:}.
% If YLBL not given, use FLDY as x-axis label. Ditto XLBLS and FLDS{:}.
%
% FIXY and/or FIXS{:} may also be a FUNCTION_HANDLE (qv.). In this case, each
% function referred to should accept a "time series" (.DATE/.DATA struct) as
% its only required argument and return a vector of indices into it. Function
% return values then index the actual values regressed, and the returned FIXY
% and FIXS{:} are the intersection of the indices RETURNED by these functions.
%
% Returns regression coefficients B, statistics struct STATS (see REGSTATS),
% updated STN struct (e.g., in case it was necessary to *create* either FLDY
% or any of FLDS{:} using VERIFY_VARIABLE, qv.), and fuzzy intersection of
% time-series indices in FIXY and FIXS{:}. (See INTERSECT_DATES for default
% fuzzy matching tolerance for time series timestamps.)
%
% If STATTOPLOT is given and non-empty, also plot that residual or score
% field against the corresponding timestamps given in FLDY.date.
%
%
% EXAMPLE:
%  >> smkf1 = load_station_data('smkf1');
%  >> rstool_station(smkf1,'sea_t',{'ct_shallow_seatemp','wind1_speed'},...
%      @ts_jas,{@ts_jas,@ts_jas},[],[],{'r','cookd'});
%
% The above example finds an optimal mixed regression between NDBC Sea
% Temperature and coincident CT Sea Temperature and Wind Speed at the SEAKEYS
% station "SMKF1" (Sombrero), restricting regression to only those indices
% return by the function TS_JAS (Jul-Aug-Sep); also do two separate figures,
% of non-normalized residuals, and Cook's distances (see REGSTATS), resp.
%
% (Note that due to memory constraints, certain REGSTATS are NOT returned.)
%
% SEE ALSO: TS_JFM, TS_AMJ, TS_JAS, TS_OND, TS_BOREAL_WARM, TS_BOREAL_COOL,
%  RSTOOL, REGSTATS, INTERSECT_DATES, DATETICK, PLOT, FUNCTION_HANDLE.
%
% Last Saved Time-stamp: <Thu 2012-12-06 15:32:21 Eastern Standard Time lew.gramer>

  B = [];
  Stats = [];

  if ( ~iscell(flds) )
    flds = { flds };
  end;

  nx = length(flds);

  stn = verify_variable(stn,fldy);
  for ix = 1:nx
    stn = verify_variable(stn,flds{ix});
  end;

  if ( ~exist('ylbl','var') || isempty(ylbl) )
    ylbl = fldy;
  end;
  if ( ~exist('xlbls','var') || isempty(xlbls) )
    xlbls = flds;
  end;

  if ( ~exist('fixy','var') || isempty(fixy) )
    fixy = 1:numel(stn.(fldy).data);
  elseif ( isa(fixy,'function_handle') )
    ylbl = [ylbl ' (' upper(func2str(fixy)) ')'];
    fixy = fixy(stn.(fldy));
  end;
  if ( ~exist('fixs','var') || isempty(fixs) )
    for ix = 1:nx
      fixs{ix} = 1:numel(stn.(flds{ix}).data);
    end;
  else
    for ix = 1:nx
      if ( isa(fixs{ix},'function_handle') )
        xlbls{ix} = [xlbls{ix} ' (' upper(func2str(fixs{ix})) ')'];
        fixs{ix} = fixs{ix}(stn.(flds{ix}));
      end;
    end;
  end;

  if ( ~exist('statToPlot','var') )
    statToPlot = {};
  elseif ( ~iscell(statToPlot) )
    statToPlot = { statToPlot };
  end;

  for ix = 1:nx
    [ixy,ixx] = intersect_dates(stn.(fldy).date(fixy),stn.(flds{ix}).date(fixs{ix}));
    fixy = fixy(ixy);
    fixs{ix} = fixs{ix}(ixx);
  end;
  for ix = 1:nx
    [ixy,ixx] = intersect_dates(stn.(fldy).date(fixy),stn.(flds{ix}).date(fixs{ix}));
    fixs{ix} = fixs{ix}(ixx);
  end;

  ny = length(fixy);
  Y(1:ny,1) = stn.(fldy).data(fixy);
  for ix = 1:nx
    X(1:ny,ix) = stn.(flds{ix}).data(fixs{ix});
  end;

  % All the REGSTATS statistics we're likely to need or want?
  regstats_stats = ...
      { ...
          'beta', ...
          'covb', ...
          'r', ...
          'mse', ...
          'rsquare', ...
          'adjrsquare', ...
          'standres', ...
          'studres', ...
          'dffit', ...
          'dffits', ...
          'covratio', ...
          'cookd', ...
          'tstat', ...
          'fstat', ...
          'dwstat', ...
      };
  Stats = regstats(Y,X,'quadratic',regstats_stats);
  Stats.rmse = sqrt(Stats.mse);
  B = Stats.beta;

  rstool(X,Y,'quadratic',0.05,xlbls,ylbl);

  if ( ~isempty(statToPlot) )
    if ( ~isfield(Stats,statToPlot) )
      warning('No field "%s" found in ROBUSTFIT Stats!', statToPlot);
    else
      for ix = 1:length(statToPlot)
        st = statToPlot{ix};
        figure;
        plot(stn.(fldy).date(fixy),Stats.(st),'.');
        maxigraph;
        datetick3;
        titlename(sprintf('ROBUSTFIT Stats.%s: %s vs. %s', st, ylbl, sprintf('%s,', xlbls{:})));
      end;
    end;
  end;

return;
