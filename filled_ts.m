function ts = filled_ts(varargin)
%function ts = filled_ts(varargin)
%
% Create a patchwork *hourly* time series TS, using every time series struct
% specified in arguments. (Each argument to FILLED_TS must be a struct with
% non-empty .date and .data fields, or a warning results.)
%
% Last Saved Time-stamp: <Thu 2011-05-19 15:05:49  Lew.Gramer>

  begdt = +Inf;
  enddt = -Inf;
  tses = struct('date',[],'data',[]);

  for ix=1:length(varargin)
    if ( ~is_valid_ts(varargin{ix}) )
      warning('Arg %g not a valid time series! Skipping...',ix);
    else
      tses(end+1) = varargin{ix};
      begdt = min([begdt,min(tses(end).date)]);
      enddt = max([enddt,max(tses(end).date)]);
    end;
  end;

  ts.date = [begdt:(1/24):enddt]';
  ts.data = repmat(nan,size(ts.date));

  % Fill in gaps in result time series from each of the arguments in turn
  for ix=1:length(tses)
    if ( ~isempty(tses(ix).date) )
      missingix = find(~isfinite(ts.data));
      [ix1,ix2] = intersect_dates(ts.date(missingix),tses(ix).date);
      if ( ~isempty(ix1) )
        ts.data(missingix(ix1)) = tses(ix).data(ix2);
      end;
    end;
  end;

  badix = find(~isfinite(ts.data));
  ts.date(badix) = [];
  ts.data(badix) = [];

return;
