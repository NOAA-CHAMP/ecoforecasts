function [climts,lhs] = grpplot_ts(ts,varargin)
%function [climts,lhs] = grpplot_ts(ts,varargin)
%
% PLOT cumulative statistic CUM on time series struct TS. All other arguments
% are passed through to GRP_TS (v.). Makes a guess (based on NPERYR returned
% by GRP_TS) at the correct timestamp format to add to X-axis of plot using
% DATETICK (v.) if extant, and calls SET_DATETICK_CURSOR (v.) if extant; this
% same guess is used to convert group numbers returned by GRP_TS into DATENUM
% (from Day 0) in CLIMTS.date. Also returns line handles LHS from PLOT (v.).
%
% CALLS: GRP_TS(TS.data,TS.date,VARARGIN{...}) with up to three args from
% VARARGIN, to calculate cumulative statistic. All args in VARARGIN after the
% third (if any) are given to PLOT(CLIMTS.date,CLIMTS.data,varargin{4:end}).
%
% Last Saved Time-stamp: <Wed 2011-12-07 18:58:59  Lew.Gramer>

  climts = [];
  lhs = [];

  nargs = numel(varargin);
  grpargs = varargin(1:min(nargs,3));

  [climts.data,climts.date,nPerYr,dt] = grp_ts(ts.data,ts.date,grpargs{:});

  dtstr=[];
  if ( climts.date(1) > 100000 )
    % DATENUM
    dtstr=2;
  elseif ( 8785 <= nPerYr && nPerYr < 17569 )
    % Year-half-hour
    climts.date=datenum(1,1,1)+(climts.date.*48);
    dtstr=6;
  elseif ( 367 <= nPerYr && nPerYr < 8785 )
    % Year-hour
    climts.date=datenum(1,1,1)+(climts.date.*24);
    dtstr=6;
  elseif ( 74 <= nPerYr && nPerYr < 367 )
    % Year-day or Julian day
    climts.date=datenum(1,1,1)+(climts.date-1);
    dtstr=6;
  elseif ( 53 <= nPerYr && nPerYr < 74 )
    % Pentad
    climts.date=datenum(1,1,1)+((climts.date-1).*5);
    dtstr=6;
  elseif ( 25 <= nPerYr && nPerYr < 53 )
    % Week
    climts.date=datenum(1,1,1)+((climts.date-1).*7);
    dtstr=6;
  elseif ( 13 <= nPerYr && nPerYr < 25 )
    % Hour of day
    climts.date=datenum(1,1,1)+(climts.date./24);
    dtstr=15;
  elseif ( 5 <= nPerYr && nPerYr < 13 )
    % Month
    climts.date=datenum(1,climts.date,1);
    dtstr=3;
  elseif ( 3 <= nPerYr && nPerYr < 5 )
    % Season / Quarter
    climts.date=datenum(1,(((climts.date-1).*3)+1),1);
    dtstr=18;
  elseif ( 2 <= nPerYr && nPerYr < 3 )
    % Semi-annual
    climts.date=datenum(1,(((climts.date-1).*6)+1),1);
    dtstr=3;
  elseif ( nPerYr < 2 )
    % Interannual
    climts.date=datenum(climts.date,1,1);
    dtstr=10;
  end;

  lhs = plot(climts.date,climts.data,varargin{4:end});

  dt = min(diff(unique(climts.date)))/2;
  xlim([nanmin(climts.date)-dt,nanmax(climts.date)+dt]);

  if ( exist('datetick') > 1 && ~isempty(dtstr) )
    datetick('x',dtstr,'keeplimits');
    if ( exist('set_datetick_cursor') > 1 )
      set_datetick_cursor;
    end;
  end;

  if ( nargout < 2 )
    clear lhs
    if ( nargout < 1 )
      clear climt
    end;
  end;

return;
