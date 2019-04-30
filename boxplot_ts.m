function [h,stats] = boxplot_ts(ts,grpfun,varargin)
%function [h,stats] = boxplot_ts(ts,grpfun,varargin)
%
% Create a BOXPLOT of the time series struct TS, grouping elements of TS by
% GRPFUN(TS.date). GRPFUN may be a function handle (DEFAULT: @GET_MONTH) or a
% string specifying grouping: 'hour', 'day', 'jday', 'annual' (uses grouping
% function GET_YEARDAY), 'triad', 'pentad', 'week', 'month', 'season', 'year'
% or any MATLAB function. All other args are name-value pairs: these args are
% passed through to BOXPLOT (v.) except for the following:
%
%  'allcolors' a single color to use for both box and outliers;
%  'indices' a vector of indices included in BOXPLOT or else a function
%     handle which accepts a time series and returns an index vector;
%  'mean' (DEFAULT: False) to also plot the NANMEAN for each period (arg may
%     simply be True, or a string specifying Marker per PLOT);
%  'std' (DEFAULT: False) to also plot the NANSTD about the mean for each
%     period (arg may be True, or a string specifying Marker per PLOT);
%  'title' a string to pass to TITLENAME (v.) or TITLE for the figure.
%  'usenames' (DEFAULT: False) to try to guess how to use names instead of
%     numbers for GROUPING: thus if GRPFUN is 'month' or @GET_MONTH, we would
%     guess that @GET_MONTH_NAME does the desired thing.
%
% BOXPLOT (v.) DEFAULTS: 'notch','on', 'whisker',2, 'positions',<GRPS>
%
% NOTE: TS may also be a cell array containing {DATE_VEC,DATA_VEC}
%
% Returns line handles H from BOXPLOT, and also an array STATS of the MAX;
% upper whisker, quartile, notch; MEDIAN; lower notch, quartile, whisker; and
% MIN resp., for each unique value returned by GRPFUN (DEFAULT: 12x9 array).
% (If 'notch'=='off', MEDIAN is substituted for both notch values in STATS.)
%
% SAMPLE CALL (all-black box-plot by year-month, black squares show MEAN):
%  >> boxplot_ts(TS,@(d)(datestr(d,12)),'allcolors','k','mean','s');
%
% Last Saved Time-stamp: <Wed 2018-08-29 17:18:56 Eastern Daylight Time gramer>

  h = [];
  stats = [];

  % if ( ~exist('ts','var') || ~is_valid_ts(ts) )
  %   error('First arg TS must be a valid time series struct');
  % end;
  if ( ~exist('ts','var') || (~is_ts(ts) && ~iscell(ts)) )
    error('First arg TS must be a time series STRUCT or CELL ARRAY {dts,dat}');
  end;
  if ( iscell(ts) )
    if ( numel(ts) ~= 2 || numel(ts{1}) ~= numel(ts{2}) || ~isnumeric(ts{1})  || ~isnumeric(ts{2}) )
      error('If first arg TS is a numeric CELL ARRAY {dts,dat}, numel(dts) must == numel(dat)');
    else
      tscel=ts; ts=[]; clear ts
      ts.date = tscel{1};
      ts.data = tscel{2};
      tscel=[]; clear tscel
    end;
  end;

  if ( ~exist('grpfun','var') || isempty(grpfun) )
    grpfun = @get_month;
  elseif ( ischar(grpfun) && isvector(grpfun) )
    f = grpfun(1:min(3,length(grpfun)));
    switch (lower(f)),
     case 'hou', grpfun=@get_hour;
     case 'day', grpfun=@floor;
     case 'jda', grpfun=@get_jday;
     case 'ann', grpfun=@get_yearday;
     case 'tri', grpfun=@get_triad;
     case 'pen', grpfun=@get_pentad;
     case 'wee', grpfun=@get_week;
     case 'mon', grpfun=@get_month;
     case 'sea', grpfun=@get_season;
     case 'yea', grpfun=@get_year;
     case 'yy',  grpfun=@get_year_2digit;
     otherwise,  grpfun=str2func(grpfun);
    end;
  end;
  if ( ~isa(grpfun,'function_handle') )
    error('Optional 2nd arg GRPFUN must be a grouping-period name or function name/handle');
  end;

  [ttl,varargin] = getarg(varargin,'TITL');

  [idxs,varargin] = getarg(varargin,'IND','default',1:numel(ts.date));
  if ( ischar(idxs) || isa(idxs,'function_handle') )
    idxs = feval(idxs,ts);
  end;
  if ( isempty(idxs) )
    error('Use of ''indices'' arg left no useable data!');
  end;
  ts.date = ts.date(idxs);
  ts.data = ts.data(idxs);

  [meanMarker,varargin] = getarg(varargin,'MEAN','default',false);
  [stdMarker,varargin] = getarg(varargin,'STD','default',false);
  [boxclr,varargin] = getarg(varargin,'ALLCOL');
  [useNames,varargin] = getarg(varargin,'USENAMES','default',false);


  dt = grpfun(ts.date);

  if ( isempty(boxclr) )
    % h = boxplot(ts.data,dt, 'notch','on', 'whisker',2, varargin{:});
    h = boxplot(ts.data,dt, 'notch','on', 'whisker',2, 'positions',dt, varargin{:});
  else
    % h = boxplot(ts.data,dt, 'notch','on', 'whisker',2, ...
    %             'symbol',[boxclr,'x'],'colors',[boxclr;boxclr;boxclr], varargin{:});
    h = boxplot(ts.data,dt, 'notch','on', 'whisker',2, 'positions',dt, ...
                'symbol',[boxclr,'x'],'colors',[boxclr;boxclr;boxclr], varargin{:});
  end;

  %ax = get(h(1),'Parent');
  ax = ancestor(h(1),'axes');

  if ( meanMarker || stdMarker )
    [cum,tid]=grp_ts(ts.data,ts.date,grpfun,@nanmean,[]);
    if ( stdMarker )
      [cums,tids]=grp_ts(ts.data,ts.date,grpfun,@nanstd,[]);
    end;
    % % Works in 2007a but not 2010a - BOXPLOT was rewritten!
    % x = get(ax,'XTick');
    holdStatus = get(ax,'NextPlot');
    set(ax,'NextPlot','add');
    if ( meanMarker && ~ischar(meanMarker) )
      %meanMarker = 's';
      %% X is already the marker for outliers!
      %meanMarker = 'x';
      meanMarker = 'd';
    end;
    if ( stdMarker && ~ischar(stdMarker) )
      stdMarker = '-';
    end;
    if ( meanMarker )
      if ( isempty(boxclr) )
        %plot(cum,['b',meanMarker],'MarkerSize',10);
        % Specify X-data ('tid') for this and all succeeding PLOT calls: this
        % is now required because we are (finally, mercifully) passing the
        % argument-pair 'positions',dt into BOXPLOT above.
        plot(tid,cum,['b',meanMarker],'MarkerSize',10);
      else
        plot(tid,cum,meanMarker,'Color',boxclr,'MarkerSize',10);
      end;
    end;
    if ( stdMarker )
      if ( isempty(boxclr) )
        plot(tid,cum+cums,['b',stdMarker],'MarkerSize',10);
        plot(tid,cum-cums,['b',stdMarker],'MarkerSize',10);
      else
        plot(tid,cum+cums,stdMarker,'Color',boxclr,'MarkerSize',10);
        plot(tid,cum-cums,stdMarker,'Color',boxclr,'MarkerSize',10);
      end;
    end;
    set(ax,'NextPlot',holdStatus);
  end;

  if ( useNames )
    xtks = xticks(ax);
    switch (upper(char(grpfun))),
     case 'GET_TRIAD',		xtkls = string(datestr(((xtks-1)*3)+1,'mmdd'));
     case 'GET_PENTAD',		xtkls = string(datestr(((xtks-1)*5)+1,'mmdd'));
     case 'GET_WEEK',		xtkls = string(datestr(((xtks-1)*7)+1,'mmdd'));
     case 'GET_MONTH',		xtkls = string(datestr(datenum(0,xtks,1),'mmm'));
     case 'GET_SEASON',		xtkls = strip(string(get_season_name(xtks)));
     otherwise,
      error('Do not know how to process "usenames" with GRPFUN %s',char(grpfun));
    end;
    xticklabels(ax,xtkls);
    if ( numel(unique(xtkls))*max(strlength(xtkls)) > 15 )
      set(ax,'XTickLabelRotation',90);
    end;
  end;

  if ( ~isempty(ttl) )
    if ( exist('titlename') > 1 )
      titlename(ttl);
    else
      title(ttl);
    end;
  end;

  if ( nargout < 1 )
    clear h;
  elseif ( nargout > 1 )
    % Returns line handles H from BOXPLOT, and also an array STATS of the MAX;
    % upper whisker, quartile, notch; MEDIAN; lower notch, quartile, whisker; and
    % MIN resp., for each unique value returned by GRPFUN (DEFAULT: 12x9 array).
    % (From HELP BOXPLOT:
    %  H has one column per box, consisting of the handles for the various
    %  parts of the box.  Each column contains 7 handles for the upper
    %  whisker, lower whisker, upper adjacent value, lower adjacent value,
    %  box, median, and outliers.)
    for perix=1:size(h,2)
      ph = h(:,perix);
      vs = get(ph(~isnan(ph)),'ydata');
      stats(perix,2) = max(vs{1});      % Upper whisker
      stats(perix,8) = min(vs{2});      % Lower whisker
      stats(perix,3) = max(vs{5});      % Upper quartile (box top)
      stats(perix,7) = min(vs{5});      % Lower quartile (box btm)
      stats(perix,5) = min(vs{6});      % Median (center line)
      if ( numel(vs) >= 7 )
        stats(perix,1) = max(vs{7});      % Max (top outlier)
        stats(perix,9) = min(vs{7});      % Min (btm outlier)
      else
        stats(perix,1) = max(vs{1});      % Upper whisker
        stats(perix,9) = min(vs{2});      % Lower whisker
      end;
      if ( numel(vs{5}) > 5 )
        stats(perix,4) = vs{5}(2);      % Upper notch (box line 2)
        stats(perix,6) = vs{5}(10);     % Lower notch (box line 10)
      else
        stats(perix,4) = min(vs{6});    % No notches
        stats(perix,6) = min(vs{6});    % No notches
      end
    end;
  end;

return;
