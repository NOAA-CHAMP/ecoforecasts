function newts = find_pct_good_periods(ts,periods,pct,goodfun,dt,maxgap)
%function newts = find_pct_good_periods(ts,periods,pct[,goodfun,[dt[,maxgap]]])
%
% Return a subset of time series TS (i.e., STRUCT with .date,.data) with only
% those points, for which there are at least PCT (0<PCT<1) good values in
% each time period: time periods are determined by PERIODS, a CHAR that must
% be one of (for AVERAGES) 'year', 'yearseason', 'yearmonth', 'yearweek',
% 'yearday', 'yearhour', or  (for CLIMATOLOGIES) 'season', 'month', 'week',
% 'jday', or 'hour'. (TO DO: 'yearpentad', 'yeartriad'; 'pentad', 'triad').
%
% "Good" indices are those returned by GOODFUN (DEFAULT: TS_ISFINITE). To run
% this function on ALL timestamps in a time series regardless of value, e.g.,
% even those whose value is NaN or Inf, specify "@(x)(1:numel(x.date))".  By
% default, however, NEWTS.data will not contain any non-finite values.
%
% Optional DT is time series resolution (DEFAULT: MEDIAN(DIFF(TS.date)).
%
% Optional MAXGAP (DEFAULT: +Inf) is the longest allowable individual time
% stamp gap in order for a period to be included in NEWTS. This is applied in
% ADDITION (increases restriction) to the percent-good test described above.
%
% CALLS: For CLIMATOLOGIES, calls RUN_LENGTH_ENCODE (v.) to determine the
% number of data points within each climatological bin, both for "all dates"
% between the first period and last period, and for the actual time series.
%
%
% SAMPLE CALL #1 - returns a new time series containing only non-NaN values,
%  and only for those weeks of the original time series that contained at
%  least 85% of their possible values, e.g., for an hourly time series, at
%  least 0.85*24*7 = 143 data points per weekly period:
%
% >> newts = find_pct_good_periods(ts,'yearseason',0.85);
%
%
% SAMPLE CALL #2 - similar to above, but allows up to 50% of each possible
% year-MONTH to be missing so long as no gap is longer than 7 days:
%
% >> newts = find_pct_good_periods(ts,'yearmonth',0.50,[],[],7);
%
%
% Last Saved Time-stamp: <Wed 2017-05-31 16:34:55 Eastern Daylight Time lew.gramer>

  if ( ~exist('ts','var') || ~is_valid_ts(ts) )
    error('TS must be a valid time series STRUCT');
  end;
  if ( ~exist('periods','var') || ~ischar(periods) || length(periods) < 2 )
    error('PERIODS must be a CHAR string at least two characters long');
  end;
  if ( ~exist('pct','var') || ~isnumeric(pct) || ~isscalar(pct) || 0 > pct || pct > 1 )
    error('PCT must be numeric scalar between 0 and 1');
  end;
  if ( ~exist('goodfun','var') || isempty(goodfun) )
    goodfun = @ts_isfinite;
  end;
  if ( ischar(goodfun) )
    goodfun = str2func(goodfun);
  end;
  if ( ~isa(goodfun,'function_handle') || exist(char(goodfun))<2 )
    error('GOODFUN if specified must be convertible to a valid FUNCTION_HANDLE');
  end;


  % Process PERIODS string
  switch ( periods ),


   % AVERAGES over each period

   case {'year'},
    grpfun = @(x)(datenum(get_year(x),1,1));
    selfun = @unique;
    begdt = grpfun(ts.date(1));
    enddt = datenum(get_year(ts.date(end))+1,1,1);

   case {'yearseason'},
    grpfun = @get_yearseason;
    selfun = @unique;
    begdt = grpfun(ts.date(1));
    [se,yr] = get_season(ts.date(end));
    % GET_SEASON returns SE between 1 and 4, so DATENUM(start of year,SE*3,1)
    % will always be the first day of the season *after* TS.DATE(end)
    enddt = datenum(yr,(se*3)+1,1);

   case {'yearmonth'},
    grpfun = @get_yearmonth;
    selfun = @unique;
    begdt = grpfun(ts.date(1));
    [mo,yr] = get_month(ts.date(end));
    enddt = datenum(yr,mo+1,1);

   case {'yearweek'},
    grpfun = @get_yearweek;
    selfun = @unique;
    begdt = grpfun(ts.date(1));
    [wk,yr] = get_week(ts.date(end));
    % GET_WEEK always returns WK between 1 and 52, so (start of year + WK*7)
    % will always be the first day of the week *after* TS.DATE(end)
    enddt = datenum(yr,1,1) + (wk*7);

   case {'yearday'},
    grpfun = @floor;
    selfun = @unique;
    begdt = grpfun(ts.date(1));
    enddt = ceil(ts.date(end));

   case {'yearhour'},
    grpfun = @get_yearhour;
    selfun = @unique;
    begdt = grpfun(ts.date(1));
    enddt = floor(ts.date(end)) + ((get_hour(ts.date(end)) + 1)/24);


   % CLIMATOLOGIES over all periods

   case {'season'},
    grpfun = @get_season;
    selfun = @(x)(run_length_encode(sort(x)));
    begdt = get_yearseason(ts.date(1));
    [yr,se] = get_yearseason(ts.date(end));
    enddt = datenum(yr,(se*3)+1,1);

   case {'month'},
    grpfun = @get_month;
    selfun = @(x)(run_length_encode(sort(x)));
    begdt = get_yearmonth(ts.date(1));
    [yr,mo] = get_yearmonth(ts.date(end))
    enddt = datenum(yr,mo+1,1);

   case {'day','jday'},
    grpfun = @get_jday;
    selfun = @(x)(run_length_encode(sort(x)));
    begdt = floor(ts.date(1));
    enddt = ceil(ts.date(end));

   case {'hour','jhour'},
    grpfun = @get_hour;
    selfun = @(x)(run_length_encode(sort(x)));
    begdt = get_yearhour(ts.date(1));
    enddt = floor(ts.date(end)) + ((get_hour(ts.date(end)) + 1)/24);

   otherwise,
    error('Period "%s" not yet implemented!',periods);
  end;

  % Time-series resolution
  if ( ~exist('dt','var') || isempty(dt) )
    dt = median(diff(ts.date));
  end;

  % Maximum allowable individual gap within a period
  if ( ~exist('maxgap','var') || isempty(maxgap) )
    maxgap = +inf;
  end;
  if ( maxgap < dt )
    error('Time series resolution DT=%g is less than MAXGAP=%g!',dt,maxgap);
  end;


  alldts = begdt:dt:enddt;
  % Make sure we don't include the beginning of the next period in our array
  alldts(alldts == enddt) = [];

  % Use GOODFUN to subset TS
  ts = subset_ts(ts,goodfun);
  if ( ts.date(1) < alldts(1) || alldts(end) < ts.date(end) )
    error('Somehow PERIODS processing returned values that do not encompass our dates?!');
  end;

  d = grpfun(ts.date);
  ud = selfun(d);

  alld = grpfun(alldts);
  ualld = selfun(alld);

  idx = [];
  % NOTE: For e.g., 'yearhour' on a big dataset, this loop could be slow-ish! 
  for uix = 1:numel(ud)
    dix = find(ismember(d,ud(uix)));
    if ( ~isempty(dix) )
      allix = find(ismember(alld,ud(uix)));
      % NOTE: MAXGAP is checked vs. actual time series time stamps!
      permaxgap = max(diff(ts.date(dix)));
      if ( ((numel(dix)/numel(allix)) >= pct) && (permaxgap <= maxgap) )
        % Would be nice to figure out a reliable way to preallocate array IDX
        idx(end+1:end+numel(dix)) = dix;
      end;
    end;
  end;

  if ( isempty(idx) )
    newts = struct('date',[],'data',[]);
  else
    newts = subset_ts(ts,idx);
  end;

return;
