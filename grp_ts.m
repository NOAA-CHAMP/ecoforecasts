function [cum,tid,nPerYr,dt,cumfun,per,dat,dts] = grp_ts(dat,dts,per_or_cumfun,sumfun,minN)
%function [cum,tid,nPerYr,dt,cumfun,per,dat,dts] = grp_ts(dat,dts,per_or_cumfun,sumfun,minN)
%
% Calculate cumulative statistic CUM on time series in vectors DAT and DTS,
% over the grouping PER_OR_CUMFUN which must either be a period string (e.g.,
% 'pentad','weekly','monthly','seasonal'; DEFAULT 'daily') or the handle of a
% cumulator function, e.g., @GET_JDAY (v.) SUMFUN specifies the statistic to
% return (DEFAULT: @NANMEDIAN), and MINN minimum count of samples/period; any
% "year-period" with fewer than MINN samples is is excluded from CUM. If MINN
% is 'none', *all* data in DAT is used; otherwise, DEFAULT is MINN equal to
% the maximum number of samples in a full PER. Also returns time index TID of
% statistic, so NUMEL(CUM)==NUMEL(TID)==NUMEL(UNIQUE(CUMFUN(DTS))). NPERYR is
% the unique elements returned by CUMFUN for a full leap-year of DT periods:
% NPERYR=NUMEL(UNIQUE(CUMFUN(DATENUM(2004,1,1):DT:DATENUM(2004,12,31,23,59,59))))
% where DT=MEDIAN(DIFF(UNIQUE(DTS))) estimates time "granularity" of DTS.
%
% Calls GRPSTATS(DAT,CUMFUN(DTS),SUMFUN) to calculate cumulative statistic.
%
% Last Saved Time-stamp: <Fri 2016-07-22 15:10:45 Eastern Daylight Time gramer>

  if ( ~exist('dat','var') || ~exist('dts','var') || numel(dat)~=numel(dts) )
    error('First and second arg must be numeric vectors of the same length!');
  end;

  if ( ~exist('per_or_cumfun','var') || isempty(per_or_cumfun) )
    per_or_cumfun = 'daily';
  end;
  if ( ~exist('sumfun','var') || isempty(sumfun) )
    if ( exist('nanmedian') > 1 )
      sumfun = @nanmedian;
    else
      sumfun = @median;
    end;
  end;
  if ( ~exist('minN','var') || isempty(minN) || strcmpi(minN,'default') )
    % DEFAULT is applied further down
    minN = [];
  end;

  if ( isa(per_or_cumfun,'function_handle') )
    cumfun = per_or_cumfun;
    per = char(per_or_cumfun); % Not currently used outside this IF block
    default_minN = 0;
  elseif ( ischar(per_or_cumfun) )
    per = per_or_cumfun;
    % Allowing leap-days might cause problems when comparing cum stats on
    % *different time series*, e.g., a time series which includes data for
    % one or more 29th's of Feb, and one which does not. If caller does not
    % like this behavior, they can simply specify custom CUMFUN and MINN.
    switch ( lower(per(1)) ),
     case {'h','hourly'},   cumfun=@get_jhour_no_leap;	default_minN=1;
     case {'d','daily'},    cumfun=@get_jday_no_leap;	default_minN=24;
     case {'p','pentad'},   cumfun=@get_pentad;		default_minN=5*24;
     case {'w','weekly'},   cumfun=@get_week;		default_minN=7*24;
     case {'m','monthly'},  cumfun=@get_month;		default_minN=28*24; %February
     case {'s','seasonal'}, cumfun=@get_season;		default_minN=90*24;
     case {'a','annual','y','year'}, cumfun=@get_year;	default_minN=365*24;
     otherwise, error('Unsupported cumulative period string "%s"!',char(per));
    end;
  else
    error('Third arg if given must be a char string or function handle!');
  end;
  clear per_or_cumfun;

  if ( nargout > 2 )
    dt = median(diff(unique(dts)));
    nPerYr = numel(unique(cumfun(datenum(2004,1,1):dt:datenum(2004,12,31,23,59,59))));
  end;

  % If caller did not pass minN of 'none', remove incomplete periods
  if ( isempty(minN) )
    minN = default_minN;
  end;
  if ( isnumeric(minN) && minN > 0 )
    % Eliminate unuseable values
    dts(~isfinite(dat)) = [];
    dat(~isfinite(dat)) = [];

    % Eliminate incomplete periods
    tiyr = get_year(dts);
    tiper = cumfun(dts);
    ti = datenum(tiyr,1,0) + tiper;
    uti = unique(ti);
    n = grpstats(dat,ti,'numel');
    badix = find(ismember(ti,uti(n < minN)));
    dts(badix) = [];
    dat(badix) = [];
  end;

  if ( isempty(dat) )
    error('No useable data to accumulate!');
  end;


  % Calculate cumulative statistic
  ti = cumfun(dts);
  if ( nargout > 1 )
    tid = unique(ti);
  end;
  cum = grpstats(dat,ti,sumfun);

return;
