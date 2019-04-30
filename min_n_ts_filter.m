function minN = min_n_ts_filter(cumfun,dt,sumfun,useleapyear)
% Over a one-year period, sampled once-per-DT (DEFAULT: DT=1/24), what is the
% minimum (or SUMFUN, e.g., MEDIAN) number of data points accumulated by the
% time-series filter function CUMFUN, e.g., MIN_N_TS_FILTER(@GET_MONTH)
% returns 28 days * 24 hours = 672 data points, or for 5-minute sampling and
% yearly accumulation, MIN_N_TS_FILTER(@GET_YEAR,(5/24/60)) returns 105120.
% This function ignores leap years, unless USELEAPYEAR is True, so that, for
% example, MIN_N_TS_FILTER(@GET_MONTH,1/24,@MIN,TRUE) returns 672.
%
% Last Saved Time-stamp: <Wed 2016-09-07 13:33:07 Eastern Daylight Time lew.gramer>

  if ( ~exist('dt','var') || isempty(dt) )
    dt = 1/24;
  end;
  if ( ~exist('sumfun','var') || isempty(sumfun) )
    sumfun = @min;
  end;
  if ( ~exist('useleapyear','var') || isempty(useleapyear) )
    useleapyear = false;
  end;

  if ( useleapyear )
    example_year = 2004;
  else
    example_year = 2005;
  end;

  dts = datenum(example_year,1,1):dt:(datenum(example_year+1,1,1)-(dt/2));
  udts = unique(cumfun(dts));
  for dtix=1:numel(udts);
    n(dtix) = numel(find(cumfun(dts) == udts(dtix)));
  end;

  minN = feval(sumfun,n);

return;
