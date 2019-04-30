function [climatology, sd, lopct, uppct] = climatologize_time_series(dts,time_series,ts_name,pct)
%function [climatology, sd, lopct, uppct] = climatologize_time_series(dts,time_series,ts_name,pct)
%
% FOR A MORE GENERAL FUNCTION, see: ANOMALIZE_TS.m.
%
% Create one or more daily climatologies from the given time stamped time
% series. DTS is a vector of MATLAB DATENUMs, TIME_SERIES is a vector of time
% series values: their sizes must match exactly. A 1x365 vector with the DAILY
% climatology is returned, as well as a 1x365 vector of standard deviations.
% (NOTE: Estimate 'sd' is NOT validated by this code: if t.s. distribution is
% far from normal, then 'lopct'/'uppct' should be used instead. See below.)
% 
% Optional third arg 'ts_name' is a name string for the variable, used in
% warnings and errors. If not specified, messages are less specific.
% 
% Optional fourth arg 'pct' is a percentage (0 < pct < 100) to use in
% returning lower and upper percentile ranges 'lopct' and 'uppct' for each
% day of the climatology. If 'pct' not specified, then 3% (the equivalent of
% 3 standard deviations for a Gaussian distribution) is assumed.
% 
% Last Saved Time-stamp: <Sun 2016-05-01 17:01:05 Eastern Daylight Time gramer>

  if ( ~exist('ts_name', 'var') || isempty(ts_name) )
    ts_name = 'unknown_var';
  end;
  if ( ~exist('pct', 'var') || isempty(pct) )
    pct = 3;
  end;

  climatology = repmat(nan, [1 365]);
  sd = repmat(nan, [1 365]);
  lopct = repmat(nan, [1 365]);
  uppct = repmat(nan, [1 365]);

  [Y M D] = datevec(dts);

  jdays = floor(dts) - datenum(Y, 1, 1) + 1;
  % Stupid leap years
  jdays(jdays == 366) = 365;

  for dayidx = 1:365
    intermediate_data = time_series(jdays == dayidx);
    intermediate_data = intermediate_data(~isnan(intermediate_data));

    % Warn about periods with too little data!
    N = length(intermediate_data);
    if ( N < 3 )
      error('Climatologize:InsufficientSample', ...
            'Insufficient data (%s)! Only %d data points for JDay %d', ...
            ts_name, N, dayidx);
    end;
    if ( N < 24 )
      warning('Climatologize:SmallSample', ...
              'Only %d hours of valid data (%s) for JDay %d', ...
              N, ts_name, dayidx);
    end;

    climatology(dayidx) = mean(intermediate_data);
    sd(dayidx) = std(intermediate_data);
    lopct(dayidx) = prctile(intermediate_data, pct);
    uppct(dayidx) = prctile(intermediate_data, (100-pct));
  end;

return;
