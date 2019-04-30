function newts = despike_ts(ts,tol,maxgap,datrng)
%function newts = despike_ts(ts,tol,maxgap,datrng)
%
% Process time series TS (struct with fields .date and .data) to remove
% spikes (e.g., from retrieving a sensor to rest on a hot ship deck). A
% "spike" is defined as any day where the one-day standard deviation is
% greater than TOL. If TOL is not specified, 1.50 times the peak one-day STD
% for the middle section (middle 80% of total record by index) is used. If TS
% has any gap longer than MAXGAP (DEFAULT: 1d), return an error. If DATRNG
% is given, calculate default TOL only from values between min/max(DATRNG).
%
% Last Saved Time-stamp: <Wed 2011-08-24 16:29:31  Lew.Gramer>

  if (~is_valid_ts(ts))
    error('Input TS does not seem to be a valid time series struct!');
  end;
  if ( ~exist('maxgap','var') || isempty(maxgap) )
    maxgap = 1;
  end;
  if ( ~exist('datrng','var') || isempty(datrng) )
    datrng = [-Inf,+Inf];
  end;
  if ( max(diff(ts.date)) > maxgap )
    error('Time series TS has a gap longer than MAXGAP=%g d',maxgap);
  end;

  begmidix = floor(0.10 * length(ts.date));
  endmidix = ceil(0.90 * length(ts.date));

  % If the caller did not specify TOL, calculate default tolerance as 150% of
  % the peak one-day Standard Deviation of the middle segment or our data.
  if ( ~exist('tol','var') || isempty(tol) )
    dts = ts.date(begmidix:endmidix);
    dat = ts.data(begmidix:endmidix);
    dts(min(datrng)>dat | dat>max(datrng)) = [];
    dat(min(datrng)>dat | dat>max(datrng)) = [];
    tol = 1.50*nanmax(grpstats(dat,floor(dts),@nanstd));
  end;

  if ( ~isnumeric(tol) || ~isscalar(tol) || tol <= 0 )
    error('Input TOL must be a numeric positive scalar!');
  end;

  % spikeix = find(abs(diff(ts.data))>tol)+1;
  dy = floor(ts.date);
  udy = unique(dy);
  std1d = grpstats(ts.data,dy,@nanstd);
  spikedy = udy(std1d>tol);
  spikeix = find(ismember(dy,spikedy));

  % Find first of tail-spikes
  endspikeix = find(spikeix>endmidix,1,'first');
  if ( ~isempty(endspikeix) )
    %endix = spikeix(endspikeix)-1;
    % Get last record for day immediately prior to first tail-spike
    endix = find(floor(ts.date)<floor(ts.date(spikeix(endspikeix))),1,'last');
    if ( isempty(endix) )
      endix = spikeix(endspikeix)-1;
    end;
  else
    endix = length(ts.date);
  end;

  % Find last of head-spikes
  begspikeix = find(spikeix<=begmidix,1,'last');
  if ( ~isempty(begspikeix) )
    %begix = spikeix(begspikeix)+1;
    % Get first record for day immediately after last head-spike
    begix = find(floor(ts.date)>floor(ts.date(spikeix(begspikeix))),1,'first');
    if ( isempty(begix) )
      begix = spikeix(begspikeix)+1;
    end;
  else
    begix = 1;
  end;

  newts.date = ts.date(begix:endix);
  newts.data = ts.data(begix:endix);

  %DEBUG:  fmg; plot_ts(ts,newts);

return;
