function newts = interp_ts(ts,sper,newsper,method,maxgap,ringlen,varargin)
%function newts = interp_ts(ts,[sper],[newsper],[method],[maxgap],[ringlen],[extrap])
%
% Interpolate on time series struct TS with expected sample period SPER, to
% produce a new regular time series struct NEWTS with sample period NEWSPER
% (DEFAULT: 1/24). Calls INTERP1 (v.) with arg METHOD (DEFAULT: 'spline') to
% interpolate; then calls FILTER_GAPS (v.) with MAXGAP (DEFAULT: greater of
% 3*SPER or 3*NEWSPR) and Gibbs ringing length RINGLEN (DEFAULT: 0 for METHOD
% 'nearest' or 'linear', SPER for all others) to remove any gaps in original
% TS from NEWTS. Optional arg EXTRAP is any of extrapolation args accepted by
% INTERP1 (v.). NOTE: If EXTRAP non-empty, then FILTER_GAPS is *not* called.
%
% NOTE: If SPER is not specified, it is calculated as MEDIAN(DIFF(TS.date)).
% If NEWSPER is the string 'same', it is set to be the same as SPER.
%
% Last Saved Time-stamp: <Mon 2018-02-19 15:22:06 Eastern Standard Time gramer>

  if ( ~isfield(ts,'date') || ~isfield(ts,'data') )
    error('First arg TS must be a time series struct!');
  end;
  if ( ~exist('sper','var') || isempty(sper) )
    %sper = min(diff(ts.date));
    sper = median(diff(ts.date));
  end;
  if ( ~exist('newsper','var') || isempty(newsper) )
    newsper = (1/24);
  elseif ( strncmpi(newsper,'same',4) )
    newsper = sper;
  end;
  if ( ~exist('method','var') || isempty(method) )
    method = 'spline';
  end;
  if ( ~exist('maxgap','var') || isempty(maxgap) )
    maxgap = max([3*sper,3*newsper]);
  end;
  if ( ~exist('ringlen','var') || isempty(ringlen) )
    if ( strcmpi(method,'nearest') || strcmpi(method,'linear') )
      ringlen = [];
    else
      ringlen = sper;
    end;
  end;

  % Only interpolate on 'proper' values (and flatten complex numbers)
  goodix = find(isfinite(ts.data));
  if ( isempty(goodix) )
    error('No valid (finite) values found!');
  end;
  ts.date = ts.date(goodix);
  ts.data = real(ts.data(goodix));
  if ( isfield(ts,'prof') )
    ts.prof = real(ts.prof(goodix,:));
  end;

  % Time series usually contain Nx1 vectors
  newts.date = [ts.date(1):newsper:ts.date(end)]';
  newts.data = interp1(ts.date,ts.data,newts.date,method,varargin{:});
  if ( isfield(ts,'prof') )
    newts.prof = interp1(ts.date,ts.prof,newts.date,method,varargin{:});
  end;

  % If no extrapolation was requested, remove all overruns and gaps
  if ( isempty(varargin) )
    x.ts = ts;
    x.newts = newts;
    if ( ~isfield(ts,'prof') )
      x = filter_gaps(x,'ts','newts',maxgap,ringlen);
    else
      x = filter_gaps_prof(x,'ts','newts',maxgap,ringlen);
    end;
    newts = x.newts;
    x = []; clear x;
  end;

return;
