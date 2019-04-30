function idx = find_pct_good_dates(dts,grpfun,pct,goodfun)
%function idx = find_pct_good_dates(dts,grpfun,pct[,goodfun])
%
% Return all indices of DATENUM vector DTS, for which there are at least PCT
% (0<PCT<1) values present in each time period: time periods are deterined by
% GRPFUN, as UNIQUE(GRPFUN(DTS)). Returns indices of all members of DTS which
% are in periods that have at least PCT good dates.
%
% NOTE: GRPFUN must be a function-handle that accepts and returns DATENUM.
%
% Last Saved Time-stamp: <Thu 2017-05-25 13:44:15 Eastern Daylight Time lew.gramer>
%
%error('UGH! This turns out to be very complex to code in a generic way... For now, see FIND_PCT_GOOD_PERIODS.m instead!');

error('UGH! This turns out to be very complex to code in a generic way... For now, see FIND_PCT_GOOD_PERIODS.m instead!');

  if ( ~exist('dts','var') || ~isnumeric(dts) || ~isvector(dts) )
    error('DTS must be a vector of DATENUM');
  end;
  if ( ~exist('grpfun','var') || ~isa(grpfun,'function_handle') )
    error('GRPFUN must be a FUNCTION_HANDLE accepting AND returning DATENUM vectors');
  end;
  if ( ~exist('pct','var') || 0 > pct || pct > 1 )
    error('PCT must be decimal between 0 and 1');
  end;

  d = grpfun(dts);
  ud = unique(d);

  % Median time series resolution 
  dt = median(diff(dts));
  % Maximum period resolution 
  DT = max(diff(ud));
  ud(end+1) = grpfun(ud(end) + (DT*1.001));

  alldts = ud(1):dt:ud(end);
  if ( dts(1) < alldts(1) || alldts(end) < dts(end) )
    error('GRPFUN returned values that do not encompass our dates!');
  end;
  alld = grpfun(alldts);
  ualld = unique(alld);

  % NOTE: The last period is unreliable - always include it for now
  idx = 1:numel(d);
  for uix = 1:numel(ud)-1
    dix = find(ismember(d,ud(uix)));
    allix = find(ismember(alld,ud(uix)));
    if ( (numel(dix)/numel(allix)) >= pct )
      idx(end+1:end+numel(dix)) = dix;
    end;
  end;
  dix = find(ismember(d,ud(end)));
  idx(end+1:end+numel(dix)) = dix;

return;
