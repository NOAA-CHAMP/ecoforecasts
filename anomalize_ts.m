function [anomts,clim,tid,asd,lopct,hipct,cumfun,per] = anomalize_ts(ts,varargin)
%function [anomts,clim,tid,asd,lopct,hipct,cumfun,per] = anomalize_ts(ts[,pct][,GRP_TS args])
%
% Calculate a year-day (or other) climatology CLIM at time values TID, for
% time series struct TS; calls GRP_TS (v.). Subtract the climatology from
% each value of TS to form new time series ANOMTS. Optionally returns vector
% of anomaly standard deviations ASD, size identical with CLIM and TID.
%
% Optional argument PCT is a percentage (0 < pct < 100) to use in returning
% lower and upper percentile ranges LOPCT and HIPCT for each period of the
% climatology. If PCT not specified, then 3% (the equivalent of 3 standard
% deviations for a Gaussian distribution) is assumed. NOTE: Estimate ASD is
% NOT validated by this function: if t.s. distribution is far from Gaussian
% normality, then LOPCT,UPPCT should be used for dispersion instead.
%
% Optional outputs CUMFUN and PER are passed through from GRP_TS (v.)
%
% SAMPLE CALL, producing a 'year-hour' (diurnal-annual) climatology:
%  [anomts,clim,tid,asd] = anomalize_ts(stn.ndbc_sea_t,@get_jhour_no_leap);
%
% Last Saved Time-stamp: <Fri 2017-02-03 13:14:37 Eastern Standard Time lew.gramer>

  args = varargin;

  % This is safe because first optional arg to GRP_TS is always non-numeric
  if ( numel(args) > 0 && isnumeric(args{1}) )
    pct = args{1};
    args(1) = [];
  else
    pct = 3;
  end;

  if ( numel(pct) == 1 )
    pctlo = pct;
    pcthi = 100 - pct;
  elseif ( numel(pct) == 2 )
    pctlo = pct(1);
    pcthi = pct(2);
  else
    error('If given, optional arg PCT must be scalar or 2-vector');
  end;
  if ( pctlo > pcthi )
    error('If optional PCT is a 2-vector, PCT(1) must be <= PCT(2)');
  end;
  if ( 0 >= pctlo || pcthi >= 100 )
    error('If given, optional arg(s) PCT must be between 0 and 100 exclusive');
  end;

  [clim,tid,nPerYr,dt,cumfun,per] = grp_ts(ts.data,ts.date,args{:});

  anomts = ts;
  [ig,ix] = ismember(cumfun(anomts.date),tid);
  anomts.data = anomts.data - clim(ix);

  if ( nargout > 3 )
    asd = grp_ts(anomts.data,anomts.date,cumfun,@nanstd);

    if ( nargout > 4 )
      lopct = grp_ts(anomts.data,anomts.date,cumfun,@(x)(prctile(x,pctlo)));
      hipct = grp_ts(anomts.data,anomts.date,cumfun,@(x)(prctile(x,pcthi)));
    end;
  end;

return;
