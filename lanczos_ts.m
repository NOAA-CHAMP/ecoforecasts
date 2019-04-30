function [newdts,newdat] = lanczos_ts(dts,dat,period)
%function [newdts,newdat] = lanczos_ts(dts,dat,[period])
%
% Apply Lanczos sine filter over PERIOD points (DEFAULT: 3) of the time
% series contained in vectors DTS (of DATENUMs) and DAT (of DOUBLE data).
% New time series NEWDTS,NEWDAT has all the gaps (and NaNs) of original data,
% and a CEIL(PERIOD/2) edge trimmed around each gap to remove Gibbs ringing.
%
% SEE: LANCZOSFILTER (MATLAB Exchange)
%
% Last Saved Time-stamp: <Tue 2010-11-16 14:03:54  Lew.Gramer>

  if ( ~exist('period','var') || isempty(period) )
    period = 3;
  end;
  Cf = 1/period;
  GibbsGap = ceil(period/2);

  newdts = [];
  newdat = [];

  mint = min(dts(:));
  gapix = [ 0 find(diff(dts) > (1.1.*mint)) length(dts) ];

  for ix = 2:length(gapix)
    ixes = (gapix(ix-1)+1):gapix(ix);
    nixes = length(ixes);
    if ( nixes > (GibbsGap*2) )
      newdts(end+1:end+nixes-(GibbsGap*2),1) = dts(ixes(1+GibbsGap:end-GibbsGap));
      y = lanczosfilter(dat(ixes),1,Cf);
      newdat(end+1:end+nixes-(GibbsGap*2),1) = y(1+GibbsGap:end-GibbsGap);
    end;
  end;

  % Filter out missing (NaN) values appropriately
  origbadix = find(~isfinite(dat(:)));
  badix = find(ismember(newdts,dts(origbadix)));
  for ix=1:GibbsGap
    badix = unique([badix(:)-1 ; badix(:) ; badix(:)+1]);
  end;
  badix(1 > badix | badix > length(newdts)) = [];
  % newdts(badix) = [];
  % newdat(badix) = [];
  newdat(badix) = nan;

return;
