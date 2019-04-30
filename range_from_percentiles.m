function ranges = range_from_percentiles(vals, pcts)
%FUNCTION ranges = range_from_percentiles(vals, pcts)
%
% Calculate seasonal ranges based on percentiles of the distribution.
%
% (NOTE: The MATLAB function PRCTILE is not available to us right now,
%  due to licensing... So approximate using the HIST function instead!)

    if ( nargin < 2 )
      %        ul  dl  vl  lo  sl  av  sh  hi  vh  dh   uh
      pcts = [ 00  01  03  15  30  70  85  97  99  100 ];
      %pcts = [ 00  02  07  15  30  70  85  93  98  100 ];
      %pcts = [ 00  05  15  25  30  70  75  85  95  100 ];
    end

    N = length(vals);
    nbins = min(100, N);
    [a h] = hist(vals, nbins);

    ranges(1) = -Inf;
    ranges(2) = min(vals);
    ranges(11) = max(vals);
    ranges(12) = +Inf;

    binsum = 0;
    curidx = 1;

    for idx = 2:length(pcts)-1

      pct = pcts(idx);
      for binidx = curidx:nbins
        binsum = binsum + a(binidx);
        if ( (binsum / N) >= (pct / 100) )
          ranges(idx+1) = h(binidx);
          break;
        end
      end
      curidx = binidx + 1;

    end

return;
