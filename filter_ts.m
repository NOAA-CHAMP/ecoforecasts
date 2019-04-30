function fts = filter_ts(ts, filterIndices, filterOrder)
%function fts = filter_ts(ts, filterIndices, filterOrder)
%
% Apply Butterworth band-pass filter to time series 'TS'.

    if ( ~exist('filterOrder', 'var') )
      filterOrder = 10;
      %filterOrder = 5;
    end

    if ( filterIndices == 0 )
        fts = ts;
    else
        HPcut = (2 / length(ts));
        LPcut = (2 / filterIndices);

        % Build Butterworth band-pass filter, with high cutoff frequency
        % (1/filterIndices), and low cutoff frequency 1/2 record length
%         [ butterB, butterA ] = ...
%                 butter(filterOrder, [HPcut LPcut]);
        [ butterB, butterA ] = ...
                butter(filterOrder, LPcut);
        fts = filtfilt(butterB, butterA, ts);
    end

return;
