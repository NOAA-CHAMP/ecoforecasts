function result = derive_prop_ranges(vals)
%FUNCTION result = derive_prop_ranges(vals)
%
% Calculate appropriate value ranges for each of the CLIPS "fuzzy value"
% categories: unbelievably high, drastic high, very high, somewhat high,
% high, average, low, somewhat low, very low, drastic low, unbelievably low
% based on a statistical treatment of the vector of values passed in...
%
% Below is a sample from the "1-m sea temp." at CMRC station, for Winter
% 2004 through 2006, showing the layout of our 'result' return structure.
%
%       vals: [1918x1 double]
%       mean: 24.8507
%         sd: 1.1729
%       incs: [-Inf -3 -2.5000 -2 -1.5000 -1 1 1.5000 2 2.5000 3 Inf]
%     ranges: [1x12 double]
%         ul: [-Inf 21.3000]
%         dl: [21.3000 21.9000]
%         vl: [21.9000 22.5000]
%         lo: [22.5000 23.1000]
%         sl: [23.1000 23.7000]
%         av: [23.7000 26]
%         sh: [26 26.6000]
%         hi: [26.6000 27.2000]
%         vh: [27.2000 27.8000]
%         dh: [27.8000 28.4000]
%         uh: [28.4000 Inf]
%
% Last Saved Time-stamp: <Thu 2008-09-11 12:27:54 Eastern Daylight Time gramer>

    result.rangeNames = { ...
        'ul', 'dl', 'vl', 'lo', 'sl', ...
        'av', ...
        'sh', 'hi', 'vh', 'dh', 'uh', ...
    };

    result.vals = [ vals ];

    % Construct an array of eleven "fuzzy ranges", i.e., unbelievably-low,
    % drastic-low, very-low, low, somewhat-low, average, somewhat-high, high,
    % very-high, drastic-high, unbelievably-high... Calculate a vector of 12
    % value boundaries between and around each pair of our 11 "fuzzy" ranges.

    result.median = median(result.vals);
    result.mean = mean(result.vals);
    result.sd = std(result.vals);


    % Use multiples of a standard deviation to calculate seasonal ranges...
    %%%% ??? (TRY USING percentile ranges INSTEAD FOR NOW!)
    %%%% ??? result.ranges = range_from_stddevs(result.vals, result.mean, result.sd);

    % Use percentiles instead of standard-deviations to calculate ranges...
    result.ranges = range_from_percentiles(result.vals);


    % Round ranges to nearest 1/100th (of a oC, psu, millisiemens, whatever...)
    %%%%result.ranges = ( round( (result.ranges .* 100.0) ) ./ 100.0 );
    result.ranges = roundn( result.ranges, -2 );

    result.ul = [ result.ranges(1) result.ranges(2) ];
    result.dl = [ result.ranges(2) result.ranges(3) ];
    result.vl = [ result.ranges(3) result.ranges(4) ];
    result.lo = [ result.ranges(4) result.ranges(5) ];
    result.sl = [ result.ranges(5) result.ranges(6) ];
    result.av = [ result.ranges(6) result.ranges(7) ];
    result.sh = [ result.ranges(7) result.ranges(8) ];
    result.hi = [ result.ranges(8) result.ranges(9) ];
    result.vh = [ result.ranges(9) result.ranges(10) ];
    result.dh = [ result.ranges(10) result.ranges(11) ];
    result.uh = [ result.ranges(11) result.ranges(12) ];

return;
