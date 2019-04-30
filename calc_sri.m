function vars = calculate_sri(stn, events)
%function vars = calculate_sri(stn, events)
%
% WORK IN PROGRESS
%
% Last Saved Time-stamp: <Wed 2008-08-13 08:24:54 Eastern Daylight Time gramer>
%

    % Now calculate a "simplified S/RI" based on number of HITS per day
    % (SUPER-periods from contiguous periods of same fuzzy maybe later)
    [Y M D h m s] = datevec(events.datenum(1:9));
    [ign begidx] = min(h);
    endidx = length(facts.datenum) - begidx + 1;
    endidx = endidx - mod(endidx, 8);
    ndays = floor((endidx - begidx - 1)/8.0);
    days = reshape(facts.datenum(begidx:endidx), [ndays 8]);
    histc(days, 0:length(fuzzies)+1, 2);

return;
