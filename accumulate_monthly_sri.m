function [ uniqmos, sri ] = accumulate_monthly_sri(stn, events)
%function [ uniqmos, sri ] = accumulate_monthly_sri(stn, events)
%
% Return a monthly cumulative S/RI for *all events* in the struct 'events'
%
% Last Saved Time-stamp: <Tue 2008-08-19 16:02:57 Eastern Daylight Time gramer>
%

    [Y M D h m s] = datevec(events.jday);
    uniqmos = datenum(Y, M, 1, 0, 0, 0);
    sri = [];
    if ( ~isempty(uniqmos) && ~isempty(events.sri) )
      [xc, yc] = consolidator(uniqmos, events.sri, 'sum');
      uniqmos = xc';
      sri = yc';
    end;


return;
