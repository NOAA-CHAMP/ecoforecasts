function newstn = unique_ts(stn)

    newstn = stn;

    repidx = find(diff(stn.datenum) == 0);
    hdrs = fieldnames(stn);
    for idx = 1:length(hdrs)
        newstn.(hdrs{idx})(repidx) = [];
    end;

return;
