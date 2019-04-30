function [ newdts, newts ] = window_func_v0(dts, ts, func, nhrs, samplehrs, plotResults)
%function [ newdts, newts ] = window_func_v0(dts, ts, func, nhrs, samplehrs, plotResults)
%
% Apply function 'func' to each 'nhrs' hour-long subperiod of time series
% 'ts'. Arg 'dts' must be a vector of monotonically increasing datenums of
% the same length as 'ts'. One value is calculated each 'samplehrs' hours.
%

    newts = [];
    newdts = [];

    % By default just return one value per nhrs-sized window
    if ( ~exist('samplehrs', 'var') || isempty(samplehrs) )
        samplehrs = nhrs;
    end;
    if ( ~exist('plotResults', 'var') || isempty(plotResults) )
        plotResults = 0;
    end;

    if ( isempty(ts) || all(isnan(ts)) )
        error('Empty or all-NaN time series "ts" specified!');
    end;
    if ( length(ts) ~= length(dts) )
        error('Date vector and time series must be same length!');
    end;
    if ( min(diff(dts)) <= 0 )
        error('Dates in "dts" must be monotonically increasing!');
    end;
    if ( length(ts) < nhrs || length(ts) < samplehrs )
        error('Time series too short (%d elts.) for your request!', length(ts));
    end;

    if ( nhrs < samplehrs )
        warning('WINDOW_FUNC:Aliasing', ...
                'Window size "nhrs" (%d) should be >= hours-per-sample (%d): Result may be severely aliased!', ...
                nhrs, samplehrs);
    end;


    % Make one value per hour - no matter how gappy our real time series is
    tdts = dts(1):(1.0/24.0):dts(end);

    % Start with a complete, regularly-spaced time series - and ignore NaNs
    tsrc = interp1(dts(~isnan(ts)), ts(~isnan(ts)), tdts);

    % How many windows of size 'nhrs' can our time series accommodate?
    nwnds = floor(length(tdts) / nhrs);


    switch ( lower(func) )
        case {'average', 'avg'},
            func = 'mean';
        case {'stdev', 's.d.', 'deviation'},
            func = 'std';
        case {'summation'},
            func = 'sum';
        case {'maximum'},
            func = 'max';
        case {'minimum'},
            func = 'min';
        case {'integral'},
            func = 'trapz';
        case {'derivative'},
            func = 'diff';
        case {'current', 'as-of'},
            func = 'subset';
    end

    % SPECIAL CASE: Much faster, if we just want one sample per window
    if ( nhrs == samplehrs )

        % Truncate time series based on even multiples of window-size
        tsrc = tsrc( (end - (nwnds * nhrs) + 1):end );

        % Rearrange our vector into nhrs-sized columns
        tsrc = reshape(tsrc, [nhrs nwnds]);

        newdts = tdts(nhrs:nhrs:end);
        newts = feval(func, tsrc);

    % For other sample frequencies, we have to do it the hard way
    else

        newidx = 0;
        for idx = nhrs:samplehrs:length(tsrc)
            newidx = newidx+1;
            smpl = tsrc((idx-nhrs+1):idx);
            if ( any(~isnan(smpl)) )
              newts(newidx) = feval(func, smpl(~isnan(smpl)));
              newdts(newidx) = tdts(idx);
            end;
        end;

    end;


    if ( plotResults )
        figure;
        plot(newdts, newts);
        datetick;
        if ( length(inputname(2)) > 0 )
            title([inputname(2) ' vs. time']);
        end;
    end;

return;
