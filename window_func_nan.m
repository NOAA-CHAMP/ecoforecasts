function [ newdts, newts ] = window_func_nan(dts, ts, func, nhrs, samplehrs, wintyp, plotResults)
%function [ newdts, newts ] = window_func_nan(dts, ts, func, nhrs, samplehrs, wintyp, plotResults)
%
% Apply function 'func' to each 'nhrs' hour-long subperiod of time series
% 'ts'. Arg 'dts' must be a vector of monotonically increasing datenums of
% the same length as 'ts'. One value is calculated each 'samplehrs' hours.
% The arg 'wintyp' will allow low-pass filtering windows of various types to
% be applied during processing: currently this argument is IGNORED. (But see
% older M-funcs window_mean and window_std for some sample implementations.)
%
% Last Saved Time-stamp: <Wed 2010-01-13 09:31:16 Eastern Standard Time Lew.Gramer>
%
%error('This function was NEVER COMPLETED!');

error('This function was NEVER COMPLETED!');


    newts = [];
    newdts = [];

    % By default just return one value per nhrs-sized window
    if ( ~exist('samplehrs', 'var') || isempty(samplehrs) )
        samplehrs = nhrs;
    end;
    if ( ~exist('wintyp', 'var') || isempty(wintyp) )
        wintyp = 'tophat';
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

    % If reading frame doesn't evenly divide sample period, aliasing likely
    if ( gcd(samplehrs, nhrs) ~= samplehrs )
        warning( 'WINDOW_FUNC:Aliasing', ...
                 'ALIASING: "nhrs" (%d) not divisible by "samplehrs" (%d)!', ...
                 nhrs, samplehrs);
    end;


    % HACK (for the current Data Parser and Campbell Scientific "brain" apps)
    ts(ts == -9.0 | ts == -999 | ts == -999.999) = nan;


    % Make one value per hour - no matter how gappy our real time series is
    tdts = dts(1):(1.0/24.0):dts(end);

    % Start with a complete, regularly-spaced time series - and ignore NaNs
    tsrc = interp1(dts(~isnan(ts)), ts(~isnan(ts)), tdts);

    filterOrder = 5;
    switch ( lower(func) )

      % Also support simple signal analysis operations (not supported in G2)
      % of low- and high-pass filtering, to allow intermediate t.s. analyses
     case {'lowpass', 'lopass', 'lp'},
      [B, A] = butter(filterOrder, (2/nhrs), 'low');
      newts = filtfilt(B, A, tsrc);
      cutoffidx = ceil(nhrs/2);
      newdts = tdts((cutoffidx+1):samplehrs:(end-cutoffidx));
      newts = newts((cutoffidx+1):samplehrs:(end-cutoffidx));

     case {'highpass', 'hipass', 'hp'},
      [B, A] = butter(filterOrder, (2/nhrs), 'high');
      val = filtfilt(B, A, tsrc);
      cutoffidx = ceil(nhrs/2);
      newdts = tdts((cutoffidx+1):samplehrs:(end-cutoffidx));
      newts = newts((cutoffidx+1):samplehrs:(end-cutoffidx));

     otherwise,

      % From ICON/G2:
      % operand is a symbol, has values average, count, sum, maximum, minimum, integral, derivative, deviation, current, as-of, or none, initially is none;
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
        case {'current', 'as-of', 'asof'},
            func = 'lastcol';
      end

      % If 'nhrs' is divisible by 'samplehrs', use a memory-intensive but FAST
      % MATLAB-optimized approach (Later perhaps: use only if it doesn't need
      % TOO much memory, i.e., number of array copies doesn't exhaust memory.)
      % if ( gcd(samplehrs, nhrs) == samplehrs && ...
      %      ((lcm(samplehrs, nhrs)/samplehrs)*length(tsrc)) < 1e8 )

      if ( gcd(samplehrs, nhrs) == samplehrs )

        newdts = repmat(nan, [1 floor(length(tdts)/samplehrs)]);
        newts = repmat(nan, [1 floor(length(tdts)/samplehrs)]);

        % "Reshape" a copy of time series, once for each sample-per-period.
        % The reason we build a vector, then "reshape" it into a rectangular
        % matrix is so that the highly optimized row-wise versions of Matlab
        % funs 'mean', 'std', 'min, 'max', 'trapz', 'diff', ... can be used.
        ncopies = lcm(samplehrs, nhrs) / samplehrs;
        for idx = 1:ncopies

            % End each copy based on a "samplehrs" reading frame: so the
            % first copy will end at the true end of the time series, the
            % second will end 'samplehrs' earlier, and so forth.
            endidx = length(tsrc) - ((idx-1)*samplehrs);
            % How many windows of size 'nhrs' can *this copy* accommodate?
            nwnds = floor(endidx / nhrs);
            % Start at the beginning of the earliest possible reading frame
            begidx = endidx - (nwnds * nhrs) + 1;

            newendidx = length(newts) - (idx-1);
            newnwnds = floor(newendidx / ncopies);
            newbegidx = newendidx - (newnwnds * ncopies) + 1;

            % Assign new time-stamps based on the END of each "nhrs" period
            newdts(newendidx:-ncopies:newbegidx) = tdts(endidx:-nhrs:begidx);

            % Rearrange vector into nwnds-sized columns, starting from END
            newtsrc = reshape(tsrc(begidx:endidx), [nhrs nwnds]);
            % Take a row-wise mean/std/min/max/etc. - one value per column!
            val = feval(func, newtsrc);
            newts(newendidx:-ncopies:newbegidx) = val(end:-1:1);

        end;
        newdts = newdts(~isnan(newts));
        newts = newts(~isnan(newts));

      % For any other window/sample-period combo, we have to do it the SLOW way
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
        newdts = newdts(~isnan(newts));
        newts = newts(~isnan(newts));

      end;

    end; % switch(lower(func)), otherwise

    if ( plotResults )
        figure;
        plot(newdts, newts);
        datetick;
        if ( length(inputname(2)) > 0 )
            title([inputname(2) ' vs. time']);
        end;
    end;

return;
