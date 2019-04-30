function [newdts,newts] = window_func(dts,ts,func,nhrs,samplehrs,wintyp,plotResults)
%function [newdts,newts] = window_func(dts,ts,func,nhrs,samplehrs,wintyp,plotResults)
%
% Apply function FUNC to each NHRS hour-long subperiod of time series TS.
% Arg DTS must be a vector of monotonically increasing DATENUMs the same
% length as TS. One value is calculated each SAMPLEHRS hours. The arg WINTYP
% will allow low-pass filtering windows of various types to be applied during
% processing: currently this argument is IGNORED. (But see older M-funcs
% WINDOW_MEAN and WINDOW_STD, for some sample implementations of this.)
%
% Last Saved Time-stamp: <Wed 2012-07-04 15:50:46  lew.gramer>


    newts = [];
    newdts = [];

    % Assume a fixed fundamental sampling rate - just for sanity checks
    hrs_per_sample = min(diff(dts))*24;

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
    if ( hrs_per_sample <= 0 )
        error('Dates in "dts" must be monotonically increasing!');
    end;
    if ( (length(ts)*hrs_per_sample) < nhrs || length(ts) < samplehrs )
        error('Time series too short (%d elts., >=%g hrs./elt.) for your request!',...
              length(ts),hrs_per_sample);
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
      badix = find(~isfinite(tsrc));
      if ( ~isempty(badix) )
        %warning('LOWPASS: Had to replace NANs/INFs with mean value!');
        tsrc(badix) = nanmean(tsrc(:));
      end;

      [B, A] = butter(filterOrder, (2/nhrs), 'low');
      newts = filtfilt(B, A, tsrc);

      cutoffidx = ceil(nhrs/2);

      % NaN-ify any gaps - and any Gibbs oscillations from those gaps
      nanix = [];
      for ix = (-cutoffidx):cutoffidx
        nanix = unique( [ nanix(:) ; (badix(:)+ix) ] );
      end;
      nanix(1 > nanix | nanix > length(newts)) = [];
      newts(nanix) = nan;

      newdts = tdts((cutoffidx+1):samplehrs:(end-cutoffidx))';
      newts = newts((cutoffidx+1):samplehrs:(end-cutoffidx))';


     case {'highpass', 'hipass', 'hp'},
      badix = find(~isfinite(tsrc));
      if ( ~isempty(badix) )
        %warning('HIGHPASS: Had to replace NANs/INFs with mean value!');
        tsrc(badix) = nanmean(tsrc(:));
      end;

      [B, A] = butter(filterOrder, (2/nhrs), 'high');
      newts = filtfilt(B, A, tsrc);

      cutoffidx = ceil(nhrs/2);

      % NaN-ify any gaps - and any Gibbs oscillations from those gaps
      nanix = [];
      for ix = (-cutoffidx):cutoffidx
        nanix = unique( [ nanix(:) , (badix(:)+ix) ] );
      end;
      nanix(1 > nanix | nanix > length(newts)) = [];
      newts(nanix) = nan;

      newdts = tdts((cutoffidx+1):samplehrs:(end-cutoffidx))';
      newts = newts((cutoffidx+1):samplehrs:(end-cutoffidx))';

     case {'decimate','downsample','subsample'},
      badix = find(~isfinite(tsrc));
      if ( ~isempty(badix) )
        %warning('DECIMATE: Had to replace NANs/INFs with mean value!');
        tsrc(badix) = nanmean(tsrc(:));
      end;

      newts = decimate(tsrc, nhrs)';
      newdts = linspace(tdts(1),tdts(end),length(newts))';


     case {'interpolate','interp','spline'},
      badix = find(~isfinite(tsrc));
      if ( ~isempty(badix) )
        %warning('INTERPOLATE: Had to replace NANs/INFs with mean value!');
        tsrc(badix) = nanmean(tsrc(:));
      end;

      newdts = [tdts(1):(nhrs/24):tdts(end)]';
      if ( strcmpi('spline',func) )
        newts = spline(tdts,tsrc,newdts);
      else
        newts = interp1(tdts,tsrc,newdts);
      end;


     otherwise,

      % From ICON/G2:
      % operand is a symbol, has values average, count, sum, maximum, minimum, integral, derivative, deviation, current, as-of, or none, initially is none;
      switch ( lower(func) )
        case {'average', 'avg', 'mean'},
            func = 'mean';
        case {'stdev', 's.d.', 'deviation', 'std', 'dev'},
            func = 'std';
        case {'variance', 'variation', 'var'},
            func = 'var';
        case {'summation', 'sum'},
            func = 'sum';
        case {'maximum', 'max'},
            func = 'max';
        case {'minimum', 'min'},
            func = 'min';
        case {'integral'},
            func = 'trapz';
        case {'derivative', 'diff'},
            func = 'diff';
        case {'finitedifference', 'findiff'},
            func = 'findiff';
        case {'current', 'as-of', 'asof'},
            func = 'lastrow';
            %func = 'lastcol';
      end

      % If 'nhrs' is divisible by 'samplehrs', use a memory-intensive but FAST
      % MATLAB-optimized approach (Later perhaps: use only if it doesn't need
      % TOO much memory, i.e., number of array copies doesn't exhaust memory.)
      % if ( gcd(samplehrs, nhrs) == samplehrs && ...
      %      ((lcm(samplehrs, nhrs)/samplehrs)*length(tsrc)) < 1e8 )

      if ( gcd(samplehrs, nhrs) == samplehrs )

        newdts = repmat(nan, [floor(length(tdts)/samplehrs) 1]);
        newts = repmat(nan, [floor(length(tdts)/samplehrs) 1]);

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

            % SPECIAL CASE: Calculate a finite-difference derivative
            if ( strcmpi(func,'findiff') )
              switch ( round(nhrs) )
               case 3,
                val = (-newtsrc(1,:) + newtsrc(3,:)) / 2;
               case 5,
                val = (newtsrc(1,:) - 8*newtsrc(2,:) + ...
                       8*newtsrc(4,:) - newtsrc(5,:)) / 12;
               case 7,
                val = (-newtsrc(1,:) + 9*newtsrc(2,:) - 45*newtsrc(3,:) + ...
                       45*newtsrc(5,:) - 9*newtsrc(6,:) + newtsrc(7,:)) / 60;
               otherwise,
                error('Do not know how to finite difference over %d hours!', nhrs);
              end;

            % Take a row-wise mean/std/min/max/etc. - one value per column!
            else
              val = feval(func, newtsrc);
            end;

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
        dataname = 'data';
        if ( length(inputname(2)) > 0 )
            dataname = inputname(2);
        end;
        title([mfilename ': ' dataname ' vs. time']);
    end;

return;
