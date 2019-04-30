function m = window_mean(ts, N, wintyp, realizable)
%function m = window_mean(ts, N, wintyp, realizable)
%
% Calculate a "running mean" (or Moving Average) over each 'N' points of time
% series 'ts', using moving window of type 'wintyp' (DEFAULT blackman-harris).
% If arg 'realizable' is present and evaluates to true, the running window is
% applied from the current point BACKWARD, rather than CENTERED on that point.
%
% NOTE: The resulting vector 'm' has the SAME LENGTH AS the input vector 'ts'.
% (A partial-sample mean is calculated for each front and back tail slot.)
%
% EXAMPLES: Specifying a wintyp of 'rect' does a Centered Moving Average.
% Specifying wintyp of 'rect' with realizable==true instead calculates a
% Simple Moving Average, i.e., the mean at index 'j' uses a rectangular
% window that ENDS at 'j' - thus matching the algorithm used in ICON/G2. A
% wintyp of 'default' would approximate an N-point low-pass filtered mean.
%
% Last Saved Time-stamp: <Wed 2008-08-13 07:55:59 Eastern Daylight Time gramer>
%

    if ( ~exist('wintyp', 'var') )
      wintyp = 'bh';
    end
    if ( ~exist('realizable', 'var') )
      realizable = false;
    end

    if ( realizable )
      N = N*2;
    end

    switch (wintyp)
     case {'tophat','rectangle','rect'}
      % True rectangular window (i.e., a running average)
      win = window(@rectwin, N);
     case {'bh','blackman-harris','default'}
      % Blackman-Harris window: reduces frequency transients.
      win = window(@blackmanharris, N);
     otherwise,
      % SEE matlab command "help window" for other options
      win = eval(['window(@' wintyp ', N)']);
    end

    if ( realizable )
      N = N/2;
      win = win(1:N);
    end


    m = repmat(nan, size(ts));

    % For tophat, n == N. But for some windows, n << N.
    n = sum(win);

    % Do a "running" mean - using the N points previous to each data point
    if ( realizable )
      for idx = N:length(m)
        m(idx) = sum( win .* ts((idx-N+1):idx) ) / n;
      end
      % Fake front tail of the time series to avoid NaNs
      for idx = 1:(N-1)
        m(idx) = ts(idx);
      end

    % Centered mean - symmetric window with "bulge" at each data point
    else
      halfN = round(N/2);
      for idx = N:(length(m)-N)
        m(idx) = sum( win .* ts((idx-halfN):(idx+halfN-1)) ) / n;
      end

      % Fake front and back tails of the time series to avoid NaNs
      for idx = 1:N
        m(idx) = sum( win .* ts(idx:(idx+N-1)) ) / n;
      end
      for idx = (length(m)-N+1):length(m)
        m(idx) = sum( win .* ts((idx-N+1):idx) ) / n;
      end
    end

return;
