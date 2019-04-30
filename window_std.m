function sd = window_std(ts, N, wintyp, realizable)
%function sd = window_std(ts, N, wintyp, realizable)
%
% Calculate a "running standard deviation" over each 'N' point sample of time
% series 'ts', using moving window of type 'wintyp' (DEFAULT tophat/rectangle).
% If arg 'realizable' is present and evaluates to true, the running window is
% applied from the current point BACKWARD, rather than CENTERED on that point.
%
% NOTE: The resulting vector 'sd' has the SAME LENGTH AS the input vector 'ts'.
% (A partial-sample s.d. is calculated for each front and back tail slot.)
%
% EXAMPLES: Specifying a wintyp of 'rect' does a Centered Subset s.d.
% Specifying wintyp of 'rect' with realizable==true instead calculates a
% Simple Moving s.d., i.e., the s.d. at index 'j' uses a rectangular
% window that ENDS at 'j' - thus matching the algorithm used in ICON/G2. A
% wintyp of 'default' would approximate an N-point low-pass filtered s.d.
%
% Last Saved Time-stamp: <Wed 2008-08-13 07:56:11 Eastern Daylight Time gramer>
%

    if ( ~exist('wintyp', 'var') )
      wintyp = 'tophat';
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


    sd = repmat(nan, size(ts));

    % For tophat, n == N. But for some windows, n << N.
    n = sum(win);

    % (N - 1)^(-0.5) is the normalizer for normal (tophat) s.d.
    % This adjusts result so s.d. is normalized with (n-1)^(-0.5).
    adj = sqrt(N - 1) / sqrt(n - 1);

    % Do a "running" mean - using the N points previous to each data point
    if ( realizable )
      for idx = N:length(sd)
        sd(idx) = adj * std( win .* ts((idx-N+1):idx) );
      end
      % Fake front tail of the time series to avoid NaNs
      for idx = 1:(N-1)
        sd(idx) = adj * std( win .* ts(1:N) );
      end

    % Centered mean - symmetric window with "bulge" at each data point
    else
      halfN = round(N/2);
      for idx = N:(length(sd)-N)
        sd(idx) = adj * std( win .* ts((idx-halfN):(idx+halfN-1)) );
      end

      % Fake front and back tails of the time series to avoid NaNs
      for idx = 1:N
        sd(idx) = adj * std( win .* ts(idx:(idx+N-1)) );
      end
      for idx = (length(sd)-N+1):length(sd)
        sd(idx) = adj * std( win .* ts((idx-N+1):idx) );
      end
    end

return;
