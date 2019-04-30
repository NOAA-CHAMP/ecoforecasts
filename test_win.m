function test_win

    dv = datevec(now);
    % End test timestamp series as of midnight last night
    endtime = datenum(dv(1), dv(2), dv(3));

    % How much data should we run our tests with?
    x = 0:((10*8000)+1);
    d = (endtime-((1.0/24.0)*(length(x)-1))):(1.0/24.0):endtime;

    % How big is our moving window?
    for w = [3 24 72]

      % How many values per window in this test?
      n = w;
      run_test('baseline', d, x, w, n);
      n = 1;
      run_test('#1', d, x, w, n);
      n = 2;
      run_test('#2', d, x, w, n);
      n = 3;
      run_test('#3', d, x, w, n);
      n = 12;
      run_test('#4', d, x, w, n);
      n = 24;
      run_test('#5', d, x, w, n);

    end;

return;

function [dts, result] = run_test(testid, d, x, w, n)

    func = 'avg';

    olddts = [0]; newdts = [0];
    oldresult = [0]; newresult = [0];

    fprintf(1, 'Test %s: %d hr-window, every %d hrs:\n', testid, w, n);
    fprintf(1, '  OLD: ');
    warnstat = warning('off', 'WINDOW_FUNC:Aliasing');
    tic;
      [olddts, oldresult] = window_func_v0(d, x, func, w, n);
    toc;
    fprintf(1, '  NEW: ');
    tic;
      [newdts, newresult] = window_func(d, x, func, w, n);
    toc;
    warning(warnstat.state, 'WINDOW_FUNC:Aliasing');

    if ( ~isequal(olddts, newdts) || ~isequal(oldresult, newresult) )
        fprintf(1, 'Old/new results do not match?? len=%d,%d dts=?%d val=?%d\n', ...
                length(oldresult), length(newresult), ...
                isequal(olddts, newdts), isequal(oldresult, newresult));
        disp([ olddts(end), newdts(end) ; ...
               oldresult(end), newresult(end) ]);
    end;

return;
