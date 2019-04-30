1;
%% SCRIPT alexanders_mit_problem.m
%% Flip N light switches in sequence: when the IX'th switch is flipped,
%% toggle all lights which are muliples of IX. At the end, how many are lit?

NN = [];
clear NN

% How many lights (N) are in the trial: try a sequence of trials
for N = 1:1000
%for N = 121
  % Pre-allocate the number of binary digits
  light = repmat(0,[1,N]);

  %tic,  % Start timing a trial
  for ix=1:N
    for ixix=ix:N;
      if ( isint(ixix/ix) )
        light(ixix) = not(light(ixix));
      end;
    end;
    %DEBUG:    disp(light);
  end;
  %toc,  % Report time consumed by a trial

  NN(N) = numel(find(light));
  %DEBUG:    disp(light);
  %DEBUG:    disp([N,NN(N)]);
end;

% How long (how many successive N's) between new values?
[rl,val]=run_length_encode(NN),
