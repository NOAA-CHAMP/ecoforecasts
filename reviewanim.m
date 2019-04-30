function reviewanim(fhs,pausesecs,restartsecs,waitsecs)
%function reviewanim(fhs,pausesecs,restartsecs,waitsecs)
%
% "Review" all figures in vector of figure handles FHS: bring each figure to
% focus in turn, with a PAUSESECS (DEFAULT 0.5) sec delay between focus
% changes to each successive figure in FHS. Describes user interface with a
% msg dialog box for WAITSECS (DEFAULT 2) sec before beginning review.
%
% NOTE: *REPEATS INDEFINITELY* until the user hits Ctrl-C or an error occurs.
% Before each complete review, pauses RESTARTSECS (DEFAULT 2*PAUSESECS) sec.
% However, if PAUSESECS is 0, will just review figures *ONE TIME*.
%
% Last Saved Time-stamp: <Sat 2016-11-05 16:42:50 Eastern Daylight Time gramer>

  if ( ~exist('fhs', 'var') )
    fhs = [];
  end;
  if ( ~exist('pausesecs', 'var') || isempty(pausesecs) )
    pausesecs = 0.5;
  end;
  if ( ~exist('restartsecs', 'var') || isempty(restartsecs) )
    restartsecs = 2*pausesecs;
  end;
  if ( ~exist('waitsecs', 'var') || isempty(waitsecs) )
    waitsecs = 2;
  end;

  if ( isempty(fhs) || all(~ishandle(fhs(:))) )
    warning('Figure handle list was empty! Reviewing ALL figures...');
    fhs = get(0, 'Children');
    if ( verLessThan('matlab','7.5') )
      [ig,ix] = sort(fhs);
    else
      fnum = get(fhs,'Number');
      [ig,ix] = sort([fnum{:}]);
    end;
    fhs = fhs(ix);
  end;
  %if ( isempty(fhs) || all(~isfinite(fhs(:))) )
  if ( isempty(fhs) || all(~ishandle(fhs(:))) )
    error('No valid figures found!');
  end;

  if ( ~isfinite(pausesecs) )
    msg = sprintf('Hit enter to change from current figure to the next\n');
  elseif ( pausesecs > 0 )
    msg = sprintf('Review will pause %g seconds on each figure\n',pausesecs);
  else
    msg = sprintf('Review will not pause between figures\n');
  end;
  if ( ~isfinite(waitsecs) )
    msg = [msg,'Hit enter to begin review: Ctrl-C to stop...'];
    disp(msg);
    msgh = msgbox(msg,upper(mfilename),'non-modal');
    pause;
    close(msgh);
  elseif ( waitsecs > 0 )
    msg = [msg,'Review will begin in ' num2str(waitsecs) ' secs: Ctrl-C to stop...'];
    disp(msg);
    msgh = msgbox(msg,upper(mfilename),'non-modal');
    pause(waitsecs);
    close(msgh);
  end;

  try
    while (1)
      %validfhs = numel(find(isfinite(fhs)));
      validfhs = numel(find(ishandle(fhs)));
      for ix = 1:numel(fhs)
        %if ( isfinite(fhs(ix)) )
        if ( ishandle(fhs(ix)) )
          if ( ~ishandle(fhs(ix)) )
            warning('Handle #%d (%f) is not a valid handle!', ix, fhs(ix));
            validfhs = validfhs - 1;
          else
            figure(fhs(ix));
            if ( ~isfinite(pausesecs) )
              pause;
            else
              pause(pausesecs);
            end;
          end;
        end;
      end;
      % If PAUSESECS==0, just run through figures ONE TIME
      if ( pausesecs == 0 )
        break;
      end;
      if ( ~isfinite(restartsecs) )
        pause;
      else
        pause(restartsecs);
      end;
      if ( ~validfhs )
        disp('No more valid figure handles! Exiting....');
        break;
      end;
    end;
  catch
    disp([upper(mfilename) ': Signal caught! Exiting...']);
  end;

return;
