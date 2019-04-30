function set_more(onOff)
%function set_more(onOff)
%
% Save system "more" state, and change to onOff (one of 'on' or 'off'). If no
% argument is given, reset to previously saved "more" state. FOR EXAMPLE:
%
%  >> set_more off   % Disable 'more' so our script does not halt at each page!
%  >> script_with_copious_screen_output;
%  >> set_more   % Reset 'more' to previous state
%
% Last Saved Time-stamp: <Fri 2010-04-23 07:31:44 Eastern Daylight Time Lew.Gramer>

  persistent persistent_more_state;
  mlock; % Keep persistent variable around even during a CLEAR ALL: in case
         % of emergency, MUNLOCK('set_more_state') and CLEAR ALL will work

  % Save current MORE state (if *not* already saved), then change it
  if ( nargin == 1 )
    switch (lower(onOff)),
     case 'off', onOff = 'off';
     otherwise, onOff = 'on';
    end;
    % Note only save state on *first* call to SET_MORE in a calling stack...
    if (isempty(persistent_more_state))
      persistent_more_state = get(0,'more');
    end;
    more(onOff);

  % No argument - *reset* MORE state, and clear saved state
  else
    % If multiple functions in the calling stack call SET_MORE with no
    % argument, just ignore all of them after the first one...
    if (isempty(persistent_more_state))
      % Ignore
    else
      more(persistent_more_state);
      persistent_more_state = [];
    end;
  end;

return;
