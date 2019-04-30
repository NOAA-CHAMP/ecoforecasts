function hrd_togo(tilwhat)
%function hrd_togo(tilwhat)
% Time left until transition to HRD

  if ( ~exist('tilwhat','var') || isempty(tilwhat) )
    tilwhat = 'trip';
  end;
  switch (lower(tilwhat)),
   case 'trip',		tgt=datenum(2019, 5,28,17, 0, 0); disp('Do not forget computer request!');
   case 'start',	tgt=datenum(2019, 6,17,10, 0, 0);
   otherwise,		error('Unknown event "%s"',tilwhat);
  end;
  x=tgt-now;
  disp(sprintf('%gd %gh',floor(x),floor((x-floor(x))*24)));
  wks = floor(x/7);
  dys = floor(x-(wks*7));
  hrs = floor((x-floor(x))*24);
  disp(sprintf('%gw %gd %gh',wks,dys,hrs));
return;
