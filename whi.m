function result = whi(lookforStr)
% whi - run 'lookfor -all' on the given string
%
% Author: Lew.Gramer@noaa.gov

  disp('ENV:');
  whos(lookforStr);
  disp('WHICH:');
  which(lookforStr, '-all');
  disp('LOOKFOR:');
  lookfor('-all', lookforStr);
return;
