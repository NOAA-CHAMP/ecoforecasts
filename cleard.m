function cleard
% cleard - clear MATLAB 'directory stack' (see 'man csh', built-in pushd)
%
% Author: LGramer@upstream.net

  global myDirStack;

  myDirStack = {};

  % dirsfname = 'c:/Documents and Settings/lew.gramer/My Documents/MATLAB/savedirs.mat';
  %[pathroot, ig, ig, ig] = fileparts(mfilename('fullpath'));
  [pathroot, ig, ig] = fileparts(mfilename('fullpath'));
  dirsfname = fullfile(pathroot,'savedirs.mat');

  if ( exist(dirsfname, 'file') )
    delete(dirsfname);
  end;

return;
