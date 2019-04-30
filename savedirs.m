function savedirs
% Save directory stack from global myDirStack into file myDirStackFilename
%
% Last Saved Time-stamp: <Tue 2016-08-16 17:12:35 Eastern Daylight Time lew.gramer>

  global myDirStackFilename;
  global myDirStack;

  % % If there's nothing in there, load our last save first...
  % if ( ~exist('myDirStack', 'var') || isempty(myDirStack) )
  %   loaddirs;
  % end;

  % This may be set in the call to INITDIRS (v.), but if not, default to
  % saving the directory stack in a .MAT file in the directory immediately
  % above where we keep this code, e.g., in the MATLAB "home" directory...
  if ( ~exist('myDirStackFilename','var') || isempty(myDirStackFilename) )
    % myDirStackFilename = 'c:/Documents and Settings/lew.gramer/My Documents/MATLAB/savedirs.mat';
    %[pathroot, ig, ig, ig] = fileparts(mfilename('fullpath'));
    [pathroot, ig, ig] = fileparts(mfilename('fullpath'));
    [pathroot, ig, ig] = fileparts(pathroot);
    myDirStackFilename = fullfile(pathroot,'savedirs.mat');
    %%%% DEBUG:    disp(['Dir stack file defaulted to ',myDirStackFilename]);
  end;

  %save(myDirStackFilename, 'myDirStack');
  %save(myDirStackFilename, 'myDirStack','-v7.3');
  save(myDirStackFilename, 'myDirStack','-v7');

return;
