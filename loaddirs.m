function loaddirs
%function loaddirs
% Load directory stack file global myDirStackFilename into global myDirStack
%
% Last Saved Time-stamp: <Tue 2016-08-16 17:24:29 Eastern Daylight Time lew.gramer>

  global myDirStackFilename;
  global myDirStack;

  % NEWmyDirStack = {};
  % if ( exist('myDirStack', 'var') && ~isempty(myDirStack) )
  %   NEWmyDirStack = myDirStack;
  %   myDirStack = {};
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

  if ( exist(myDirStackFilename, 'file') )
    %%%% DEBUG:    disp(['Loading ',myDirStackFilename]);
    load(myDirStackFilename, 'myDirStack');
  else
    %%%% DEBUG:    warning('No savedirs.mat file!');
    myDirStack = {};
  end;

  % myDirStack = { myDirStack{:} NEWmyDirStack{:} };

return;
