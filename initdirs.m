function initdirs(aDirStackFilename)
% initdirs - initialize C-Shell-like directory stack
%
% Last Saved Time-stamp: <Tue 2016-08-16 17:11:37 Eastern Daylight Time lew.gramer>

  global myDirStackFilename;
  global myDirStack;

  if ( ~exist('aDirStackFilename','var') )
    aDirStackFilename = [];
  end;

  if ( ~isempty(aDirStackFilename) && isempty(regexp(aDirStackFilename,'.mat$')) )
    % Maybe they just passed us a pathname? Maybe not!
    if ( ~exist(aDirStackFilename,'dir') )
      warning('Invalid first argument: using default directory stack file path');
      aDirStackFilename = [];
    end;
  end;

  if ( isempty(aDirStackFilename) )
    % If caller did not pass a directory, default to saving the directory
    % stack in a .MAT file in the directory immediately above where we keep
    % this code, e.g., in the MATLAB "home" directory...

    % myDirStackFilename = 'c:/Documents and Settings/lew.gramer/My Documents/MATLAB/savedirs.mat';
    %[pathroot, ig, ig, ig] = fileparts(mfilename('fullpath'));
    [pathroot, ig, ig] = fileparts(mfilename('fullpath'));
    [pathroot, ig, ig] = fileparts(pathroot);
    myDirStackFilename = fullfile(pathroot,'savedirs.mat');
    %%%% DEBUG:    disp(['Defaulted dir stack file to ',myDirStackFilename]);
  elseif ( ~isempty(regexp(aDirStackFilename,'.mat$')) )
    myDirStackFilename = aDirStackFilename;
    %%%% DEBUG:    disp(['Dir stack file set to ',myDirStackFilename]);
  elseif ( exist(aDirStackFilename,'dir') )
    % They just passed us a pathname!
    myDirStackFilename = fullfile(aDirStackFilename,'savedirs.mat');
    %%%% DEBUG:    disp(['Dir stack file set to ',myDirStackFilename]);
  else
    error('Logic error: Pass in a valid directory stack path, filename, or nothing!');
  end;
  clear aDirStackFilename

  if ( ~exist('myDirStack', 'var') || isempty(myDirStack) )
    % myDirStack = {};
    loaddirs;
  end

  if ( ~isempty(myDirStack) )
    newdir = myDirStack{1};
    try
      cd(newdir);
    catch
      warning('No such dir "%s": trying cdtomfile.', newdir);
      try
        cdtomfile(newdir);
      catch
        try
          % Last ditch effort - try the next dir on the stack, if any...
          popd;
        catch
          try
            % Maybe third time is the charm? If not, give up...
            popd;
          catch
          end
        end
      end
    end
  end;

  savedirs;

return;
