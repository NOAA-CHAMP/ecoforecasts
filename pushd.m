function result = pushd(newdir)
% pushd - like C-Shell "pushd" builtin, but for MATLAB environment
%
% Author: LGramer@upstream.net

  global myDirStack;

  if ( nargout > 0 )
    result = [];
  end

  if ( ~exist('myDirStack', 'var') || isempty(myDirStack) )
    % myDirStack = {};
    loaddirs;
  end;

  if ( isempty(myDirStack) || ~strcmpi(myDirStack{1}, pwd) )
    myDirStack = { pwd myDirStack{:} };
  end;


  % Handle just "pd" - with no arguments
  if ( nargin < 1 )
    switch numel(myDirStack);
     case {0, 1},
      myDirStack = { pwd };
      savedirs;
      return;
     case 2,
      newdir = myDirStack{2};
      myDirStack = { newdir myDirStack{1} };
     otherwise,
      newdir = myDirStack{2};
      myDirStack = { newdir myDirStack{[1 3:end]} };
    end;

  % Handle, e.g., "pushd +2" (see Linux 'man csh')
  else
    % Convert "pushd +N" to be the same as "pushd(+N)"
    newdir_num = str2double(newdir);
    if ( ~isnumeric(newdir) && ~isnan(newdir_num) )
      newdir = newdir_num;
    end;

    if ( isnumeric(newdir) )
      if ( (newdir < 0) || ((newdir+1) > numel(myDirStack)) )
        savedirs;
        error('Not that many directories on stack! PWD not changed...');
      end;
      newdir = myDirStack{newdir+1};
    end;
  end

  olddir = pwd;

  if ( strcmpi(olddir, newdir) )
    oldix = find( strcmpi(myDirStack, olddir) );
    myDirStack(oldix) = [];
    myDirStack = { olddir myDirStack{1:end} };
    savedirs;
    return;
  end;


  try
    cd(newdir);
  catch
    warning('No such dir "%s": trying cdtomfile.', newdir);
    try
      cdtomfile(newdir);
    catch
      savedirs;
      error('No such dir OR m-file "%s"! CWD not changed...', newdir);
    end
  end

  newdir = pwd;

  oldix = find( strcmpi(myDirStack, olddir) );
  myDirStack(oldix) = [];
  newix = find( strcmpi(myDirStack, newdir) );
  myDirStack(newix) = [];

  myDirStack = { newdir olddir myDirStack{1:end} };
  savedirs;

  % disp( sprintf('%s\n', myDirStack{:}) );
  disp( sprintf('   \t\t %s', myDirStack{1}) );
  for ix = 2:numel(myDirStack)
    disp( sprintf('%2d \t\t %s', (ix-1), myDirStack{ix}) );
  end;

  if ( nargout > 0 )
    result = pwd;
  end

return;
