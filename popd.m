function result = popd(nPops)
% popd - like C-Shell 'popd' builtin, but for MATLAB environment
%
% Author: LGramer@upstream.net

  global myDirStack;

  if ( nargin < 1 )
    nPops = 0;
  end;

  if ( nargout > 0 )
    result = [];
  end;

  if ( ~exist('myDirStack', 'var') || isempty(myDirStack) )
    loaddirs;
  end;

  if ( isempty(myDirStack) || ~strcmpi(myDirStack{1}, pwd) )
    myDirStack = { pwd myDirStack{:} };
  end;

  % Convert "popd +N" to be the same as "popd(+N)"
  nPops_num = str2double(nPops);
  if ( ~isnumeric(nPops) && ~isnan(nPops_num) )
    nPops = nPops_num;
  end;
  if ( ~isnumeric(nPops) )
    error('USAGE: popd [#pops]. Optional arg #pops must evaluate numeric!');
  end;

  if ( nPops < 0 )
    error('USAGE: popd [#pops]. Optional arg #pops must be >= 0!');
  end;

  if ( (nPops+1) > numel(myDirStack) )
    error('Not enough directories on stack to pop %d!', nPops);
  end


  if ( nPops > 0 )
    myDirStack(nPops+1) = [];
  else % I.e., nPops == 0
    newdir = myDirStack{2};
    try
      cd(newdir);
      myDirStack(1) = [];
    catch
      warning('No such dir: "%s", PWD not changed, removing dir...', newdir);
      myDirStack(2) = [];
    end
  end;

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
