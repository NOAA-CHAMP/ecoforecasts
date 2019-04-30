function dirs
% dirs - like C-Shell 'dirs' builtin, but for MATLAB environment
%
% Author: LGramer@upstream.net

  global myDirStack;

  if ( ~exist('myDirStack', 'var') || isempty(myDirStack) )
    % myDirStack = {};
    loaddirs;
  end

  if ( isempty(myDirStack) || ~strcmpi(myDirStack{1}, pwd) )
    myDirStack = { pwd myDirStack{:} };
  end;

  % disp( sprintf('%s\n', myDirStack{:}) );
  disp( sprintf('   \t\t %s', myDirStack{1}) );
  for ix = 2:numel(myDirStack)
    disp( sprintf('%2d \t\t %s', (ix-1), myDirStack{ix}) );
  end;

return;
