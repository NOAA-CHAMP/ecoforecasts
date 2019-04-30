function str = textize(str)
%function str = textize(str): Return string suitable for use in TITLE, XLABEL, etc.
  try,
    str = strrep(str,'_','\_');
  catch,
    % For "character vector" weirdness in more recent MATLAB versions
    str = strrep(string(str),'_','\_');
  end;
return;
