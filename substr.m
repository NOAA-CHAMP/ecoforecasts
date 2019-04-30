function str = substr(str,n)
%function str = substr(str,n)
  if ( isstring(str) || iscellstr(str) )
    for ix=1:numel(str)
      aStr = str{ix};
      str(ix) = aStr(1:n);
    end;
  elseif ( ischar(str) )
    str = str(1:n);
  else
    error('First argument has unknown class');
  end;
return;
