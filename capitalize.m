function res = capitalize(str)
%function res = capitalize(str)
%
% Capitalize string: First char, and any characters following ' ' or '_'.
% For non-vector character array, capitalize each row; for cell array,
% capitalize each cell element which is a character array. Non-strings are
% passed through unchanged (NOTE: this behavior differs from UPPER/LOWER.)
%
% Last Saved Time-stamp: <Tue 2016-04-12 12:21:14 Eastern Daylight Time gramer>

  res = str;

  if ( iscell(str) )
    for ix = 1:numel(str);
      res{ix} = capitalize(str{ix});
    end;
  elseif ( ischar(str) && ~isempty(str) )
    if ( ~isvector(str) )
      for ix = 1:size(str,1);
        res(ix,:) = capitalize(str(ix,:));
      end;
    else
      res(1) = upper(str(1));
      uscs = union(strfind(str,' '),strfind(str,'_'));
      for ix=uscs(:)';
        if ( length(str) < ix+1 )
          break;
        end;
        res(ix+1) = upper(str(ix+1));
      end;
    end;
  end;

return;
