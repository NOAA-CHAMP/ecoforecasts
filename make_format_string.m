function fmt = make_format_string(cellfmt,delim,eol)
%function fmt = make_format_string(cellfmt,delim,eol)
%
% Make a format string, e.g., for TEXTSCAN (v.), using the cell array format
% specifier CELLFMT. CELLFMT consists of a cell array of 2-tuples, first a
% counter, then a format specifier, e.g., {2,'%*[^,]',2,'%[^,]'} makes a
% format string to ignore the first two comma-delimited fields and return the
% next two. Optional DELIM (DEFAULT: ',') is the separator between format
% fields in the final string, and EOL (DEFAULT: '\n') the terminator.
%
% Last Saved Time-stamp: <Wed 2016-08-17 12:28:36 Eastern Daylight Time gramer>

  if ( ~exist('delim','var') || ~ischar(delim) )
    delim = ',';
  end;
  if ( ~exist('eol','var') || ~ischar(eol) )
    eol = '\n';
  end;

  fmt = '';
  curdelim = '';
  for ix = 1:2:numel(cellfmt)
    str = cellfmt{ix+1};
    if ( isempty(str) )
      str = '%*[^,]';
    elseif ( strcmp(str,'%s') )
      str = '%[^,]';
    end;
    for nix = 1:cellfmt{ix}
      fmt = [fmt,curdelim,str];
      curdelim = delim;
    end;
  end;
  fmt = [fmt,eol];

return;
