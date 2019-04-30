function flds = grepstruct(strct,rgxp)
%function flds = grepstruct(strct,rgxp)
%
% Find all field names in STRCT matching regular expression RGXP. If no
% return arg FLDS specified, list of matching fields is output with DISP.
%
% Last Saved Time-stamp: <Wed 2011-03-30 09:08:02  Lew.Gramer>

  flds = {};

  allflds = fieldnames(strct);
  for fldix = 1:length(allflds)
    fld = allflds{fldix};
    if ( ~isempty(regexp(fld,rgxp)) )
      flds = { flds{:} fld };
    end;
  end;  

  % FIELDNAMES returns Nx1 cell array - do the same here
  flds = flds';

  if ( nargout < 1 )
    disp(flds(:));
    clear flds;
  end;

return;
