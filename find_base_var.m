function basevar = find_base_var(stn, varname)

  flds = fieldnames(stn);

  uscrs = strfind(varname, '_');

  for idx = length(uscrs):-1:1
    findx = find(strcmp(flds, varname(1:uscrs(idx)-1)));
    if ( ~isempty(findx) )
      basevar = flds{findx(1)};
      break;
    end;
  end;

return;
