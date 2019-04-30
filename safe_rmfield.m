function s = safe_rmfield(s,flds)
%function s = safe_rmfield(s,flds)
%
% Identical to RMFIELD (v.), except that when any field in FLDS does not
% exist in struct (or matrix of structs) S, no error is returned: instead, a
% WARNING with ID 'rmfield:NoField' (which may be toggled OFF) is displayed.
%
% Last Saved Time-stamp: <Fri 2012-03-23 12:26:41  Lew.Gramer>

  if ( ~iscellstr(flds) && ~ischar(flds) )
    error('Second arg must be a string array or cell array of strings!');
  end;

  if ( ischar(flds) )
    flds = cellstr(flds);
  end;

  all_flds = fieldnames(s);
  badix = (~ismember(flds,all_flds));

  s = rmfield(s,flds(~badix));
  ixes = find(badix)';
  if ( ~isempty(ixes) )
    for ix=ixes(:)';
      warning('rmfield:NoField','Missing field "%s"',flds{ix});
    end;
  end;

return;
