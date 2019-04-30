function txt = interpmethod_to_text(method)
%function txt = interpmethod_to_text(method)
% Convert an interpolation METHOD specifier (e.g., the METHOD argument to the
% function INTERP_FIELD, v.) into a character string.

  txt = '';
  if ( ischar(method) )
    txt = method;
  elseif ( isa(method,'function_handle') )
    txt = char(method);
  elseif ( iscell(method) )
    txt = char(method{1});
    for ix=2:numel(method)
      if ( isnumeric(method{ix}) )
        txt = [txt,'_',num2str(method{ix})];
      else
        txt = [txt,'_',char(method{ix})];
      end;
    end;
  end;
  txt = upper(txt);

return;
