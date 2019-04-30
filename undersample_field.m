function newfld = undersample_field(fld,kernel,method)
%function newfld = undersample_field(fld,kernel,method)
% Downsample an NxM or PxNxM  matrix FLD, using KERNEL (DEFAULT: [5,5]) and
% METHOD (DEFAULT: @nanmean): Accumulates METHOD on each kernel(2)xkernel(1)
% rectangle of FLD: result NEWFLD is a matrix whose last two dimensions are
% FLOOR(N/kernel(2)) and FLOOR(M/kernel(1)), resp.

  if ( ~exist('kernel','var') || isempty(kernel) )
    kr = 5;
    kc = 5;
  elseif ( isnumeric(kernel) && numel(kernel) == 2 )
    kr = kernel(2);
    kc = kernel(1);
  elseif ( isnumeric(kernel) && numel(kernel) == 1 )
    kr = kernel;
    kc = kernel;
  else
    error('KERNEL must be empty (DEFAULT), or a numeric 1 or 2-vector');
  end;
  if ( kr == 1 &&  kc == 1 )
    newfld = fld;
    %%%% EARLY RETURN
    return;
  end;
  if ( ~exist('method','var') || isempty(method) )
    method = @nanmean;
  elseif ( ischar(method) )
    method = str2func(method);
  end;

  if ( ndims(fld) > 2 )
    error('3D matrices not implemented yet!');
  end;
  [n,m] = size(fld);
  newn = floor(n/kr);
  newm = floor(m/kc);
  newfld = repmat(nan,[newn,newm]);

  if ( newn*newm > 1e6 )
    warning('This will be slow!');
  end;
  for rowix=kr:kr:n;
    for colix=kc:kc:m;
      tmpl = fld(rowix-kr+1:rowix,colix-kc+1:colix);
      newfld(rowix/kr,colix/kc) = method(tmpl(:));
    end;
%     tmpl = fld(rowix-kr+1:rowix,:);
% keyboard;
%     newfld(rowix/kr,:) = method(tmpl(:));
  end;

return;
