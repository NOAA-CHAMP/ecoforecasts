function toomany = checknansforspline(V)
%function toomany = checknansforspline(V)
  toomany = false;

  [inan, jnan] = find(isnan(V));
  ncolnan = length(unique(jnan));
  nrownan = length(unique(inan));
  if ( ncolnan > 0 || nrownan > 0 )
    % Minimize loss of data. Strip rows instead of cols if there are less rows
    if ncolnan > nrownan
      toomany = checknansforspline_helper(V.');
    else
      toomany = checknansforspline_helper(V);
    end;
  end;
return;


function toomany = checknansforspline_helper(V)
%function toomany = checknansforspline_helper(V)
  numd = ndims(V);
  sizev = size(V);
  if numd > 2
    V = reshape(V,prod(sizev(1:(numd-1))),sizev(numd));
  end;
  nanv = find(sum(isnan(V),1));
  if ~isempty(nanv)
    V(:,nanv) = [];
  end;
  if numd > 2
    ncol = size(V,2);
    sizev(numd) = ncol;
    V = reshape(V,sizev);
  end;

  toomany = ( isempty(V) || isvector(V) );
return;
