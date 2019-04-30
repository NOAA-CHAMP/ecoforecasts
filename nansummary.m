function smry = nansummary(vals)
%function smry = nansummary(vals)
%
% Report MIN,25th,50th,75th PRCTILE,MAX, MEAN and STD of finite values in vector VALS

  vals = vals(isfinite(vals));
  vals = vals(:);

  if ( isempty(vals) )
    error('VALS has no valid finite values!');
  end;

  smry(1,1) = min(vals);
  smry(2,1) = prctile(vals,25);
  smry(3,1) = median(vals);
  smry(4,1) = prctile(vals,75);
  smry(5,1) = max(vals);
  smry(6,1) = mean(vals);
  smry(7,1) = std(vals);

  nm = inputname(1);
  if ( isempty(nm) )
    nm = 'Summary';
  end;

  if ( nargout < 1 )
    disp([nm ': ',num2str(smry')])
    clear smry;
  end;

return;
