function nums = month_str2num(strs)
%function nums = month_str2num(strs)
% 
% Return a vector containing the month number (1-12) for each month name in
% the character array or cellstr STRS. This is generally faster than, e.g.,
% GET_MONTH(DATENUM(STRS,'mmm')) (v.) especially for long arrays.
%
% Last Saved Time-stamp: <Fri 2014-10-24 13:26:25 Eastern Daylight Time gramer>

  if ( iscellstr(strs) )
    nstrs = numel(strs);
  elseif ( ischar(strs) )
    nstrs = size(strs,1);
  else
    error('First arg must be a char array or CELLSTR');
  end;

  nums = repmat(nan,[nstrs,1]);
  nms = get_monthnames;
  for mo=1:12
    ix = strmatch(nms(mo,:),strs);
    nums(ix) = mo;
  end;

return;
