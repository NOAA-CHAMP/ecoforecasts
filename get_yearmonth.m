function yrmos = get_yearmonth(dts,inc)
%function yrmos = get_yearmonth(dts,inc)
% Return DATENUM of start of the calendar month, for each DATENUM in DTS.
%
% If optional INC (DEFAULT: 0) is +1, return start of SUCCEEDING month; if
% -1, start of PRECEDING month; INC may be any integer scalar or vector; if a
% vector, then NUMEL(INC) must equal NUMEL(DTS).
% 
% 
% Last Saved Time-stamp: <Sun 2017-06-04 16:10:11 Eastern Daylight Time gramer>

  if ( ~exist('inc','var') || isempty(inc) )
    inc=0;
  end;

  [yrs,mos,dys] = datevec(dts);
  yrmos = datenum(yrs,mos+inc,1);

return;
