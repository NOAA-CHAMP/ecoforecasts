function yrse = get_yearseason(dts,inc)
%function yrse = get_yearseason(dts,inc)
% Return DATENUM of start of the season (v. GET_SEASON), for each DATENUM in DTS.
%
% If optional INC (DEFAULT: 0) is +1, return start of SUCCEEDING season; if
% -1, start of PRECEDING season; INC may be any integer scalar or vector; if
% a vector, then NUMEL(INC) must equal NUMEL(DTS).
% 
% Last Saved Time-stamp: <Sun 2017-06-04 16:10:06 Eastern Daylight Time gramer>

  if ( ~exist('inc','var') || isempty(inc) )
    inc=0;
  end;

  [s,yrs] = get_season(dts);
  s = s+inc;
  yrse = datenum(yrs,((s-1)*3)+1,1);

return;
