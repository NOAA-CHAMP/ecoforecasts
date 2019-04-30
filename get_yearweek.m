function yrwks = get_yearweek(dts,inc)
%function yrwks = get_yearweek(dts,inc)
% Return the DATENUM of the start of week, for each DATENUM in DTS. Note this
% is NOT based on day of week (Mon-Sun); rather any date from, e.g., midnight
% Jan 1 to midnight Jan 7 of year YR will return 00:00 Jan 1, and any date
% from julian day 358 to New Years night will return DATENUM(YR,1,1)+358-1.
%
% If optional INC (DEFAULT: 0) is +1, return start of SUCCEEDING week; if
% -1, start of PRECEDING week; INC may be any integer scalar or vector; if a
% vector, then NUMEL(INC) must equal NUMEL(DTS).
% 
% Last Saved Time-stamp: <Sun 2017-06-04 16:09:54 Eastern Daylight Time gramer>

  if ( ~exist('inc','var') || isempty(inc) )
    inc=0;
  end;

  [wk,yr] = get_week(dts);

  % GET_WEEK always returns WK between 1 and 52
  wk = wk+inc;
  yrwks = datenum(yr,1,1) + ((wk-1)*7);

return;
