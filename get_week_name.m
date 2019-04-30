function wknms = get_week_name(dts,longnames)
%function wknms = get_week_name(dts,longnames)
% 
% For each DATENUM in DTS, return short string for first date of "week",
% i.e., septad (always starts on Jan. 1). If LONGNAMES (DEFAULT: False), use
% monthname for current Locale, e.g., 'Jan01' through 'Dec22' for US-English.
% Otherwise, return month numbers. Calls DATESTR with 'mmdd' or 'mmmdd'.
% 
% Last Saved Time-stamp: <Wed 2018-08-29 15:43:50 Eastern Daylight Time gramer>

  if ( ~exist('longnames','var') || isempty(longnames) )
    longnames = false;
  end;

  if ( longnames )
    wknms = datestr(((get_week(dts)-1)*7)+1,'mmmdd');
  else
    wknms = datestr(((get_week(dts)-1)*7)+1,'mmdd');
  end;

return;
