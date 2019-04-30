function jday = ymd2jday(yr, mo, dy)
%function jday = ymd2jday(yr, mo, dy)
% Convert Year, Month, and Day to "Julian" day (year-day)

    jday = datenum(yr, mo, dy) - datenum(yr, 1, 1) + 1;

    if ( nargout < 1 )
      disp(jday);
    end;

return;
