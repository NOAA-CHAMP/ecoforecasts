function [yr, mo, dy] = jday2ymd(yr, jday)
%function [yr, mo, dy] = jday2ymd(yr, jday)
% Convert "Julian" day (year-day) to Year, Month, and Day

    dt = datenum(yr, 1, 1) + jday - 1;
    dvec = datevec(dt);
    yr = dvec(1);
    mo = dvec(2);
    dy = dvec(3);

    if ( nargout < 1 )
      disp(datestr(dt));
    end;

return;