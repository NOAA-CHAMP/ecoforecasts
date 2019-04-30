function stn = station_tides_to_idepths(stn)
%function stn = station_tides_to_idepths(stn)
%
% Calculate site-depth time series from tides ('tide' and/or 'ndbc_tide').
% Also adds new "_m" tide fields in [m] instead of [FEET].
%
% Last Saved Time-stamp: <Tue 2010-05-25 11:13:12 Eastern Daylight Time lew.gramer>

  try
    [ig,ig,stz] = get_station_coords(stn.station_name);
  catch
    stz = 2;
    warning('Using default mean site depth %g',stz);
  end;

  %DEBUG:  disp(stz);

  tfld = 'tide';
  dfld = 'tide_m';
  ifld = 'tide_i_depth';
  if ( isfield(stn,tfld) )
    stn.(dfld) = stn.(tfld);
    stn.(dfld).data = stn.(dfld).data .* 0.3048;

    goodix = find(5 < stn.(tfld).data & stn.(tfld).data < 15);
    stn.(ifld).date = stn.(tfld).date(goodix);
    stn.(ifld).data = stn.(tfld).data(goodix) - 10;
    % Convert from feet to [m]
    stn.(ifld).data = stn.(ifld).data .* 0.3048;
    stn.(ifld).data = stn.(ifld).data + stz;
  end;

  tfld = 'ndbc_tide';
  dfld = 'ndbc_tide_m';
  ifld = 'ndbc_tide_i_depth';
  if ( isfield(stn,tfld) )
    stn.(dfld) = stn.(tfld);
    stn.(dfld).data = stn.(dfld).data .* 0.3048;

    goodix = find(-5 < stn.(tfld).data & stn.(tfld).data < 5);
    stn.(ifld).date = stn.(tfld).date(goodix);
    stn.(ifld).data = stn.(tfld).data(goodix);
    % Convert from feet to [m]
    stn.(ifld).data = stn.(ifld).data .* 0.3048;
    stn.(ifld).data = stn.(ifld).data + stz;
  end;

return;
