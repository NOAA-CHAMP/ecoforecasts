function stn = station_dewp_to_relhumid(stn, afld, dfld, qfld)
%function stn = station_dewp_to_relhumid(stn, afld, dfld, qfld)
%
% Call DEWP_TO_RELHUMID (v.) on data in two fields 'afld' (air temperature)
% and 'dfld' (dewpoint temperature - assumed to be measured at the same
% elevation), to create a new field 'qfld' containing Relative Humidity.
%
% Last Saved Time-stamp: <Sat 2009-07-04 11:14:18 Eastern Daylight Time gramer>

  ix = find(~isnan(stn.(dfld).data));

  [aix,dix] = intersect_dates(stn.(afld).date, stn.(dfld).date(ix));

  if ( isfield(stn, qfld) )
    stn = rmfield(stn, qfld);
  end;
  stn.(qfld).date = stn.(dfld).date(ix(dix));
  stn.(qfld).data = dewp_to_relhumid(stn.(afld).data(aix), stn.(dfld).data(ix(dix)));

  stn.(qfld).date(isnan(stn.(qfld).data)) = [];
  stn.(qfld).data(isnan(stn.(qfld).data)) = [];

return;
