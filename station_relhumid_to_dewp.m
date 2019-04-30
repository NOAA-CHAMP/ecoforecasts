function stn = station_relhumid_to_dewp(stn,afld,qfld,dfld)
%function stn = station_relhumid_to_dewp(stn,afld,qfld,dfld)
%
% Call RELHUMID_TO_DEWP (v.) on data in fields name AFLD (air temperature)
% and QFLD (Relative Humidity 0-100%, assumed to be measured at the same
% elevation), to create a new field DFLD containing dewpoint temperature.
% QFLD may also be a constant value (e.g., 100 for saturated dew temp).
%
% Last Saved Time-stamp: <Fri 2011-07-01 08:29:41  Lew.Gramer>

  if ( ischar(qfld) )
    ix = find(~isnan(stn.(qfld).data));
    [aix,qix] = intersect_dates(stn.(afld).date, stn.(qfld).date(ix));
    q = stn.(qfld).data(ix(qix));
  elseif ( isnumeric(qfld) && isscalar(qfld) )
    aix = 1:length(stn.(afld).date);
    q = repmat(qfld,size(stn.(afld).data));
  else
    error('QFLD must either be a fieldname string or a numeric scalar (constant)!');
  end;

  if ( isfield(stn, dfld) )
    stn = rmfield(stn, dfld);
  end;
  stn.(dfld).date = stn.(afld).date(aix);
  stn.(dfld).data = relhumid_to_dewp(stn.(afld).data(aix),q);

  stn.(dfld).date(isnan(stn.(dfld).data)) = [];
  stn.(dfld).data(isnan(stn.(dfld).data)) = [];

return;
