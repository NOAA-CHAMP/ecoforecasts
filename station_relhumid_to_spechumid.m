function stn = station_relhumid_to_spechumid(stn,afld,qfld,sfld)
%function stn = station_relhumid_to_spechumid(stn,afld,qfld,sfld)
%
% Call RELHUMID_TO_SPECHUMID (v.) on data in the two fields STN.(AFLD) (air
% temperature) and STN.(QFLD) (Relative Humidity - assumed to be measured or
% calculated for the same elevation), to create a new field STN.(SFLD)
% containing Specific Humidity [kg/kg]. (NOTE: STN.SFLD is removed first!)
% QFLD may also be a constant value (e.g., 100 for saturated q).
%
% Last Saved Time-stamp: <Fri 2011-07-01 08:26:04  Lew.Gramer>

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

  if ( isfield(stn, sfld) )
    stn = rmfield(stn, sfld);
  end;
  stn.(sfld).date = stn.(afld).date(aix);
  stn.(sfld).data = relhumid_to_spechumid(stn.(afld).data(aix),q);

  stn.(sfld).date(isnan(stn.(sfld).data)) = [];
  stn.(sfld).data(isnan(stn.(sfld).data)) = [];

return;
