function stn = station_spechumid_to_relhumid(stn, afld, sfld, qfld)
%function stn = station_spechumid_to_relhumid(stn, afld, sfld, qfld)
%
% Call SPECHUMID_TO_RELHUMID (v.) on data in the two fields STN.(AFLD) (air
% temperature) and STN.(SFLD) (Specific Humidity [kg/kg] - assumed measured
% or calculated for the same elevation), to create a new field STN.(QFLD)
% containing Relative Humidity (0-100%). (NOTE: STN.QFLD is removed first!)
%
% Last Saved Time-stamp: <Tue 2010-03-23 13:01:29 Eastern Daylight Time Lew.Gramer>

  ix = find(~isnan(stn.(sfld).data));

  [aix,cix] = intersect_dates(stn.(afld).date, stn.(sfld).date(ix));

  if ( isfield(stn, qfld) )
    stn = rmfield(stn, qfld);
  end;
  stn.(qfld).date = stn.(sfld).date(ix(cix));
  stn.(qfld).data = spechumid_to_relhumid(stn.(afld).data(aix), stn.(sfld).data(ix(cix)));

  stn.(qfld).date(isnan(stn.(qfld).data)) = [];
  stn.(qfld).data(isnan(stn.(qfld).data)) = [];

return;
