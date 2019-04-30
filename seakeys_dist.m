function d = seakeys_dist(stnm1, stnm2)

error('This function has been replaced by STATION_DIST!');

  d = nan;

  SEAKEYS = get_seakeys;
  ix1 = find(strcmp(SEAKEYS.codes, upper(stnm1)));
  ix2 = find(strcmp(SEAKEYS.codes, upper(stnm2)));

  if ( isempty(ix1) || isempty(ix2) )
    return;
  end;

  lat1 = SEAKEYS.lats(ix1);
  lat2 = SEAKEYS.lats(ix2);
  lon1 = SEAKEYS.lons(ix1);
  lon2 = SEAKEYS.lons(ix2);

  d = sw_dist([lat1 lat2], [lon1 lon2], 'km');

return;
