1;

stn1.station_name='aoat_adcp_1'; stn1.lon=-80.618633; stn1.lat=24.898016;
stn1 = get_ngdc_bathy_station(stn1);
stn1.ngdc_92m_depth = interpn(stn1.ngdc_92m_bathy.lat,stn1.ngdc_92m_bathy.lon,stn1.ngdc_92m_bathy.field,stn1.lat,stn1.lon);
stn1.ngdc_92m_depth,
%-4.1318

stn2.station_name='aoat_adcp_2'; stn2.lon=-80.618257; stn2.lat=24.897121;
stn2 = get_ngdc_bathy_station(stn2);
stn2.ngdc_92m_depth = interpn(stn2.ngdc_92m_bathy.lat,stn2.ngdc_92m_bathy.lon,stn2.ngdc_92m_bathy.field,stn2.lat,stn2.lon);
stn2.ngdc_92m_depth,
%-3.5491

stn3.station_name='aoat_adcp_3'; stn3.lon=-80.617453; stn3.lat=24.897930;
stn3 = get_ngdc_bathy_station(stn3);
stn3.ngdc_92m_depth = interpn(stn3.ngdc_92m_bathy.lat,stn3.ngdc_92m_bathy.lon,stn3.ngdc_92m_bathy.field,stn3.lat,stn3.lon);
stn3.ngdc_92m_depth,
%-3.6975

stnNW.station_name='NW'; stnNW.lat=24.898019; stnNW.lon=-80.618631;
% Divers reported 15' depth = 4.6m
stnNW.depth=4.6;
stnNW.ngdc_92m_depth = interpn(stn3.ngdc_92m_bathy.lat,stn3.ngdc_92m_bathy.lon,stn3.ngdc_92m_bathy.field,stnNW.lat,stnNW.lon);
stnNW.ngdc_92m_depth,
%

stnNE.station_name='NE'; stnNE.lat=24.897888; stnNE.lon=-80.616147;
% Diver (Fajans) reports 16' depth = 4.9m
stnNE.depth=4.9;
stnNE.ngdc_92m_depth = interpn(stn3.ngdc_92m_bathy.lat,stn3.ngdc_92m_bathy.lon,stn3.ngdc_92m_bathy.field,stnNE.lat,stnNE.lon);
stnNE.ngdc_92m_depth,
%

stnS.station_name='S'; stnS.lat=24.895926; stnS.lon=-80.617765;
% Divers reported 20' depth = 6.1m
stnS.depth=6.1;
stnS.ngdc_92m_depth = interpn(stn3.ngdc_92m_bathy.lat,stn3.ngdc_92m_bathy.lon,stn3.ngdc_92m_bathy.field,stnS.lat,stnS.lon);
stnS.ngdc_92m_depth,
%

