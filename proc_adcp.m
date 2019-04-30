1;

% stn.station_name='aoat_adcp_3'; stn.lon=-80.617453; stn.lat=24.897930;
% stn = get_ngdc_bathy_station(stn3;

stnNW.station_name='NW'; stnNW.lat=24.898019; stnNW.lon=-80.618631;
% Divers reported 15' depth = 4.6m
stnNW.depth=4.6;
% stnNW.ngdc_92m_depth = interpn(stn.ngdc_92m_bathy.lat,stn.ngdc_92m_bathy.lon,stn.ngdc_92m_bathy.field,stnNW.lat,stnNW.lon);
% stnNW.ngdc_92m_depth,

stnNE.station_name='NE'; stnNE.lat=24.897888; stnNE.lon=-80.616147;
% Diver (Fajans) reports 16' depth = 4.9m
stnNE.depth=4.9;
% stnNE.ngdc_92m_depth = interpn(stn.ngdc_92m_bathy.lat,stn.ngdc_92m_bathy.lon,stn.ngdc_92m_bathy.field,stnNE.lat,stnNE.lon);
% stnNE.ngdc_92m_depth,

stnS.station_name='S'; stnS.lat=24.895926; stnS.lon=-80.617765;
% Divers reported 20' depth = 6.1m
stnS.depth=6.1;
% stnS.ngdc_92m_depth = interpn(stn.ngdc_92m_bathy.lat,stn.ngdc_92m_bathy.lon,stn.ngdc_92m_bathy.field,stnS.lat,stnS.lon);
% stnS.ngdc_92m_depth,


% Adjust for atmospheric pressure at pressure-calibration times:
%  E3660 zeroed at 2012 Jan 09 16:53 Local (750KHz w/waves)
%   ICON rooftop barometric pressure: 1018.3
%  E3773 zeroed at 2012 Jan 09 18:40 Local (1.5MHz no waves)
%   ICON rooftop barometric pressure: 1018.2
%  E3772 zeroed at 2012 Jan 09 19:45 Local (1.5MHz w/waves)
%   ICON rooftop barometric pressure: 1018.3

% Adjust for local MAGNETIC VARIANCE at deployment site (should be ~6oW)
% 

