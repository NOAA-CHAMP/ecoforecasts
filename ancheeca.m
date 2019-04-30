1;

stn.station_name = 'cheeca'; stn.lon=-80.61475; stn.lat=24.904083;
stn = get_ngdc_bathy_station(stn);


rix=430:500; cix=400:470; lons=stn.ngdc_92m_bathy.lon(rix,cix); lats=stn.ngdc_92m_bathy.lat(rix,cix); fmg; map_freef([min(lons(:)),max(lons(:)),min(lats(:)),max(lats(:))],'none'); [cs,ch]=contour(lons,lats,stn.ngdc_92m_bathy.field(rix,cix),[-2:-1:-20]); clabel(cs,ch,'LabelSpacing',288); plot(stn.lon,stn.lat,'kp'); titlename('Cheeca Rocks bathymetry');
