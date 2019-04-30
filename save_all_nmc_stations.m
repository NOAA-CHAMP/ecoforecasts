1;

datapath = get_ecoforecasts_path('data');

STATIONS = get_all_station_metadata();

flds = {'nmc_wind_u','nmc_wind_v'};
for stix=1:length(STATIONS.codes)
  stations(stix).station_name = STATIONS.codes{stix};
  stations(stix).lon = STATIONS.lons(stix);
  stations(stix).lat = STATIONS.lats(stix);
  stations(stix).depth = STATIONS.depths(stix);
  for ix = 1:length(flds)
    stations(stix).(flds{ix}) = [];
    stations(stix).(['native_' flds{ix}]) = [];
  end;
end;

stations = get_nmc_winds(stations);

for stix=1:length(stations)
  station = stations(stix);
  stnm = lower(station.station_name);
  matfname = fullfile(datapath,[stnm '_nmc_winds.mat']);
  if ( exist(matfname,'file') )
    warning('MAT already exists - will NOT overwrite! "%s"',matfname);
  else
    %DEBUG:
    disp(matfname);
    save(matfname,'station');
  end;
  station = []; clear station;
end;

%%%% ??? stations = []; clear stations;
