1;

more off;

stnms = {'fwyf1','mlrf1','lonf1','smkf1','sanf1','dryf1'};
stnms = {'smkf1','sanf1','dryf1'};

for ix = 1:length(stnms)
  stnm = stnms{ix};
  dts{ix}.station_name = stnm;

  station = []; clear station;
  station = load_station_data(stnm);

  dts{ix}.pyrs = [];

  station.station_name, fields(station),

  if ( isfield(station, 'licor_surf_par') );
    dvec = datevec(station.licor_surf_par.date);
    dts{ix}.pyrs = unique(dvec(:,1));
  end;
end;

station = []; clear station;
clear datapath stnms stnm ix dvec;

fprintf( 1, '\n\nREPORT:\n' );
for ix = 1:length(dts)
  fprintf( 1, '\n %s: P:%d-%d\n', ...
           dts{ix}.station_name, ...
           min(dts{ix}.pyrs), max(dts{ix}.pyrs) );
end;

more on;
