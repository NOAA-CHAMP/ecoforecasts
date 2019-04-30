1;

error('This fix-it script should be run ONLY ONCE EVER');

for stnc = { ...
    '42003', ...
    'cryf1', ...
    'dryf1', ...
    'fwyf1', ...
    'lkwf1', ...
    'looe1', ...
    'lonf1', ...
    'mlrf1', ...
    'sanf1', ...
    'smkf1', ...
    'tnrf1', ...
           };

  stnm = stnc{:};
  disp(stnm);

  station = []; clear station
  load(['data/' stnm '_tmd_tide.mat'],'station');

  station.tmd_tide.date = station.tmd_tide.date(:);
  station.tmd_tide_u.date = station.tmd_tide.date;
  station.tmd_tide_v.date = station.tmd_tide.date;
  station.tmd_tide_i_depth.date = station.tmd_tide_i_depth.date(:);

  station.tpxo_tide.date = station.tpxo_tide.date(:);
  station.tpxo_tide_u.date = station.tpxo_tide.date;
  station.tpxo_tide_v.date = station.tpxo_tide.date;
  station.tpxo_tide_i_depth.date = station.tpxo_tide_i_depth.date(:);

  save(['data/' stnm '_tmd_tide.mat'],'station');
end;
