1;

datapath = 'data';

regions = {'asam','freef','ecarib','gbr'};

for ix = 1:length(regions)
  region = regions{ix},
  tlon = load(fullfile(datapath, [region '.lons.dat']));
  tlat = load(fullfile(datapath, [region '.lats.dat']));
  topo = load(fullfile(datapath, [region '.topo.dat']));

  save(fullfile(datapath, [region '.topo.mat']), 'tlon', 'tlat', 'topo');
end;
