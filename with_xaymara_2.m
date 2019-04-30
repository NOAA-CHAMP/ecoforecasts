1;

datpath = get_ecoforecasts_path('data');

if ( ~exist('deps','var') )

  rawin = importdata(fullfile(datpath,'Pastreoides individuals by site FL Keys.xlsx'));

  snms = rawin.textdata.Sheet1(2:end,2);
  lats = rawin.data.Sheet1(1:end,1);
  lons = rawin.data.Sheet1(1:end,2);
  deps = rawin.data.Sheet1(1:end,4);

  usnms = unique(snms);
end;


for ix=1:numel(usnms);
  six = find(strcmp(snms,usnms{ix}));
  disp({usnms{ix},numel(six),median(deps(six)),roundn(range(deps(six)),-2),roundn(range(deps(six)),-2)/100,});
end;
