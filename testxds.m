1;

f{1} = { anfld, fcfld{1} };
f{2} = { 'sea_t', anfld };
f{3} = { 'wind1_speed', fcfld{1} };

tol = 3;
maxix = 40000;

for ix = 1:length(f)
  disp(f{ix}{1}); disp(f{ix}{2});
  tic;
    [anix,fcix1] = intersect_dates(station.(f{ix}{1}).date(1:maxix), station.((f{ix}{2})).date(1:maxix), tol);
  toc;
  disp(numel(anix));
  tic;
    [anix,fcix1] = intersect_dates_slow(station.(f{ix}{1}).date(1:maxix), station.((f{ix}{2})).date(1:maxix), tol);
  toc;
  disp(numel(anix));

  tic;
    [anix,fcix1] = intersect_dates(station.(f{ix}{1}).date(1:3:maxix), station.((f{ix}{2})).date(1:maxix), tol);
  toc;
  disp(numel(anix));
  tic;
    [anix,fcix1] = intersect_dates_slow(station.(f{ix}{1}).date(1:3:maxix), station.((f{ix}{2})).date(1:maxix), tol);
  toc;
  disp(numel(anix));

  tic;
    [anix,fcix1] = intersect_dates(station.(f{ix}{1}).date(1:maxix), station.((f{ix}{2})).date(1:3:maxix), tol);
  toc;
  disp(numel(anix));
  tic;
    [anix,fcix1] = intersect_dates_slow(station.(f{ix}{1}).date(1:maxix), station.((f{ix}{2})).date(1:3:maxix), tol);
  toc;
  disp(numel(anix));
end;
