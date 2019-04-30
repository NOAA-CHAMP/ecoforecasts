function stn = load_ncdc_winds(stn,fname)

% m = csvread('data/ncdc_pago_pago_6653724097667dat.txt',2,0);
% rawdat = importdata('data/ncdc_pago_pago_6653724097667dat.txt');

  rawdat = importdata(fname);

  flds = rawdat.textdata(:,5);

  rwndix = find(strcmp(flds,'RWND'));
  rwndix = rwndix - 2; % Two header lines (all text)
  rwndmtx = rawdat.data(rwndix,3:2:end);
  rwnd = reshape(rwndmtx', [1 numel(rwndmtx)])';

  rdirix = find(strcmp(flds,'RDIR'));
  rdirix = rdirix - 2; % Two header lines (all text)
  rdirmtx = rawdat.data(rdirix,3:2:end);
  rdir = reshape(rdirmtx', [1 numel(rdirmtx)])';

  badix = find(abs(rwnd) == 99999 | abs(rdir) == 99999);
  rwnd(badix) = [];
  rdir(badix) = [];

  rwnd = rwnd ./ 10;

  % Convert from miles/hr to kts
  rwnd = rwnd ./ 1.15155;

  yr = floor(rawdat.data(rwndix,1)/100);
  mo = rawdat.data(rwndix,1) - (yr*100);
  dts = [datenum(yr(1),mo(1),1):(datenum(yr(1),mo(1),1)+numel(rwnd)-1)]';


  stn.ncdc_wind_speed.date = dts;
  stn.ncdc_wind_speed.data = rwnd;

  stn.ncdc_wind_dir.date = dts;
  stn.ncdc_wind_dir.data = rdir;

  stn.ncdc_wind_u.date = dts;
  stn.ncdc_wind_u.data = rwnd .* (-sind(rdir));

  stn.ncdc_wind_v.date = dts;
  stn.ncdc_wind_v.data = rwnd .* (-cosd(rdir));

  % save('data/asam_stns_0km.mat', '-append', 'stns');

return;
