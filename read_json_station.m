function stn = read_json_station(stn,fname)

  if ( ~exist('stn','var') ); stn=[]; end;
  res = read_json(fname);
  dts = datenum(1970,1,1)+ (res.index/3600/24/1e3);
  stn = merge_station_data(newstn);

return;
