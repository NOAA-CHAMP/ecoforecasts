function disp_langdon_stn(stn)
  disp(sprintf('%s,%s,%g-%g,%g-%g,%.2f',stn.station_name,stn.long_name,...
               degrees2dm(stn.lat),degrees2dm(stn.lon),unitsratio('ft','m')*stn.depth));
return;
