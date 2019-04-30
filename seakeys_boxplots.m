1;

for stn = {smkf1,lonf1,mlrf1,fwyf1}

  station_boxplots(stn{:},'ndbc_air_t','Air Temperature (^oC)',[5 35],1987:2009,[],0.75,[true,false,false,false]);
  station_boxplots(stn{:},'ndbc_air_t','Air Temperature (^oC)',[5 35],1993:2010,[],[],[false,true,false,true]);

  station_boxplots(stn{:},'ndbc_sea_t','Sea Temperature (^oC)',[5 35],1987:2009,[],0.75,[true,false,false,false]);
  station_boxplots(stn{:},'ndbc_sea_t','Sea Temperature (^oC)',[5 35],1993:2010,[],[],[false,true,false,true]);

  station_boxplots(stn{:},'ndbc_wind1_speed','Wind Speed (kts)',[0 50],1987:2009,[],0.75,[true,false,false,false]);
  station_boxplots(stn{:},'ndbc_wind1_speed','Wind Speed (kts)',[0 50],1993:2010,[],[],[false,true,false,true]);

  stn = []; clear stn;

end;
