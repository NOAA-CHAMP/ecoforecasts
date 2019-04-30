1;

for cstnm={'fwyf1','mlrf1','lonf1','smkf1','sanf1','dryf1',};

  stn=[]; clear stn

  stnm = cstnm{:};
  stn = get_station_from_station_name(stnm); stn = load_all_ndbc_data(stn);
  [ta,ts] = intersect_tses(stn.ndbc_air_t,stn.ndbc_sea_t);

  [acl,asd] = climatologize_time_series(ts.date,ta.data,'Ta');
  [scl,ssd] = climatologize_time_series(ts.date,ts.data,'Ts');
  %{
  fmg; plot(1:365,acl,1:365,scl);
  fmg; plot(1:365,acl,'b-',1:365,acl-asd,'b:',1:365,acl+asd,'b:',1:365,scl,'r-',1:365,scl-ssd,'r:',1:365,scl+ssd,'r:');
  scatter_fit_ts_seasons(ta,ts,[],[],'Ta','Ts');
  corr2(ta.data,ts.data)
  jfmix = find(get_season(ta.date)==1);
  amjix = find(get_season(ta.date)==2);
  jasix = find(get_season(ta.date)==3);
  ondix = find(get_season(ta.date)==4);
  corr2(ta.data(jfmix),ts.data(jfmix))
  corr2(ta.data(amjix),ts.data(amjix))
  corr2(ta.data(jasix),ts.data(jasix))
  corr2(ta.data(ondix),ts.data(ondix))
  %}

  tsa.date=ts.date; tsa.data=ts.data-scl(get_jday_no_leap(ts.date))';
  taa.date=ta.date; taa.data=ta.data-scl(get_jday_no_leap(ta.date))';
  corr2(taa.data,tsa.data)
  %scatter_fit_ts_seasons(taa,tsa,[],[],'Ta','Ts');

  for mo=1:12;
    moix = find(get_month(taa.date)==mo);
    N(mo) = numel(moix);
    R(mo) = corr2(taa.data(moix),tsa.data(moix));
  end;
  N,
  R,

  fmg;
  bar(1:12,R.^2);
  axis([0.5,12.5,0,1]);
  titlename([upper(stnm),' Ta - Ts R^2']);

end;
