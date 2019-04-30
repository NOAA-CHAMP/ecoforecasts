1;

if ( ~exist('field','var') )
  field=[]; clear field
  read_nwps_gribs;
end;
if ( ~exist('stns','var') )
  %stns(1) = get_station_from_station_name('mlrf1');
  stns = {};
  for stnm=string({'sanf1','smkf1','lonf1','mlrf1','fwyf1','lkwf1',});
    %disp(stnm);
    stn = get_station_from_station_name(stnm{:});
    stn = load_all_ndbc_data(stn);
    stn = station_bulk_windstress(stn,'ndbc_wind_stress','ndbc_wind1_speed',10,...
                                  'ndbc_wind1_speed','ndbc_wind1_speed',[],'ndbc_barom');
    
    stn = verify_variable(stn,{'ndbc_wind1_speed_3_d_lp','ndbc_air_t_3_d_lp','ndbc_sea_t_3_d_lp','ndbc_wind1_stress_3_d_lp',});
    stns{end+1} = stn;
    clear stn stnm
  end;
end;

for stix=1:numel(stns)
  [lonerr,lonix] = min(abs(field.lon-stns{stix}.lon));
  [laterr,latix] = min(abs(field.lat-stns{stix}.lat));  
  stns{stix}.nwps_sigwavehgt.date = field.date;
  if ( lonerr <= 0.10 && laterr <= 0.10 )
    stns{stix}.nwps_sigwavehgt.data = field.nwps_sigwavehgt.field(:,latix,lonix);
  else
    stns{stix}.nwps_sigwavehgt.data = repmat(nan,size(field.date));
  end;
end;

for stix=1:numel(stns)
  ax=[];
  xl = stns{stix}.nwps_sigwavehgt.date([1,end]);
  if ( isfield(stns{stix},'ndbc_sea_t_3_d_lp') && is_valid_ts(stns{stix}.ndbc_sea_t_3_d_lp) )
    fmg; ax(1)=spt(2,1,1); plot_ts(stns{stix}.nwps_sigwavehgt,stns{stix}.ndbc_wind_stress); xlim(xl); ylim([0.0,3.5]); datetick3; grid on; grid minor; legend('Hs','\tau');
    ax(2)=spt(2,1,2); plot_ts(stns{stix}.ndbc_air_t,stns{stix}.ndbc_sea_t,stns{stix}.ndbc_sea_t_3_d_lp,'k-'); xlim(xl); ylim([21,34]); grid on; grid minor; legend('T_A','T_S','T_S^3^d');
  else
    fmg; plot_ts(stns{stix}.nwps_sigwavehgt); xlim(xl); ylim([0.0,3.5]); ax(1)=gca; legend('Hs');
  end;
  title(ax(1),[upper(stns{stix}.station_name),' ']);
end;
