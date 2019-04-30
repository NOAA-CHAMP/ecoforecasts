function stn = test_station_heat_flux

  x=load('data\test3_0.txt');%read file with 1hr-average data; set your local directory 

  jdy=x(:,1);%time in the form YYYYMMDDHHSS.SS
  dts = datenum(num2str(jdy),'yyyymmddHHMM.SS')

  w=x(:,2); %true wind speed, m/s; etl sonic anemometer
  stn.ndbc_wind1_speed.data = w;

  st=x(:,3);%sea snake temperature, C (0.05 m depth)
  stn.seatemp.data = st;

  a=x(:,4);%air temperature, C (z=14.5 m)
  stn.ndbc_air_t.data = a;

  qa=x(:,5)./1000;%g/kg->kg/kg%air specific humidity, g/kg (z=14.5  m)
  stn.ndbc_spechumid.data = qa;

  dsrf=x(:,6);%downward solar flux, W/m^2 (ETL units)
  stn.erai_dsrf.data = dsrf;

  dlrf=x(:,7);%downward IR flux, W/m^2 (ETL units)
  stn.erai_dlrf.data = dlrf;

  prcp=x(:,8);%rainrate, mm/hr (ETL STI optical rain gauge, uncorrected)
  stn.erai_precip.data = prcp;

  t=x(:,11);%6-m deotg T from MSP, C    
  stn.ndbc_sea_t.data = t;

  p=repmat(1008,size(dts));                     %air pressure
  stn.ndbc_barom.data = p;

  pblz = 600;

  % Assume surface current projects onto 2% of wind velocity [e.g., Ardhuin et al. 2009]
  warning('Estimating projected ocean currents from wind speed');
  ou = 0.020 .* w;

  flds = fieldnames(stn);
  for ix=1:length(flds)
    fld = flds{ix};
    stn.(fld).date = dts;
  end;

  stn.lat=mean(x(:,9));%latitude, deg  (SCS pcode)
  stn.lon=mean(x(:,10));%longitude, deg (SCS pcode)

  wz=15;%anemometer ht
  az=15;%air T height
  qz=15;%humidity height
  stz=6;%bulk water temperature sensor depth, U/APL MSP&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

  jcool=1;
  jwarm=1;
  jwave=0;

  result = real( cor30a_warm(dts,w,ou,t,a,ss,sa,dsrf,dlrf,prcp,pblz,p,stn.lon,stn.lat,wz,az,az,stz,jwarm,jcool,jwave) );

  tau = result(:,3);
  shf = -result(:,1);
  lhf = -result(:,2);
  rhf = -result(:,14);

  diags = [PFX 'cordiags'];
  disp(['Updating ' diags]);
  stn.(diags).date = dts;
  stn.(diags).ustar = result(:,8);
  stn.(diags).tstar = result(:,9);
  stn.(diags).qstar = result(:,10);
  stn.(diags).Cd = result(:,16);
  stn.(diags).Ch = result(:,17);
  stn.(diags).Ce = result(:,18);
  stn.(diags).Ug = result(:,22);
  if ( size(result,2) >= 24 )
    stn.(diags).dtwarm = result(:,23);
    stn.(diags).dxwarm = result(:,24);
  end;

  wind_stress = [PFX 'wind_stress'];
  if ( isfield(stn,wind_stress) ); stn = rmfield(stn,wind_stress); end;
  stn.(wind_stress).date = dts;
  stn.(wind_stress).data = tau;

  sensible_heat_flux = [PFX 'sensible_heat_flux'];
  if ( isfield(stn,sensible_heat_flux) ); stn = rmfield(stn,sensible_heat_flux); end;
  stn.(sensible_heat_flux).date = dts;
  stn.(sensible_heat_flux).data = shf;

  latent_heat_flux = [PFX 'latent_heat_flux'];
  if ( isfield(stn,latent_heat_flux) ); stn = rmfield(stn,latent_heat_flux); end;
  stn.(latent_heat_flux).date = dts;
  stn.(latent_heat_flux).data = lhf;

  net_heat_flux = [PFX 'net_heat_flux'];
  if ( isfield(stn,net_heat_flux) ); stn = rmfield(stn,net_heat_flux); end;
  stn.(net_heat_flux).date = dts;
  stn.(net_heat_flux).data = shf + lhf;

  % Add in rain heat flux, if we were given rain data
  if (~isempty(rhf))
    rain_heat_flux = [PFX 'rain_heat_flux'];
    if ( isfield(stn,rain_heat_flux) ); stn = rmfield(stn,rain_heat_flux); end;
    stn.(rain_heat_flux).date = dts;
    stn.(rain_heat_flux).data = rhf;

    stn.(net_heat_flux).data = stn.(net_heat_flux).data + stn.(rain_heat_flux).data;
  end;

return;
