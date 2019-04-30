1;

more off;
warnstat = warning('OFF','GFS:BadURL');
diary('factize_gfs.log');

timenow;

args.bbpath = get_ecoforecasts_path('bb');
args.model = string('gfs_0p5');

%https://www.ncei.noaa.gov/thredds/dodsC/gfs-g4-anl-files/201007/20100730/gfsanl_4_20100730_0600_000.grb2.html
%https://www.ncei.noaa.gov/thredds/dodsC/gfs-g4-anl-files/201607/20160730/gfsanl_4_20160730_0600_000.grb2.html

if ( args.model == 'gfs_0p5' )
  BASEURL = 'https://www.ncei.noaa.gov/thredds/dodsC/gfs-g4-anl-files/';
elseif ( args.model == 'gfs_0p25' )
  BASEURL = 'https://rda.ucar.edu/thredds/dodsC/files/g/ds084.1/';
else
  error('Unrecognized model code %s',args.model);
end;

alldts = datenum(2007,1,1):floor(now-1);
%alldts = datenum(2007,1,1);
%alldts = datenum(2017,7,20);

% Run ("reference") hours
allhrs = 0:6:18;
%allhrs = 12;

% Forecast hours
%allfhrs = 0:3:6;
allfhrs = 0:3:3;
%allfhrs = 6;

STATIONS = get_all_station_metadata;
STATIONS.lons = mod(360+STATIONS.lons,360);

for stix=1:numel(STATIONS.codes)
  stnm = lower(string(STATIONS.codes{stix}));
  stns(stix).station_name = stnm;
  stns(stix).lon = STATIONS.lons(stix);
  stns(stix).lat = STATIONS.lats(stix);
  stns(stix).depth = STATIONS.depths(stix);

  stns(stix).u.date = repmat(nan,[numel(alldts)*numel(allhrs)*numel(allfhrs),1]);
  stns(stix).u.data = stns(stix).u.date;
  stns(stix).v = stns(stix).u;
  stns(stix).speed = stns(stix).u;
  stns(stix).dir = stns(stix).u;
  stns(stix).airtemp = stns(stix).u;
end;

% >> url = 'https://www.ncei.noaa.gov/thredds/dodsC/gfs-g4-anl-files/201707/20170720/gfsanl_4_20170720_1200_006.grb2';
% >>       url = sprintf('%s%04d%02d/%04d%02d%02d/gfsanl_4_%04d%02d%02d_%02d00_%03d.grb2',BASEURL,yr,mo,yr,mo,dy,yr,mo,dy,hr,fhr);
% >> url
% url =
% https://www.ncei.noaa.gov/thredds/dodsC/gfs-g4-anl-files/201707/20170720/gfsanl_4_20170720_1200_006.grb2

tsix = 1;

for dtix = 1:numel(alldts)
  dt = alldts(dtix);
  %DEBUG:
  disp(datestr(dt));

  yr = get_year(dt);
  mo = get_month(dt);
  dy = get_dom(dt);
  jd = get_jday(dt);

  %fname = fullfile(args.bbpath,sprintf('gfs_0p5_%04d%02d%02d.bb',yr,mo,dy));
  fname = fullfile(args.bbpath,sprintf('%s_wind_%04d%03d.bb',args.model,yr,jd));
  fid = fopen(fname,'w+');
  if ( isempty(fid) || fid<2 )
    error('Unable to open output file "%s"',fname);
  end;

  for hrix = 1:numel(allhrs)
    hr = allhrs(hrix);
    for fhrix = 1:numel(allfhrs)
      fhr = allfhrs(fhrix);

      if ( args.model == 'gfs_0p5' )
        %dp = pydap.client.open_url('https://www.ncei.noaa.gov/thredds/dodsC/gfs-g4-anl-files/201807/20180720/gfsanl_4_20180720_1200_006.grb2')
        url = sprintf('%s%04d%02d/%04d%02d%02d/gfsanl_4_%04d%02d%02d_%02d00_%03d.grb2',BASEURL,yr,mo,yr,mo,dy,yr,mo,dy,hr,fhr);
      elseif ( args.model == 'gfs_0p25' )
        error('%s MODEL NOT YET IMPLEMENTED',args.model);
        % SAMPLE DATAURL: https://rda.ucar.edu/thredds/dodsC/files/g/ds084.1/2015/20150115/gfs.0p25.2015011500.f000.grib2?Temperature_height_above_ground[0:1:0][0:1:0][258][1120],Dewpoint_temperature_height_above_ground[0:1:0][0:1:0][258][1120],u-component_of_wind_height_above_ground[0:1:0][0:1:0][258][1120],v-component_of_wind_height_above_ground[0:1:0][0:1:0][258][1120]
        % pyses = pydap.cas.urs.setup_session(username, password, check_url=url)
        % dp = pydap.client.open_url(url, session=session);
      end;

      try,
        nc = mDataset(url);
        if ( ~exist('lons','var') || isempty(lons) )
          lons = cast(nc{'lon'}(:,:),'double');
          lats = cast(nc{'lat'}(:,:),'double');
        end;
        % T = cast(nc{'Temperature_height_above_ground'}(1,1,:,:),'double');
        T = cast(nc{'Temperature_surface'}(1,:,:),'double');
        u = cast(nc{'u-component_of_wind_height_above_ground'}(1,1,:,:),'double');
        v = cast(nc{'v-component_of_wind_height_above_ground'}(1,1,:,:),'double');
        % u = cast(nc{'u-component_of_wind_planetary_boundary'}(1,:,:),'double');
        % v = cast(nc{'v-component_of_wind_planetary_boundary'}(1,:,:),'double');
        close(nc); clear nc

        [sp,dr] = uv_to_spddir(u,v);
        T = T - 273.15;
      catch ME,
        try,
          close(nc); clear nc
          %fclose(fid);
        end;
        warning('GFS:BadURL','Unable to open "%s"',url);
        %throw(ME);
      end; %try,

      if ( ~exist('u','var') || isempty(u) || isempty(v) || isempty(T) )
        warning('Missing or empty data: "%s"',url);

      else

        % (mlrf1 07/09/2008 07 191 1754 amodis seadas_sst 27.75)
        % (pvgf1 07/24/2018 18 205 0000 gfs speed 13.45)
        %for stix=33 % MLRF1
        for stix=1:numel(STATIONS.codes)
          stnm = lower(string(STATIONS.codes{stix}));
          [lonerr,lonix] = min(abs(lons-STATIONS.lons(stix)));
          [laterr,latix] = min(abs(lats-STATIONS.lats(stix)));
          % (pvgf1 07/24/2018 18 205 0000 gfs u 1.52)
          fprintf(fid,'(%s %02d/%02d/%04d %02d %03d %02d00 gfs u %g)\n',...
                  stnm,mo,dy,yr,mod(yr,100),jd,hr+fhr,u(latix,lonix));
          fprintf(fid,'(%s %02d/%02d/%04d %02d %03d %02d00 gfs v %g)\n',...
                  stnm,mo,dy,yr,mod(yr,100),jd,hr+fhr,v(latix,lonix));
          fprintf(fid,'(%s %02d/%02d/%04d %02d %03d %02d00 gfs speed %g)\n',...
                  stnm,mo,dy,yr,mod(yr,100),jd,hr+fhr,sp(latix,lonix));
          fprintf(fid,'(%s %02d/%02d/%04d %02d %03d %02d00 gfs dir %g)\n',...
                  stnm,mo,dy,yr,mod(yr,100),jd,hr+fhr,dr(latix,lonix));
          fprintf(fid,'(%s %02d/%02d/%04d %02d %03d %02d00 gfs airtemp %g)\n',...
                  stnm,mo,dy,yr,mod(yr,100),jd,hr+fhr,T(latix,lonix));

          stns(stix).u.date(tsix) = dt + ((hr+fhr)/24);
          stns(stix).u.data(tsix) = u(latix,lonix);
          stns(stix).v.date(tsix) = dt + ((hr+fhr)/24);
          stns(stix).v.data(tsix) = v(latix,lonix);
          stns(stix).speed.date(tsix) = dt + ((hr+fhr)/24);
          stns(stix).speed.data(tsix) = sp(latix,lonix);
          stns(stix).dir.date(tsix) = dt + ((hr+fhr)/24);
          stns(stix).dir.data(tsix) = dr(latix,lonix);
          stns(stix).airtemp.date(tsix) = dt + ((hr+fhr)/24);
          stns(stix).airtemp.data(tsix) = T(latix,lonix);
        end; %for stix=1:numel(STATIONS.codes)
        tsix = tsix + 1;

        u=[]; v=[]; sp=[]; dr=[]; T=[]; clear u v sp dr T lonerr lonix laterr latix stix stnm
      end; %if ( exist('u','var') )

    end; %for fhrix = 1:numel(allfhrs)
  end; %for hrix = 1:numel(allhrs)

  fclose(fid);
  clear fid fname hrix jd 

end; %for dtix = 1:numel(alldts)

timenow;

save(fullfile(get_ecoforecasts_path('data'),'factize_gfs.mat'),'stns','-v7.3');

stns(1).u_3_d_avg=[];
stns(1).v_3_d_avg=[];
stns(1).speed_3_d_avg=[];
stns(1).dir_3_d_avg=[];
stns(1).airtemp_3_d_avg=[];
for stix=1:numel(STATIONS.codes)
  stnm = lower(string(STATIONS.codes{stix}));
  stns(stix).lat = STATIONS.lats(stix);
  stns(stix).u.date(tsix:end) = [];
  stns(stix).u.data(tsix:end) = [];
  stns(stix).v.date(tsix:end) = [];
  stns(stix).v.data(tsix:end) = [];
  stns(stix).speed.date(tsix:end) = [];
  stns(stix).speed.data(tsix:end) = [];
  stns(stix).dir.date(tsix:end) = [];
  stns(stix).dir.data(tsix:end) = [];
  stns(stix).airtemp.date(tsix:end) = [];
  stns(stix).airtemp.data(tsix:end) = [];

  stns(stix)=verify_variable(stns(stix),{'u_3_d_avg','v_3_d_avg','speed_3_d_avg','airtemp_3_d_avg'},true);
  stns(stix).dir_3_d_avg.date = stns(stix).u_3_d_avg.date;
  stns(stix).dir_3_d_avg.data = uv_to_dir(stns(stix).u_3_d_avg.data,stns(stix).v_3_d_avg.data);

  station=stns(stix);
  disp(stnm);
  save(fullfile(get_ecoforecasts_path('data'),['factize_gfs_',char(stnm)]),'station','-v7.3');
  station=[]; clear station;
end; %for stix=1:numel(STATIONS.codes)

save(fullfile(get_ecoforecasts_path('data'),'factize_gfs.mat'),'stns','lats','lons','-v7.3');

clear ans args BASEURL dt dtix dy fhr fhrix hr laterr latix lonerr lonix mo url yr
%clear alldts allfhrs lats lons STATIONS

timenow;

diary off;
warning(warnstat); clear warnstat
more on;


if 0;
  for ix=1:numel(stns);
    fh=fmg; boxplot_ts(stns(ix).speed_3_d_avg);
    titlename(textize(stns(ix).station_name));
    waitfor(fh); try, close(fh); clear fh; end;
  end;
end;


% (pvgf1 07/24/2018 18 205 0000 gfs count 1.00)
% (pvgf1 07/24/2018 18 205 0000 gfs u 1.52)
% (pvgf1 07/24/2018 18 205 0000 gfs v 6.75)
% (pvgf1 07/24/2018 18 205 0000 gfs speed 13.45)
% (pvgf1 07/24/2018 18 205 0000 gfs dir 192.70)

% # grep 'pvgf1.*speed' bb/gfs_wind_2018_205.bb
% (pvgf1 07/24/2018 18 205 0000 gfs speed 13.45)
% (pvgf1 07/24/2018 18 205 0300 gfs speed 14.00)
% (pvgf1 07/24/2018 18 205 0600 gfs speed 13.21)
% (pvgf1 07/24/2018 18 205 0900 gfs speed 16.28)
% (pvgf1 07/24/2018 18 205 1200 gfs speed 11.34)
% (pvgf1 07/24/2018 18 205 1500 gfs speed 14.25)
% (pvgf1 07/24/2018 18 205 1800 gfs speed 11.65)
% (pvgf1 07/24/2018 18 205 2100 gfs speed 14.43)


%<bound method DatasetType.keys of
% {'LatLon_Projection', 'lat', 'lon', 'reftime', 'time', 'time_bounds', 'time1', 'isobaric', 'height_above_ground_layer', 'height_above_ground_layer_bounds', 'height_above_ground', 'sigma', 'depth_below_surface_layer', 'depth_below_surface_layer_bounds', 'pressure_difference_layer', 'pressure_difference_layer_bounds', 'isobaric1', 'isobaric2', 'height_above_ground1', 'height_above_ground2', 'isobaric3', 'altitude_above_msl', 'height_above_ground3', 'height_above_ground_layer1', 'height_above_ground_layer1_bounds', 'pressure_difference_layer1', 'pressure_difference_layer1_bounds', 'isobaric4', 'pressure_difference_layer2', 'pressure_difference_layer2_bounds', 'sigma_layer', 'sigma_layer_bounds', 'height_above_ground4', 'isobaric5', 'potential_vorticity_surface', 'Absolute_vorticity_isobaric': {'Absolute_vorticity_isobaric', 'time1', 'isobaric4', 'lat', 'lon'},
% 'Albedo_surface_6_Hour_Average': {'time', 'lat', 'lon'},
% 'Apparent_temperature_height_above_ground': {'time1', 'height_above_ground', 'lat', 'lon'},
% 'Cloud_mixing_ratio_isobaric': {'time1', 'isobaric3', 'lat', 'lon'},
% 'Cloud_water_entire_atmosphere_single_layer': {'time1', 'lat', 'lon'},
% 'Convective_available_potential_energy_surface': {'time1', 'lat', 'lon'},
% 'Convective_available_potential_energy_pressure_difference_layer': {'time1', 'pressure_difference_layer2', 'lat', 'lon'},
% 'Convective_inhibition_surface': {'time1', 'lat', 'lon'},
% 'Convective_inhibition_pressure_difference_layer': {'time1', 'pressure_difference_layer2', 'lat', 'lon'},
% 'Convective_precipitation_surface_6_Hour_Accumulation': {'time', 'lat', 'lon'},
% 'Dewpoint_temperature_height_above_ground': {'time1', 'height_above_ground', 'lat', 'lon'},
% 'Geopotential_height_potential_vorticity_surface': {'time1', 'potential_vorticity_surface', 'lat', 'lon'},
% 'Geopotential_height_surface': {'time1', 'lat', 'lon'},
% 'Geopotential_height_isobaric': {'time1', 'isobaric', 'lat', 'lon'},
% 'Geopotential_height_zeroDegC_isotherm': {'time1', 'lat', 'lon'},
% 'Geopotential_height_maximum_wind': {'time1', 'lat', 'lon'},
% 'Geopotential_height_tropopause': {'time1', 'lat', 'lon'},
% 'Geopotential_height_highest_tropospheric_freezing': {'time1', 'lat', 'lon'},
% 'Haines_index_surface': {'time1', 'lat', 'lon'},
% 'ICAO_Standard_Atmosphere_Reference_Height_maximum_wind': {'time1', 'lat', 'lon'},
% 'ICAO_Standard_Atmosphere_Reference_Height_tropopause': {'time1', 'lat', 'lon'},
% 'Ice_cover_surface': {'time1', 'lat', 'lon'},
% 'Land_cover_0__sea_1__land_surface': {'time1', 'lat', 'lon'},
% 'Latent_heat_net_flux_surface_6_Hour_Average': {'time', 'lat', 'lon'},
% 'Maximum_temperature_height_above_ground_6_Hour_Maximum': {'time', 'height_above_ground', 'lat', 'lon'},
% 'Minimum_temperature_height_above_ground_6_Hour_Minimum': {'time', 'height_above_ground', 'lat', 'lon'},
% 'Momentum_flux_u-component_surface_6_Hour_Average': {'time', 'lat', 'lon'},
% 'Momentum_flux_v-component_surface_6_Hour_Average': {'time', 'lat', 'lon'},
% 'Per_cent_frozen_precipitation_surface': {'time1', 'lat', 'lon'},
% 'Potential_temperature_sigma': {'time1', 'sigma', 'lat', 'lon'},
% 'Precipitable_water_entire_atmosphere_single_layer': {'time1', 'lat', 'lon'},
% 'Precipitation_rate_surface_6_Hour_Average': {'time', 'lat', 'lon'},
% 'Pressure_potential_vorticity_surface': {'time1', 'potential_vorticity_surface', 'lat', 'lon'},
% 'Pressure_convective_cloud_bottom': {'time1', 'lat', 'lon'},
% 'Pressure_low_cloud_bottom_6_Hour_Average': {'time', 'lat', 'lon'},
% 'Pressure_middle_cloud_bottom_6_Hour_Average': {'time', 'lat', 'lon'},
% 'Pressure_high_cloud_bottom_6_Hour_Average': {'time', 'lat', 'lon'},
% 'Pressure_high_cloud_top_6_Hour_Average': {'time', 'lat', 'lon'},
% 'Pressure_convective_cloud_top': {'time1', 'lat', 'lon'},
% 'Pressure_low_cloud_top_6_Hour_Average': {'time', 'lat', 'lon'},
% 'Pressure_middle_cloud_top_6_Hour_Average': {'time', 'lat', 'lon'},
% 'Pressure_surface': {'time1', 'lat', 'lon'},
% 'Pressure_maximum_wind': {'time1', 'lat', 'lon'},
% 'Pressure_tropopause': {'time1', 'lat', 'lon'},
% 'Pressure_height_above_ground': {'time1', 'height_above_ground2', 'lat', 'lon'},
% 'Pressure_reduced_to_MSL_msl': {'time1', 'lat', 'lon'},
% 'Relative_humidity_entire_atmosphere_single_layer': {'time1', 'lat', 'lon'},
% 'Relative_humidity_highest_tropospheric_freezing': {'time1', 'lat', 'lon'},
% 'Relative_humidity_pressure_difference_layer': {'time1', 'pressure_difference_layer', 'lat', 'lon'},
% 'Relative_humidity_sigma_layer': {'time1', 'sigma_layer', 'lat', 'lon'},
% 'Relative_humidity_isobaric': {'time1', 'isobaric', 'lat', 'lon'},
% 'Relative_humidity_zeroDegC_isotherm': {'time1', 'lat', 'lon'},
% 'Relative_humidity_height_above_ground': {'time1', 'height_above_ground', 'lat', 'lon'},
% 'Relative_humidity_sigma': {'time1', 'sigma', 'lat', 'lon'},
% 'Sensible_heat_net_flux_surface_6_Hour_Average': {'time', 'lat', 'lon'},
% 'Snow_depth_surface': {'time1', 'lat', 'lon'},
% 'Soil_temperature_depth_below_surface_layer': {'time1', 'depth_below_surface_layer', 'lat', 'lon'},
% 'Specific_humidity_pressure_difference_layer': {'time1', 'pressure_difference_layer', 'lat', 'lon'},
% 'Specific_humidity_height_above_ground': {'time1', 'height_above_ground1', 'lat', 'lon'},
% 'Storm_relative_helicity_height_above_ground_layer': {'time1', 'height_above_ground_layer', 'lat', 'lon'},
% 'Temperature_potential_vorticity_surface': {'time1', 'potential_vorticity_surface', 'lat', 'lon'},
% 'Temperature_low_cloud_top_6_Hour_Average': {'time', 'lat', 'lon'},
% 'Temperature_middle_cloud_top_6_Hour_Average': {'time', 'lat', 'lon'},
% 'Temperature_high_cloud_top_6_Hour_Average': {'time', 'lat', 'lon'},
% 'Temperature_surface': {'time1', 'lat', 'lon'},
% 'Temperature_pressure_difference_layer': {'time1', 'pressure_difference_layer', 'lat', 'lon'},
% 'Temperature_isobaric': {'time1', 'isobaric', 'lat', 'lon'},
% 'Temperature_maximum_wind': {'time1', 'lat', 'lon'},
% 'Temperature_altitude_above_msl': {'time1', 'altitude_above_msl', 'lat', 'lon'},
% 'Temperature_height_above_ground': {'time1', 'height_above_ground4', 'lat', 'lon'},
% 'Temperature_tropopause': {'time1', 'lat', 'lon'},
% 'Temperature_sigma': {'time1', 'sigma', 'lat', 'lon'},
% 'Total_cloud_cover_entire_atmosphere_6_Hour_Average': {'time', 'lat', 'lon'},
% 'Total_cloud_cover_low_cloud_6_Hour_Average': {'time', 'lat', 'lon'},
% 'Total_cloud_cover_middle_cloud_6_Hour_Average': {'time', 'lat', 'lon'},
% 'Total_cloud_cover_high_cloud_6_Hour_Average': {'time', 'lat', 'lon'},
% 'Total_cloud_cover_convective_cloud': {'time1', 'lat', 'lon'},
% 'Total_cloud_cover_boundary_layer_cloud_6_Hour_Average': {'time', 'lat', 'lon'},
% 'Total_ozone_entire_atmosphere_single_layer': {'time1', 'lat', 'lon'},
% 'Total_precipitation_surface_6_Hour_Accumulation': {'time', 'lat', 'lon'},
% 'Categorical_Rain_surface_6_Hour_Average': {'time', 'lat', 'lon'},
% 'Categorical_Freezing_Rain_surface_6_Hour_Average': {'time', 'lat', 'lon'},
% 'Categorical_Ice_Pellets_surface_6_Hour_Average': {'time', 'lat', 'lon'},
% 'Categorical_Snow_surface_6_Hour_Average': {'time', 'lat', 'lon'},
% 'Convective_Precipitation_Rate_surface_6_Hour_Average': {'time', 'lat', 'lon'},
% 'Potential_Evaporation_Rate_surface': {'time1', 'lat', 'lon'},
% 'Ozone_Mixing_Ratio_isobaric': {'time1', 'isobaric5', 'lat', 'lon'},
% 'Icing_severity_isobaric': {'time1', 'isobaric1', 'lat', 'lon'},
% 'Vertical_Speed_Shear_tropopause': {'time1', 'lat', 'lon'},
% 'Vertical_Speed_Shear_potential_vorticity_surface': {'time1', 'potential_vorticity_surface', 'lat', 'lon'},
% 'U-Component_Storm_Motion_height_above_ground_layer': {'time1', 'height_above_ground_layer1', 'lat', 'lon'},
% 'V-Component_Storm_Motion_height_above_ground_layer': {'time1', 'height_above_ground_layer1', 'lat', 'lon'},
% 'Ventilation_Rate_planetary_boundary': {'time1', 'lat', 'lon'},
% 'MSLP_Eta_model_reduction_msl': {'time1', 'lat', 'lon'},
% '5-Wave_Geopotential_Height_isobaric': {'time1', 'isobaric2', 'lat', 'lon'},
% 'Zonal_Flux_of_Gravity_Wave_Stress_surface_6_Hour_Average': {'time', 'lat', 'lon'},
% 'Meridional_Flux_of_Gravity_Wave_Stress_surface_6_Hour_Average': {'time', 'lat', 'lon'},
% 'Planetary_Boundary_Layer_Height_surface': {'time1', 'lat', 'lon'},
% 'Pressure_of_level_from_which_parcel_was_lifted_pressure_difference_layer': {'time1', 'pressure_difference_layer1', 'lat', 'lon'},
% 'Downward_Short-Wave_Radiation_Flux_surface_6_Hour_Average': {'time', 'lat', 'lon'},
% 'Upward_Short-Wave_Radiation_Flux_surface_6_Hour_Average': {'time', 'lat', 'lon'},
% 'Upward_Short-Wave_Radiation_Flux_atmosphere_top_6_Hour_Average': {'time', 'lat', 'lon'},
% 'Downward_Long-Wave_Radp_Flux_surface_6_Hour_Average': {'time', 'lat', 'lon'},
% 'Upward_Long-Wave_Radp_Flux_surface_6_Hour_Average': {'time', 'lat', 'lon'},
% 'Upward_Long-Wave_Radp_Flux_atmosphere_top_6_Hour_Average': {'time', 'lat', 'lon'},
% 'Cloud_Work_Function_entire_atmosphere_single_layer_6_Hour_Average': {'time', 'lat', 'lon'},
% 'Sunshine_Duration_surface': {'time1', 'lat', 'lon'},
% 'Surface_Lifted_Index_surface': {'time1', 'lat', 'lon'},
% 'Best_4_layer_Lifted_Index_surface': {'time1', 'lat', 'lon'},
% 'Volumetric_Soil_Moisture_Content_depth_below_surface_layer': {'time1', 'depth_below_surface_layer', 'lat', 'lon'},
% 'Ground_Heat_Flux_surface_6_Hour_Average': {'time', 'lat', 'lon'},
% 'Wilting_Point_surface': {'time1', 'lat', 'lon'},
% 'Land-sea_coverage_nearest_neighbor_land1sea0_surface': {'time1', 'lat', 'lon'},
% 'Field_Capacity_surface': {'time1', 'lat', 'lon'},
% 'Vertical_velocity_pressure_isobaric': {'time1', 'isobaric3', 'lat', 'lon'},
% 'Vertical_velocity_pressure_sigma': {'time1', 'sigma', 'lat', 'lon'},
% 'Visibility_surface': {'time1', 'lat', 'lon'},
% 'Water_equivalent_of_accumulated_snow_depth_surface': {'time1', 'lat', 'lon'},
% 'Water_runoff_surface_6_Hour_Accumulation': {'time', 'lat', 'lon'},
% 'Wind_speed_gust_surface': {'time1', 'lat', 'lon'},
% 'u-component_of_wind_planetary_boundary': {'time1', 'lat', 'lon'},
% 'u-component_of_wind_potential_vorticity_surface': {'time1', 'potential_vorticity_surface', 'lat', 'lon'},
% 'u-component_of_wind_pressure_difference_layer': {'time1', 'pressure_difference_layer', 'lat', 'lon'},
% 'u-component_of_wind_isobaric': {'time1', 'isobaric', 'lat', 'lon'},
% 'u-component_of_wind_maximum_wind': {'time1', 'lat', 'lon'},
% 'u-component_of_wind_altitude_above_msl': {'time1', 'altitude_above_msl', 'lat', 'lon'},
% 'u-component_of_wind_height_above_ground': {'time1', 'height_above_ground3', 'lat', 'lon'},
% 'u-component_of_wind_tropopause': {'time1', 'lat', 'lon'},
% 'u-component_of_wind_sigma': {'time1', 'sigma', 'lat', 'lon'},
% 'v-component_of_wind_planetary_boundary': {'time1', 'lat', 'lon'},
% 'v-component_of_wind_potential_vorticity_surface': {'time1', 'potential_vorticity_surface', 'lat', 'lon'},
% 'v-component_of_wind_pressure_difference_layer': {'time1', 'pressure_difference_layer', 'lat', 'lon'},
% 'v-component_of_wind_isobaric': {'time1', 'isobaric', 'lat', 'lon'},
% 'v-component_of_wind_maximum_wind': {'time1', 'lat', 'lon'},
% 'v-component_of_wind_altitude_above_msl': {'time1', 'altitude_above_msl', 'lat', 'lon'},
% 'v-component_of_wind_height_above_ground': {'time1', 'height_above_ground3', 'lat', 'lon'},
% 'v-component_of_wind_tropopause': {'time1', 'lat', 'lon'},
% 'v-component_of_wind_sigma': {'time1', 'sigma', 'lat', 'lon'}}>
