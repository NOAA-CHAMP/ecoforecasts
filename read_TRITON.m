
WIND_DIRECTION
WIND_SPEED
AIR_TEMPERATURE
RELATIVE_HUMIDITY
AIR_PRESSURE
PRECIPITATION
SHORTWAVE_RADIATION
NORTHWARD_CURRENT
EASTWARD_CURRENT
VERTICAL_CURRENT
SURFACE_SALINITY
SURFACE_TEMPERATURE
VERTICAL_MEAN_TEMPERATURE
TWENTY_DEGREES_DEPTH
DYNAMIC_HEIGHT
NORTHWARD_WIND
EASTWARD_WIND
CONDUCTIVITY
WATER_TEMPERATURE
WATER_PRESSURE
SALINITY
DENSITY

nc = mDataset('d:/ecoforecasts/data/TRITON_H_10_200109-201807.nc');
t = cast(nc{'TIME'}(:),'double');
dt = datenum(1970,1,1,0,0,0) + (t/3600/24);
dr = cast(nc{'WIND_DIRECTION'}(:,:,:),'double');
sp = cast(nc{'WIND_SPEED'}(:,:,:),'double');
close(nc); clear nc ans
