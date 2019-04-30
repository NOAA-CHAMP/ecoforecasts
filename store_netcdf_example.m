1;

%Open the file
ncid = netcdf.create(['./thefile.nc'],'NC_WRITE')
 
%Define the dimensions
dimidt = netcdf.defDim(ncid,'time',mytimesize);
dimidp = netcdf.defDim(ncid,'pressure',mypressuresize);
dimidlat = netcdf.defDim(ncid,'latitude',mylatitudesize);
dimidlon = netcdf.defDim(ncid,'longitude',mylongitudesize);
 
%Define IDs for the dimension variables (pressure,time,latitude,...)
time_ID=netcdf.defVar(ncid,'time','double',[dimidt]);
pressure_ID=netcdf.defVar(ncid,'pressure','double',[dimidp]);
latitude_ID=netcdf.defVar(ncid,'latitude','double',[dimidlat]);
longitude_ID=netcdf.defVar(ncid,'longitude','double',[dimidlon]);
 
%Define the main variable ()
temperature_ID = netcdf.defVar(ncid,'temperature','double',[dimidt dimidp dimidlat dimidlon]);
 
%We are done defining the NetCdf
netcdf.endDef(ncid);
 
%Then store the dimension variables in
netcdf.putVar(ncid,time_ID,mytimearray);
netcdf.putVar(ncid,pressure_ID,mypressurearray);
netcdf.putVar(ncid,latitude_ID,mylatitudearray);
netcdf.putVar(ncid,longitude_ID,mylongitudearray);
 
%Then store my main variable
netcdf.putVar(ncid,temperature_ID,mytemperaturearray);
 
%We're done, close the netcdf
netcdf.close(ncid)
