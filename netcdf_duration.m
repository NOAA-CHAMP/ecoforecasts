function strs = netcdf_duration(dts)
%function strs = netcdf_duration(dts)
% Returns array of Strings representing each (presumably relative) DATENUM in
% DTS as a standard netCDF "duration", e.g., 'P15Y8M26DT23H25M0S'.
  vecs = datevec(dts);
  strs = sprintf('P%dY%dM%dDT%dH%dM%dS',vecs(1),vecs(2),vecs(3),vecs(4),vecs(5),round(vecs(6)));
return;
