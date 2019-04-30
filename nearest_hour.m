function hrs = nearest_hour(dts)
%function hrs = nearest_hour(dts)
% Return DATENUM of nearest whole hour for each timestamp DTS

  hrs = round(dts*24)/24;

return;
