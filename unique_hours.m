function hrs = unique_hours(dts)
%function hrs = unique_hours(dts)
% Return sorted DATENUM of unique whole hours nearest each element of DTS

  hrs = unique(round(dts*24)/24);

return;
