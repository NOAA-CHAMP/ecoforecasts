function [brkix,brkln] = get_section_breaks(bath,stationCoords)
%function [brkix,brkln] = get_section_breaks(bath,stationCoords)
%
% Divide cruise data into cross-shore sections
%
% Last Saved Time-stamp: <Fri 2016-04-15 10:18:49 Eastern Daylight Time gramer>

  coastCoords = get_isobath(bath,0);
  lat = stationCoords(2,:);
  dx = distance_wgs84(lat,lon,coastCoords(2,idx)',coastCoords(1,idx)');

  brkix = [0;find(diff(dx)<0);length(dx)]+1;
  brkln = diff(brkix)-1;

return;
