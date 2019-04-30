function [lat,lon,res,bbox] = get_ool_region(rgn)
%function [lat,lon,res,bbox] = get_ool_region(rgn)
%
% Return vectors of latitudes and longitudes associated with University of
% South Florida (USF) Optical Oceanography Lab (OOL) dataset region RGN,
% which may be, e.g., one of 'SE_FL', 'CNMI', 'SAMOA', 'ECARIB', 'GCOOS',
% 'BERMUDA', etc. Also returns a scalar of nominal dataset resolution (in
% [m]) and the dataset's bounding box (as used, e.g., by BBOXINSIDE, v.)
%
% CALLERS: coral/CRCP/Sediment/{extract_ool_250m,local_extract_ool_250m}.m
%
% Last Saved Time-stamp: <Sun 2018-08-19 00:28:10 Eastern Daylight Time gramer>

  if ( ~exist('rgn','var') || ~ischar(rgn) )
    error('First arg should be a REGION character string! See http://optics.marine.usf.edu');
  end;

  switch (upper(rgn)),
   case 'SE_FL',
    % 880x572
    % SE_FL: The South East Florida Shelf region is designed to show the lower South
    % East part of the Florida Shelf and is bounded by these coordinates: 27.5°N
    % 25.5°N 79°W and 80.3°W.
    lat = linspace(25.5,27.5,880);
    lon = linspace(-80.3,-79.0,572);
   case 'CNMI',
    % 704x528
    % CNMI: The Commonwealth of the Northern Mariana Islands region is bounded by
    % these coordinates: 13.9°N - 15.5°N, 144.9°E - 146.1°E.
    lat = linspace(13.90,15.50,704);
    lon = linspace(144.9,146.1,528);
   case 'SAMOA',
    % 440x792
    % SAMOA: Samoa Islands region is bounded by these coordinates: 14.7°S -
    % 13.7°S, 171°W - 169.2°W.
    lat = linspace(-14.70,-13.70,440);
    lon = linspace(-171.0,-169.2,792);
   case 'ECARIB',
    % 1430x1650 - 23°N 10°N 60°W and 75°W
    lat = linspace(10,23,1430);
    lon = linspace(-75,-60,1650);
   case 'GCOOS',
    % 1430x2090 - 31°N 18°N 79°W and 98°W
    lat = linspace(18,31,1430);
    lon = linspace(-98,-79,2090);
   case 'BERMUDA',
    %Bermuda is an area bounded within these coordinates: 37.0°N 27.0°N 59.0°W and 69.0°W.
    % 1100x1100 - 31°N 18°N 79°W and 98°W
    lat = linspace(27,37,1100);
    lon = linspace(-69,-59,1100);
   case 'FLKEYS',
    %The Florida Keys region is designed to show the lower West part of Florida including
    % Florida Bay and is bounded by these coordinates: 880x1320: 26°N 24°N 80°W and 83°W.
    lat = linspace(24,26,880);
    lon = linspace(-83,-80,1320);
   otherwise,
    error('Do not yet know coordinates for region "%s"! See http://optics.marine.usf.edu',rgn);
  end;
  
  % All images use cylindrical equidistant projection, 250/1000 m resolution. 

  %res = [min(diff(lon)),min(diff(lat))]*111e3;
  res = distance_wgs84(lat(1),lon(1),lat(2),lon(1))*1e3;

  bbox = [min(lon),max(lon),min(lat),max(lat)];

return;
