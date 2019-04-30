1;

fmg;
% Universal Transverse Mercator projection
axesm utm;
% UTM zone for areas of Florida
setm(gca,'zone','17r');
% Draw the map grid and label it
setm(gca,'grid','on','meridianlabel','on','parallellabel','on')
load coast
plotm(lat,long);

%states = shaperead('usastatehi', 'UseGeoCoords', true);
%framem
%faceColors = makesymbolspec('Polygon', {'INDEX', [1 numel(states)], 'FaceColor', polcmap(numel(states))});
%geoshow(states,'DisplayType', 'polygon', 'SymbolSpec', faceColors)

[latlim, lonlim] = utmzone('17r')
%[x,y] = mfwdtran(latlim, lonlim)

%z1 = utmzone(p1)
% Obtain the suggested ellipsoid vector and name for this zone
[ellipsoid,estr] = utmgeoid('17r');
% Set up the UTM coordinate system based on this information
utmstruct = defaultm('utm'); 
utmstruct.zone = '17r'; 
utmstruct.geoid = ellipsoid; 
utmstruct = defaultm(utmstruct);
[x,y] = mfwdtran(utmstruct,latlim,lonlim)

[lat,lon] = utm2deg(949622.38,589331.52,'17 R')
[lat,lon] = minvtran(utmstruct,949622.38,589331.52)

%% Call UTMZONE recursively to obtain the limits of the UTM zone within which a point location falls
%[zonelats zonelons] = utmzone(utmzone(40.7, -74.0))

