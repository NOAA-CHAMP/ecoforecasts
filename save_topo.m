1;

load ETOPO1_Bed_g_gmt4.mat;

latix = find(-29 <= V_lat & V_lat <= -6);
lonix = find(-180 <= V_lon & V_lon <= -159);
lats = V_lat(latix);
lons = V_lon(lonix);
topo = V_z(latix,lonix);
region = 'asam';
% .MAT on phodnet MATLAB7 is not readable on PC MATLAB7! Idiots...
% save([region '.topo.phodnet.mat'],'lats','lons','topo');
save([region '.lats.dat'],'-ascii','lats');
save([region '.lons.dat'],'-ascii','lons');
save([region '.topo.dat'],'-ascii','topo');
size(topo),
% ans =
%         1381        1261


latix = find(8 <= V_lat & V_lat <= 31);
lonix = find(-97 <= V_lon & V_lon <= -74);
lats = V_lat(latix);
lons = V_lon(lonix);
topo = V_z(latix,lonix);
region = 'freef';
% .MAT on phodnet MATLAB7 is not readable on PC MATLAB7! Idiots...
% save([region '.topo.phodnet.mat'],'lats','lons','topo');
save([region '.lats.dat'],'-ascii','lats');
save([region '.lons.dat'],'-ascii','lons');
save([region '.topo.dat'],'-ascii','topo');
size(topo),
% ans =
%         1381        1381

latix = find(1 <= V_lat & V_lat <= 24);
lonix = find(-79 <= V_lon & V_lon <= -56);
lats = V_lat(latix);
lons = V_lon(lonix);
topo = V_z(latix,lonix);
region = 'ecarib';
% .MAT on phodnet MATLAB7 is not readable on PC MATLAB7! Idiots...
% save([region '.topo.phodnet.mat'],'lats','lons','topo');
save([region '.lats.dat'],'-ascii','lats');
save([region '.lons.dat'],'-ascii','lons');
save([region '.topo.dat'],'-ascii','topo');
size(topo),
% ans =
%         1381        1381

latix = find(-32 <= V_lat & V_lat <= -9);
lonix = find(135 <= V_lon & V_lon <= 158);
lats = V_lat(latix);
lons = V_lon(lonix);
topo = V_z(latix,lonix);
region = 'gbr';
% .MAT on phodnet MATLAB7 is not readable on PC MATLAB7! Idiots...
% save([region '.topo.phodnet.mat'],'lats','lons','topo');
save([region '.lats.dat'],'-ascii','lats');
save([region '.lons.dat'],'-ascii','lons');
save([region '.topo.dat'],'-ascii','topo');
size(topo),
% ans =
%         1381        1381
