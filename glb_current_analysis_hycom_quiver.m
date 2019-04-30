1;

%nc = mDataset('d:/thesis/data/hycom/GOM2/archv.2014_001_12.nc');
%nc = mDataset('d:/thesis/data/hycom/eFKEYS_01-1-7-2012/eFKEYS_archv.2012_001_06_3zuvwts.nc');

%nc = mDataset('http://tds.hycom.org/thredds/dodsC/glb_current_analysis');
nc = mDataset('http://tds.hycom.org/thredds/dodsC/GLBa0.08/expt_91.2');

lats = cast(nc{'Latitude'}(:),'double');
lons = cast(nc{'Longitude'}(:),'double');
deps = cast(nc{'Depth'}(:),'double');
u = cast(nc{'u'}(1,1:3,:,:),'double');
v = cast(nc{'v'}(1,1:3,:,:),'double');
close(nc); clear nc

jx=1059:1084; 
ix=1042:1061; 
fmg;
quiver(squeeze(u(1,jx,ix)),squeeze(v(1,jx,ix)));
quiver(squeeze(u(2,jx,ix)),squeeze(v(2,jx,ix)));

fmg;
quiver([0,0,0]',-deps(1:3),squeeze(u(1:3,1000,1000)),squeeze(v(1:3,1000,1000)));
view(0,90)

u1 = u(1,1:1000,1000:2001); u2=u(2,1:1000,1000:2001);
scatter_fit(u1(:),u2(:),[],[],[],[],[],true)
v1 = v(1,1:1000,1000:2001); v2=v(2,1:1000,1000:2001);
scatter_fit(v1(:),v2(:),[],[],[],[],[],true)
