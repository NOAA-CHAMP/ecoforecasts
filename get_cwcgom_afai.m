1;

nc = mDataset('http://cwcgom.aoml.noaa.gov/thredds/dodsC/AFAI/USFAFAI3D.nc');
nj_info(nc),
t
t = nc{'time'}(:);
xt = (t/3600/24) + datenum(1970,01,01);
datestr(xt)
nj_info(nc),
AFAI = cast(nc{'AFAI'}(1,:,:),'double');
fmg; contourf(AFAI); colorbar;
close(nc); clear nc
