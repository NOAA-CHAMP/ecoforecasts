1;

bda1.lon = -64.69243; bda1.lat = 32.3598833;
bda2.lon = -64.65738; bda2.lat = 32.3367167;
bda1 = plot_hires_bathymetry(bda1,-[0:2:80],[7e3,7e3],[],[],[],[],true);

x=[]; x.bda1=bda1; x.bda2=bda2;
[x,bda1.ngdc_hires_bathy] = find_ngdc_slope_sites(x,bda1.ngdc_hires_bathy,5,{@nanmean,4,4,5});
x.depths(:)',x.betas(:)',
%ans =
%  -10.7929  -16.7136
%ans =
%    0.0053    0.1234
