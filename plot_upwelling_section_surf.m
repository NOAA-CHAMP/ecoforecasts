1;
%% SCRIPT plot_upwelling_section_surf.m:
%

if ( ~exist('doPrint','var') )
  doPrint = false;
end;
if ( ~exist('fontsz','var') )
  fontsz = 24;
end;

if ( ~exist('fld','var') )
  %fld = 'seatemp';	fldstr = 'Sea temperature [\circC]';
  fld = 'N';		fldstr = '[NO3+NO2] [\mumol]';
end;

%fmg; plot(sfomc.ne_buoy.N.date,sfomc.ne_buoy.N.prof(:,end-3:end),sfomc.sw_buoy.N.date,sfomc.sw_buoy.N.prof(:,end),sfomc.nw_w_btm.N.date,sfomc.nw_w_btm.N.prof(:,end)); xlim(datenum(1999,[7,7],[27,29])); datetick3; ylim([-0.5,13]); legend(string([sfomc.ne_buoy.N.depths(end-3:end),sfomc.sw_buoy.N.depths(end),sfomc.nw_w_btm.N.depths(end)])); set(get(gca,'child'),'linew',2);

if ( strcmp(fld,'seatemp') )
  dt=datenum(1999,7,28,14,46,55);
  swdt = datenum(1999,7,27,16,37,04);
  [ig,fwydix]=min(abs(fwyf1.ndbc_air_t.date-dt));
  endpt = fwyf1.ndbc_air_t.data(fwydix);
elseif ( strcmp(fld,'N') )
  dt = datenum(1999,7,28,10,40,00);
  swdt = datenum(1999,7,27,17,50,00);
  endpt = 0;
end;
[ig,nedix]=min(abs(sfomc.ne_buoy.(fld).date-dt));
[ig,swdix]=min(abs(sfomc.sw_buoy.(fld).date-swdt));
[ig,nwdix]=min(abs(sfomc.nw_w_btm.(fld).date-dt));

nezix = [1,6,7,8];

[nw_dx,nw_az,nw_nearestCoords] = get_contour_distance(sfomc.c_buoy.ngdc_hires_bathy,sfomc.nw_w_btm);

[c_dx,c_az,c_nearestCoords] = get_contour_distance(sfomc.c_buoy.ngdc_hires_bathy,sfomc.c_buoy);

[ne_dx,ne_az,ne_nearestCoords] = get_contour_distance(sfomc.c_buoy.ngdc_hires_bathy,sfomc.ne_buoy);
ne_dx = 3.20;

ds = [...
    0,...
    nw_dx, repmat(nw_dx,[1,numel(sfomc.nw_w_btm.(fld).depths)]), nw_dx,...
    c_dx, repmat(c_dx,[1,numel(sfomc.sw_buoy.(fld).depths)]), c_dx,...
    ne_dx, repmat(ne_dx,[1,numel(sfomc.ne_buoy.(fld).depths(nezix))]), ne_dx,...
     ];
zs = [ ...
    0;...
    0; sfomc.nw_w_btm.(fld).depths'; -20;...
    0; sfomc.sw_buoy.(fld).depths'; -50;...
    0; sfomc.ne_buoy.(fld).depths(nezix)'; -50; ...
     ];
Ts = [ ...
    endpt,...
    endpt, sfomc.nw_w_btm.(fld).prof(nwdix,:), sfomc.sw_buoy.(fld).prof(swdix,end),...
    endpt, sfomc.sw_buoy.(fld).prof(swdix,:), sfomc.ne_buoy.(fld).prof(nedix,nezix(end)),...
    endpt, sfomc.ne_buoy.(fld).prof(nedix,nezix), sfomc.ne_buoy.(fld).prof(nedix,nezix(end)),...
     ];

xq = 0:0.01:ne_dx;
yq = -50.15:0.05:0;
[Xq,Yq] = meshgrid(xq,yq);

Vq = griddata(ds,zs,Ts,Xq,Yq);

brown = [0.4,0.0,0.0];

fmg;

[res,sfomc.c_buoy,lh]=plot_bathy_transect({sfomc.c_buoy,c_nearestCoords},ne_dx+0.100,c_az+180);
badix = find(isnan(res.depths));
res.lon(badix) = []; res.lat(badix) = []; res.field(badix) = []; res.range(badix) = []; res.depths(badix) = []; 
set(lh,'LineWidth',9,'Color',brown);
legend off
%xlim([nw_dx-0.100,ne_dx+0.100]); ylim([-53,3]);
xlim([1.00,3.49]); ylim([-53,3]);
axis(axis);

BWx = [res.range,res.range(end),res.range(1),res.range(1)];
BWy = [res.depths,min(res.depths),min(res.depths),res.depths(1)];

%BWix=[];
if ( ~exist('BWix','var') )
  tic, BWix = find(inpolygon(Xq,Yq,BWx,BWy)); toc,
end;
Vq(BWix) = nan;
sh = surf(Xq,Yq,Vq);
ws = warning('OFF','MATLAB:handle_graphics:Patch:NumColorsNotEqualNumVertsException');
shading interp
warning(ws); clear ws
set(sh,'FaceAlpha',0.8);
flh = fill(BWx,BWy,brown);

%cbh=colorbar('Location','EastOutside','FontSize',fontsz);
cbh=colorbar('Location','East');
%cbpos=get(cbh,'pos'); set(cbh,'pos',[cbpos(1),0.3,cbpos(3),0.6]);
xlabel(cbh,fldstr,'FontSize',fontsz,'Color','k');
plot3(nw_dx,-sfomc.nw_w_btm.depth,1,'rp','MarkerSize',24,'LineWidth',2);
plot3(c_dx,-sfomc.c_buoy.depth,1,'rp','MarkerSize',24,'LineWidth',2);
plot3(ne_dx,-sfomc.ne_buoy.depth,1,'rp','MarkerSize',24,'LineWidth',2);

xlabel('Distance from shore [km]');
set(gca,'FontSize',fontsz);

set(gca,'CLim',[0,3]);

titlename(datestr(dt));
if ( doPrint )
  print('-dpng',fullfile(get_coral_path,'CRCP','Upwelling','CoRIS',[mfilename,'_',fld,'.png']));
end;
