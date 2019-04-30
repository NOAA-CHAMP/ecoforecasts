1;

% Fix H as needed
if ( size(h,1) ~= numel(lat) ); disp('Transposing fields...'); h = h'; end;
if ( nanmin(h(:)) > -eps ); disp('Inverting H sign...'); h = -h; end;

clear DBG
% DBG.point = [1295,3586]; % BNPON
% DBG.point = [1297,3574]; % DBG.coords = [-80.1764,25.3641]; %BNPMI
if ( ~exist('stnm','var') )
  stnm = 'MLRF1';
  % stnm = 'CCAF1';
  % stnm = 'CNDF1';
  % stnm = 'CGRF1';
  % stnm = 'CPIF1';
end;
DBG = get_station_from_station_name(stnm);
[DBG.laterr,DBG.latix] = min(abs(lat-DBG.lat));
[DBG.lonerr,DBG.lonix] = min(abs(lon-DBG.lon));
if (DBG.laterr>0.01||DBG.lonerr>0.01); error('Coords outside field: %s',DBG.station_name); end;

dx = min(diff(lat))*111e3;

DBG.latrng = 0.020;
DBG.latixen = find((DBG.lat-DBG.latrng)<=lat & lat<=(DBG.lat+DBG.latrng));

DBG.lonrng = 0.020;
DBG.lonixen = find((DBG.lon-DBG.lonrng)<=lon & lon<=(DBG.lon+DBG.lonrng));

DBG.lons = lon(DBG.lonixen);
DBG.lats = lat(DBG.latixen);
DBG.h = h(DBG.latixen,DBG.lonixen);
[DBG.dzdx,DBG.dzdy] = gradient(-DBG.h,dx);
DBG.bet = uv_to_spd(DBG.dzdx,DBG.dzdy);
DBG.hch = hch(DBG.latixen,DBG.lonixen);
DBG.hcrng = hcrng(DBG.latixen,DBG.lonixen);
flatix = find(DBG.bet<0.0010); % Minimum slope for Horizontal Convection?
DBG.dzdx(flatix) = 0; DBG.dzdy(flatix) = 0; DBG.bet(flatix) = 0;

if ( ~exist('fh','var') || ~ishandle(fh) )
  fh=fmg;
else
  fmg(fh);
end;
contourf(DBG.lons,DBG.lats,DBG.h,-[0:1:30]);
caxis([-30,0]);
colormap(viridis);
colorbar('Location','East');
quiver(DBG.lons,DBG.lats,DBG.dzdx,DBG.dzdy,'w','LineWidth',1.5);
% [LONS,LATS]=meshgrid(DBG.lons,DBG.lats); text(LONS(:),LATS(:),num2str(DBG.bet(:),'%.3f')); clear LONS LATS
plot(DBG.lon,DBG.lat,'rp','MarkerSize',16);
titlename([textize(DBG.station_name),': h and \beta']);


fmg;
contourf(DBG.lons,DBG.lats,DBG.hcrng,[0:30:300]);
caxis([0,300]);
colormap(viridis);
colorbar('Location','East');
quiver(DBG.lons,DBG.lats,DBG.dzdx,DBG.dzdy,'w','LineWidth',1.5);
% [LONS,LATS]=meshgrid(DBG.lons,DBG.lats); text(LONS(:),LATS(:),num2str(DBG.bet(:),'%.3f')); clear LONS LATS
plot(DBG.lon,DBG.lat,'rp','MarkerSize',16);
titlename([textize(DBG.station_name),': HC rng [m] and \beta']);
