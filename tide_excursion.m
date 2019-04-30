1;
%% SCRIPT TIDE_EXCURSION.m: Plot bathymetry and a "tidal excursion" at site
%
% Plots motion of a water parcel under tide solution from STATION_TMD_TIDE,
% using tidal velocity *at* station STNM (DEFAULT 'looe1'). Uses TIDEMODEL
% (DEFAULT 'tmd' for Gulf of Mexico) to plot a tide-day (point every DHR h,
% DEFAULT 1) every DDY days (DEFAULT 4) for a period of NDYS (DEFAULT 29)
% from day DY0 (DEFAULT DATENUM(2010,8,1)). Print if DOPRINT (DEFAULT False). 
% 
% Last Saved Time-stamp: <Fri 2016-11-11 15:07:36 Eastern Standard Time gramer>

if ( ~exist('doPrint','var') || isempty(doPrint) )
  doPrint = false;
end;
if ( ~exist('stnm','var') || isempty(stnm) )
  % stnm = 'fwyf1';
  stnm = 'mlrf1';
  % stnm = 'lonf1';
  % stnm = 'smkf1';
  % stnm = 'looe1';
  % stnm = 'sanf1';
  % stnm = 'dryf1';
end;
if ( ~exist('tidemodel','var') || isempty(tidemodel) )
  tidemodel = 'tmd';
  % tidemodel = 'tpxo';
  % % tidemodel = 'ao';
end;
if ( ~exist('ddy','var') || isempty(ddy) )
  %ddy = 1;
  ddy = 4;
end;
if ( ~exist('ndys','var') || isempty(ndys) )
  % ndys = 4;
  % ndys = 2;
  % ndys = 4;
  ndys = 29;
end;
if ( ~exist('dhr','var') || isempty(dhr) )
  dhr = 1;
end;
if ( ~exist('dy0','var') || isempty(dy0) )
  % dy0 = -inf;
  dy0 = datenum(2010,8,1);
end;


%% Get data

stn = get_station_from_station_name(stnm);
stn = station_tmd_tide(stn);

ufld = [tidemodel,'_tide_u'];
vfld = [tidemodel,'_tide_v'];

dy0ix = find(stn.(ufld).date >= dy0,1) - 1;
ndys = ndys - mod(ndys,ddy) + 1;
nhrs = 24*ndys;

dlon = distance_wgs84(stn.lat,stn.lon,stn.lat,stn.lon+1)*1e3;
dlat = distance_wgs84(stn.lat,stn.lon,stn.lat+1,stn.lon)*1e3;

stn = plot_hires_bathymetry(stn,-[0:4:80],[40e3,40e3],true,@contour,[],[],false);
plot(stn.lon,stn.lat,'rp','MarkerSize',16);
titlename([upper(tidemodel),' ',num2str(ndys),'d tidal ellipse for ',upper(stnm)]);

clrs = {'k','b',[0,.5,0],'r','m','c','y',[.5,.5,.5]};
nclrs = numel(clrs);
for dix = 1:ddy:ndys
  x = stn.lon;
  y = stn.lat;
  for hix = 1:dhr:24
    ix = ((dix-1)*24) + hix;
    dx = stn.(ufld).data(dy0ix+ix) * 3600 / dlon;
    dy = stn.(vfld).data(dy0ix+ix) * 3600 / dlat;
    x(hix+1) = x(hix) + dx;
    y(hix+1) = y(hix) + dy;
  end;
  plot( x,y,'-','Color',clrs{mod(dix-1,nclrs)+1} );
end;

if ( doPrint )
  print('-dpng',...
        fullfile(get_ecoforecasts_path('figs'),...
                 [lower(stnm),'_',tidemodel,'_',datestr(dy0,'yyyymmdd'),'.png']));
end;
