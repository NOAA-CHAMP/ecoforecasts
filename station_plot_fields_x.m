function fh = station_plot_fields_x(stn,fldsnm,dts)
%function fh = station_plot_fields_x(stn,fldsnm,dts)
%
% Plot a succession of SURF (v.) subplots (v.) from time series of 2-d fields
% STN.(FLDSNM). The following sub-fields are expected, or used if available:
%   STN.(FLDSNM).date :  REQUIRED Tx1 vector of DATENUMs for 2-d fields
%   STN.(FLDSNM).field : REQUIRED TxNxM matrix of successive 2-d fields
%   STN.(FLDSNM).lon :   OPTIONAL Mx1 vector of longitudes
%   STN.(FLDSNM).lat :   OPTIONAL Nx1 vector of latitudes
%   STN.lon :            OPTIONAL scalar longitude of station site
%   STN.lat :            OPTIONAL scalar latitude of station site
%   STN.station_name :   OPTIONAL string identifying the station
%
% NOTE: Surface-plots field values *relative to the station position*, if given.
% CALLS: INTERSECT_DATES,INTERPN,SUBPLOT,SURF,PLOT3,MAXIGRAPH,SUPTITLE.
%
% EXAMPLE:
%  >> % Plot all 1km AVHRR SST fields around Molasses Reef between May 15 and 30, 2007
%  >> stn.station_name = 'mlrf1';
%  >> [stn.lon,stn.lat] = get_station_coords(stn.station_name);
%  >> stn = get_avhrr_weekly_field(stn);
%  >> station_plot_fields(stn,'avhrr_weekly_sst',datenum(2007,5,15):(1/24):datenum(2007,5,30));
%
% Last Saved Time-stamp: <Wed 2011-07-06 06:24:37  lew.gramer>

  fh = [];

  [ig,dix] = intersect_dates(dts,stn.(fldsnm).date);

  flddts = stn.(fldsnm).date(dix);
  ndts = length(flddts);
  %DEBUG:  ndts,
  if ( ndts > 25 )
    error('Too many intersecting dates (%g > 25)!', ndts);
  end;

  nrows = floor(sqrt(ndts));
  ncols = ceil(ndts/nrows);
  %DEBUG:  nrows, ncols,

  nx = size(stn.(fldsnm).field,3);
  ny = size(stn.(fldsnm).field,2);

  stnlon = [];
  stnlat = [];
  if ( isfield(stn,'lon') )
    stnlon = stn.lon;
    stnlat = stn.lat;
  end;
  %DEBUG:  stnlon, stnlat,

  if ( ~isfield(stn.(fldsnm),'lon') )
    lons = 1:nx;
    lats = 1:ny;
  else
    lons = stn.(fldsnm).lon;
    lats = stn.(fldsnm).lat;
  end;

  %DEBUG:  size(lons), size(lats),
  fld = stn.(fldsnm).field(dix,:,:);
  stnsst = repmat(0,size(flddts));
  if ( ~isempty(stnlon) )
    % stnsst = interpn(flddts,lons,lats,fld,flddts,stnlon,stnlat,'spline');
    % stnsst = squeeze(fld(:,round(nx/2),round(ny/2)));
    stnsst = interpn(flddts,lats,lons,fld,flddts,stnlat,stnlon,'spline');
  end;
  for ix = 1:length(dix)
    fld(ix,:,:) = fld(ix,:,:) - stnsst(ix);
  end;
  xmin = min(lons(:));  xmax = max(lons(:));
  ymin = min(lats(:));  ymax = max(lats(:));
  cmin = min(fld(isfinite(fld(:))));
  cmax = max(fld(isfinite(fld(:))));
  cmin = -max(abs(cmin),abs(cmax));
  cmax = max(abs(cmin),abs(cmax));
  %DEBUG:  cmin, cmax,
  cmin = -0.5; cmax = 0.5;


  centerx = interp1(lons,1:nx,stnlon);
  centery = interp1(lats,1:ny,stnlat);


  fh = figure;
  for ix = 1:length(dix)
    subplot(nrows,ncols,ix);
    surf(1:nx,1:ny,squeeze(fld(ix,:,:)));
    set(gca,'ydir','reverse');
    % surf(lons,lats,squeeze(fld(ix,:,:)));
    % contourf(lons,lats,squeeze(fld(ix,:,:)));
    % xlim([xmin xmax]);
    xlim([1 nx]);
    % ylim([ymin ymax]);
    ylim([1 ny]);
    caxis([cmin cmax]);
    view(2);
    title(datestr(flddts(ix)));
    if ( ~isempty(stnlon) )
      hold on;
      plot3(centerx,centery,cmax+10,'k*');
      [ig,uix] = min(abs(stn.quasi_eulerian_u.date - flddts(ix)));

      euu = stn.quasi_eulerian_u.data(uix) .* (3600/1e3);
      euv = stn.quasi_eulerian_v.data(uix) .* (3600/1e3);

      % U and V are in target direction: subtract for source point!
      sourcex = centerx - euu;
      %sourcey = centery - euv;
      % BUT rows are normally reversed in SST fields, so *add* V! Grumble...
      sourcey = centery + euv;

      % quiver3(stnlon,stnlat,cmax+10,-stn.quasi_eulerian_u.data(uix)*6*3600/111e3,stn.quasi_eulerian_v.data(uix)*6*3600/111e3,0);
      plot3(sourcex,sourcey,cmax+10,'r.');
    end;
  end;
  maxigraph;

  % Create colorbar and super-title
  axes('Units','norm', 'Position',[0.10 0.10 0.88 0.83], 'Visible','off');
  colorbar;
  caxis([cmin cmax]);
  stnm = 'Station';
  if ( isfield(stn,'station_name') )
    stnm = stn.station_name;
  end;
  suptitle(sprintf('%s.%s %s-%s',...
                   stnm,strrep(fldsnm,'_','\_'),...
                   datestr(flddts(1)),datestr(flddts(end))));

return;
