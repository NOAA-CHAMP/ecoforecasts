function fh = station_plot_fields(stn,fldnm,dts,isAnom)
%function fh = station_plot_fields(stn,fldnm,dts,isAnom)
%
% Plot a succession of SURF (v.) subplots (v.) from time series of 2-d fields
% STN.(FLDNM). The following sub-fields are expected, or used if available:
%   STN.(FLDNM).date :  REQUIRED Tx1 vector of DATENUMs for 2-d fields
%   STN.(FLDNM).field : REQUIRED TxNxM matrix of successive 2-d fields
%   STN.(FLDNM).lon :   OPTIONAL Mx1 vector of longitudes
%   STN.(FLDNM).lat :   OPTIONAL Nx1 vector of latitudes
%   STN.lon :           OPTIONAL scalar longitude of station site
%   STN.lat :           OPTIONAL scalar latitude of station site
%   STN.station_name :  OPTIONAL string identifying the station
%
% If otional ISANOM is true, sets CAXIS limits appropriate for an anomaly field.
%
% NOTE: Surface-plots field values *relative to the station position*, if given.
%
% CALLS: INTERSECT_DATES,INTERPN,SUBPLOT,SURF,PLOT3,MAXIGRAPH,SUPTITLE.
%
% EXAMPLE:
%  >> % Plot all 1km AVHRR SST fields around Molasses Reef between May 15 and 30, 2007
%  >> stn.station_name = 'mlrf1';
%  >> [stn.lon,stn.lat] = get_station_coords(stn.station_name);
%  >> stn = get_avhrr_weekly_field(stn);
%  >> station_plot_fields(stn,'avhrr_weekly_sst_field',datenum(2007,5,15):(1/24):datenum(2007,5,30));
%
% Last Saved Time-stamp: <Wed 2011-07-06 06:24:50  lew.gramer>

  if ( ~exist('isAnom','var') || isempty(isAnom) )
      isAnom = false;
  end;

  fh = [];

  [ig,dix] = intersect_dates(dts,stn.(fldnm).date);

  flddts = stn.(fldnm).date(dix);
  ndts = length(flddts);
  %DEBUG:  ndts,
  if ( ndts > 25 )
    error('Too many intersecting dates (%g > 25)!', ndts);
  end;

  nrows = floor(sqrt(ndts));
  ncols = ceil(ndts/nrows);
  %DEBUG:  nrows, ncols,

  nx = size(stn.(fldnm).field,3);
  ny = size(stn.(fldnm).field,2);

  stnlon = [];
  stnlat = [];
  if ( isfield(stn,'lon') )
    stnlon = stn.lon;
    stnlat = stn.lat;
  end;
  %DEBUG:  stnlon, stnlat,

  if ( ~isfield(stn.(fldnm),'lon') )
    lons = 1:nx;
    lats = 1:ny;
  else
    lons = stn.(fldnm).lon;
    lats = stn.(fldnm).lat;
  end;

  %DEBUG:  size(lons), size(lats),
  fld = stn.(fldnm).field(dix,:,:);
  stnsst = repmat(0,size(flddts));
  if ( isAnom )
      if ( ~isempty(stnlon) )
          % stnsst = interpn(flddts,lons,lats,fld,flddts,stnlon,stnlat,'spline');
          stnsst = squeeze(fld(:,round(nx/2),round(ny/2)));
      end;
      for ix = 1:length(dix)
          fld(ix,:,:) = fld(ix,:,:) - stnsst(ix);
      end;
      cmin = -max(abs(cmin),abs(cmax));
      cmax = max(abs(cmin),abs(cmax));
      %cmin = -0.5; cmax = 0.5;
  else
      cmin = min(fld(isfinite(fld(:))));
      cmax = max(fld(isfinite(fld(:))));
  end;
  %DEBUG:  cmin, cmax,
  xmin = min(lons(:));  xmax = max(lons(:));
  ymin = min(lats(:));  ymax = max(lats(:));


  fh = figure;
  for ix = 1:length(dix)
    subplot(nrows,ncols,ix);
    surf(lons,lats,squeeze(fld(ix,:,:)));
    xlim([xmin xmax]);
    ylim([ymin ymax]);
    caxis([cmin cmax]);
    view(2);
    title(datestr(flddts(ix)));
    if ( ~isempty(stnlon) )
      hold on;
      plot3(stnlon,stnlat,cmax+10,'k*');
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
                   stnm,strrep(fldnm,'_','\_'),...
                   datestr(flddts(1)),datestr(flddts(end))));

return;
