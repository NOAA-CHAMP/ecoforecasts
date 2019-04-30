function quiver_field(stn,dtrng,ufld,vfld,bathfld,bathcntrs)
%function quiver_field(stn,dtrng,ufld,vfld,bathfld,bathcntrs)
%
% CONTOUR (v.) plot bathymetry, overlaid with QUIVER (v.) plot of NANMEAN of
% U and V vector components in fields STN.(UFLD) and STN.(VFLD).  All fields
% STN.(UFLD),STN.(VFLD),STN.(BATHFLD) must have subfields .lon,.lat,.field.
% If BATHFLD is 'none', no bathymetry is plotted.  BATHCNTRS lists isobaths
% to contour (DEFAULT: [-2,-5,-10,-20,-30,-80]).  If DTRNG is a long vector
% of DATENUM, plot result of INTERSECT_DATES(DTRNG,STN.(UFLD).date); but if
% NUMEL(DTRNG)==1, plot field nearest that date, and if NUMEL(DTRNG)==2, plot
% all dates between those two dates, incl.; if ISA(DTRNG,'function_handle'),
% plot result of DTRNG(STN.(UFLD)); if 'default', plot mean of *all* dates.
%
% Last Saved Time-stamp: <Wed 2011-09-07 14:45:36  lew.gramer>

  if ( ~exist('dtrng','var') || strcmpi(dtrng,'default') )
    dtrng = [-Inf,+Inf];
  end;
  if ( ~exist('ufld','var') || isempty(ufld) )
    ufld = 'fkeys_hycom_u_field';
  end;
  if ( ~exist('vfld','var') || isempty(vfld) )
    vfld = 'fkeys_hycom_v_field';
  end;
  if ( ~exist('bathfld','var') || isempty(bathfld) )
    bathfld = 'ngdc_92m_bathy';
  end;
  if ( ~exist('bathcntrs','var') || isempty(bathcntrs) )
    bathcntrs = [-2,-5,-10,-20,-30,-80];
  end;

  if ( ~isfield(stn,ufld) || ~isfield(stn,vfld) )
    error('First arg STN must have fields %s and %s',ufld,vfld);
  end;

  if ( isempty(dtrng) )
    error('Empty DTRNG specified');
  elseif ( isa(dtrng,'function_handle') )
    dtix = dtrng(stn.(ufld));
    dtstr = [ '@' upper(char(dtrng)) ];
  elseif ( numel(dtrng) == 1 )
    [ig,dtix]=min(abs(dtrng-stn.(ufld).date));
    dtstr = datestr(stn.(ufld).date(dtix));
  elseif ( numel(dtrng) == 2 )
    dtix=find(min(dtrng)<=stn.(ufld).date & stn.(ufld).date<=max(dtrng));
    dtstr = [datestr(stn.(ufld).date(dtix(1))) ' - ' datestr(stn.(ufld).date(dtix(end)))];
  else
    [ig,dtix]=intersect_dates(dtrng,stn.(ufld).date);
    dtstr = [datestr(stn.(ufld).date(dtix(1))) ' - ' datestr(stn.(ufld).date(dtix(end)))];
  end;

  if ( isempty(dtix) )
    error('No dates matching DTRNG found in STN.(UFLD).date');
  end;

  dtstr = [dtstr ' (N=' num2str(numel(dtix)) ')'];


  minlon = min(stn.(ufld).lon(:));
  maxlon = max(stn.(ufld).lon(:));
  dlon = min(diff(unique(stn.(ufld).lon(:))));
  minlat = min(stn.(ufld).lat(:));
  maxlat = max(stn.(ufld).lat(:));
  dlat = min(diff(unique(stn.(ufld).lat(:))));

  [LON,LAT] = meshgrid(stn.(ufld).lon,stn.(ufld).lat);
  u = squeeze(nanmean(stn.(ufld).field(dtix,:,:),1));
  v = squeeze(nanmean(stn.(vfld).field(dtix,:,:),1));


  % Do plotting

  if ( strcmpi(bathfld,'default') )
    map_freef([minlon-2*dlon,maxlon+2*dlon,minlat-2*dlat,maxlat+2*dlat], bathcntrs);
  elseif ( isfield(stn,bathfld) )
    if ( ndims(stn.(bathfld).field) == 2 )
      [cs,ch] = contour(stn.(bathfld).lon,stn.(bathfld).lat,stn.(bathfld).field,bathcntrs);
    else
      [cs,ch] = contour(stn.(bathfld).lon,stn.(bathfld).lat,squeeze(nanmean(stn.(bathfld).field,1)),bathcntrs);
    end;
    % legend(ch,num2str(bathcntrs(:)),'Location','SouthEast');
    clabel(cs,ch);
  end;

  legendlon = minlon-dlon;
  legendlat = maxlat+dlat;

  qh = quiver( [LON(:);legendlon],...
               [LAT(:);legendlat],...
               [u(:);1.0],[v(:);0.0] );

  qh_children=get(qh,'Children');
  xdat=get(qh_children(3),'XData'); qh_legend_x = xdat(end);
  ydat=get(qh_children(3),'YData'); qh_legend_y = ydat(end);
  text(qh_legend_x,qh_legend_y-dlat,'1.0 m/s', ...
       'HorizontalAlignment','Left', ...
       'VerticalAlignment','Bottom', ...
       'BackgroundColor',[1,1,1]);

  yasp = 1;
  if ( isfield(stn,'lon') && isfield(stn,'lat') )
    plot(stn.lon,stn.lat,'kp', 'MarkerSize',18);
    yasp = cosd(stn.lat);
  end;

  daspect([1,yasp,1]);
  axis([minlon-2*dlon,maxlon+2*dlon,minlat-2*dlat,maxlat+2*dlat]);
  ttlstr = [upper(stn.station_name) ' ' ufld '/v ' dtstr];
  titlename(strrep(ttlstr,'_','\_'));

return;
