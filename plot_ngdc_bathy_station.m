function stn = plot_ngdc_bathy_station(stn,my_contours,rad,doCoast,contourFun)
%function stn = plot_ngdc_bathy_station(stn,my_contours,rad,doCoast,contourFun)
%
% Extract NGDC 3-arcsecond bathymetry for RAD km (DEFAULT: 40) surrounding
% STN (struct or name string): plot filled contours at depths in optional
% arg MY_CONTOURS (DEFAULT: -[0:4:160]). If optional DOCOAST True (DEFAULT),
% also plot high-resolution coast-line (if STN in SOFLA_COAST.DAT region).
%
% GLOBALS: SOFLA_COAST
% CALLS: GET_NGDC_BATHY_STATION (Ecoforecasts); CONTOURF, PLOT, FILL. Caller
% may specify an alternative countour-function (e.g., CONTOUR) as last arg.
%
% Last Saved Time-stamp: <Sat 2013-10-19 23:28:43 Eastern Daylight Time gramer>


  global sofla_coast;

  if ( ~exist('my_contours','var') || isempty(my_contours) )
    my_contours = -[ 0:4:160 ];
  end;
  if ( ~exist('rad','var') || isempty(rad) )
    rad = [];
  end;
  if ( ~exist('doCoast','var') || isempty(doCoast) )
    doCoast = true;
  end;
  if ( ~exist('contourFun','var') || isempty(contourFun) )
    contourFun = @contourf;
  end;

  if ( ~isfield(stn,'ngdc_92m_bathy') )
    disp('Extracting station bathymetry');
    stn = get_ngdc_bathy_station(stn,rad,[],false);
  end;

  figure;
  hold on;

  % For most coral reef areas, "bathymetry" above MHHW is just a nuisance
  %contourf(stn.ngdc_92m_bathy.lon,stn.ngdc_92m_bathy.lat,stn.ngdc_92m_bathy.field,my_contours);
  feval(contourFun,stn.ngdc_92m_bathy.lon,stn.ngdc_92m_bathy.lat,stn.ngdc_92m_bathy.field,my_contours);
  maxigraph;
  colorbar;
  % Plot station location as a white pentagram
  if ( isfield(stn,'lon') && isfield(stn,'lat') )
    stnlon = stn.lon;
    stnlat = stn.lat;
    plot(stnlon,stnlat,'wp', 'MarkerEdgeColor','black', 'MarkerFaceColor','white');
  else
    stnlon = mean(stn.ngdc_92m_bathy.lon(:));
    stnlat = mean(stn.ngdc_92m_bathy.lat(:));
    plot(stnlon,stnlat,'wo', 'MarkerEdgeColor','black', 'MarkerFaceColor','white');
  end;

  if ( doCoast )
    if ( ~exist('sofla_coast', 'var') || isempty(sofla_coast) )
      disp('Reloading coastline');
      % sofla_coast = load('sofla_coast_low.dat');
      % sofla_coast = load('sofla_coast_medium.dat');
      sofla_coast = load('sofla_coast.dat');
    end;
    %ch=fill(sofla_coast(:,1), sofla_coast(:,2), [0.4 0.3 0.2], 'LineWidth',5);
    ch=fill(sofla_coast(:,1), sofla_coast(:,2), [0.0 0.0 0.0], 'LineWidth',2);
    % In datamode, user is interested in bathymetry, not coastline
    set(ch,'HitTest','off')
  end;

  grid on;
  titlename(['NGDC bathymetry surrounding ' upper(stn.station_name)]);

  if ( nargout < 1 )
    stn = []; clear stn;
  end;

return;
