function [stn,cs,h,clh] = plot_ngdc_bathymetry(stn_or_stnm_or_loc,my_contours,rad,doCoast,contourFun,doLabels,bathfile)
%function [stn,cs,h,clh] = plot_ngdc_bathymetry(stn_or_stnm_or_loc,my_contours,rad,doCoast,contourFun,doLabels,bathfile)
%
% Extract NGDC 3-arcsecond bathymetry for RAD m (DEFAULT: 40e3) surrounding
% STN (struct), STNM (name string) or LOC (lon/lat 2-vector): plot filled
% contours at depths in optional arg MY_CONTOURS (DEFAULT: -[0:4:160]). If
% optional DOCOAST True (DEFAULT), also plot high-res coast-line (if STN in
% SOFLA_COAST.DAT region). SEE older function: PLOT_NGDC_BATHY_STATION.
%
% GLOBALS: SOFLA_COAST
% CALLS: READ_NGDC_BATHYMETRY (Ecoforecasts); CONTOURF, PLOT, FILL. Caller
% may specify CONTOURFUN-ction (DEFAULT: @CONTOURF). If DOLABELS (DEFAULT:
% empty) is logical TRUE, pass *two* return values from @CONTOURFUN to CLABEL
% (v.); if DOLABELS is a numeric vector, then only label those contour levels;
% if nonempty, call also returns [CS,H]=@CONTOURFUN(...) and CLH=CLABEL(...).
% Finally, if DOLABELS is a cell, one element should be a logical for option
% INLINELABELS (DEFAULT: false); if true, call CLABEL(CS)/CLABEL(CS,DOLABELS).
% Optional BATHFILE if present is passed through to READ_NGDC_BATYMETRY.
%
% Last Saved Time-stamp: <Wed 2018-08-15 12:08:09 Eastern Daylight Time gramer>

  global sofla_coast;

  warning('ool_250m:usehighestres',...
          'Higher-resolution bathymetry may be available... See PLOT_HIRES_BATHYMETRY.M');

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
  if ( ~exist('doLabels','var') || isempty(doLabels) )
    doLabels = [];
  end;
  if ( ~exist('bathfile','var') || isempty(bathfile) )
    bathfile = [];
  end;

  inlineLabels = false;
  if ( iscell(doLabels) )
    if ( islogical(doLabels{1}) )
      inlineLabels = doLabels{1};
      doLabels = doLabels{2};
    else
      error('If DOLABELS is a cell, the first element should be a logical INLINELABELS');
    end;
  end;

  if ( ~isfield(stn_or_stnm_or_loc,'ngdc_92m_bathy') )
    disp('Extracting NGDC bathymetry');
    stn = read_ngdc_bathymetry(stn_or_stnm_or_loc,rad,bathfile);
  else
    stn = stn_or_stnm_or_loc;
  end;
  clear stn_or_stnm_or_loc;

  % figure;
  % hold on;
  % maxigraph;
  fmg;

  % For most coral reef areas, "bathymetry" above MHHW is just a nuisance
  if ( isempty(doLabels) || (islogical(doLabels) && ~doLabels) )
    %contourf(stn.ngdc_92m_bathy.lon,stn.ngdc_92m_bathy.lat,stn.ngdc_92m_bathy.field,my_contours);
    feval(contourFun,stn.ngdc_92m_bathy.lon,stn.ngdc_92m_bathy.lat,stn.ngdc_92m_bathy.field,my_contours);
  elseif ( ~inlineLabels )
    %contourf(stn.ngdc_92m_bathy.lon,stn.ngdc_92m_bathy.lat,stn.ngdc_92m_bathy.field,my_contours);
    [cs,h]=feval(contourFun,stn.ngdc_92m_bathy.lon,stn.ngdc_92m_bathy.lat,stn.ngdc_92m_bathy.field,my_contours);
    if ( islogical(doLabels) )
      clh = clabel(cs);
    else
      clh = clabel(cs,doLabels);
    end;
  else
    %contourf(stn.ngdc_92m_bathy.lon,stn.ngdc_92m_bathy.lat,stn.ngdc_92m_bathy.field,my_contours);
    [cs,h]=feval(contourFun,stn.ngdc_92m_bathy.lon,stn.ngdc_92m_bathy.lat,stn.ngdc_92m_bathy.field,my_contours);
    if ( islogical(doLabels) )
      clh = clabel(cs,h);
    else
      clh = clabel(cs,h,doLabels);
    end;
  end;

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
    %fill(sofla_coast(:,1), sofla_coast(:,2), [.4 .3 .2]);
    ch=fill(sofla_coast(:,1), sofla_coast(:,2), [0.0 0.0 0.0], 'LineWidth',2);
    % In datamode, user is interested in bathymetry, not coastline
    set(ch,'HitTest','off')
  end;

  grid on;
  if ( isfield(stn,'station_name') )
    titlename(['NGDC bathymetry surrounding ' strrep(upper(stn.station_name),'_','\_')]);
  else
    titlename(['NGDC bathymetry surrounding site ']);
  end;

  daspect([1,cosd(stnlat),1]);

  if ( nargout < 1 )
    stn = []; clear stn;
  end;

return;
