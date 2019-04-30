function [stn,cs,h,clh,fh,cbh] = plot_hires_bathymetry(stn_or_stnm_or_loc,my_contours,rad,doCoast,contourFun,doLabels,bathfiles,useHighestRes,fh,fld)
%function [stn,cs,h,clh,fh,cbh] = plot_hires_bathymetry(stn_or_stnm_or_loc,my_contours,rad,doCoast,contourFun,doLabels,bathfiles,useHighestRes,fh,fld)
%
% Plot high-resolution bathymetric contour map around a specified location.
%
% Returns optional graphic handles from plotting (described below), and STN
% STRUCT populated with field STN.ngdc_hires_bathy (if not already present).
%
% If STN is a STNM (name string) or LOC ([lon,lat] 2-vector), or a STRUCT
% that is lacking field STN.ngdc_hires_bathy, call READ_HIRES_BATHYMETRY (v.)
% to extract high-resolution bathymetry for RAD m (DEFAULT: 40e3) around STN.
%
% If optional CHAR FLD is some other field name than 'ngdc_hires_bathy', then
% absence of that field in the STRUCT STN passed in results in an ERROR.
%
% Plots filled contours by default (but see CONTOURFUN) at depths in optional
% arg MY_CONTOURS (DEFAULT: -[0:4:160]). If RAD is a 2-vector, plot rectangle
% of RAD(2) (E-W) by RAD(1) (N-S) around site. If optional DOCOAST is True
% (DEFAULT), also plot a heavy coastline contour (0-isobath). Calls FEVAL
% on CONTOURFUN (DEFAULT: @CONTOURF). Places a star at STN coords with PLOT.
%
% If DOLABELS (DEFAULT: []) is logical TRUE, pass *two* return values from
% @CONTOURFUN to CLABEL (v.): by DEFAULT, labels are *inline*. If DOLABELS
% is a numeric vector, then only label those contour levels. If DOLABELS is a
% cell, DOLABELS{1} must be logical (DEFAULT: True) whether to do INLINE
% labels; if DOLABELS{1}, call CLABEL(CS,DOLABELS{2:end}). If DOLABELS is
% nonempty, call also returns [CS,H]=@CONTOURFUN(...) and CLH=CLABEL(...).
% CLABEL default 'FontSize' is 20, 'LabelSpacing' is 4*72 (4 inches). Caller
% modifies these or other options by passing DOLABELS={doInline,{args...}}.
% Also returns CBH=COLORBAR.
%
% Optional BATHFILES is passed through to READ_HIRES_BATYMETRY; optional
% USEHIGHESTRES is also passed to READ_HIRES_BATYMETRY. If FH is a valid
% Figure handle, plot in it. If empty (DEFAULT), return FH = FMG (see).
%
% (See also older functions PLOT_NGDC_BATHYMETRY and READ_NGDC_BATHYMETRY.)
%
% Last Saved Time-stamp: <Wed 2018-10-17 10:56:45 EDT lew.gramer>

  stn=[];
  cs=[];
  h=[];
  clh=[];

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
  if ( ~exist('bathfiles','var') || isempty(bathfiles) )
    bathfiles = [];
  end;
  if ( ~exist('useHighestRes','var') || isempty(useHighestRes) )
    useHighestRes = [];
  end;
  if ( ~exist('fh','var') || isempty(fh) )
    fh = [];
  end;

  default_fld = 'ngdc_hires_bathy';
  if ( ~exist('fld','var') || isempty(fld) )
    fld = default_fld;
  end;

  inlineLabels = true;
  if ( iscell(doLabels) )
    if ( islogical(doLabels{1}) )
      inlineLabels = doLabels{1};
      doLabels = doLabels{2};
    else
      error('If DOLABELS is a cell, the first element should be a logical INLINELABELS');
    end;
  end;

  % Region Of Interest - defaults to whole bathymetry extract
  roi = [];
  % Special case for STRUCT arrays: the first bathymetry we get will fill in
  % empty fields for all of the remaining STRUCTs in the array, so...
  if ( isfield(stn_or_stnm_or_loc,fld) && isempty(stn_or_stnm_or_loc.(fld)) )
    stn_or_stnm_or_loc = rmfield(stn_or_stnm_or_loc,fld);
  end;
  if ( ~isfield(stn_or_stnm_or_loc,fld) )
    if ( ~strcmp(fld,default_fld) )
      error('Field STN.%s not found!',fld);
    else
      disp('Extracting NGDC high-resolution bathymetry');
      [stn,rad] = read_hires_bathymetry(stn_or_stnm_or_loc,rad,bathfiles,useHighestRes);
    end;
  else
    stn = stn_or_stnm_or_loc;
    if ( ~isempty(rad) && isnumeric(rad) )
      roi = true;
    end;
  end;
  clear stn_or_stnm_or_loc;

  if ( isfield(stn,'lon') && isfield(stn,'lat') )
    stnlon = stn.lon;
    stnlat = stn.lat;
  else
    stnlon = mean(stn.(fld).lon(:));
    stnlat = mean(stn.(fld).lat(:));
  end;

  % If we may have been passed a narrower (or wider) radius than bathymetry
  % was initially loaded at, use that as the limit on our map plot (below)
  if ( roi )
    if ( isscalar(rad) )
      rad = [rad,rad];
    end;
    [lats,lons] = reckon_wgs84(stnlat,stnlon,rad([2,1,2,1])./1e3,[0,90,180,270]);
    roi = [min(lons),max(lons),min(lats),max(lats)];
  end;


  if ( ishandle(fh) )
    figure(fh);
  else
    % fh = figure;
    % hold on;
    % maxigraph;
    fh = fmg;
  end;

  % For most coral reef areas, "bathymetry" above MHHW is just a nuisance
  if ( isempty(doLabels) || (islogical(doLabels) && ~doLabels) )
    %contourf(stn.(fld).lon,stn.(fld).lat,stn.(fld).field,my_contours);
    [cs,h]=feval(contourFun,stn.(fld).lon,stn.(fld).lat,stn.(fld).field,my_contours);
  elseif ( ~inlineLabels )
    %contourf(stn.(fld).lon,stn.(fld).lat,stn.(fld).field,my_contours);
    [cs,h]=feval(contourFun,stn.(fld).lon,stn.(fld).lat,stn.(fld).field,my_contours);
    if ( islogical(doLabels) )
      clh = clabel(cs,'FontSize',20);
    elseif ( iscell(doLabels) )
      clh = clabel(cs,doLabels{:});
    else
      clh = clabel(cs,doLabels,'FontSize',20);
    end;
  else
    %contourf(stn.(fld).lon,stn.(fld).lat,stn.(fld).field,my_contours);
    [cs,h]=feval(contourFun,stn.(fld).lon,stn.(fld).lat,stn.(fld).field,my_contours);
    if ( islogical(doLabels) )
      clabel(cs,h,'FontSize',20,'LabelSpacing',4*72);
    elseif ( iscell(doLabels) )
      clabel(cs,h,doLabels{:});
    else
      clabel(cs,h,doLabels,'FontSize',20,'LabelSpacing',4*72);
    end;
    clh = h;
  end;

  cbh = colorbar;
  set(get(cbh,'Label'),'String','Depth [m]');

  if ( exist('set_surf_cursor','file') )
    set_surf_cursor;
  end;

  % Plot station location as a white pentagram bordered in black
  plot(stnlon,stnlat,'wp', 'MarkerEdgeColor','black', 'MarkerFaceColor','white');

  if ( doCoast )
    plot_hires_coastline(stn.(fld));
  end;

  % Limit plot to Region Of Interest (if that isn't whole field: see above)
  if ( numel(roi) == 4 )
    axis(roi);
  end;

  grid on;
  if ( isfield(stn,'station_name') )
    titlename(['Bathymetry surrounding ' char(upper(textize(stn.station_name)))]);
  else
    titlename(['Bathymetry surrounding site ']);
  end;

  daspect([1,cosd(stnlat),1]);

  if ( nargout < 1 )
    stn = []; clear stn;
  end;

return;
