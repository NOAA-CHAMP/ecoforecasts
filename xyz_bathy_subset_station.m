function stn = xyz_bathy_subset_station(rawxyz,stn,rad,doPlot,doOverwrite,doReload)
%function stn = xyz_bathy_subset_station(rawxyz,stn,rad,doPlot,doOverwrite,doReload)
%
% Load NGDC 3-arcsecond (~92m) resolution Coastal Relief Model bathymetric
% data within RAD km the site specified by STN struct or STNM station name.
% RAD may be scalar, or vector [XRAD YRAD] (DEFAULT: 40 km). Return result
% in a STN.ngdc_92m_bathy field struct (with fields .lon, .lat, and .field),
% also saved in FULLFILE(DATAPATH,[STNM '_ngdc_bathy_<XRAD>x<YRAD>.mat']).
%
% Arg RAWXYZ must be an Nx3 matrix, or struct with fields .lon,.lat,.depth.
%
% If optional arg DOPLOT is TRUE, CONTOURF plot loaded sea floor topography.
% If optional arg DORELOAD is TRUE, try to load data from station MAT file.
%
% Last Saved Time-stamp: <Tue 2010-11-09 10:34:28  Lew.Gramer>

  datapath = get_ecoforecasts_path('data');
  coastpath = get_ecoforecasts_path('coast');

  if ( ~isstruct(rawxyz) )
    if ( size(rawxyz,2) >= 3 && isnumeric(rawxyz) )
      lon = rawxyz(:,1);
      lat = rawxyz(:,2);
      depth = rawxyz(:,3);
    else
      error('First arg RAWXYZ was non-numeric or had too few columns!');
    end;
  elseif ( all(isfield(rawxyz,{'lon','lat','depth'})) )
    lon = rawxyz.lon;
    lat = rawxyz.lat;
    depth = rawxyz.depth;
  else
    error('First arg RAWXYZ was a struct missing needed fields!');
  end;
  rawxyz = []; clear rawxyz;

  if ( isfield(stn,'station_name') )
    stnm = lower(stn.station_name);
  elseif ( ischar(stn) )
    stnm = lower(stn);
    clear stn;
    stn.station_name = upper(stnm);
  end;
  if ( ~all(isfield(stn,{'lat','lon'})) )
    [stn.lon,stn.lat,stn.depth] = get_station_coords(stn.station_name);
  end;

  if ( ~exist('rad','var') || isempty(rad) )
    rad = 40;
  end;
  if ( isscalar(rad) )
    rad = [ rad rad ];
  end;

  if ( ~exist('doPlot','var') || isempty(doPlot) )
    doPlot = false;
  end;

  if ( ~exist('doOverwrite','var') || isempty(doOverwrite) )
    doOverwrite = false;
  end;
  if ( ~exist('doReload','var') || isempty(doReload) )
    doReload = false;
  end;

  result = [];

  matfname = fullfile(datapath,sprintf('%s_ngdc_bathy_%gx%g.mat',stnm,rad(1),rad(2)));

  if ( exist(matfname,'file') )
    if (doOverwrite)
      warning('Deleting old file %s', matfname);
      delete(matfname);
    elseif (doReload)
      load(matfname,'result');
    else
      error('File already exists and doOVERWRITE,doRELOAD both FALSE: "%s"',matfname);
    end;
  end;

  if ( ~isfield(result,'ngdc_92m_bathy') )

    disp(['Subsetting raw XYZ data for ' upper(stnm)]);

    % Find all points "inside" our bounding rectangle
    inix = find( (abs(stn.lon-lon)<=(rad(1)/(111.12*cosd(stn.lat)))) ...
                 & (abs(stn.lat-lat)<=(rad(2)/(111.12))) );

    if ( numel(inix) < 4 )
      error('Insufficient data near %g,%g!',stn.lon,stn.lat);
    end;

    lon = lon(inix);
    lat = lat(inix);
    depth = depth(inix);

    [LON,LAT] = meshgrid(unique(lon),unique(lat));
    result.ngdc_92m_bathy.lon = LON';
    result.ngdc_92m_bathy.lat = LAT';
    result.ngdc_92m_bathy.field = griddata(lon,lat,depth,LON',LAT');

    disp(['Saving result to ' matfname]);
    save(matfname,'result');

  end;

  stn.ngdc_92m_bathy = result.ngdc_92m_bathy;
  result = []; clear result;

  if ( doPlot )
    figure;
    hold on;
    % contourf(stn.ngdc_92m_bathy.lon,stn.ngdc_92m_bathy.lat,stn.ngdc_92m_bathy.field);
    contourf(stn.ngdc_92m_bathy.lon,stn.ngdc_92m_bathy.lat,stn.ngdc_92m_bathy.field,-[0:5:80]);
    maxigraph;
    colorbar;
    plot(stn.lon,stn.lat,'wp');
    grid on;
    titlename(['NDGC bathymetry surrounding ' upper(stnm)]);
  end;

return;
