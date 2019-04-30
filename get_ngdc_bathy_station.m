function stn = get_ngdc_bathy_station(stn_or_stnm,rad,rawxyz,doPlot)
%function stn = get_ngdc_bathy_station(stn_or_stnm,rad,rawxyz,doPlot)
%
% Load NGDC 3-arcsecond (~92m) resolution Coastal Relief Model bathymetric
% data within RAD km the site specified by STN struct or STNM station name.
% RAD may be scalar, or vector [XRAD YRAD] (DEFAULT: 40 km). Return result
% in a STN.ngdc_92m_bathy field struct (with fields .lon, .lat, and .field),
% also saved in FULLFILE(DATAPATH,[STNM '_ngdc_bathy_<XRAD>x<YRAD>.mat']).
%
% If optional arg RAWXYZ is a struct with .lon,.lat,.depth fields, extract
% topography from that struct directly; don't (re)load NGDC bathymetry data,
% and don't bother to try reloading station's NGDC 92m MAT file (see above).
%
% If optional arg DOPLOT is TRUE, CONTOURF plot loaded sea floor topography.
%
% Last Saved Time-stamp: <Wed 2012-07-25 17:17:21  lew.gramer>

  datapath = get_ecoforecasts_path('data');
  coastpath = get_ecoforecasts_path('coast');

  stn = get_station_from_station_name(stn_or_stnm);
  stnm = lower(stn.station_name);
  clear stn_or_stnm;

  if ( ~exist('rad','var') || isempty(rad) )
    rad = 40;
  end;
  if ( isscalar(rad) )
    rad = [ rad rad ];
  end;

  if ( ~exist('rawxyz','var') || isempty(rawxyz) )
    rawxyz = [];
  end;

  if ( ~exist('doPlot','var') || isempty(doPlot) )
    doPlot = false;
  end;

  matfname = fullfile(datapath,sprintf('%s_ngdc_bathy_%gx%g.mat',stnm,rad(1),rad(2)));

  if ( ~exist(matfname,'file') && ~all(isfield(stn,{'lat','lon'})) )
    % NOTE: Only try to get stn coords by name if we definitely need them!
    [stn.lon,stn.lat,stn.depth] = get_station_coords(stn.station_name);
  end;

  % If user passed us a RAWXYZ, (re-)subset our bathymetry from it
  if ( all(isfield(rawxyz,{'lon','lat','depth'})) )
    disp('Subsetting raw data from RAWXYZ struct');
    result = xyz_bathy_subset_station(rawxyz,stn,rad);

  % Otherwise, if we previously subsetted, just reload it from MAT file
  elseif ( exist(matfname,'file') )
    load(matfname,'result');

  % Otherwise, subset from our default (very large) raw XYZ MAT file
  else
    rawmatfname = fullfile(coastpath,'LGramer1-80.mat');
    disp(['Subsetting raw data from ' rawmatfname]);
    if ( ~exist(rawmatfname,'file') )
      error('Cannot find raw MAT file "%s"!',rawmatfname);
    end;
    rawxyz = load(rawmatfname,'lon','lat','depth');
    result = xyz_bathy_subset_station(rawxyz,stn,rad);
  end;

  rawxyz = []; clear rawxyz;


  stn.ngdc_92m_bathy = result.ngdc_92m_bathy;
  result = []; clear result;

  if ( doPlot )
    plot_ngdc_bathy_station(stn);
  end;

  if ( nargout < 1 )
    stn = []; clear stn;
  end;

return;
