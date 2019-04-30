function result = station_field_transect(stn,fldnm,rng,az,org,dts,method)
%function result = station_field_transect(stn,fldnm,rng,az,org,dts,method)
%
% Extract a transect from field STN.(FLDNM) along azimuth AZ at distances
% in RNG (in [km]), beginning at the location of STN (STN.lon,STN.lat). If
% STN.(FLDNM) is a time series of fields, extract transect for each datenum
% in DTS (DEFAULT: all times in time series) using INTERP3 (v.). Otherwise,
% use INTERP2 (v.). In either case, call interpolator using METHOD (DEFAULT:
% for 3-D, 'spline'; for 2-D, 'nearest'). CALLS: TRANSECT_WGS84.
%
% Optional arg ORG may be a two-vector ([LON,LAT]) specifying the origin
% point of the transect(s): DEFAULT is [stn.lon,stn.lat].
%
% Last Saved Time-stamp: <Wed 2018-08-15 15:00:23 Eastern Daylight Time gramer>

  if ( any(~isfield(stn,{'lon','lat'})) )
    error('First arg STN must be a STRUCT with .lon and .lat fields');
  end;

  if ( ~exist('org','var') || isempty(org) )
    org = [stn.lon,stn.lat];
  end;
  if ( ~exist('dts','var') || isempty(dts) )
    if ( isfield(stn.(fldnm),'date') )
      dts = stn.(fldnm).date;
    else
      dts = [];
    end;
  end;

  [result.lon,result.lat] = transect_wgs84(org(1),org(2),rng,az);

  if ( isfield(stn.(fldnm),'date') )
    if ( ~exist('method','var') || isempty(method) )
      method = 'spline';
    end;
    [X,Y,Z] = meshgrid(stn.(fldnm).lon,stn.(fldnm).lat,stn.(fldnm).date);
    result.date = dts;
    result.field = interp3(X,Y,Z,...
                           permute(stn.(fldnm).field,[2 3 1]), ...
                           result.lon,result.lat,result.date,method);
  else
    if ( ~exist('method','var') || isempty(method) )
      method = 'nearest';
    end;
    % result.field = interp2(stn.(fldnm).lon,stn.(fldnm).lat,stn.(fldnm).field, ...
    [LONS,LATS] = meshgrid(unique(stn.(fldnm).lon),unique(stn.(fldnm).lat));
    %result.field = interp2(LONS,LATS,stn.(fldnm).field', ...
    ws = warning('OFF','MATLAB:interp2:NaNstrip');
    result.field = interp2(LONS,LATS,stn.(fldnm).field, ...
                           result.lon,result.lat,method);
    warning(ws);
  end;

  result.dx = rng;

return;
