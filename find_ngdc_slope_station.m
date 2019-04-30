function [bet,ang,iso,dep,stn] = find_ngdc_slope_station(stn,rad_m,avg_m,fdif_m,method,extrapval)
%function [bet,ang,iso,dep,stn] = find_ngdc_slope_station(stn[,rad_m[,avg_m[,fdif_m[,method[,extrapval]]]]])
%
% Find seafloor slope BET from high-resolution bathymetry for station STN. If
% not present, calls READ_HIRES_BATHYMETRY to create STN.ngdc_hires_bathy,
% and calls FIND_NGDC_SLOPE to add fields .beta_deg and .beta to it.
%
% If AVG_M is given, calls INTERP_FIELD (v.) to downscale bathymetry to that
% resolution.
%
% ANG, angle of seafloor slope in degrees clockwise from True North, and
% local isobath angle ISO can also be returned (this is MOD(ANG-90,360)).
%
% Last Saved Time-stamp: <Mon 2018-07-02 20:06:27 Eastern Daylight Time gramer>

  if ( ~exist('stn','var') || ~isfield(stn,'lon') || ~isfield(stn,'lat') )
    error('First argument must be a station STRUCT');
  end;
  if ( ~exist('rad_m','var') || isempty(rad_m) || all(isnan(rad_m)) )
    rad_m = 1e3;
  end;
  if ( ~exist('avg_m','var') )
    avg_m = [];
  end;
  if ( ~exist('fdif_m','var') || isempty(fdif_m) )
    fdif_m = nan;
  end;
  if ( ~exist('method','var') )
    method = [];
  end;
  if ( ~exist('extrapval','var') )
    extrapval = [];
  end;

  minrad = nanmin([rad_m(:)',avg_m,fdif_m]);

  if ( numel(rad_m) == 1 )
    radx = nanmax([rad_m,fdif_m]);
    rady = nanmax([rad_m,fdif_m]);
  elseif ( numel(rad_m) == 2 )
    radx = nanmax([rad_m(1),fdif_m]);
    rady = nanmax([rad_m(2),fdif_m]);
  else
    error('Radius RAD_M must be a numeric 1- or 2-vector in [m]');
  end;

  if ( ~isfield(stn,'ngdc_hires_bathy') )
    stn = read_hires_bathymetry(stn,[radx,rady]);
  end;
  res_m = min([stn.ngdc_hires_bathy.xres(:)',stn.ngdc_hires_bathy.yres(:)']);

  avgpts = ceil(avg_m / res_m);
  if ( avgpts > 1 )
    % We now know AVG_M is not NaN
    interpMethod = {@nanmean,avgpts};
    lons = stn.ngdc_hires_bathy.lon(1:avgpts:end);
    lats = stn.ngdc_hires_bathy.lat(1:avgpts:end);
    [LONS,LATS] = meshgrid(lons,lats);
    fld = interp_field(stn.ngdc_hires_bathy.lat,stn.ngdc_hires_bathy.lon,stn.ngdc_hires_bathy.field,LATS,LONS,interpMethod);
    stn.ngdc_hires_bathy.lon = lons;
    stn.ngdc_hires_bathy.lat = lats;
    lat_res_deg = min(diff(unique(LATS(:))));
    lon_res_deg = min(diff(unique(LONS(:))));
    stn.ngdc_hires_bathy.xres = distance_wgs84(stn.lat,stn.lon,stn.lat,stn.lon+lon_res_deg)*1e3;
    stn.ngdc_hires_bathy.yres = distance_wgs84(stn.lat,stn.lon,stn.lat+lat_res_deg,stn.lon)*1e3;
    stn.ngdc_hires_bathy.files = {'@nanmean',num2str(avgpts)};
    stn.ngdc_hires_bathy.field = reshape(fld,[numel(lats),numel(lons)]);
    res_m = min([stn.ngdc_hires_bathy.xres,stn.ngdc_hires_bathy.yres]);
    LONS=[]; LATS=[]; clear LONS LATS
    fld=[]; clear fld
  end;
  if ( isnan(fdif_m) )
    npts = 11;
  else
    npts = ceil(fdif_m / res_m);
    npts = interp1([2,3:2:11],[2,3:2:11],npts,'nearest','extrap');
  end;

  [bet,ang,iso,stn.ngdc_hires_bathy] = find_ngdc_slope(stn.ngdc_hires_bathy,stn.lon,stn.lat,npts,method,extrapval);
  if ( nargout > 3 )
    dep = interp2(stn.ngdc_hires_bathy.lon,stn.ngdc_hires_bathy.lat,stn.ngdc_hires_bathy.field,stn.lon,stn.lat,'nearest');
  end;

return;
