function newfld = downsample_field(fld,avgpts,avgmethod)
%function newfld = downsample_field(fld[,avgpts[,avgmethod]])
%
% Call INTERP_FIELD (v.) to average (DEFAULT AVGMETHOD=@nanmean) over every
% AVGPTS (DEFAULT: 5) geographic points in 2D or 3D-field STRUCT FLD, which
% must have fields .lon, .lat, .field.
%
% Last Saved Time-stamp: <Fri 2018-08-17 10:16:31 Eastern Daylight Time gramer>
  
  if ( ~exist('avgpts','var') || isempty(avgpts) )
    avgpts = 5;
  end;
  if ( ~exist('avgmethod','var') || isempty(avgmethod) )
    avgmethod = @nanmean;
  end;
  
  interpMethod = {avgmethod,avgpts};
  newfld.lon = fld.lon(1:avgpts:end);
  newfld.lat = fld.lat(1:avgpts:end);
  [LONS,LATS] = meshgrid(newfld.lon,newfld.lat);
  dat = interp_field(fld.lat,fld.lon,fld.field,LATS,LONS,interpMethod);
  lat_res_deg = min(diff(unique(LATS(:))));
  lon_res_deg = min(diff(unique(LONS(:))));
  newfld.xres = distance_wgs84(LATS(1),LONS(1),LATS(1),LONS(1)+lon_res_deg)*1e3;
  newfld.yres = distance_wgs84(LATS(1),LONS(1),LATS(1)+lat_res_deg,LONS(1))*1e3;
  newfld.files = {'@nanmean',num2str(avgpts)};
  newfld.field = reshape(dat,[numel(lats),numel(lons)]);
  
return;
