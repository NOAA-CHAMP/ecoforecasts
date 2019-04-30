function stn = read_ngdc_bathymetry(stn_or_stnm_or_locn,rad,bathfile,useHighestRes)
%function stn = read_ngdc_bathymetry(stn_or_stnm_or_locn,rad,bathfile,useHighestRes)
%
% Retrieve bathymetry surrounding station struct STN, named site STNM, or
% lon/lat 2-vector LOCN, from netCDF bathymetry file BATHFILE (DEFAULT:
% "ECOFORECASTS/coast/fl_east_gom_crm_v1.nc") using njToolbox (v.). Unlike
% GET_NGDC_BATHY_STATION (v.), this function always reads the netCDF
% bathymetry file: result is never stored in intermediate MAT file so station
% name STNM or STN.station_name string is *not* required. Subsets bathymetry
% to a RAD(2)xRAD(1) m (DEFAULT: [40e3,40e3], i.e., 40 km) rectangle.
% If USEHIGHESTRES (DEFAULT: true), and bathymetry is available on the local
% disk at better than 92 m res, use it instead. (Field names do not change.)
%
% Last Saved Time-stamp: <Wed 2015-06-17 15:36:15 Eastern Daylight Time gramer>

  if ( isnumeric(stn_or_stnm_or_locn) )
    stn.lon = stn_or_stnm_or_locn(1);
    stn.lat = stn_or_stnm_or_locn(2);
    stn.station_name = num2str([stn.lon,stn.lat]);
  else
    stn = get_station_from_station_name(stn_or_stnm_or_locn);
  end;
  clear stn_or_stnm_or_locn;

  if ( ~exist('rad','var') || isempty(rad) )
    rad = 40e3;
  end;
  if ( numel(rad) == 1 )
    rad = [rad,rad];
  end;

  if ( ~exist('useHighestRes','var') || isempty(useHighestRes) )
    useHighestRes = true;
  end;

  coastpath = get_ecoforecasts_path('coast');
  if ( ~exist('bathfile','var') || isempty(bathfile) )
    % Use the finest available bathymetry
    if ( useHighestRes && ...
         24<=stn.lat-(rad(2)/111e3) && stn.lat+(rad(2)/111e3)<=38 ...
         && -90.75<=stn.lon-(rad(1)/111e3) && stn.lon+(rad(1)/111e3)<=-85.0 )
      bathfile = fullfile(coastpath,'northern_gulf_coast_mhw.grd');
    elseif ( 24<=stn.lat && stn.lat<=35 && -87<=stn.lon && stn.lon<=-78 )
      bathfile = fullfile(coastpath,'fl_east_gom_crm_v1.nc');
    elseif ( 24<=stn.lat && stn.lat<=38 && -108<=stn.lon && stn.lon<=-94 )
      bathfile = fullfile(coastpath,'western_gom_crm_v1.nc');
    elseif ( 24<=stn.lat && stn.lat<=36 &&  -94<=stn.lon && stn.lon<=-87 )
      bathfile = fullfile(coastpath,'central_gom_crm_v1.nc');
    elseif ( 16<=stn.lat && stn.lat<=20 &&  -68<=stn.lon && stn.lon<=-64 )
      bathfile = fullfile(coastpath,'puerto_rico_crm_v1.nc');
    elseif ( 18<=stn.lat && stn.lat<=24 && -162<=stn.lon && stn.lon<=-152 )
      bathfile = fullfile(coastpath,'hawaii_crm_v1.nc');
    else
      error('Please specify an NGDC CRM netCDF file including those coordinates??');
    end;
    disp(['Using ',bathfile]);
  end;

  nc = mDataset(bathfile);
  if ( isempty(nc) )
    error('MDATASET unable to open %s',bathfile);
  end;
  try,
    all_lats = cast(nc{'y'}(:),'double');
    all_lons = cast(nc{'x'}(:),'double');
    nlats = numel(all_lats);
    nlons = numel(all_lons);
    lat_res_deg = min(diff(unique(all_lats(:))));
    lon_res_deg = min(diff(unique(all_lons(:))));
    lat_res = distance_wgs84(stn.lat,stn.lon,stn.lat+lat_res_deg,stn.lon)*1e3;
    lon_res = distance_wgs84(stn.lat,stn.lon,stn.lat,stn.lon+lon_res_deg)*1e3;
    % dlat=rad(2)/92.5;
    % dlon=rad(1)/cosd(25)/92.5;
    dlat = ceil(rad(2)/lat_res);
    dlon = ceil(rad(1)/cosd(stn.lat)/lon_res);

    [lonerr,lonix]=min(abs(all_lons-stn.lon));
    [laterr,latix]=min(abs(all_lats-stn.lat));
    latixen = latix-dlat:latix+dlat;
    lonixen = lonix-dlon:lonix+dlon;
    latixen(1 > latixen | latixen > nlats) = [];
    lonixen(1 > lonixen | lonixen > nlons) = [];
    if ( ~isempty(latixen) && ~isempty(lonixen) )
      z = cast(nc{'z'}(latixen,lonixen),'double')';
    end;
  catch,
    catchwarn(['Failure reading ',bathfile]);
  end;
  close(nc);

  if ( lonerr > (2*lon_res) || laterr > (2*lat_res) || isempty(latixen) || isempty(lonixen) )
    error('Location %f,%f lies outside domain of "%s"',stn.lon,stn.lat,bathfile);
  end;

  %stn.ngdc_92m_bathy.bathfile = bathfile;
  stn.ngdc_92m_bathy.lonix = lonix;
  stn.ngdc_92m_bathy.latix = latix;
  % stn.ngdc_92m_bathy.lon = all_lons(lonix-dlon:lonix+dlon);
  % stn.ngdc_92m_bathy.lat = all_lats(latix-dlat:latix+dlat);
  stn.ngdc_92m_bathy.lon = all_lons(lonixen);
  stn.ngdc_92m_bathy.lat = all_lats(latixen);
  stn.ngdc_92m_bathy.field = z';

return;
