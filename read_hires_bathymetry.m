function [stn,rad] = read_hires_bathymetry(stn_or_stnm_or_locn,rad,bathfiles,useHighestRes,fld)
%function [stn,rad] = read_hires_bathymetry(stn_or_stnm_or_locn,rad,bathfiles,useHighestRes,fld)
%
% Retrieve bathymetry surrounding station struct STN, named site STNM, or
% lon/lat 2-vector LOCN, from bathymetry files BATHFILES in either netCDF/GRD
% (using njToolbox, v.) or ASC format. Unlike GET_NGDC_BATHY_STATION (v.),
% this function always reads the bathymetry file: result is never stored in
% intermediate MAT file so station name STNM or STN.station_name string is
% *not* required. Subsets to a RAD(2)xRAD(1) m (DEFAULT: [40e3,40e3], i.e.,
% 40 km N-S by 40 km E-W) rectangle. If USEHIGHESTRES (DEFAULT: true), and
% bathymetry finer than 92 m res is available on disk, use it instead. (Field
% names don't change.) If multiple files are within the requested radius,
% stitch them together to produce the resulting field. 
%
% NOTE: When stitching together, files with highest resolution and/or most
% SHELF COVERAGE will overlay lower-resolution/more-oceanic files: uses #
% valid gridpoints <=200 m depth as a proxy for resolution/shelf coverage.
%
% Last Saved Time-stamp: <Mon 2019-04-29 16:34:22 Eastern Daylight Time gramer>

  if ( isnumeric(stn_or_stnm_or_locn) )
    stn.lon = stn_or_stnm_or_locn(1);
    stn.lat = stn_or_stnm_or_locn(2);
    stn.station_name = num2str([stn.lon,stn.lat]);
  elseif ( iscell(stn_or_stnm_or_locn) )
    stn.lon = stn_or_stnm_or_locn{1};
    stn.lat = stn_or_stnm_or_locn{2};
    if ( numel(stn_or_stnm_or_locn) >= 3 )
      stn.station_name = stn_or_stnm_or_locn{3};
    else
      stn.station_name = num2str([stn.lon,stn.lat]);
    end;
  else
    stn = get_station_from_station_name(stn_or_stnm_or_locn);
    if ( ~isfield(stn,'lon') || ~isfield(stn,'lat') )
      if ( isfield(stn,'station_name') )
        error('Do not know coordinates for %s',stn.station_name);
      else
        error('No station name or lon/lat given');
      end;
    end;
  end;
  clear stn_or_stnm_or_locn;
  if ( ~isscalar(stn.lon) || ~isscalar(stn.lat) )
    error('First arg should specify scalar LON and LAT');
  end;

  if ( ~exist('rad','var') || isempty(rad) )
    rad = 40e3;
  end;
  if ( numel(rad) == 1 )
    rad = [rad,rad];
  end;
  rad = rad(:)';
  if ( ~exist('useHighestRes','var') || isempty(useHighestRes) )
    useHighestRes = true;
  end;
  if ( ~exist('fld','var') || isempty(fld) )
    fld = 'ngdc_hires_bathy';
  end;

  coastpath = get_ecoforecasts_path('coast');

  [lats,lons] = reckon_wgs84(stn.lat,stn.lon,rad([2,1,2,1])./1e3,[0,90,180,270]);
  bbox = [min(lons),max(lons),min(lats),max(lats)];

  if ( ~exist('bathfiles','var') || isempty(bathfiles) )
    bathfiles = {};

    % Do any of the files we've accumulated so far CONTAIN OUR LOCATION?
    stn_covered = false;

    % Use the finest available bathymetry - but don't mix resolutions
    if ( useHighestRes )

      % [bathfiles,stn_covered] = ...
      %     read_hires_bathymetry_check_file(stn,bbox,filebbox,file,bathfiles,stn_covered);

      % % 5 m products (UH/NOAA-CRED: Pacific Islands Benthic Habitat Mapping Center)
      % % NOTE: Coverage of these files can be sketchy! Use with care...
      % [bathfiles,stn_covered] = ...
      %     read_hires_bathymetry_check_file(stn,bbox,[145.5485,145.8534,15.0654,15.3132],...
      %                                      'sai_mb_li_db',bathfiles,stn_covered); % Saipan - CNMI
      % [bathfiles,stn_covered] = ...
      %     read_hires_bathymetry_check_file(stn,bbox,[145.1050,145.3150,14.0434,14.2292],...
      %                                      'rot_dbmb_5m',bathfiles,stn_covered); % Rota - CNMI
      % [bathfiles,stn_covered] = ...
      %     read_hires_bathymetry_check_file(stn,bbox,[145.5016,145.7559,14.8029,15.1510],...
      %                                      'tinmblidbmos',bathfiles,stn_covered); % Tinian - CNMI

      % [bathfiles,stn_covered] = ...
      %     read_hires_bathymetry_check_file(stn,bbox,[-170.9282,-170.4849,-14.3839,-14.1952],...
      %                                      'tut_dbmb',bathfiles,stn_covered); % Tutuila - AmSamoa
      % [bathfiles,stn_covered] = ...
      %     read_hires_bathymetry_check_file(stn,bbox,[-169.7054,-169.5838,-14.2032,-14.1334],...
      %                                      'oo_dbmb_mos4',bathfiles,stn_covered); % Ofu-Olosega - AmSamoa
      % [bathfiles,stn_covered] = ...
      %     read_hires_bathymetry_check_file(stn,bbox,[-169.5975,-169.4116,-14.2858,-14.1858],...
      %                                      'tau_dbmb_mos',bathfiles,stn_covered); % Ta'u - AmSamoa

      % 30 m INTERPOLATED products (UH/NOAA-CRED: Pacific Islands Benthic Habitat Mapping Center)
      % NOTE: Digital Elevation Model/bathymetry data interleaved, downsampled
      % and spline-interpolated by Lew.Gramer@noaa.gov! Use with care...
      % Saipan - CNMI
      [bathfiles,stn_covered] = ...
          read_hires_bathymetry_check_file(stn,bbox,[  145.5486  145.8532   15.0657   15.3130],...
                                           'pibhmc_bathy_5m_saipan.dods.30m.usgs_dem_10m_saipan.dods',...
                                           bathfiles,stn_covered);
      % Rota - CNMI
      [bathfiles,stn_covered] = ...
          read_hires_bathymetry_check_file(stn,bbox,[  145.1024  145.3165   14.0871   14.2309],...
                                           'pibhmc_bathy_5m_rota.dods.30m.usgs_dem_10m_rota.dods',...
                                           bathfiles,stn_covered);
      % Tinian - CNMI
      [bathfiles,stn_covered] = ...
          read_hires_bathymetry_check_file(stn,bbox,[  145.5220  145.7368   14.7984   15.1270],...
                                           'pibhmc_bathy_5m_tinian.PARTIAL.dods.30m.usgs_dem_10m_tinian.dods',...
                                           bathfiles,stn_covered);

      % Tutuila - AmSamoa
      [bathfiles,stn_covered] = ...
          read_hires_bathymetry_check_file(stn,bbox,[ -170.9280 -170.4850  -14.3839  -14.1957],...
                                           'pibhmc_bathy_5m_tutuila.dods.30m.usgs_dem_10m_tutuila.dods',...
                                           bathfiles,stn_covered);
      % Ofu-Olosega - AmSamoa
      [bathfiles,stn_covered] = ...
          read_hires_bathymetry_check_file(stn,bbox,[ -169.7051 -169.5840  -14.2032  -14.1339],...
                                           'pibhmc_bathy_5m_ofuolosega.dods.30m.usgs_dem_10m_ofuolosega.dods',...
                                           bathfiles,stn_covered);
      % Ta'u - AmSamoa
      [bathfiles,stn_covered] = ...
          read_hires_bathymetry_check_file(stn,bbox,[ -169.5973 -169.4117  -14.2858  -14.1863],...
                                           'pibhmc_bathy_5m_tau.dods.30m.usgs_dem_10m_tau.dods',...
                                           bathfiles,stn_covered);


      if ( ~stn_covered )
        %DEBUG:        disp(['Station not covered']);
        bathfiles = {};
      end;


      % 1/3-arcsecond products (~10 m)
      if ( isempty(bathfiles) )
        [bathfiles,stn_covered] = ...
            read_hires_bathymetry_check_file(stn,bbox,[-88.30,-87.65,30.00,31.00],...
                                             fullfile(coastpath,'mobile_al_mhw.grd'),bathfiles,stn_covered);
        [bathfiles,stn_covered] = ...
            read_hires_bathymetry_check_file(stn,bbox,[-91.60,-88.80,28.60,29.70],...
                                             fullfile(coastpath,'southern_louisiana_mhw.grd'),bathfiles,stn_covered);
        % [bathfiles,stn_covered] = ...
        %     read_hires_bathymetry_check_file(stn,bbox,[-80.3600,-79.4201,26.2900,27.3099],...
        %                                      'pb_mhw',bathfiles,stn_covered); % Palm Beach, Florida
        [bathfiles,stn_covered] = ...
            read_hires_bathymetry_check_file(stn,bbox,[-80.3600,-79.8901,26.7999,27.3099],...
                                             'pb_mhw_NW',bathfiles,stn_covered); % Palm Beach, Florida
        [bathfiles,stn_covered] = ...
            read_hires_bathymetry_check_file(stn,bbox,[-79.8901,-79.4201,26.7999,27.3099],...
                                             'pb_mhw_NE',bathfiles,stn_covered); % Palm Beach, Florida
        [bathfiles,stn_covered] = ...
            read_hires_bathymetry_check_file(stn,bbox,[-80.3600,-79.8901,26.2900,26.7999],...
                                             'pb_mhw_SW',bathfiles,stn_covered); % Palm Beach, Florida
        [bathfiles,stn_covered] = ...
            read_hires_bathymetry_check_file(stn,bbox,[-79.8901,-79.4201,26.2900,26.7999],...
                                             'pb_mhw_SE',bathfiles,stn_covered); % Palm Beach, Florida
        [bathfiles,stn_covered] = ...
            read_hires_bathymetry_check_file(stn,bbox,[-89.30,-88.30,29.70,30.60],...
                                             'biloxi_ms',bathfiles,stn_covered);
        [bathfiles,stn_covered] = ...
            read_hires_bathymetry_check_file(stn,bbox,[-86.10,-85.20,29.55,30.50],...
                                             fullfile(coastpath,'panama_city_fl_1-3_arc-second_mhw_netcdf.grd'),bathfiles,stn_covered);
        % [bathfiles,stn_covered] = ...
        %     read_hires_bathymetry_check_file(stn,bbox,[-80.4100,-79.4000,25.2500,26.3200],...
        %                                      'miami_fl_mhw',bathfiles,stn_covered);
        [bathfiles,stn_covered] = ...
            read_hires_bathymetry_check_file(stn,bbox,[-80.4100,-79.9050,25.2500,25.7850],...
                                             'miami_fl_mhw_SW',bathfiles,stn_covered);
        [bathfiles,stn_covered] = ...
            read_hires_bathymetry_check_file(stn,bbox,[-79.9051,-79.4000,25.2500,25.7850],...
                                             'miami_fl_mhw_SE',bathfiles,stn_covered);
        [bathfiles,stn_covered] = ...
            read_hires_bathymetry_check_file(stn,bbox,[-80.4100,-79.9050,25.7849,26.3200],...
                                             'miami_fl_mhw_NW',bathfiles,stn_covered);
        [bathfiles,stn_covered] = ...
            read_hires_bathymetry_check_file(stn,bbox,[-79.9051,-79.4000,25.7849,26.3200],...
                                             'miami_fl_mhw_NE',bathfiles,stn_covered);
        % [bathfiles,stn_covered] = ...
        %     read_hires_bathymetry_check_file(stn,bbox,[-82.18,-81.27,23.98,25.02],...
        %                                      'key_west_fl_mhw',bathfiles,stn_covered);
        [bathfiles,stn_covered] = ...
            read_hires_bathymetry_check_file(stn,bbox,[-82.1800,-81.7250,23.9800,24.5000],...
                                             'key_west_fl_mhw_SW',bathfiles,stn_covered);
        [bathfiles,stn_covered] = ...
            read_hires_bathymetry_check_file(stn,bbox,[-81.7251,-81.2700,23.9800,24.5000],...
                                             'key_west_fl_mhw_SE',bathfiles,stn_covered);
        [bathfiles,stn_covered] = ...
            read_hires_bathymetry_check_file(stn,bbox,[-82.1800,-81.7250,24.4999,25.0200],...
                                             'key_west_fl_mhw_NW',bathfiles,stn_covered);
        [bathfiles,stn_covered] = ...
            read_hires_bathymetry_check_file(stn,bbox,[-81.7251,-81.2700,24.4999,25.0200],...
                                             'key_west_fl_mhw_NE',bathfiles,stn_covered);
        [bathfiles,stn_covered] = ...
            read_hires_bathymetry_check_file(stn,bbox,[-65.04,-64.40,17.61,17.88],...
                                             'st_croix_mhw_update',bathfiles,stn_covered);
        [bathfiles,stn_covered] = ...
            read_hires_bathymetry_check_file(stn,bbox,[-65.15,-64.64,18.17,18.48],...
                                             'sts_thoms_john_mhw_update',bathfiles,stn_covered);
        [bathfiles,stn_covered] = ...
            read_hires_bathymetry_check_file(stn,bbox,[-65.90,-65.35,18.05,18.60],...
                                             'fajardo',bathfiles,stn_covered);
        [bathfiles,stn_covered] = ...
            read_hires_bathymetry_check_file(stn,bbox,[-67.60,-67.10,17.90,18.60],...
                                             'mayaguez',bathfiles,stn_covered);
        [bathfiles,stn_covered] = ...
            read_hires_bathymetry_check_file(stn,bbox,[-67.10,-66.40,17.70,18.05],...
                                             'ponce',bathfiles,stn_covered);
        [bathfiles,stn_covered] = ...
            read_hires_bathymetry_check_file(stn,bbox,[145.66995,145.84793,15.07495,15.30237],...
                                             'cnmi-saipan_c3-ssl-REVERSE',bathfiles,stn_covered);
        [bathfiles,stn_covered] = ...
            read_hires_bathymetry_check_file(stn,bbox,[144.30,145.20,13.00,13.90],...
                                             'guam_1_3s_20081014',bathfiles,stn_covered);

        if ( ~stn_covered )
          %DEBUG:        disp(['Station not covered']);
          bathfiles = {};
        end;
      end; %if ( isempty(bathfiles) )

      % 1-arcsecond (or 30 m) products
      if ( isempty(bathfiles) )
        % % Stupid USGS DEM format took forever to decipher for South Florida 30 m file
        % [bathfiles,stn_covered] = ...
        %   read_hires_bathymetry_check_file(stn,bbox,[-82.1389,-79.9876,24.3906,26.2517],...
        %                                    'F010_25081C1_BIG.dem',bathfiles,stn_covered);
        [bathfiles,stn_covered] = ...
          read_hires_bathymetry_check_file(stn,bbox,[-82.1389,-81.0633,24.3906,25.3212],...
                                           'F010_25081C1_BIG.dem_SW',bathfiles,stn_covered);
        [bathfiles,stn_covered] = ...
          read_hires_bathymetry_check_file(stn,bbox,[-81.0636,-79.9876,24.3906,25.3212],...
                                           'F010_25081C1_BIG.dem_SE',bathfiles,stn_covered);
        [bathfiles,stn_covered] = ...
          read_hires_bathymetry_check_file(stn,bbox,[-82.1389,-81.0633,25.3209,26.2517],...
                                           'F010_25081C1_BIG.dem_NW',bathfiles,stn_covered);
        [bathfiles,stn_covered] = ...
          read_hires_bathymetry_check_file(stn,bbox,[-81.0636,-79.9876,25.3209,26.2517],...
                                           'F010_25081C1_BIG.dem_NE',bathfiles,stn_covered);

        % And HEY PRESTO! NOAA finally makes these data available via netCDF
        % just when I probably can't make good use of them any more. :(
        % Lew.Gramer@noaa.gov, 2019 Apr 29
        %% https://www.ngdc.noaa.gov/thredds/catalog/regional/catalog.html?dataset=regionalDatasetScan/southern_florida_F010_2018.nc

        % Old file was: 'us_virgin_islands_vi_1s_mhw';
        [bathfiles,stn_covered] = ...
          read_hires_bathymetry_check_file(stn,bbox,[-65.15,-64.00,17.00,19.00],...
                                           'us_virgin_islands_vi_1s_mhw_update',bathfiles,stn_covered);
        [bathfiles,stn_covered] = ...
          read_hires_bathymetry_check_file(stn,bbox,[-68.00,-65.00,17.00,19.00],...
                                           'pr_1s',bathfiles,stn_covered);
        [bathfiles,stn_covered] = ...
          read_hires_bathymetry_check_file(stn,bbox,[-90.75,-85.00,24.00,38.00],...
                                           fullfile(coastpath,'northern_gulf_coast_mhw.grd'),bathfiles,stn_covered);
        [bathfiles,stn_covered] = ...
          read_hires_bathymetry_check_file(stn,bbox,[-65.06,-64.51,32.17,32.55],...
                                           'Bermuda_DEM_5012_1_msl',bathfiles,stn_covered);
      end; %if ( isempty(bathfiles) )

    end; %if ( useHighestRes )

    if ( ~stn_covered )
      %DEBUG:      disp(['Station not covered']);
      bathfiles = {};
    end;

    % 3-arcsecond products (6-arcsecond for Mariana region)
    if ( isempty(bathfiles) ) 
      [bathfiles,stn_covered] = ...
          read_hires_bathymetry_check_file(stn,bbox,[-87,-78,24,35],...
                                           fullfile(coastpath,'fl_east_gom_crm_v1.nc'),bathfiles,stn_covered);
      [bathfiles,stn_covered] = ...
          read_hires_bathymetry_check_file(stn,bbox,[-108,-94,24,38],...
                                           fullfile(coastpath,'western_gom_crm_v1.nc'),bathfiles,stn_covered);
      [bathfiles,stn_covered] = ...
          read_hires_bathymetry_check_file(stn,bbox,[-94,-87,24,36],...
                                           fullfile(coastpath,'central_gom_crm_v1.nc'),bathfiles,stn_covered);
      [bathfiles,stn_covered] = ...
          read_hires_bathymetry_check_file(stn,bbox,[-65.5,-64.0,31.5,33.0],...
                                           'Bermuda_DEM_5013_3_msl',bathfiles,stn_covered);
      [bathfiles,stn_covered] = ...
          read_hires_bathymetry_check_file(stn,bbox,[-162,-152,18,24],...
                                           fullfile(coastpath,'hawaii_crm_v1.nc'),bathfiles,stn_covered);
      [bathfiles,stn_covered] = ...
          read_hires_bathymetry_check_file(stn,bbox,[-171.14,-169.30,-14.73,-13.83],...
                                           'pago_pago_3s',bathfiles,stn_covered);
      [bathfiles,stn_covered] = ...
          read_hires_bathymetry_check_file(stn,bbox,[139.00,150.00,10.00,22.00],...
                                           'MarianaTrench_6as_v1',bathfiles,stn_covered);
      [bathfiles,stn_covered] = ...
          read_hires_bathymetry_check_file(stn,bbox,[-85,-68,31,40],...
                                           fullfile(coastpath,'se_atl_crm_v1.nc'),bathfiles,stn_covered);
      [bathfiles,stn_covered] = ...
          read_hires_bathymetry_check_file(stn,bbox,[-128,-117,37,44],...
                                           fullfile(coastpath,'central_pacific_crm_v1.nc'),bathfiles,stn_covered);
      [bathfiles,stn_covered] = ...
          read_hires_bathymetry_check_file(stn,bbox,[-80,-64,40,48],...
                                           fullfile(coastpath,'ne_atl_crm_v1.nc'),bathfiles,stn_covered);
      [bathfiles,stn_covered] = ...
          read_hires_bathymetry_check_file(stn,bbox,[-128,-116,44,49],...
                                           fullfile(coastpath,'nw_pacific_crm_v1.nc.nc'),bathfiles,stn_covered);
    end; %if ( isempty(bathfiles) ) % 3-arcsecond products (6 for Marianas)


    % Regional Fallbacks: 30-arcsecond GEBCO_8 data:
    %  GEBCO_2014_2D_-180.0_-30.0_-100.0_35.0.nc
    %  GEBCO_2014_2D_-100.0_-30.0_-10.0_35.0.nc
    %  GEBCO_2014_2D_30.0_-30.0_120.0_35.0.nc
    %  GEBCO_2014_2D_120.0_-30.0_180.0_35.0.nc
    if ( isempty(bathfiles) ) 
      [bathfiles,stn_covered] = ...
          read_hires_bathymetry_check_file(stn,bbox,[-180.0,-100.0,-30.0,+35.0],...
                                           fullfile(coastpath,'GEBCO_2014_2D_-180.0_-30.0_-100.0_35.0.nc'),...
                                           bathfiles,stn_covered);
    end;
    if ( isempty(bathfiles) ) 
      [bathfiles,stn_covered] = ...
          read_hires_bathymetry_check_file(stn,bbox,[-100.0,-10.0,-30.0,+35.0],...
                                           fullfile(coastpath,'GEBCO_2014_2D_-100.0_-30.0_-10.0_35.0.nc'),...
                                           bathfiles,stn_covered);
    end;
    if ( isempty(bathfiles) ) 
      [bathfiles,stn_covered] = ...
          read_hires_bathymetry_check_file(stn,bbox,[-30.0,120.0,-30.0,+35.0],...
                                           fullfile(coastpath,'GEBCO_2014_2D_30.0_-30.0_120.0_35.0.nc'),...
                                           bathfiles,stn_covered);
    end;
    if ( isempty(bathfiles) ) 
      [bathfiles,stn_covered] = ...
          read_hires_bathymetry_check_file(stn,bbox,[120.0,180.0,-30.0,+35.0],...
                                           fullfile(coastpath,'GEBCO_2014_2D_120.0_-30.0_180.0_35.0.nc'),...
                                           bathfiles,stn_covered);
    end;


    % GLOBAL FALLBACK: 1-arcminute ETOPO data
    if ( isempty(bathfiles) ) 
      [bathfiles,stn_covered] = ...
          read_hires_bathymetry_check_file(stn,bbox,[-180.0,+180.0,-90.0,+90.0],...
                                           fullfile(coastpath,'ETOPO1_Bed_g_gmt4.grd'),bathfiles,stn_covered);
    end;

    if ( ~stn_covered )
      warning('Location not in any file found! Try USEHIGHESTRES=False?');
    end;

  end; %if ( ~exist('bathfiles','var') || isempty(bathfiles) )


  % Don't force caller to wrap a single filename string in '{}'
  if ( ischar(bathfiles) )
    bathfiles = {bathfiles};
  end;

  if ( isempty(bathfiles) )
    error('Please specify a netCDF, GRD or ASC bathymetry file containing those coordinates??');
  end;

  minlon = +inf;
  maxlon = -inf;
  minlat = +inf;
  maxlat = -inf;
  for bix = 1:numel(bathfiles)
    bathfile = bathfiles{bix};
    %DEBUG:    disp(['Reading ',bathfile]);

    [pth,fnm,fex]=fileparts(bathfile);
    doAsc = ~( strcmpi(fex,'.nc') || strcmpi(fex,'.grd') );
    if ( doAsc )
      [all_lons,all_lats,all_zs] = asc2mat(bathfile);
    else
      %DEBUG:      disp(['Loading ',bathfile]);
      nc = mDataset(bathfile);
      if ( isempty(nc) )
        error('MDATASET unable to open %s',bathfile);
      end;
      % Non-standard variable names from GEBCO
      if ( ~isempty(strfind(bathfile,'GEBCO')) )
        all_lats = cast(nc{'lat'}(:),'double');
        all_lons = cast(nc{'lon'}(:),'double');
      else
        all_lats = cast(nc{'y'}(:),'double');
        all_lons = cast(nc{'x'}(:),'double');
      end;
    end; %if ( doAsc ) else

    nlats = numel(all_lats);
    nlons = numel(all_lons);
    lat_res_deg = min(diff(unique(all_lats(:))));
    lon_res_deg = min(diff(unique(all_lons(:))));
    lat_res = distance_wgs84(stn.lat,stn.lon,stn.lat+lat_res_deg,stn.lon)*1e3;
    lon_res = distance_wgs84(stn.lat,stn.lon,stn.lat,stn.lon+lon_res_deg)*1e3;
    dlat = ceil(rad(2)/lat_res);
    dlon = ceil(rad(1)/cosd(stn.lat)/lon_res);

    % NOTE: If we're stitching together multiple files, be careful around
    % file boundaries to adjust the radius of the rectangle we subset to!
    latix = round(interp1(all_lats,1:numel(all_lats),stn.lat,'linear','extrap'));
    lonix = round(interp1(all_lons,1:numel(all_lons),stn.lon,'linear','extrap'));

    latixen = latix-dlat:latix+dlat;
    lonixen = lonix-dlon:lonix+dlon;
    latixen(1 > latixen | latixen > nlats) = [];
    lonixen(1 > lonixen | lonixen > nlons) = [];
    if ( isempty(latixen) || isempty(lonixen) )
      error('Location %f,%f lies outside %fx%f neighborhood of "%s"',stn.lon,stn.lat,rad(1),rad(2),bathfile);
    end;

    if ( doAsc )
      z = all_zs(latixen,lonixen)';
    else
      % Non-standard variable names from GEBCO
      if ( ~isempty(strfind(bathfile,'GEBCO')) )
        z = cast(nc{'elevation'}(latixen,lonixen),'double')';
      else
        z = cast(nc{'z'}(latixen,lonixen),'double')';
      end; %if ( ~isempty(strfind(bathfile,'GEBCO')) )
      close(nc);
    end; %if ( doAsc )

    x(bix).doAsc = doAsc;
    x(bix).lon = all_lons(lonixen);
    x(bix).lat = all_lats(latixen);
    x(bix).field = z';
    x(bix).dlon = lon_res_deg;
    x(bix).dlat = lat_res_deg;
    x(bix).dlon_km = lon_res;
    x(bix).dlat_km = lat_res;
    x(bix).shelf_n = numel(find(z>-200));

    minlon = min([minlon,x(bix).lon(:)']);
    maxlon = max([maxlon,x(bix).lon(:)']);
    minlat = min([minlat,x(bix).lat(:)']);
    maxlat = max([maxlat,x(bix).lat(:)']);

  end; %for bix = 1:numel(bathfiles)

  dlon = median([x.dlon]);
  dlat = median([x.dlat]);
  lon = minlon:dlon:maxlon;
  lat = minlat:dlat:maxlat;
  stn.(fld).lon = lon;
  stn.(fld).lat = lat;
  stn.(fld).requested_rad = rad;

  % Quilt rectangles from different datasets onto one result rectangle
  if ( numel(x) == 1 )
    stn.(fld).field = x.field;
    stn.(fld).xres = x.dlon_km;
    stn.(fld).yres = x.dlat_km;
  else
    if ( any(abs(diff([x.dlon])) >= dlon/2) || any(abs(diff([x.dlat])) >= dlat/2) )
      %DEBUG:      disp('RESOLUTION MISMATCH!'); keyboard;
      error('Resolution mismatch between datasets!');
    end;

    % Preallocate result rectangle with NaNs
    stn.(fld).field = repmat(nan,[numel(lat),numel(lon)]);

    % Apply datasets in order from smallest to largest area of intersection in
    % shelf waters (<200 m depth). This way, less salient data is overwritten.
    [ig,ixs] = sort([x.shelf_n]);
    for bix=ixs(:)';
      lonix = interp1(lon,1:numel(lon),x(bix).lon,'nearest','extrap');
      latix = interp1(lat,1:numel(lat),x(bix).lat,'nearest','extrap');
      stn.(fld).lon(lonix) = x(bix).lon;
      stn.(fld).lat(latix) = x(bix).lat;
      stn.(fld).field(latix,lonix) = x(bix).field;
      stn.(fld).xres(bix) = x(bix).dlon_km;
      stn.(fld).yres(bix) = x(bix).dlat_km;
    end; %for bix=ixs(:)';
  end; %if ( numel(x) == 1 ) else

  stn.(fld).files = bathfiles;

  %DEBUG:  disp([mfilename,' DONE']);

return;


%%%% 
%%%% PRIVATE FUNCTIONS
%%%% 

function [bathfiles,stn_covered] = read_hires_bathymetry_check_file(stn,bbox,filebbox,file,bathfiles,stn_covered)
  if ( bboxint(bbox,filebbox) > eps )
    bathfiles{end+1} = file;
  end;
  if ( ~isempty(bboxinside(stn.lon,stn.lat,filebbox,false,0)) )
    stn_covered = true;
  end;
return;
