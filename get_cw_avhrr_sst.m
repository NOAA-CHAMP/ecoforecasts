function stn = get_cw_avhrr_sst(stn_or_loc,rad,forceExtract)
%function stn = get_cw_avhrr_sst(stn_or_loc,rad,forceExtract)
%
% Get AVHRR SST product from CoastWatch/OceanWatch Caribbean/Gulf of Mexico
% node. Gets RAD(1)xRAD(2) (DEFAULT: 10x10) pixel field around STN location.
% If FORCEEXTRACT (DEFAULT: false), ignore MAT file and get fresh NC data.
%
% Last Saved Time-stamp: <Wed 2016-01-27 14:28:35 Eastern Standard Time lew.gramer>

  % Retrieve data from a path relative to this M-file's local directory
  if ( ~exist('datapath', 'var') || isempty(datapath) )
    datapath = get_ecoforecasts_path('data');
  end;
  if ( ~exist(datapath, 'dir') )
      error('CW:NoDataPath', ...
            'The data source path "%s" does not exist!', ...
            datapath);
  end;

  if ( isnumeric(stn_or_loc) && numel(stn_or_loc) == 2 )
    % STN_OR_LOC is a vector with [LON,LAT]
    stn.lon = stn_or_loc(1);
    stn.lat = stn_or_loc(2);
  elseif ( ischar(stn_or_loc) )
    stn = get_station_from_station_name(stn_or_loc);
  else
    stn = stn_or_loc;
  end;
  clear stn_or_loc

  if ( ~exist('rad','var') || isempty(rad) )
    rad = [10,10];
  end;
  if ( isnumeric(rad) && isscalar(rad) )
    rad = [rad,rad];
  end;
  if ( ~isnumeric(rad) || numel(rad) ~= 2 )
    error('RAD arg must be a numeric one- or two-vector');
  end;
  if ( ~exist('forceExtract','var') || isempty(forceExtract) )
    forceExtract = false;
  end;

  matfname = [];
  if ( isfield(stn,'station_name') )
    matfname = fullfile(datapath,[lower(stn.station_name),'-cw_avhrr_sst.mat']);
  end;

  % Reload previously-extracted raw CW AVHRR SST from a MAT file
  if ( ~forceExtract && exist(matfname,'file') )
    disp(['Loading ',matfname]);
    stn.cw_avhrr_sst = load(matfname);

  % Extract raw CW AVHRR SST from online netCDF data
  else
    stn.avhrr_ncfname = 'http://cwcgom.aoml.noaa.gov/thredds/dodsC/AVHRRSST/SST.nc';
    nc = mDataset(stn.avhrr_ncfname);
    if ( isempty(nc) )
      error('Unable to open %s',stn.avhrr_ncfname);
    end;
    %nj_info(nc),
    t = (cast(nc{'time'}(:),'double')/3600/24)+datenum(1970,1,1,0,0,0);
    lon = cast(nc{'x'}(:),'double');
    lat = cast(nc{'y'}(:),'double');
    [lonerr,lonix] = min(abs(stn.lon-lon))
    [laterr,latix] = min(abs(stn.lat-lat))
    if ( lonerr > 1 || laterr > 1 )
      close(nc); clear nc
      error('Station location not within product domain!');
    end;
    lonixen=lonix-rad(1):lonix+rad(1);
    latixen=latix-rad(2):latix+rad(2);
    sst = cast(nc{'sst'}(:,latixen,lonixen),'double');
    close(nc); clear nc

    stn.cw_avhrr_sst.date = t;
    stn.cw_avhrr_sst.data = sst(:,10,10);
    stn.cw_avhrr_sst.lon = lon(lonixen);
    stn.cw_avhrr_sst.lat = lat(latixen);
    stn.cw_avhrr_sst.field = sst;
    clear lat latixen lon lonixen sst t

    if ( isempty(matfname) )
      disp('Unable to save MAT file - no Station Name');
    else
      cw_avhrr_sst = stn.cw_avhrr_sst;
      disp(['Saving ',matfname]);
      save(matfname,'-struct','cw_avhrr_sst');
      clear cw_avhrr_sst
    end;
  end;

return;
