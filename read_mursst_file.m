function stn = read_mursst_file(stn_or_stnm,date_or_fname)
%function stn = read_mursst_file(stn_or_stnm,date_or_fname)

  datapath = get_ecoforecasts_path('data');
  murpath = fullfile(datapath,'mursst');

  stn = get_station_from_station_name(stn_or_stnm);
  if ( ~exist('date_or_fname','var') || isempty(date_or_fname) )
    dt = datenum(2011,1,1);
    dtstr = '20110101';
    fname = fullfile(murpath,'20110101-JPL-L4UHfnd-GLOB-v01-fv03-MUR.nc');
  elseif ( ischar(date_or_fname) )
    fname = date_or_fname;
    dtstr = 'UNKNOWN'; % for now...
    dt = NaN; % for now...
  elseif ( isnumeric(date_or_fname) )
    dt = date_or_fname;
    dtstr = datestr(dt,'yyyymmdd');
    fname = fullfile(murpath,[dtstr '-JPL-L4UHfnd-GLOB-v01-fv03-MUR.nc']);
  else
    error('Second arg DATE_OR_FNAME must be a DATENUM or file name string');
  end;

  if ( ~exist(fname,'file') )
    error('Cannot find netCDF file "%s"',fname);
  end;


  lonrad=6; nlon = (2*lonrad)+1;
  latrad=5; nlat = (2*latrad)+1;

  nc = mDataset(fname);
  if ( isempty(nc) )
    error('Unable to open netCDF file "%s"',fname);
  end;
  if ( ~isfield(stn,'mursst_lonix') || ~isfield(stn,'mursst_latix') )
    lons=cast(nc{'lon'}(:),'double');
    lats=cast(nc{'lat'}(:),'double');
    [stn.mursst_lonerr,stn.mursst_lonix]=min(abs(stn.lon-lons));
    [stn.mursst_laterr,stn.mursst_latix]=min(abs(stn.lat-lats));
    stn.mursst_field.latix=stn.mursst_latix-latrad:stn.mursst_latix+latrad;
    stn.mursst_field.lonix=stn.mursst_lonix-lonrad:stn.mursst_lonix+lonrad;
    stn.mursst_field.lon=lons(stn.mursst_field.lonix);
    stn.mursst_field.lat=lats(stn.mursst_field.latix);
    stn.mursst_field.date=[];
    stn.mursst_field.field=repmat(nan,[0,nlat,nlon]);
  end;
  stn.mursst_field.date(end+1,1) = dt;
  stn.mursst_field.field(end+1,:,:) = ...
      cast(nc{'analysed_sst'}(1,stn.mursst_field.latix,stn.mursst_field.lonix),'double') - 273.14;
  close(nc); clear nc

return;
