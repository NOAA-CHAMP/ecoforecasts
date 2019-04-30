1;
% Sample URL:
%  ftp://ftp.aoml.noaa.gov/pub/phod/matthieu.lehenaff/archv.2014_001_12.nc
% Sample file:
%  c:/Users/gramer/Documents/rsmas/Coastal/thesis/data/hycom/GOM2/archv.2014_001_12.nc

if ( ~exist('doFTP','var') || isempty(doFTP) )
  doFTP = false;
end;

datpath = get_ecoforecasts_path('data');
hycpath = get_thesis_path('../data/hycom/GOM2');

%%%% HACK %%%
% dts = datenum(2014,1,[1:365],12,0,0);
% dts = datenum(2014,1,1,12,0,0);
% dts = datenum(2014,1,[1:99],12,0,0);
dts = datenum(2014,1,[1:31],12,0,0);
% dts = datenum(2014,1,[1:200],12,0,0);

matfname = fullfile(datpath,sprintf('GOM2-%04d-%03d-%03d.mat',...
                                    get_year(dts(1)),get_jday(dts(1)),get_jday(dts(end))));

if ( exist(matfname,'file') )
  disp(['Loading ',matfname]);
  load(matfname);

else
  disp(['Extracting ',datestr(dts(1)),' - ',datestr(dts(end))]);
  mdl=[]; clear mdl

  % For this analysis, restrict ourselves to just the FRT reefs
  % mdl.lon_range = [-83.0,-79.5];
  % mdl.lat_range = [+24.0,+27.5];
  mdl.lon_range = [-84.0,-79.0];
  mdl.lat_range = [+24.0,+28.0];
  mdl.z_range = [0,50];

  ftph = [];
  if ( doFTP )
    try,
      ftph = ftp('ftp.aoml.noaa.gov','anonymous','lew.gramer@noaa.gov');
      cd(ftph,'/pub/phod/matthieu.lehenaff');
    catch,
      catchwarn;
      try, close(ftph); catch, end;
      ftph = [];
    end;
  end;

  for dtix=1:numel(dts)
    dt = dts(dtix);
    yr = get_year(dt);
    jd = get_jday(dt);
    hr = get_hour(dt);
    % disp(datestr(dt));
    if ( mod(jd,7) == 0 )
      disp(datestr(dt));
    end;

    ncname = sprintf('archv.%04d_%03d_12.nc',yr,jd);
    ncpath = fullfile(hycpath,ncname);

    if ( ~isempty(ftph) && ~exist(ncpath,'file') )
      try,
        mget(ftph,ncname,hycpath);
      catch,
        catchwarn;
      end;
    end;

    if ( ~exist(ncpath,'file') )
      warning('Skipping %s',ncpath);

    else
      mdl.d(dtix,1) = dt;
      nc = mDataset(ncpath);
      try,
        if ( ~isfield(mdl,'lon') )
          mdl.true_lon = cast(nc{'Longitude'}(:),'double');
          mdl.true_lat = cast(nc{'Latitude'}(:),'double');
          mdl.true_z = cast(nc{'Depth'}(:),'double');
          mdl.true_nlon = numel(mdl.true_lon);
          mdl.true_nlat = numel(mdl.true_lat);
          mdl.true_nz = numel(mdl.true_z);

          mdl.lonix = find(min(mdl.lon_range)<=mdl.true_lon & mdl.true_lon<=max(mdl.lon_range));
          mdl.latix = find(min(mdl.lat_range)<=mdl.true_lat & mdl.true_lat<=max(mdl.lat_range));
          mdl.zix = find(min(mdl.z_range)<=mdl.true_z & mdl.true_z<=max(mdl.z_range));

          mdl.lon = mdl.true_lon(mdl.lonix);
          mdl.lat = mdl.true_lat(mdl.latix);
          mdl.z = mdl.true_z(mdl.zix);

          % Preallocate matrices for ocean model variables
          mdl.nlon = numel(mdl.lon);
          mdl.nlat = numel(mdl.lat);
          mdl.nz = numel(mdl.z);

          maxsz = numel(dts)*mdl.nlat*mdl.nlon*mdl.nz*4*8;
          if ( maxsz > 1.1e9 )
            disp(['PLEASE CONFIRM: STRUCT may occupy ',num2str(maxsz/(2^30)),' Gb']);
            keyboard;
          end;

          mdl.u = repmat(nan,[numel(dts),mdl.nz,mdl.nlat,mdl.nlon]);
          mdl.v = repmat(nan,[numel(dts),mdl.nz,mdl.nlat,mdl.nlon]);
          %mdl.w = repmat(nan,[numel(dts),mdl.nz,mdl.nlat,mdl.nlon]);
          mdl.t = repmat(nan,[numel(dts),mdl.nz,mdl.nlat,mdl.nlon]);
          mdl.s = repmat(nan,[numel(dts),mdl.nz,mdl.nlat,mdl.nlon]);
          mdl.mld = repmat(nan,[numel(dts),mdl.nlat,mdl.nlon]);
        end;

        mdl.u(dtix,:,:,:) = squeeze(cast(nc{'u'}(:,mdl.zix,mdl.latix,mdl.lonix),'double'));
        mdl.v(dtix,:,:,:) = squeeze(cast(nc{'v'}(:,mdl.zix,mdl.latix,mdl.lonix),'double'));
        %mdl.w(dtix,:,:,:) = squeeze(cast(nc{'w_velocity'}(:,mdl.zix,mdl.latix,mdl.lonix),'double'));
        mdl.t(dtix,:,:,:) = squeeze(cast(nc{'temperature'}(:,mdl.zix,mdl.latix,mdl.lonix),'double'));
        mdl.s(dtix,:,:,:) = squeeze(cast(nc{'salinity'}(:,mdl.zix,mdl.latix,mdl.lonix),'double'));
        mdl.mld(dtix,:,:) = -squeeze(cast(nc{'mld'}(:,mdl.latix,mdl.lonix),'double'));
      catch,
        catchwarn;
      end;
      close(nc); clear nc;
    end;
  end; %for dtix=1:numel(dts)

  if ( ~isempty(ftph) )
    try,
      close(ftph); clear ftph;
    catch,
      catchwarn;
    end;
  end;


  % disp(['*NOT* Saving ',matfname]);
  disp(['Saving ',matfname]);
  save(matfname,'mdl');

end; %if ( exist(matfname,'file') ) else
