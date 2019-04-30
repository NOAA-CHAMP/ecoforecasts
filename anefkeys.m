1;
% Sample file:
%  c:/Users/gramer/Documents/rsmas/Coastal/thesis/data/hycom/eFKEYS_01-1-7-2012/eFKEYS_archv.2012_001_12_3zuvwts.nc

datpath = get_ecoforecasts_path('data');
hycpath = get_thesis_path('../data/hycom/eFKEYS_01-1-7-2012');

matfname = fullfile(datpath,'eFKEYS_01-1-7-2012.mat');

if ( exist(matfname,'file') )
  disp(['Loading ',matfname]);
  load(matfname);

else
  mdl=[]; clear mdl

  % For this analysis, restrict ourselves to just the FRT reefs
  mdl.lon_range = [-83.0,-79.5];
  mdl.lat_range = [+24.0,+27.5];
  % "Surface" (0 m) eFKEYS layer is almost identical to 1 m layer
  mdl.z_range = [0,140];

  dts = datenum(2012,1,1,[6:6:(7*24)-1],0,0);

  for dtix=1:numel(dts)
    dt = dts(dtix);
    yr = get_year(dt);
    jd = get_jday(dt);
    hr = get_hour(dt);
    disp(datestr(dt));
    % if ( hr == 0 )
    %   disp(datestr(dt));
    % end;

    ncname = sprintf('eFKEYS_archv.%04d_%03d_%02d_3zuvwts.nc',yr,jd,hr);
    ncpath = fullfile(hycpath,ncname);
    if ( ~exist(ncpath,'file') )
      warning('Skipping %s',ncpath);
    else
      mdl.d(dtix,1) = dt;
      nc = mDataset(ncpath);
      try,
        if ( ~isfield(mdl,'lon') )
          % HACK ALERT: Can't fit 17 levels in memory for even ONE variable!
          % This hack keeps file size down below 27*653*589*4*8*3 = 9.97 Gb.
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
        end;

        mdl.u(dtix,:,:,:) = squeeze(cast(nc{'u'}(:,mdl.zix,mdl.latix,mdl.lonix),'double'));
        mdl.v(dtix,:,:,:) = squeeze(cast(nc{'v'}(:,mdl.zix,mdl.latix,mdl.lonix),'double'));
        %mdl.w(dtix,:,:,:) = squeeze(cast(nc{'w_velocity'}(:,mdl.zix,mdl.latix,mdl.lonix),'double'));
        mdl.t(dtix,:,:,:) = squeeze(cast(nc{'temperature'}(:,mdl.zix,mdl.latix,mdl.lonix),'double'));
        mdl.s(dtix,:,:,:) = squeeze(cast(nc{'salinity'}(:,mdl.zix,mdl.latix,mdl.lonix),'double'));
      catch,
        catchwarn;
      end;
      close(nc); clear nc;
    end;
  end; %for dtix=1:numel(dts)

  disp(['Saving ',matfname]);
  save(matfname,'mdl');

end; %if ( exist(matfname,'file') ) else
