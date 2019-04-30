function stn = calc_field_terms(stn,fnm,tnm,interpMethod,site_lat,site_lon,grdtmplt,doSecondOrder,doHourly,varargin)
%function stn = calc_field_terms(stn,fnm[,tnm[,interpMethod[,site_lat,site_lon[,grdtmplt[,doSecondOrder[,doHourly...]]]]]])
%
% Calculate two-dimensional gradient and Laplacian of 2-D field time series
% STN.(FNM). Expects STN.(FNM) to have sub-fields .date,.lon,.lat,.field, and
% STN.(FNM).field be a TxMxN matrix with T time steps. Four TxMxN fields are
% added to STN.(FNM): .gradient_x,.gradient_y,.gradient_t,.laplacian. If time
% series field name TNM is given (DEFAULT none), interpolate field terms onto
% station location, returning results in new fields [TNM '_x'], [TNM '_y'],
% [TNM '_l']. If optional arg INTERPMETHOD is 'midpoint' (DEFAULT), choose
% point in each field nearest to center; otherwise, pass to INTERP_FIELD (v.)
% with coords SITE_LAT and SITE_LON if present, else STN.lat and STN.lon.
% Optional GRDTMPLT, an integer [3,5,7], specifies centered finite-difference
% template size for GRADIENTN(v.). DEFAULT: 3 points (identical to GRADIENT).
% NOTE: This argument still does not affect Laplacian result in any way...
% If DOSECONDORDER (DEFAULT: true), then also calculate 2nd order gradients,
% and Direct estimate of the Laplacian (as XX+YY, *not* from calling DEL2).
% If DOHOURLY (DEFAULT: true) and string TNM is present, also interpolate all
% gradients and Laplacians we calculate into hourly time series (INTERP_TS).
% Any additional args after DOHOURLY are passed through to INTERP_TS.
%
% Sample Calling Sequences:
%
%  >> % Load FKEYS HYCOM data for Molasses Reef SEAKEYS site MLRF1
%  >> stn = get_fkeys_hycom('mlrf1');
%  >> % Calculate del and del^2 fields for sea surface temperature
%  >> stn = calc_field_terms(stn,'fkeys_hycom_seatemp_field');
%
%  >> % Load USF weekly 1km AVHRR SST for Sombrero SEAKEYS site SMKF1
%  >> stn = get_avhrr_weekly_field('smkf1');
%  >> % Calc field terms using 7-point finite difference for SST gradients
%  >> stn =  calc_field_terms(stn,'avhrr_weekly_sst_field',[],[],[],[],7);
%  >> % Recalculate using 5-point finite difference for SST gradients, and
%  >> % linear-triangularly interpolate to time series at station location
%  >> stn =  calc_field_terms(stn,'avhrr_weekly_sst_field',...
%  >>                         'avhrr_weekly_sst','tri,linear',...
%  >>                         stn.lon,stn.lat,5);
%
% Last Saved Time-stamp: <Tue 2012-05-08 15:22:06  Lew.Gramer>

  if ( ~exist('tnm','var') || isempty(tnm) )
    tnm = [];
  end;
  if ( ~exist('interpMethod','var') || isempty(interpMethod) )
    interpMethod = [];
  end;
  if ( ~exist('site_lat','var') || isempty(site_lat) )
    site_lat = [];
  end;
  if ( ~exist('site_lon','var') || isempty(site_lon) )
    site_lon = [];
  end;
  if ( isempty(site_lat) && isfield(stn,'lat') )
    site_lat = stn.lat;
  end;
  if ( isempty(site_lon) && isfield(stn,'lon') )
    site_lon = stn.lon;
  end;

  if ( ~exist('grdtmplt','var') || isempty(grdtmplt) )
    grdtmplt = 3;
  else
    %DEBUG:
    warning('Gradient template %g points',grdtmplt);
  end;
  if ( ~exist('doSecondOrder','var') || isempty(doSecondOrder) )
    doSecondOrder = true;
  end;
  if ( ~exist('doHourly','var') || isempty(doHourly) )
    doHourly = true;
  end;

  oldflds={'date','lon','lat','field'};
  if ( ~isfield(stn,fnm) || ~all(isfield(stn.(fnm),oldflds)) || ndims(stn.(fnm).field)~=3 )
    error('STN.(FNM) must exist and be a well-formed time series of 2-D fields!');
  end;
  if ( length(unique(stn.(fnm).lon)) < grdtmplt && length(unique(stn.(fnm).lat)) < grdtmplt )
    error('STN.(FNM).field must have at least %d points to calculate a gradient!',grdtmplt);
  end;

  newflds={'gradient_x','gradient_y','gradient_t','laplacian'};
  for ix=1:length(newflds)
    try; stn.(fnm)=rmfield(stn.(fnm),newflds{ix}); catch; end;
  end;

  dts = stn.(fnm).date(:);
  lon = stn.(fnm).lon(:);
  lat = stn.(fnm).lat(:);
  fld = stn.(fnm).field;

  % Reflect U/D and/or L/R so our field has monotonically-increasing coordinates
  fliplon = false;
  fliplat = false;
  if ( any(diff(lon)<=0) )
    if ( any(diff(lon)>=0) )
      error('Longitude not monotonic!');
    else
      lon = lon(end:-1:1);
      fld = flipdim(fld,3);
      fliplon = true;
    end;
  end;
  if ( any(diff(lat)<=0) )
    if ( any(diff(lat)>=0) )
      error('Latitude not monotonic!');
    else
      lat = lat(end:-1:1);
      fld = flipdim(fld,2);
      fliplat = true;
    end;
  end;

  midlon = repmat(lon(round(end/2)),size(lat));
  midlat = repmat(lat(round(end/2)),size(lon));

  % Use MKS units [but time in index units] for gradient and Laplacian
  dt = 1;
  dx = [0 ; cumsum(sw_dist(midlat(:),lon(:),'km'))]*1e3;
  dy = [0 ; cumsum(sw_dist(lat(:),midlon(:),'km'))]*1e3;

  % We probably don't really care about 'dTdt', but this form of call to
  % GRADIENT is generally much faster than calling it once per time unit

  rotfld = permute(fld,[2 3 1]);
  [dTdx,dTdy,dTdt] = gradientn(rotfld,grdtmplt,dx,dy,dt);
  % [dTdx,dTdy,dTdt] = gradient(rotfld,dx,dy,dt);
  % [dTdx,dTdy,dTdt] = gradientn(rotfld,3,dx,dy,dt);
  % [dTdx,dTdy,dTdt] = gradientn(rotfld,5,dx,dy,dt);
  % [dTdx,dTdy,dTdt] = gradientn(rotfld,7,dx,dy,dt);

  if ( doSecondOrder )
    [dTdxx,dTdxy,dTdxt] = gradientn(dTdx,grdtmplt,dx,dy,dt);
    [dTdyx,dTdyy,dTdyt] = gradientn(dTdy,grdtmplt,dx,dy,dt);
  end;

  dTdt = permute(dTdt,[3 1 2]);
  dTdx = permute(dTdx,[3 1 2]);
  dTdy = permute(dTdy,[3 1 2]);

  if ( doSecondOrder )
    dTdxx = permute(dTdxx,[3 1 2]);
    dTdxy = permute(dTdxy,[3 1 2]);
    dTdyx = permute(dTdyx,[3 1 2]);
    dTdyy = permute(dTdyy,[3 1 2]);
  end;

  l = repmat(nan,size(fld));
  % For the Laplacian, no way that I can see to avoid this slow loop
  %DEBUG:
  disp(['Laplacian of ' fnm]);
  for ix=1:size(fld,1)
    % See HELP DEL2: must multiply result by 4
    l(ix,:,:) = 4.*del2(squeeze(fld(ix,:,:)),dx,dy);
  end;

  % If we had to reflect input field, reflect results back to match original
  if ( fliplat )
    dTdt = flipdim(dTdt,2);
    dTdx = flipdim(dTdx,2);
    dTdy = flipdim(dTdy,2);
    l = flipdim(l,2);
    if ( doSecondOrder )
      dTdxx = flipdim(dTdxx,2);
      dTdxy = flipdim(dTdxy,2);
      dTdyx = flipdim(dTdyx,2);
      dTdyy = flipdim(dTdyy,2);
    end;
  end;
  if ( fliplon )
    dTdt = flipdim(dTdt,3);
    dTdx = flipdim(dTdx,3);
    dTdy = flipdim(dTdy,3);
    l = flipdim(l,3);
    if ( doSecondOrder )
      dTdxx = flipdim(dTdxx,3);
      dTdxy = flipdim(dTdxy,3);
      dTdyx = flipdim(dTdyx,3);
      dTdyy = flipdim(dTdyy,3);
    end;
  end;

  stn.(fnm).gradient_t = dTdt;
  stn.(fnm).gradient_x = dTdx;
  stn.(fnm).gradient_y = dTdy;
  stn.(fnm).laplacian = l;

  if ( doSecondOrder )
    stn.(fnm).gradient_xx = dTdxx;
    stn.(fnm).gradient_xy = dTdxy;
    stn.(fnm).gradient_yx = dTdyx;
    stn.(fnm).gradient_yy = dTdyy;
    % Direct Laplacian (not DEL2 estimate)
    stn.(fnm).gradient_dl = dTdxx + dTdyy;
  end;

  % Is caller also requested station-based time series of field terms
  if ( ~isempty(tnm) )

    xfld = [tnm '_x'];
    yfld = [tnm '_y'];
    lfld = [tnm '_l'];
    stn.(xfld).date = dts;
    stn.(yfld).date = dts;
    stn.(lfld).date = dts;
    if ( strncmpi(interpMethod,'midpoint',3) )
      % We are only interested in the central point
      stn.(xfld).data = squeeze(dTdx(:,midx,midy));
      stn.(yfld).data = squeeze(dTdy(:,midx,midy));
      stn.(lfld).data = squeeze(l(:,midx,midy));
    else
      if ( isempty(site_lat) || isempty(site_lon) )
        error('Interpolation requested but site LAT, LON not both specified!');
      end;
      stn.(xfld).data = interp_field(lat,lon,dTdx,site_lat,site_lon,interpMethod);
      stn.(yfld).data = interp_field(lat,lon,dTdy,site_lat,site_lon,interpMethod);
      stn.(lfld).data = interp_field(lat,lon,l,site_lat,site_lon,interpMethod);
    end;

    if ( doSecondOrder )
      xxfld = [tnm '_xx'];
      xyfld = [tnm '_xy'];
      yxfld = [tnm '_yx'];
      yyfld = [tnm '_yy'];
      dlfld = [tnm '_dl'];
      stn.(xxfld).date = dts;
      stn.(xyfld).date = dts;
      stn.(yxfld).date = dts;
      stn.(yyfld).date = dts;
      stn.(dlfld).date = dts;
      if ( strncmpi(interpMethod,'midpoint',3) )
        % We are only interested in the central point
        stn.(xxfld).data = squeeze(dTdxx(:,midx,midy));
        stn.(xyfld).data = squeeze(dTdxy(:,midx,midy));
        stn.(yxfld).data = squeeze(dTdyx(:,midx,midy));
        stn.(yyfld).data = squeeze(dTdyy(:,midx,midy));
        stn.(dlfld).data = squeeze(dTdxx(:,midx,midy)+dTdyy(:,midx,midy));
      else
        stn.(xxfld).data = interp_field(lat,lon,dTdxx,site_lat,site_lon,interpMethod);
        stn.(xyfld).data = interp_field(lat,lon,dTdxy,site_lat,site_lon,interpMethod);
        stn.(yxfld).data = interp_field(lat,lon,dTdyx,site_lat,site_lon,interpMethod);
        stn.(yyfld).data = interp_field(lat,lon,dTdyy,site_lat,site_lon,interpMethod);
        stn.(dlfld).data = interp_field(lat,lon,dTdxx+dTdyy,site_lat,site_lon,interpMethod);
      end;
    end;


    % Also spline-fit native time increment into an hourly time series
    if ( doHourly )
      disp(['Hourly time series in hourly_',tnm,'_*']);
      hxfld = ['hourly_' xfld]; stn.(hxfld) = struct('date',[],'data',[]);
      hyfld = ['hourly_' yfld]; stn.(hyfld) = struct('date',[],'data',[]);
      hlfld = ['hourly_' lfld]; stn.(hlfld) = struct('date',[],'data',[]);
      if ( is_valid_ts(stn.(xfld)) ); stn.(hxfld) = interp_ts(stn.(xfld),varargin{:}); end;
      if ( is_valid_ts(stn.(yfld)) ); stn.(hyfld) = interp_ts(stn.(yfld),varargin{:}); end;
      if ( is_valid_ts(stn.(lfld)) ); stn.(hlfld) = interp_ts(stn.(lfld),varargin{:}); end;

      if ( doSecondOrder )
        hxxfld = ['hourly_' xxfld]; stn.(hxxfld) = struct('date',[],'data',[]);
        hxyfld = ['hourly_' xyfld]; stn.(hxyfld) = struct('date',[],'data',[]);
        hyxfld = ['hourly_' yxfld]; stn.(hyxfld) = struct('date',[],'data',[]);
        hyyfld = ['hourly_' yyfld]; stn.(hyyfld) = struct('date',[],'data',[]);
        hdlfld = ['hourly_' dlfld]; stn.(hdlfld) = struct('date',[],'data',[]);
        if ( is_valid_ts(stn.(xxfld)) ); stn.(hxxfld) = interp_ts(stn.(xxfld),varargin{:}); end;
        if ( is_valid_ts(stn.(xyfld)) ); stn.(hxyfld) = interp_ts(stn.(xyfld),varargin{:}); end;
        if ( is_valid_ts(stn.(yxfld)) ); stn.(hyxfld) = interp_ts(stn.(yxfld),varargin{:}); end;
        if ( is_valid_ts(stn.(yyfld)) ); stn.(hyyfld) = interp_ts(stn.(yyfld),varargin{:}); end;
        % Should we add interpolated 2nd order gradients here, instead of interpolating DL?
        if ( is_valid_ts(stn.(dlfld)) ); stn.(hdlfld) = interp_ts(stn.(dlfld),varargin{:}); end;
      end;
    end;

  end;

return;
