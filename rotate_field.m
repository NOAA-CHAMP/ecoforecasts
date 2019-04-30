function [NEWLAT,NEWLON,newfld,dx]=rotate_field(lats,lons,fld,ori,sitelat,sitelon,dx,interpMethod,extrapVal,tri)
%function [NEWLAT,NEWLON,newfld,dx]=rotate_field(lats,lons,fld,ori[,sitelat[,sitelon[,dx[,interpMethod[,extrapVal[,tri]]]]])
%
% Return *time series field* NEWFLD by rotating all gridpoints of FLD by ORI
% degrees clockwise: LATS and LONS are monotonic vectors of length N and M,
% resp.; FLD is an array of size DxNxM; SITELAT, SITELON are scalar coords of
% a "center of rotation". Optional INTERPMETHOD (DEFAULT 'linear'), EXTRAPVAL
% (DEFAULT NaN), and TRI (DEFAULT false) are passed to INTERP_FIELD (v.).
% NEWFLD is a squared (DxRxR) time series field, with gridpoints spaced DX
% [km] (DEFAULT: min grid spacing in FLD) apart, at coords NEWLAT, NEWLON.
%
% NOTE: ROTATE_FIELD also handles special cases where FLD is 1xNxM or NxM.
%
% Last Saved Time-stamp: <Tue 2011-11-01 18:31:22  lew.gramer>

  if ( nargin < 4 )
    error('First four arguments are REQUIRED');
  end;
  if ( ~isnumeric(lats) || ~isnumeric(lons) || ~isnumeric(fld) || ~isnumeric(ori) )
    error('First four args must be numeric (vector, vector, matrix, scalar, resp.)');
  end;

  % Singleton special case: convert NxM array FLD to 1xNxM 3-d array
  if ( ndims(fld) == 2 && size(fld,1) == numel(lats) && size(fld,2) == numel(lons) )
    fld = reshape(fld,[1,size(fld,1),size(fld,2)]);
  end;
  if ( ndims(fld) ~= 3 || size(fld,2) ~= numel(lats) || size(fld,3) ~= numel(lons) )
    error('Args LATS,LONS,FLD must have N elements, M elements, and be DxNxM, resp.!');
  end;

  if ( ~( isnumeric(ori) && isscalar(ori) && 0<=ori && ori<=360 ) )
    error('Fourth arg ORI must be a scalar isobath orientation 0-360 degT');
  end;

  if ( ~exist('sitelat','var') || isempty(sitelat) )
    sitelat = mean(lats);
  end;
  if ( ~exist('sitelon','var') || isempty(sitelon) )
    sitelon = mean(lons);
  end;
  if ( ~isnumeric(sitelat) || ~isnumeric(sitelon) || ~isscalar(sitelat) || ~isscalar(sitelon) )
    error('Optional SITELAT,SITELON must be numerical scalars!');
  end;

  if ( ~exist('interpMethod','var') || isempty(interpMethod) )
    interpMethod = 'linear';
  end;
  if ( ~exist('extrapVal','var') || isempty(extrapVal) )
    extrapVal = nan;
  end;
  if ( ~exist('tri','var') || isempty(tri) )
    tri = false;
  end;

  % Find minimum gridpoint distance (assume it lies at one of the corners)
  if ( ~exist('dx','var') || isempty(dx) )
    uldlon = distance_wgs84(lats(1),lons(1),lats(1),lons(2));
    uldlat = distance_wgs84(lats(1),lons(1),lats(2),lons(1));
    urdlon = distance_wgs84(lats(1),lons(end),lats(1),lons(end-1));
    urdlat = distance_wgs84(lats(1),lons(end),lats(2),lons(end));

    lrdlon = distance_wgs84(lats(end),lons(end),lats(end),lons(end-1));
    lrdlat = distance_wgs84(lats(end),lons(end),lats(end-1),lons(end));
    lldlon = distance_wgs84(lats(end),lons(1),lats(end),lons(2));
    lldlat = distance_wgs84(lats(end),lons(1),lats(end-1),lons(1));

    dx = min([uldlon,uldlat,urdlon,urdlat,lrdlon,lrdlat,lldlon,lldlat]);
    %DEBUG:    disp(dx);
  end;

  D = size(fld,1);
  R = min(numel(lats),numel(lons));
  rngs = dx .* [ceil(-R/2):floor(R/2)];

  [llats,llons]=reckon_wgs84(sitelat,sitelon,rngs,ori);
  [xlats,xlons]=reckon_wgs84(sitelat,sitelon,rngs,ori+90);

  NEWLAT = repmat(nan,[numel(llats),numel(xlons)]);
  NEWLON = repmat(nan,[numel(llats),numel(xlons)]);
  for rix=1:numel(llats)
    [NEWLAT(rix,:),NEWLON(rix,:)]=reckon_wgs84(llats(rix),llons(rix),rngs,ori+90);
  end;

  newfld = interp_field(lats,lons,fld,NEWLAT(:),NEWLON(:),interpMethod,extrapVal,tri);

  newfld = squeeze(reshape(newfld,[D,R,R]));

  % TO DO LATER: Remove any row and/or column that is ALWAYS NaN

return;
