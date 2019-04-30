function [stn,rotfld] = station_rotate_field(stn,ori,fld,rotfld,calcFieldTerms,grdtmplt,sitelat,sitelon,dx,interpMethod,extrapVal,tri)
%function [stn,rotfld] = station_rotate_field(stn,ori,fld[,rotfld[,calcFieldTerms[,grdtmplt[,sitelat[,sitelon[,dx[,interpMethod[,extrapVal[,tri]]]]]]]]])
%
% This function differs from STATION_REORIENT_VECTORS (v.), in that it
% creates a new time series field STN.(ROTFLD) whose gridpoints are rotated
% in two-dimensional space (by ORI or STN.(ORI) degrees True) relative to
% STN.(FLD), but with the same horizontal resolution (in KM) as STN.(FLD).
%
% Arguments STN, ORI, and FLD are all required. Other argument DEFAULTS:
% ROTFLD=[FLD '_rotated']; CALCFIELDTERMS=false; GRDTMPLT=3 (v. GRADIENTN);
% SITELAT=STN.lat if present in STN, else []; SITELON=STN.lon if present,
% else []; all other args have the same defaults as in ROTATE_FIELD (v.).
%
% NOTE: Function *WILL* overwrite STN.(ROTFLD) with a warning, if it exists.
%
% Last Saved Time-stamp: <Tue 2011-11-01 19:03:53  lew.gramer>

  % Validate our station struct input
  if ( ~exist('stn','var') || isempty(stn) || ~isstruct(stn) )
    error('First arg STN must be a station struct');
  end;

  % Validate our isobath orientation input
  if ( ischar(ori) && isfield(stn,ori) )
    ori = stn.(ori);
  end;
  if ( ~( isnumeric(ori) && isscalar(ori) && 0<=ori && ori<=360 ) )
    error('Second arg ORI must be a scalar or fieldname in STN specifying orientation, 0-360oT');
  end;

  % Validate our time series field to be rotated
  if ( ~exist('fld','var') || isempty(fld) || ~isfield(stn,fld) )
    error('Third arg FLD must be a fieldname in STN');
  end;
  if ( ~isfield(stn.(fld),'lon') || ~isfield(stn.(fld),'lat') || ~isfield(stn.(fld),'field') )
    error('Field STN.(FLD) must be a valid TIME SERIES FIELD (with .lon,.lat,.field)');
  end;

  if ( ~exist('rotfld','var') || isempty(rotfld) )
    rotfld = [fld '_rotated'];
  end;
  if ( isfield(stn,rotfld) )
    warning('Ecoforecasts:RotField:Overwrite',...
            'Overwriting existing STN.(%s)!',rotfld);
  end;

  if ( ~exist('calcFieldTerms','var') || isempty(calcFieldTerms) )
    calcFieldTerms = false;
  end;
  if ( ~exist('grdtmplt','var') || isempty(grdtmplt) )
    grdtmplt = 3;
  end;

  if ( ~exist('sitelat','var') )
    if ( isfield(stn,'lat') )
      sitelat=stn.lat;
    else
      sitelat = [];
    end;
  end;
  if ( ~exist('sitelon','var') )
    if ( isfield(stn,'lon') )
      sitelon=stn.lon;
    else
      sitelon = [];
    end;
  end;

  if ( ~exist('dx','var') )
    dx = [];
  end;
  if ( ~exist('interpMethod','var') )
    interpMethod = [];
  end;
  if ( ~exist('extrapVal','var') )
    extrapVal = [];
  end;
  if ( ~exist('tri','var') )
    tri = [];
  end;


  stn.(rotfld).date = stn.(fld).date;

  [stn.(rotfld).lat,stn.(rotfld).lon,stn.(rotfld).field,stn.(rotfld).dx] = ...
      rotate_field(stn.(fld).lat,stn.(fld).lon,stn.(fld).field,ori,sitelat,sitelon,interpMethod,extrapVal,tri);


  if ( calcFieldTerms )
    f = stn.(rotfld).field;
    dx = stn.(rotfld).dx.*1e3;

    dt = 1;
    rf = permute(f,[2 3 1]);
    [dTdx,dTdy,dTdt] = gradientn(rf,grdtmplt,dx,dx,dt);
    dTdx = permute(dTdx,[3 1 2]);
    dTdy = permute(dTdy,[3 1 2]);
    dTdt = permute(dTdt,[3 1 2]);
    rf=[]; clear rf

    l = repmat(nan,size(f));
    % For the Laplacian, no way that I can see to avoid this slow loop
    %DEBUG:
    disp(['Laplacian of ' rotfld]);
    for ix=1:size(f,1)
      % See HELP DEL2: must multiply result by 4
      l(ix,:,:) = 4.*del2(squeeze(f(ix,:,:)),dx,dx);
    end;

    stn.(rotfld).gradient_x = dTdx;
    stn.(rotfld).gradient_y = dTdy;
    stn.(rotfld).gradient_t = dTdt;
    stn.(rotfld).laplacian = l;

    nablatmplt = min(7,grdtmplt+2);

    rf = permute(dTdx,[2 3 1]);
    [dTdxx,dTdyx,dTdtx] = gradientn(rf,nablatmplt,dx,dx,dt);
    dTdxx = permute(dTdxx,[3 1 2]);
    stn.(rotfld).gradient_xx = dTdxx;

    rf = permute(dTdy,[2 3 1]);
    [dTdxy,dTdyy,dTdty] = gradientn(rf,nablatmplt,dx,dx,dt);
    dTdyy = permute(dTdyy,[3 1 2]);
    stn.(rotfld).gradient_yy = dTdyy;

    stn.(rotfld).nabla = stn.(rotfld).gradient_xx + stn.(rotfld).gradient_yy;
  end;

return;
