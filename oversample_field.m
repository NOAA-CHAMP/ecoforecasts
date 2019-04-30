function fld = oversample_field(lofld,hixy,fldnm,varargin)
%function fld = oversample_field(lofld,hixy[,fldnm[,tol][,interpMethod,extrapVal,triangular,...])
%
% Call INTERP_FIELD (v.) to oversample field LOFLD, a STRUCT with fields
% .lon, .lat, and .(FLDNM) (DEFAULT: .field), and optionally .date, onto the
% coordinates in the vector fields of STRUCT HIXY (i.e., non-plaid .lon,
% .lat, and optionally .date). If optional TOL is given (a numeric value), it
% is passed to INTERSECT_DATES. Any other args are passed to INTERP_FIELD.
%
% SEE ALSO oversample_attenuate_field.m
%
% Last Saved Time-stamp: <Mon 2017-07-03 21:20:47 Eastern Daylight Time gramer>

  if ( ~exist('fldnm','var') || isempty(fldnm) )
    fldnm = 'field';
  end;

  tol = [];
  args = varargin;
  if ( numel(args) >= 1 && isnumeric(args{1}) )
    tol = args{1};
    args(1) = [];
  end;

  [LON,LAT] = meshgrid(hixy.lon,hixy.lat);
  %function vals = interp_field(lats,lons,fld,sitelat,sitelon,interpMethod,extrapVal,triangular,doTiming)
  fldvec = interp_field(lofld.lat,lofld.lon,lofld.(fldnm),LAT,LON,args{:});

  if ( ~isfield(hixy,'date') || ~isfield(lofld,'date') )
    fld = reshape(fldvec(:),[numel(hixy.lat),numel(hixy.lon)]);

  else
    [loix,hiix] = intersect_dates(lofld.date,hixy.date,tol);
    fldvec = fldvec(loix,:);
    fld = reshape(fldvec(:),[numel(hiix),numel(hixy.lat),numel(hixy.lon)]);
  end;

return;
