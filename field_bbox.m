function bbox = field_bbox(fld,gutter)
%function bbox = field_bbox(fld,gutter)
%
% Return bounding box for coordinates of field FLD. FLD may be a STRUCT with
% .lon and .lat fields, 2xN or Nx2 matrix, or cell array with two vectors of
% LON and LAT coords., resp. User can specify an optional GUTTER (DEFAULT: no
% gutter): if the string 'gutter', one extra gridpoint-width; if scalar, that
% length in degrees; if a 2-vector, number of degrees longitude and latitude,
% resp.; if a 4-vector, left-, right-, btm-, top-bounds, resp.
%
% Last Saved Time-stamp: <Sat 2017-02-11 16:34:23 Eastern Standard Time gramer>

  bbox = [];
  if ( isfield(fld,'lon') && isfield(fld,'lat') )
    lons = fld.lon(:);
    lats = fld.lat(:);
  elseif ( iscell(fld) && numel(fld) == 2 && isnumeric(fld{1}) && isnumeric(fld{2}) )
    lons = fld{1}(:);
    lats = fld{2}(:);
  elseif ( isnumeric(fld) && (size(fld,1) == 2 || size(fld,2) == 2) )
    if ( size(fld,1) == 2 )
      lons = fld(1,:)';
      lats = fld(2,:)';
    else
      lons = fld(:,1);
      lats = fld(:,2);
    end;
  else
    error('First arg must be a STRUCT, 2-CELL of numeric vectors, or numeric 2xN matrix');
  end;

  if ( ~exist('gutter','var') || isempty(gutter) )
    gutter = 0;
  elseif ( ischar(gutter) && ( ...
      strncmpi(gutter,'gutter',3) && strncmpi(gutter,'margin',3) && ...
      strncmpi(gutter,'gridpoint',3) ) )
    gutter = max([max(diff(unique(lons(:)))),max(diff(unique(lats(:))))]);
  elseif ( ~isnumeric(gutter) )
    error('GUTTER must be empty, the string "gutter", or numeric');
  end;

  if ( isscalar(gutter) )
    lgutter = gutter;
    rgutter = gutter;
    bgutter = gutter;
    tgutter = gutter;
  elseif ( numel(gutter) == 2 )
    lgutter = gutter(1);
    rgutter = gutter(1);
    bgutter = gutter(2);
    tgutter = gutter(2);
  elseif ( numel(gutter) == 4 )
    lgutter = gutter(1);
    rgutter = gutter(2);
    bgutter = gutter(3);
    tgutter = gutter(4);
  else
    error('Numeric GUTTER must be a 1-, 2-, or 4-vector');
  end;
    
  bbox = [min(lons)-lgutter max(lons)+rgutter min(lats)-bgutter max(lats)+tgutter];

return;
