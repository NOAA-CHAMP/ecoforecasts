function [nearestCoords,d,az,locus,idx,locusZn,fromZn] = get_nearest_points(locus_or_bath,fromCoords,locusZn,fromZn)
%function [nearestCoords,d,az,locus,idx,locusZn,fromZn] = get_nearest_points(locus_or_bath,fromCoords,locusZn,fromZn)
%
% For each point in 2xN coordinate matrix FROMCOORDS ([LON;LAT]), find the
% point nearest it in a locus of points, e.g., the shore line. LOCUS_OR_BATH
% is either locus of points (2xN LON;LAT coordinate matrix) or a bathymetry
% field (STRUCT with .lon,.lat,.field); if the latter, the 0 m isobath is
% used as the locus of points (v. GET_ISOBATH). Returns the 2xN [LON;LAT]
% coordinates matrix, the full LOCUS of points used, index IDX of points in
% LOCUS nearest to each point in FROMCOORDS, and distance (NOTE: in KM) from
% each point in fromCoords to its nearest point in LOCUS.
%
% If LOCUS_OR_BATH is a CELL with a STRUCT first element and a numeric second 
% element, then that numeric 2nd elt. is passed to GET_ISOBATH (see above).
%
% If it is possible to convert LOCUS and FROMCOORDINATES into UTM in the same
% zone, call KNNSEARCH (v.) which is more efficient. Otherwise, we must loop
% over all points in LOCUS, which can be much slower.
%
% If LOCUSZN is a non-empty UTM Zone, assume LOCUS is already in UTM (not
% LON;LAT) for that zone. If FROMZN is non-empty, assume FROMCOORDS is UTM.
%
% In all cases, returned NEARESTCOORDS and LOCUS are *always* LON;LAT.
%
% CALLS: DISTANCE_WGS84, GET_ISOBATH, KNNSEARCH, LATLON2UTM
%
% Last Saved Time-stamp: <Thu 2019-02-07 16:44:41 Eastern Standard Time gramer>

  iso = 0;
  if ( iscell(locus_or_bath) )
    if ( isstruct(locus_or_bath{1}) && isnumeric(locus_or_bath{2}) )
      iso = locus_or_bath{2};
      locus_or_bath = locus_or_bath{1};
    else
      error('If first arg is a CELL, it must be {BATHY_STRUCT,ISOBATH}');
    end;
  end;
  if ( isstruct(locus_or_bath) )
    locus = get_isobath(locus_or_bath,iso);
  elseif ( isnumeric(locus_or_bath) && (size(locus_or_bath,1) == 2) )
    locus = locus_or_bath;
  else
    error('First arg LOCUS_OR_BATH should be a STRUCT or 2xN coordinate vector');
  end;
  clear locus_or_bath

  if ( ~exist('fromCoords','var') || ~isnumeric(fromCoords) || size(fromCoords,1) ~= 2 )
    error('Second arg FROMCOORDS should be 2xN coordinate (LON;LAT or X;Y) vector!');
  end;
  if ( ~exist('locusZn','var') || isempty(locusZn) )
    locusZn = [];
  end;
  if ( ~exist('fromZn','var') || isempty(fromZn) )
    fromZn = [];
  end;
  if ( xor(isempty(locusZn),isempty(fromZn)) )
    error('If either LOCUS or FROMCOORDS is in UTM, both must be');
  end;

  %  Special handling required for LAT,LON as (LON1-LON2)^2~=DISTANCE_WGS84!
  locusUTM = locus;
  if ( isempty(locusZn) )
    [locusUTM(1,:),locusUTM(2,:),locusZn] = latlon2utm(locus(2,:),locus(1,:));
    locusZn = string(locusZn);
  end;
  fromCoordsUTM = fromCoords;
  if ( isempty(fromZn) )
    [fromCoordsUTM(1,:),fromCoordsUTM(2,:),fromZn] = latlon2utm(fromCoords(2,:),fromCoords(1,:));
    fromZn = string(fromZn);
  end;

  if ( numel(unique(locusZn)) == 1 && numel(unique(fromZn)) == 1 && locusZn(1) == fromZn(1) )
    % Nearest-neighbor search (works on co-zonal UTM coordinates)
    [idx,ds] = knnsearch(locusUTM',fromCoordsUTM');

    % % An alternative approach from: https://www.mathworks.com/matlabcentral/answers/424346-calculate-distance-between-a-point-and-a-contour-line
    % % Calculate closest distance between trajectory and contour
    % [K,d] = dsearchn([Clon,Clat],[lon,lat]);
    % % Check if the pts are inside or outside of the contour
    % in = inpolygon(lon,lat,poly(:,1),poly(:,2));
    % % Change sign of those inside
    % d(in) = -d(in);
    
    [locus(2,:),locus(1,:)] = utm2latlon(locusUTM(1,:),locusUTM(2,:),locusZn{1});
    [fromCoords(2,:),fromCoords(1,:)] = utm2latlon(fromCoordsUTM(1,:),fromCoordsUTM(2,:),fromZn{1});
    
    [d,az] = distance_wgs84(fromCoords(2,:),fromCoords(1,:),locus(2,idx),locus(1,idx));

  else
    warning('Multiple UTM zones! This may take a while...');
    %  Brute force required across UTM zones because (LON1-LON2)^2~=DISTANCE_WGS84!
    d = repmat(nan,[1,size(fromCoords,2)]);
    idx = repmat(nan,[1,size(fromCoords,2)]);
    az = repmat(nan,[1,size(fromCoords,2)]);
    for ix = 1:size(fromCoords,2);
      [ds,azs] = distance_wgs84(fromCoords(2,ix),fromCoords(1,ix),locus(2,:),locus(1,:));
      [d(ix),idx(ix)] = min(ds);
      az(ix) = azs(idx(ix));
    end;
  end;

  nearestCoords = locus(:,idx);

return;
