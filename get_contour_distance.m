function [dx,az,nearestCoords] = get_contour_distance(contourCoords_or_bath,fromCoords)
%function [dx,az,nearestCoords] = get_contour_distance(contourCoords_or_bath,fromCoords)
%
% Get the distance in KM between each point in FROMCOORDS (a 2xN matrix
% LON;LAT, OR a matrix of STRUCT with .lon and .lat files) and its nearest
% point in contour or locus of points CONTOURCOORDS; first arg may also be a
% bathymetry STRUCT with fields .lon,.lat,.field, or a CELL array containing
% {BATHY_STRUCT,NUMERIC_ISOBATH}. In any case, the first arg is passed
% through untouched to GET_NEAREST_POINTS (v.)
%
% Optionally also returns AZ azimuths [deg True] to shore, and NEARESTCOORDS
% the nearest points on the contour from each point. DX and AZ should be of
% size 1xN, NEARESTCOORDS of size 2xN.
%
% Last Saved Time-stamp: <Wed 2019-02-06 16:50:51 EST lew.gramer>

  if ( isstruct(fromCoords) && isfield(fromCoords,'lon') && isfield(fromCoords,'lat') )
    lon = [fromCoords(:).lon];
    lat = [fromCoords(:).lat];
    clear fromCoords
    fromCoords(1,1:numel(lon)) = lon;
    fromCoords(2,1:numel(lat)) = lat;
  elseif ( ~isnumeric(fromCoords) || ~ismatrix(fromCoords) || size(fromCoords,1) ~= 2 )
    error('Second arg FROMCOORDS should be a field STRUCT or 2xN coordinate matrix');
  end;

  [nearestCoords,dx,az] = get_nearest_points(contourCoords_or_bath,fromCoords);

return;
