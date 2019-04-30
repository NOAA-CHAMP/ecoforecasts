function ixes = bbox2ind(bbox,LON,LAT)
%function ixes = bbox2ind(bbox,LON,LAT)
%
% Return a 4x2 matrix of indices IXES in LON,LAT for bounding box BBOX (must
% be in the form [MINX, MAXX, MINY,MAXY]; indices are for the lon and lat of
% lower-left, lower-right, upper-right, and upper-left corners of BBOX, resp.
% Sizes of all non-singleton dimensions of LON must equal those of LAT. NOTE:
% If LON and LAT are vectors, uses MESHGRID (v.) to put them in plaid form.
% 
% Last Saved Time-stamp: <Fri 2016-05-20 11:10:03 Eastern Daylight Time gramer>

  if ( bbox(2)<=bbox(1) || bbox(4)<=bbox(3) )
    error('Ecoforecasts:BBox:badbbox','BBOX should have form [MINX,MAXX,MINY,MAXY]');
  end;
  if ( ~isnumeric(LON) || ~isnumeric(LAT) || ~ismatrix(LON) || ~ismatrix(LAT) )
    error('Ecoforecasts:BBox:badargs','LON and LAT should be numeric matrices');
  end;

  % Remove any singleton dimensions
  LON = squeeze(LON);
  LAT = squeeze(LAT);

  if ( isvector(LON) && isvector(LAT) )
    warning('Ecoforecasts:BBox:nonplaid','Using MESHGRID(SORT(LON(:)),SORT(LAT(:)))');
    lon = sort(LON(:));
    lat = sort(LAT(:));
    [LON,LAT] = meshgrid(lon,lat);
  elseif ( isvector(LON) || isvector(LAT) )
    error('Ecoforecasts:BBox:badargs','If either LON or LAT is a vector, both must be vectors');
  else
    if ( LON(1,1) == LON(1,2) )
      warning('Ecoforecasts:BBox:transposed','Transposing LON');
      LON = LON';
    end;
    if ( LAT(1,1) == LAT(2,1) )
      warning('Ecoforecasts:BBox:transposed','Transposing LAT');
      LAT = LAT';
    end;
    lon = LON(1,:);
    lat = LAT(:,1);
  end;

  if ( ndims(LON) ~= ndims(LAT) || any(size(LON) ~= size(LAT)) )
    error('Ecoforecasts:BBox:badargs','If matrices, LON and LAT should be equal-sized');
  end;

  maxerr_ul = abs(LON(1,1)-LON(2,2)) + abs(LAT(1,1)-LAT(2,2));
  maxerr_lr = abs(LON(end-1,end-1)-LON(end,end)) + abs(LAT(end-1,end-1)-LAT(end,end));
  maxerr = max(maxerr_ul,maxerr_lr)*sqrt(2);

  [err,ix] = min( abs(LON(:)-bbox(1)) + abs(LAT(:)-bbox(3)) );
  if (err>maxerr+eps); warning('Ecoforecasts:BBox:badbbox','MINX,MINY may be outside region'); end;
  [ixes(1,2),ixes(1,1)] = ind2sub(size(LON),ix);

  [err,ix] = min( abs(LON(:)-bbox(2)) + abs(LAT(:)-bbox(3)) );
  if (err>maxerr+eps); warning('Ecoforecasts:BBox:badbbox','MAXX,MINY may be outside region'); end;
  [ixes(2,2),ixes(2,1)] = ind2sub(size(LON),ix);

  [err,ix] = min( abs(LON(:)-bbox(2)) + abs(LAT(:)-bbox(4)) );
  if (err>maxerr+eps); warning('Ecoforecasts:BBox:badbbox','MAXX,MAXY may be outside region'); end;
  [ixes(3,2),ixes(3,1)] = ind2sub(size(LON),ix);

  [err,ix] = min( abs(LON(:)-bbox(1)) + abs(LAT(:)-bbox(4)) );
  if (err>maxerr+eps); warning('Ecoforecasts:BBox:badbbox','MINX,MAXY may be outside region'); end;
  [ixes(4,2),ixes(4,1)] = ind2sub(size(LON),ix);

return;
