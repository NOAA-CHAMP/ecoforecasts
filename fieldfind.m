function [lonix,latix,rmse] = fieldfind(fldlon,fldlat,lon,lat)
%function [lonix,latix,rmse] = fieldfind(fldlon,fldlat,lon,lat)
%
% Return the longitude and latitude indices of position(s) LON,LAT within the
% 2-D field having coordinate grid FLDLON,FLDLAT (vectors or plaid matrices).
% If FLDLON,FLDLAT plaid, LON and LAT must have the same number of elements:
% in that case, the root-mean-squared error RMSE can also be returned.
%
% Last Saved Time-stamp: <Fri 2017-04-28 13:23:18 Eastern Daylight Time gramer>

  lonix=[];
  latix=[];

  if ( isvector(fldlon) )
    if ( ~isvector(fldlat) )
      error('If either of FLDLON,FLDLAT is a vector, both must be');
    end;
    lonix = interp1(fldlon,1:numel(fldlon),lon,'nearest',nan);
    latix = interp1(fldlat,1:numel(fldlat),lat,'nearest',nan);

  elseif ( ndims(fldlon) ~= ndims(fldlat) || any(size(fldlon) ~= size(fldlat)) )
    error('If matrices, FLDLON and FLDLAT should be plaid');

  elseif ( numel(lon) ~= numel(lat) )
    error('If FLDLON and FLDLAT are plaid, LON and LAT must have equal NUMEL');

  else
    % WARNING: THIS LOOP MAY BE VERY SLOOOOOOOOOW...
    for ix=1:numel(lon)
      [rmse(ix),fldix] = min(sqrt( ((fldlon(:)-lon(ix)).^2) ...
                                  + ((fldlat(:)-lat(ix)).^2) ));
      if ( 1 > fldix || fldix > numel(fldlon) )
        lonix(ix) = nan;
        latix(ix) = nan;
      else
        [lonix(ix),latix(ix)] = ind2sub(size(fldlon),fldix);
      end;
    end; %for ix=1:numel(lon)
  end; %if ( isvector(fldlon) ) elseif else

return;
