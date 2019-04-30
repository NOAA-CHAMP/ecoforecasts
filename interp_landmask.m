function fld = interp_landmask(fld,fldnm)
%function fld = interp_landmask(fld,fldnm)
%
% Take a gappy field (e.g., one extracted from a satellite product with an
% overactive landmask) and interpolate it to a gapless field using Spline
% interpolation of all of the nearest points. LANDMASK points are selected by
% finding every index in the field for which ALL values are NaN; for a 2D
% field (e.g., bathymetry), landmask is any NaN value.
%
% Last Saved Time-stamp: <Sat 2015-08-29 22:35:40 Eastern Daylight Time gramer>

  if ( ndims(fld.(fldnm)) == 2 )
    % Spline-fit the 2D field
    [xix,yix]=meshgrid(fld.lat,fld.lon);
    s = warning('OFF','MATLAB:chckxy:IgnoreNaN');  % Ignoring NaNs is the point
    newfld = interp2(fld.lat,fld.lon,fld.(fldnm),xix,yix,'spline',nan);
    warning(s.state,'MATLAB:chckxy:IgnoreNaN');

    % Fill in NaNs with spline-fit values
    [landrix,landcix] = find(~isfinite(squeeze(nanmean(fld.(fldnm)))));
    fld.(fldnm)(:,landrix,landcix) = newfld(:,landrix,landcix);

  elseif ( ndims(fld.(fldnm)) == 3 && size(fld.(fldnm),1) > 1 )
    % Spline-fit the whole 3D field. Note this also interpolates missing dates
    % at non-land pixels, which we do not want for purposes of this function.
    [xix,yix,zix]=meshgrid(fld.lat,fld.date,fld.lon);
    s = warning('OFF','MATLAB:chckxy:IgnoreNaN');  % Ignoring NaNs is the point
    newfld = interp3(fld.lat,fld.date,fld.lon,fld.(fldnm),xix,yix,zix,'spline',nan);
    warning(s.state,'MATLAB:chckxy:IgnoreNaN');

    % Find "landmask" pixels and fill in with spline-fit values
    [landrix,landcix] = find(~isfinite(squeeze(nanmean(fld.(fldnm)))));
    fld.(fldnm)(:,landrix,landcix) = newfld(:,landrix,landcix);

  else
    error('Arg FLD.(FLDNM) must be a 2-D or (non-trivially) 3-D array');
  end;

return;
