function [bath,rad] = read_hires_bathymetry_for_field(fld,useHighestRes,buff)
%function [bath,rad] = read_hires_bathymetry_for_field(fld,useHighestRes,buff)
%
% Return STRUCT BATH containing a high-resolution bathymetry field centered
% at MEAN(fld.lon(:)), MEAN(fld.lat(:)), of "radius" as close as possible to
% RAD = GEOGRAPHIC_RADIUS_M(fld.lon,fld.lat). CALLS: READ_HIRES_BATHYMETRY.
% FLD may also be a 2x cell array of N-vectors {LONS,LATS}.
%
% If optional BUFF is specified, add that much space (in meters) around the
% edges of bathymetry field. BUFF may be scalar or a 2-vector ([LON,LAT]).
%
%
% Last Saved Time-stamp: <Tue 2017-06-20 11:50:35 Eastern Daylight Time gramer>

  if ( exist('fld','var') && iscell(fld) )
    fldcel = fld; fld=[]; clear fld
    fld.lon = fldcel{1};
    fld.lat = fldcel{2};
  end;
  if ( ~exist('fld','var') || ~isfield(fld,'lon') || ~isfield(fld,'lat') )
    error('First arg must be a STRUCT with .lon,.lat fields or a CELL of {LONS,LATS}');
  end;
  if ( ~exist('useHighestRes','var') || isempty(useHighestRes) )
    useHighestRes = [];
  end;

  [rad(1),rad(2)] = geographic_radius_m(fld.lon,fld.lat);
  % Increase radius slightly
  rad = rad.*1.01;

  if ( exist('buff','var') )
    if ( ~isnumeric(buff) || numel(buff)<1 || numel(buff)>2 )
      error('BUFF if specified must be a numeric scalar or 2-vector');
    end;
    if ( isscalar(buff) )
      rad = rad+buff;
    else
      rad(1) = rad(1)+buff(1);
      rad(2) = rad(2)+buff(2);
    end;
  end;


  lon = mean(fld.lon(:));
  lat = mean(fld.lat(:));
  %[stn,rad] = read_hires_bathymetry(stn_or_stnm_or_locn,rad,bathfiles,useHighestRes,fld)
  [bath,rad] = read_hires_bathymetry([lon,lat],rad,[],useHighestRes);

return;
