function [sp,dr] = ts_uv_to_spddir(u,v,doWind)
%function [sp,dr] = ts_uv_to_spddir(u,v,doWind)
%
% Calculate speed and direction time series structs SP and DR, from U and V
% components of current vector. If optional DOWIND is True, convert speed
% from [m/s] to [Kts.] and convert direction from TARGET to SOURCE. If any of
% .prof, .rawprof, or .field are present in both U and V, also convert those.
% SEE: SPDDIR_TO_UV, SPDDIR_TO_UV_CURR.
%
% Last Saved Time-stamp: <Fri 2018-03-23 14:34:29 Eastern Daylight Time gramer>

  if ( ~is_valid_ts(u) || ~is_valid_ts(v) )
    error('First two args must be valid Time Series structs');
  end;
  if ( ~exist('doWind','var') || isempty(doWind) )
    doWind = false;
  end;

  [u,v] = intersect_tses(u,v);

  sp.date = u.date(:);
  dr.date = u.date(:);

  if ( isfield(u,'data') && isfield(v,'data') )
    sp.data = uv_to_spd(u.data(:),v.data(:));
    if ( ~doWind )
      dr.data = uv_to_dir_curr(u.data(:),v.data(:));
    else
      dr.data = uv_to_dir(u.data(:),v.data(:));
      % Convert meters/sec to knots
      sp.data = mps2kts(sp.data);
    end;
  end;

  if ( isfield(u,'depths') && isfield(v,'depths') )
    sp.depths = u.depths;
    dr.depths = u.depths;
  end;
  if ( isfield(u,'prof') && isfield(v,'prof') )
    sp.prof = uv_to_spd(u.prof,v.prof);
    if ( ~doWind )
      dr.prof = uv_to_dir_curr(u.prof,v.prof);
    else
      dr.prof = uv_to_dir(u.prof,v.prof);
      % Convert meters/sec to knots
      sp.prof = mps2kts(sp.prof);
    end;
  end;
  if ( isfield(u,'rawprof') && isfield(v,'rawprof') )
    sp.rawprof = uv_to_spd(u.rawprof,v.rawprof);
    if ( ~doWind )
      dr.rawprof = uv_to_dir_curr(u.rawprof,v.rawprof);
    else
      dr.rawprof = uv_to_dir(u.rawprof,v.rawprof);
      % Convert meters/sec to knots
      sp.rawprof = mps2kts(sp.rawprof);
    end;
  end;

  if ( isfield(u,'lon') && isfield(v,'lon') )
    sp.lon = u.lon;
    dr.lon = u.lon;
  end;
  if ( isfield(u,'lat') && isfield(v,'lat') )
    sp.lat = u.lat;
    dr.lat = u.lat;
  end;
  if ( isfield(u,'field') && isfield(v,'field') )
    sp.field = uv_to_spd(u.field,v.field);
    if ( ~doWind )
      dr.field = uv_to_dir_curr(u.field,v.field);
    else
      dr.field = uv_to_dir(u.field,v.field);
      % Convert meters/sec to knots
      sp.field = mps2kts(sp.field);
    end;
  end;

return;
