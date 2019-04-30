function [u,v] = ts_spddir_to_uv(sp,dr,doWind)
%function [u,v] = ts_spddir_to_uv(sp,dr,doWind)
%
% Calculate U and V components of current vector from speed and direction
% time series structs SP and DR. If optional DOWIND is True, convert speeds
% from [Kts.] to [m/s] and convert direction from TARGET to SOURCE. If either
% .prof or .field are present in both SP and DR, also calculate those in U,V.
% SEE ALSO: SPDDIR_TO_UV, UV_TO_SPDDIR, SPDDIR_TO_UV_CURR, UV_TO_SPDDIR_CURR.
%
% Last Saved Time-stamp: <Fri 2018-03-23 14:35:37 Eastern Daylight Time gramer>

  if ( ~is_valid_ts(sp) || ~is_valid_ts(dr) )
    error('First two args must be valid Time Series structs');
  end;
  if ( ~exist('doWind','var') || isempty(doWind) )
    doWind = false;
  end;

  [sp,dr] = intersect_tses(sp,dr);

  u.date = sp.date(:);
  v.date = sp.date(:);

  if ( isfield(sp,'data') && isfield(dr,'data') )
    if ( ~doWind )
      [u.data,v.data] = spddir_to_uv_curr(sp.data(:),dr.data(:));
    else
      [u.data,v.data] = spddir_to_uv(sp.data(:),dr.data(:));
      % Convert knots to meters/sec
      u.data = kts2mps(u.data);
      v.data = kts2mps(v.data);
    end;
  end;

  if ( isfield(sp,'depths') && isfield(dr,'depths') )
    u.depths = sp.depths;
    v.depths = sp.depths;
  end;
  if ( isfield(sp,'prof') && isfield(dr,'prof') )
    if ( ~doWind )
      [u.prof,v.prof] = spddir_to_uv_curr(sp.prof,dr.prof);
    else
      [u.prof,v.prof] = spddir_to_uv(sp.prof,dr.prof);
      % Convert knots to meters/sec
      u.prof = kts2mps(u.prof);
      v.prof = kts2mps(v.prof);
    end;
  end;
  if ( isfield(sp,'rawprof') && isfield(dr,'rawprof') )
    if ( ~doWind )
      [u.rawprof,v.rawprof] = spddir_to_uv_curr(sp.rawprof,dr.rawprof);
    else
      [u.rawprof,v.rawprof] = spddir_to_uv(sp.rawprof,dr.rawprof);
      % Convert knots to meters/sec
      u.rawprof = kts2mps(u.rawprof);
      v.rawprof = kts2mps(v.rawprof);
    end;
  end;

  if ( isfield(sp,'lon') && isfield(dr,'lon') )
    u.lon = sp.lon;
    v.lon = sp.lon;
  end;
  if ( isfield(sp,'lat') && isfield(dr,'lat') )
    u.lat = sp.lat;
    v.lat = sp.lat;
  end;
  if ( isfield(sp,'field') && isfield(dr,'field') )
    if ( ~doWind )
      [u.field,v.field] = spddir_to_uv_curr(sp.field,dr.field);
    else
      [u.field,v.field] = spddir_to_uv(sp.field,dr.field);
      % Convert knots to meters/sec
      u.field = kts2mps(u.field);
      v.field = kts2mps(v.field);
    end;
  end;

return;
