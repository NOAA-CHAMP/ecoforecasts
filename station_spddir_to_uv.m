function stn = station_spddir_to_uv(stn,sfld,dfld,ufld,vfld,doCurrents,doOverwrite)
%function stn = station_spddir_to_uv(stn,sfld,dfld,ufld,vfld,doCurrents,doOverwrite)
%
% Calculate U and V components of wind (or if optional DOCURRENTS is True,
% currents) in [m/s], from time series of speed and direction contained in
% structs STN.(SFLD) and STN.(DFLD), resp. (If DOCURRENTS, speed is assumed
% to be also in [m/s], otherwise we assume it is in KTS.) Optional UFLD and
% VFLD field names are constructed from SFLD by replacing "speed" with "u"
% and "v", resp.) If optional DOOVERWRITE is True, UFLD and VFLD may be
% identical to SFLD and DFLD, resp. Otherwise, it is an error for either pair
% of field names to be identical. SEE: SPDDIR_TO_UV, SPDDIR_TO_UV_CURR.
%
% Last Saved Time-stamp: <Fri 2012-05-18 14:55:27  Lew.Gramer>

  if ( ~exist('ufld','var') || isempty(ufld) )
    ufld = strrep(sfld,'speed','u');
  end;
  if ( ~exist('vfld','var') || isempty(vfld) )
    vfld = strrep(sfld,'speed','v');
  end;
  if ( ~exist('doCurrents','var') || isempty(doCurrents) )
    doCurrents = false;
  end;
  if ( ~exist('doOverwrite','var') || isempty(doOverwrite) )
    doOverwrite = false;
  end;

  if ( ~doOverwrite && (strcmp(ufld,sfld) || strcmp(vfld,sfld)) )
    error('Please specify UFLD and VFLD, or else set DOOVERWRITE to True');
  end;

  [s,d] = intersect_tses(stn.(sfld),stn.(dfld));

  stn.(ufld).date = s.date(:);
  stn.(vfld).date = s.date(:);
  if ( doCurrents )
    [stn.(ufld).data,stn.(vfld).data] = spddir_to_uv_curr(s.data(:),d.data(:));
  else
    [stn.(ufld).data,stn.(vfld).data] = spddir_to_uv(s.data(:),d.data(:));
    % Convert knots to meters/sec
    stn.(ufld).data = kts2mps(stn.(ufld).data(:));
    stn.(vfld).data = kts2mps(stn.(vfld).data(:));
  end;

return;
