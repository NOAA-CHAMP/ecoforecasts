function [stn,xfld,lfld] = station_reorient_vectors(stn,ori,ufld,vfld,xfld,lfld)
%function [stn,xfld,lfld] = station_reorient_vectors(stn,ori,ufld,vfld,xfld,lfld)
%
% Reorient time series STN.(UFLD) and STN.(VFLD) (x and y vector components)
% to local isobath orientation ORI (degrees T; or fieldname): add time series
% structs XFLD and LFLD for cross-shore and long-shore components, resp. If
% XFLD is absent or empty, tries to replace '_u' in UFLD with '_xshore'. And
% if LFLD is absent or empty, tries to replace '_v' in VFLD with '_lshore'.
% If this fails, then tries to replace '_x' in UFLD with '_xshore', and '_y'
% in VFLD with '_lshore': this is useful for reorienting, e.g., GRADIENTs.
%
% Note STATION_REORIENT_VECTORS *never* overwrites STN.(UFLD) or STN.(VFLD),
% unless caller explicitly respecifies those same names in XFLD and LFLD.
%
% As a convenience for ADCP and section data, if STN.(UFLD) and STN.(VFLD)
% both have .prof fields, also add .prof field to STN.(XFLD) and STN.(LFLD).
%
% Last Saved Time-stamp: <Fri 2016-08-26 14:46:16 Eastern Daylight Time gramer>

  if ( ischar(ori) && isfield(stn,ori) )
    ori = stn.(ori);
  end;
  if ( ~( isnumeric(ori) && isscalar(ori) && -360<=ori && ori<=360 ) )
    error('ORI must be a scalar orientation between -360 and +360 or a field name in STN!');
  end;

  % If caller gave us none, try to munge reasonable fieldnames. But NEVER
  % overwrite U and V fields, unless caller told us to explicitly. NOTE: The
  % UFLD~'_x' and VFLD~'_y' cases are useful for field GRADIENT time series.

  if ( ~exist('xfld','var') || isempty(xfld) )
    xfld = regexprep(ufld,'_u$','_xshore');
    if ( strcmpi(ufld,xfld) || strcmpi(vfld,xfld) )
      xfld = regexprep(ufld,'_u_','_xshore_');
      if ( strcmpi(ufld,xfld) || strcmpi(vfld,xfld) )
        xfld = regexprep(ufld,'_x$','_xshore');
        if ( strcmpi(ufld,xfld) || strcmpi(vfld,xfld) )
          error('Specify a unique XFLD name, or use UFLD again to overwrite!');
        end;
      end;
    end;
  end;
  if ( ~exist('lfld','var') || isempty(lfld) )
    lfld = regexprep(vfld,'_v$','_lshore');
    if ( strcmpi(ufld,lfld) || strcmpi(vfld,lfld) )
      lfld = regexprep(vfld,'_v_','_lshore_');
      if ( strcmpi(ufld,lfld) || strcmpi(vfld,lfld) )
        lfld = regexprep(vfld,'_y$','_lshore');
        if ( strcmpi(ufld,lfld) || strcmpi(vfld,lfld) )
          error('Specify a unique LFLD name, or use VFLD again to overwrite!');
        end;
      end;
    end;
  end;

  [uix,vix] = intersect_dates(stn.(ufld).date,stn.(vfld).date);

  u = stn.(ufld).data(uix);
  v = stn.(vfld).data(vix);

  [x,l] = reorient_vectors(ori,u,v);

  stn.(xfld).date = stn.(ufld).date(uix);
  stn.(xfld).data = x;
  stn.(lfld).date = stn.(vfld).date(vix);
  stn.(lfld).data = l;

  % For ADCP and horizontal (ship section) profiles
  if ( isfield(stn.(ufld),'prof') && isfield(stn.(vfld),'prof') )
    uprof = stn.(ufld).prof(uix,:);
    vprof = stn.(vfld).prof(vix,:);

    [xprof,lprof] = reorient_vectors(ori,uprof,vprof);

    stn.(xfld).prof = xprof;
    stn.(lfld).prof = lprof;
  end;

  if ( nargout < 2 )
    clear xfld lfld;
  end;

return;
