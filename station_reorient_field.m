function [stn,xfld,lfld] = station_reorient_field(stn,ori,fldnm,ufld,vfld,xfld,lfld)
%function [stn,xfld,lfld] = station_reorient_field(stn,ori,fldnm,ufld,vfld,xfld,lfld)
%
% Reorient field time matrices STN.(FLDNM).(UFLD) and STN.(FLDNM).(VFLD) (x
% and y field components of STN.(FLDNM), e.g., gradients) to local isobath
% orientation ORI (degrees T; or fieldname), adding time matrices XFLD and
% LFLD for cross-shore and long-shore field components, resp. If XFLD is
% absent or empty, tries to replace '_u' in UFLD with '_xshore'. And if LFLD
% is absent or empty, tries to replace '_v' in VFLD with '_lshore'. If this
% fails, then tries to replace '_x' in UFLD with '_xshore', and '_y' in VFLD
% with '_lshore': this is useful for reorienting, e.g., GRADIENTs.
%
% STN, ORI, and FLDNM are all required. Other DEFAULTS: UFLD='gradient_x',
% VFLD='gradient_y', XFLD='gradient_xshore', and LFLD='gradient_lshore'.
%
% Note function *never* overwrites STN.(FLDNM).(UFLD) or STN.(FLDNM).(VFLD),
% unless caller explicitly respecifies those same names in XFLD and LFLD.
%
% Last Saved Time-stamp: <Thu 2011-05-05 08:26:50  lew.gramer>

  if ( ischar(ori) && isfield(stn,ori) )
    ori = stn.(ori);
  end;
  if ( ~( isnumeric(ori) && isscalar(ori) && 0<=ori && ori<=360 ) )
    error('ORI must be a scalar orientation 0-360 or a field name in STN!');
  end;

  if ( ~exist('ufld','var') || isempty(ufld) )
    ufld = 'gradient_x';
  end;
  if ( ~exist('vfld','var') || isempty(vfld) )
    vfld = 'gradient_y';
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

  u = stn.(fldnm).(ufld);
  v = stn.(fldnm).(vfld);

  % NOTE: REORIENT_VECTORS preserves input matrix shape, incl. time dimension
  [x,l] = reorient_vectors(ori,u,v);

  stn.(fldnm).(xfld) = x;
  stn.(fldnm).(lfld) = l;

  if ( nargout < 2 )
    clear xfld lfld;
  end;

return;
