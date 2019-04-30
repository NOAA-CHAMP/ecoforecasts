function stn = verify_variable(stn, varname, forceRebuild)
%function stn = verify_variable(stn, varname, forceRebuild)
%
% Calls BUILD_COMBO_VAR and BUILD_DERIVED_VAR (q.v.). If the optional arg
% FORCEREBUILD is True, first removes field VARNAME from STN if it exists.
% VARNAME may be a cellstr, in which case multiple variables are verified.
%
% Last Saved Time-stamp: <Wed 2011-12-28 13:07:50  lew.gramer>

  if ( ~exist('forceRebuild','var') || isempty(forceRebuild) )
    forceRebuild = false;
  end;

  if ( iscellstr(varname) )
    % For multiple variable-name inputs, recurse...
    for cvarname=varname(:)'
      stn = verify_variable(stn,cvarname{:},forceRebuild);
    end;

  else
    if ( forceRebuild && isfield(stn,varname) )
      stn = rmfield(stn,varname);
    end;

    % First ensure overarching 'combo' variable is built, if any
    stn = build_combo_var(stn, varname);

    % Then make sure all other 'derived' variables get built too
    stn = build_derived_var(stn, varname);

    % (Note: Building a Combo Var automatically forces its second input
    % variable to be built (if it is not yet). Separate build_derived_var
    % call above ensures any vars derived FROM a combo var also get built.)
  end;

return;
