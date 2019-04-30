function stns = multi_verify_variable(stns, varname, indices)
%function stns = multi_verify_variable(stns, varname, indices)
%
% This ugly Hack-in-a-Box m-function is necessitated by MATLAB's stupid lack
% of proper support for adding new fields to vectors of structs...
%
% Add an empty field named VARNAME to every struct in vector of structs STNS.
% Then use VERIFY_VARIABLE (qv.) to populate fields STNS(INDICES).(VARNAME)
% with actual derived-variable data as well. (INDICES DEFAULT: All structs.)
%
% NOTE: Once MULTI_VERIFY_VARIABLE is called with numel(INDICES)<numel(STNS)
% for a given VARNAME, the fields STNS(:).(VARNAME) for all *other* indices
% will always *remain empty*. Amazing how ugly MATLAB struct vectors are...
% Anyway, for this reason, use of optional arg INDICES is *not recommended*.
%
% Last Saved Time-stamp: <Tue 2010-02-09 19:51:39 Eastern Standard Time gramer>

  if ( ~exist('indices','var') || isempty(indices) )
    indices = 1:length(stns);
  end;
  % Ugliness required to make simple MATLAB FOR loop (below) behave
  if ( size(indices,1) ~= 1 )
    indices = [indices(:)]';
  end;

  ss = stns;
  % Add blank field to ALL structs in STNS - not just those in INDICES.
  % Do not add .VARNAME to SS (above), as it would break VERFIY_VARIABLE.
  stns(indices(1)).(varname) = [];

  for ix = indices
    s = ss(ix);
    if ( isfield(s,varname) && isempty(s.(varname)) )
      s = rmfield(s, varname);
    end;
    s = verify_variable(s, varname);
    if ( isfield(s,varname) )
      stns(ix) = s;
    end;
    s = []; clear s;
  end;

  ss = []; clear ss;

return;
