function fix_varnamelengths(vars)
%function fix_varnamelengths(vars)
%
% *CAUTION*: Calling this m-function may have the side effect of changing the
% contents of *every string* (char vector) in the calling workspace! Only use
% if your workspace is full of long variable-name variables needing shortened.
%
% Evaluate all variables from the cellstr VARS in the 'caller' workspace, and
% fix any UNIQUE character string contents to meet MATLAB valid variable name
% criteria with GENVARNAME (v.). If omitted, VARS = EVALIN('caller','who').
%
% Last Saved Time-stamp: <Wed 2011-07-13 09:04:24  lew.gramer>

  if ( ~exist('vars','var') || isempty(vars) )
    vars = evalin('caller','who');
  end;

  charvars = {};
  charvals = {};
  for ix = 1:length(vars)
    var = vars{ix};
    val = evalin('caller',var);
    if ( ischar(val) && isvector(val) && ~iscell(val) )
      charvars{end+1} = var;
      charvals{end+1} = val;
    end;
  end;

  [ig,uniqix,nonuniqix] = unique(charvals);
  newcharvals = genvarname(charvals(uniqix));
  newcharvals = newcharvals(nonuniqix);

  for ix = 1:length(charvars)
    assignin('caller',charvars{ix},newcharvals{ix});
  end;

return;
