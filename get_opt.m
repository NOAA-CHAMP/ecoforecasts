function val = get_opt(opts,nm,dflt)
%function val = get_opt(opts,nm,dflt)
%
% Return the value VAL of the option field OPTS.(NM). If there is no field
% field named NM in OPTS, either return DFLT; or if DFLT is also a struct
% with field DFLT.(NM), return that value; or if DFLT is not given, VAL=[].
%
% Last Saved Time-stamp: <Thu 2011-03-24 10:32:06  lew.gramer>

  if ( isfield(opts,nm) )
    val = opts.(nm);
  elseif ( nargin < 3 )
    val = [];
  elseif ( isfield(dflt,nm) )
    val = dflt.(nm);
  else
    val = dflt;
  end;

return;
