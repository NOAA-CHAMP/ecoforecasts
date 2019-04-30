function stns = verify_variable_multi(stns,varname,stnix)
%function stns = verify_variable_multi(stns,varname,stnix)
%
% Calls build_combo_var and build_derived_var (q.v.)
%
% Last Saved Time-stamp: <Sat 2011-04-02 10:36:07  Lew.Gramer>

  if ( ~exist('stnix','var') )
    stnix = 1:numel(stns);
  end;

  allix = stnix(:);
  for ix = 1:numel(stns)
    if ( ismember(ix,allix) )
      stns(ix).(varname) = [];

      warnstat = warning('query','BuildDerivedVar:EmptyField');
      warning('off','BuildDerivedVar:EmptyField');

      % First ensure overarching 'combo' variable is built, if any
      stns(ix) = build_combo_var(stns(ix), varname);

      % Then make sure all other 'derived' variables get built too
      stns(ix) = build_derived_var(stns(ix), varname);

      warning(warnstat);
    end;
  end;

return;
