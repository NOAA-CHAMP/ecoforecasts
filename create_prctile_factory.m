function fac = create_prctile_factory(stn,varname,PRCTILES,FUZZIES)
%function fac = create_prctile_factory(stn,varname,PRCTILES,FUZZIES)
%
% Use percentiles (v. PRCTILE, MATLAB Statistics Toolbox) to generate simple
% fact fuzzy ranges for variable VARNAME in station STN. Optional 2nd and 3rd
% args PRCTILE, FUZZIES: DEFAULTs are per CREATE_PRCTILE_FUZZY. (qv.) 
%
% Last Saved Time-stamp: <Sun 2010-02-14 19:25:39 Eastern Standard Time gramer>

  if ( ~exist('PRCTILES','var') )
    PRCTILES = [];
  end;
  if ( ~exist('FUZZIES','var') )
    FUZZIES = [];
  end;

  fac.variable = varname;

  [cutoffs,ranges] = create_prctile_fuzzy(stn,varname,PRCTILES,FUZZIES);
  fac.fuzzies = ranges;

return;
