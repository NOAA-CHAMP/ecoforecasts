function [cutoffs,ranges,FUZZIES] = create_prctile_fuzzy(stn,varname,PRCTILES,FUZZIES)
%function [cutoffs,ranges,FUZZIES] = create_prctile_fuzzy(stn,varname,PRCTILES,FUZZIES)
%
% Use percentiles (v. PRCTILE, MATLAB Statistics Toolbox) to generate simple
% fact fuzzy ranges for variable VARNAME in station STN. PRCTILES: optional
% 2nd arg to PRCTILE, DEFAULT: [00.62 02.27 06.68 15.87 30.85 69.15 84.13 93.32 97.72 99.38]
%
% Last Saved Time-stamp: <Sun 2010-02-14 18:47:43 Eastern Standard Time gramer>

% SYMBOLIC VALUE                  LOWER BOUND PERCENTILE
% 'unbelievably-low',             0
% 'drastic-low',                  00.62
% 'very-low',                     02.27
% 'low',                          06.68
% 'somewhat-low',                 15.87
% 'average',                      30.85
% 'somewhat-high',                69.15
% 'high',                         84.13
% 'very-high',                    93.32
% 'drastic-high',                 97.72
% 'unbelievably-high',            99.38


  % Works well for normal (unimodal) distributions
  if ( ~exist('PRCTILES','var') || isempty(PRCTILES) )
    PRCTILES = [00.62 02.27 06.68 15.87 30.85 69.15 84.13 93.32 97.72 99.38];
  end;

  if ( ~exist('FUZZIES','var') || isempty(FUZZIES) )
    FUZZIES = { ...
        'drastic-low', ...
        'very-low', ...
        'low', ...
        'somewhat-low', ...
        'average', ...
        'somewhat-high', ...
        'high', ...
        'very-high', ...
        'drastic-high', ...
              };
  end;

  cutoffs = prctile(stn.(varname).data, PRCTILES);
  for cix = 1:(length(cutoffs)-1)
    ranges{cix} = { FUZZIES{cix} , cutoffs([cix cix+1]) };
  end;

return;
