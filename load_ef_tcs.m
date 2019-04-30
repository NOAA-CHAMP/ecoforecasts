function ef = load_ef_tcs()
%function ef = load_ef_tcs()
%
% Load ecoforecast model for "thermal-coral-stress"
%
% Last Saved Time-stamp: <Fri 2009-01-30 00:48:03 Eastern Standard Time gramer>
%

  ef.name = 'Thermal-Coral-Stress';
  ef.shortdesc = 'Environmental stress on coral community may be conducive to bleaching';
  ef.desc = [ 'Thermal stress on coral community is extreme' ...
              'enough that it may be conducive to' ...
              'bleaching among multiple coral species.' ...
            ];

  ef.rules.it.facts = { { 'sea1m' } };
  ef.rules.it.fuzzies = {{'drastic-high'}};

%   ef.rules.it.facts = { { 'bleaching_seatemp' } };
%   ef.rules.it.fuzzies = {{'drastic-high'}};

%   ef.rules.mmit.facts = { {'bleaching_seatemp_monthly','seandbc_monthly','sea1m_monthly'} };
%   ef.rules.mmit.fuzzies = {{'high','very-high','drastic-high'}};

return;
