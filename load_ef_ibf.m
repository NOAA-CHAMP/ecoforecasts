function ef = load_ef_ibf()
%function ef = load_ef_ibf()
%
% Load ecoforecast model for "internal-bore-flux"
%
% Last Saved Time-stamp: <Wed 2008-08-20 13:48:31 Eastern Daylight Time gramer>
%

  ef.name = 'Internal-Bore-Flux';
  ef.shortdesc = 'Internal tidal bore breaking onto the reef';
  ef.desc = ['Tide-driven internal waves breaking on the reef slope may be ' ...
             'bringing allochthonous water onto the reef crest and back-reef'];

  ef.rules.tvw.desc = 'Possible nutrient flux from internal wave breaking (high sea temp. variability + low wind)';
  ef.rules.tvw.facts = {'seandbc_variability', 'windsp_3day'};
  ef.rules.tvw.fuzzies = {{'high', 'very-high', 'drastic-high'}, ...
                    {'drastic-low','very-low','low','somewhat-low','average'}};

  ef.rules.tvek.desc = 'Possible nutrient flux from internal wave breaking (high sea temp. variability + low wind-driven Ekman flux)';
  ef.rules.tvek.facts = {'seandbc_variability', 'ekmandir_7day', 'windsp_3day'};
  ef.rules.tvek.fuzzies = {{'high', 'very-high', 'drastic-high'}, ...
                     {'downshore', 'upshore'}, ...
                     {'drastic-high','very-high','high','somewhat-high'}};

return;
