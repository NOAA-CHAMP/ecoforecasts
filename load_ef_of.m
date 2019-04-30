function ef = load_ef_of()
%function ef = load_ef_of()
%
% Load ecoforecast model for "onshore-flux"
%
% Last Saved Time-stamp: <Wed 2008-08-20 07:49:25 Eastern Daylight Time gramer>
%

  ef.name = 'Onshore-Flux';
  ef.shortdesc = 'Extremes of temperature, salinity or nutrients being brought over the reef';
  ef.desc = ['Wind- or tide-driven circulation around the reef slope may be ' ...
             'bringing allochthonous water onto the reef crest and back-reef'];

  ef.rules.ft.desc = 'Currents bringing atypical but low-nutrient water onto the reef (normal chlorophyll a + high sea temp. variability + weak wind var.)';
  ef.rules.ft.facts = {'seandbc_variability', 'windsd_7day'};
  ef.rules.ft.fuzzies = {{'high', 'very-high', 'drastic-high'}, ...
                      {'drastic-low','very-low','low'}};

  ef.rules.w.desc = 'Winds bringing atypical but low-nutrient water onto the reef (normal chlorophyll_a + high sea temp. variability + strong wind var.)';
  ef.rules.w.facts = {'seandbc_variability', 'windsd_7day'};
  ef.rules.w.fuzzies = {{'high', 'very-high', 'drastic-high'}, ...
                      {'high', 'very-high', 'drastic-high'}};

  ef.rules.sv.desc = 'Atypical but low-nutrient water near the reef - source uncertain (normal chlorophyll_a + high sea temp. variability + normal wind var.)';
  ef.rules.sv.facts = {'seandbc_variability'};
  ef.rules.sv.fuzzies = {{'high', 'very-high', 'drastic-high'}};

  ef.rules.cft.desc = 'Upwelling? Currents bringing atypical, nutrient-rich water onto the reef (high chlorophyll_a + high sea temp. variability + weak wind var.)';
  ef.rules.cft.facts = {'seandbc_variability', 'windsd_7day', 'sat_chlor_a'};
  ef.rules.cft.fuzzies = {{'high', 'very-high', 'drastic-high'}, ...
                      {'drastic-low','very-low','low'}, ...
                      {'high', 'very-high', 'drastic-high'}};

  ef.rules.ct.desc = 'Upwelling? Atypical, nutrient-rich water near the reef (high chlorophyll_a + high sea temp. variability + average wind var.)';
  ef.rules.ct.facts = {'seandbc_variability', 'sat_chlor_a'};
  ef.rules.ct.fuzzies = {{'high', 'very-high', 'drastic-high'}, ...
                      {'high', 'very-high', 'drastic-high'}};

  ef.rules.cw.desc = 'Upwelling? Winds bringing atypical, nutrient-rich water onto the reef (high chlorophyll_a + high sea temp. variability + strong wind var.)';
  ef.rules.cw.facts = {'seandbc_variability', 'windsd_7day', 'sat_chlor_a'};
  ef.rules.cw.fuzzies = {{'high', 'very-high', 'drastic-high'}, ...
                      {'high', 'very-high', 'drastic-high'}, ...
                      {'high', 'very-high', 'drastic-high'}};

  ef.rules.runoff.desc = 'Contamination? Atypical, nutrient-rich water near the reef (high chlorophyll_a + weak sea temp. variability + weak wind var.)';
  ef.rules.runoff.facts = {'seandbc_variability', 'windsd_7day', 'sat_chlor_a'};
  ef.rules.runoff.fuzzies = {{'drastic-low','very-low','low','somewhat-low','average'}, ...
                      {'drastic-low','very-low','low','somewhat-low','average'}, ...
                      {'high', 'very-high', 'drastic-high'}};

  ef.rules.ch.desc = 'Nutrient delivery of uncertain mechanism (high chlorophyll_a + normal sea temp. variability + normal wind var.)';
  ef.rules.ch.facts = {'sat_chlor_a'};
  ef.rules.ch.fuzzies = {{'high', 'very-high', 'drastic-high'}};

  ef.rules.mix.desc = 'Contamination? Winds bringing atypical, nutrient-rich water onto the reef (high chlorophyll_a + weak sea temp. variability + strong wind var.)';
  ef.rules.mix.facts = {'seandbc_variability', 'windsd_7day', 'sat_chlor_a'};
  ef.rules.mix.fuzzies = {{'drastic-low','very-low','low','somewhat-low','average'}, ...
                      {'high', 'very-high', 'drastic-high'}, ...
                      {'high', 'very-high', 'drastic-high'}};

return;
