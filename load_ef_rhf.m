function ef = load_ef_rhf()
%function ef = load_ef_rhf()
%
% Load ecoforecast model for "rapid-heat-flux"
%
% Last Saved Time-stamp: <Mon 2009-03-23 16:42:51 Eastern Daylight Time gramer>
%

  ef.name = 'Rapid-Heat-Flux';
  ef.shortdesc = 'Extremes of temperature over the reef';
  ef.desc = ['Circulation or high surface flux, causing high sea ' ...
             'temperature variance over the crest and back-reef'];

  % ---
  ef.rules.hww.desc = 'Offshore circulation affecting the reef (high sea temp. variability + weak air temp. var. + normal-to-weak wind var.)';
  ef.rules.hww.facts = {'seandbc_variability','airt_variability','windsd_7day'};
  ef.rules.hww.fuzzies = {{'high','very-high','drastic-high'}, ...
                      {'drastic-low','very-low','low'}, ...
                      {'drastic-low','very-low','low','somewhat-low','average','somewhat-high'}};

  % ---
  ef.rules.hwh.desc = 'Wind-driven circulation affecting the reef (high sea temp. variability + weak air temp. var. + high wind var.)';
  ef.rules.hwh.facts = {'seandbc_variability','airt_variability','windsd_7day'};
  ef.rules.hwh.fuzzies = {{'high','very-high','drastic-high'}, ...
                      {'drastic-low','very-low','low'}, ...
                      {'high','very-high','drastic-high'}};

  % ---
  ef.rules.hnh.desc = 'Combination of wind-driven circulation and surface flux (high sea temp. variability + high-to-normal air temp. var. + high wind var.)';
  ef.rules.hnh.facts = {'seandbc_variability','airt_variability','windsd_7day'};
  ef.rules.hnh.fuzzies = {{'high','very-high','drastic-high'}, ...
                      {'somewhat-low','average','somewhat-high','high','very-high','drastic-high'}, ...
                      {'high','very-high','drastic-high'}};

  % ---
  ef.rules.hhn.desc = 'Surface heat flux, possible offshore circulation over reef (high sea temp. variability + high air temp. var. + normal-to-weak wind var.)';
  ef.rules.hhn.facts = {'seandbc_variability','airt_variability','windsd_7day'};
  ef.rules.hhn.fuzzies = {{'high','very-high','drastic-high'}, ...
                      {'high','very-high','drastic-high'}, ...
                      {'drastic-low','very-low','low','somewhat-low','average','somewhat-high'}};

  % ---
  ef.rules.hnn.desc = 'Offshore circulation, possible surface heat flux over reef (high sea temp. variability + normal air temp. var. + normal-to-weak wind var.)';
  ef.rules.hnn.facts = {'seandbc_variability','airt_variability','windsd_7day'};
  ef.rules.hnn.fuzzies = {{'high','very-high','drastic-high'}, ...
                      {'somewhat-low','average','somewhat-high'}, ...
                      {'drastic-low','very-low','low','somewhat-low','average','somewhat-high'}};

  % ---
  ef.rules.huu.desc = 'Insufficient data to attribute cause of reef variability (high sea temp. variability + unknown air temp. var. + unknown wind var.)';
  ef.rules.huu.facts = {'seandbc_variability'};
  ef.rules.huu.fuzzies = {{'high','very-high','drastic-high'}};


return;
