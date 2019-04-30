function ef = load_ef_pas()
%function ef = load_ef_pas()
%
% Load ecoforecast model for "porites-astreoides-spawning"
%
% Last Saved Time-stamp: <Wed 2008-08-13 08:28:40 Eastern Daylight Time gramer>
%

  ef.name = 'Porites-Astreoides-Spawning';
  ef.shortdesc = 'Spawning (planula release) by mustard-hill coral';
  ef.desc = ['Moon phase and environmental conditions are conducive to ' ...
             'spawning of planulae by brooding coral Porites astreoides'];

  ef.tp.facts = {'spawning_seatemp', 'photo_accum'};
  ef.tp.fuzzies = {{'conducive'}, {'high','very-high','drastic-high'}};

return;
