function stn = build_factories(stn)
%function stn = build_factories(stn)
%
% Add some additional fact-factories to the SOMBRERO KEY (smkf1) station
% data struct 'stn' - factories that were never implemented in ICON/G2.
% NOTE: This is meant to be called IN ADDITION TO 'load_factories'.
%
% Last Saved Time-stamp: <Sun 2009-03-22 14:01:36 Eastern Daylight Time gramer>
%

  if ( isfield(stn.factories, 'parsurf') )
    stn.factories = rmfield(stn.factories, 'parsurf');
  end;
  stn.factories.parsurf.variable = 'model_surf_par';
  stn = model_surf_par('smkf1', stn, 1988, 2008);
  stn.factories.parsurf.fuzzies = ...
      { ...
          { 'drastic-low',	[ 0	50 ] }, ...
          { 'very-low',		[ 50	100 ] }, ...
          { 'low',		[ 100	300 ] }, ...
          { 'somewhat-low',	[ 300	500 ] }, ...
          { 'average',		[ 500	700 ] }, ...
          { 'somewhat-high',	[ 700	800 ] }, ...
          { 'high',		[ 900	1100 ] }, ...
          { 'very-high',	[ 1100	1300 ] }, ...
          { 'drastic-high',	[ 1300	1600 ] }, ...
      };

return;
