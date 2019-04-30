function sp = barom_to_surf(p,a,pz)
%function sp = barom_to_surf(p,a,pz)
%
% Convert barometric air pressure and air temperature to sea-level pressure
%
% Last Saved Time-stamp: <Sat 2011-01-29 12:32:21 Eastern Standard Time gramer>

  sp = p ./ exp( -pz ./ ( (a + 273.15) .* 29.263 ) );

return;
