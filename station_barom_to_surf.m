function stn = station_barom_to_surf(stn,pfld,afld,spfld,pz)
%function stn = station_barom_to_surf(stn,pfld,afld,spfld,pz)
%
% Convert station barometric air pressure to sea-level
%
% Last Saved Time-stamp: <Sat 2011-01-29 12:34:03 Eastern Standard Time gramer>

  if ( ~exist('spfld','var') || isempty(spfld) )
    spfld = [pfld '_sealevel'];
  end;
  if ( ~exist('pz','var') || isempty(pz) )
    [wz,az,pz,stz,dtz,slz,dlz] = station_instrument_heights(stn.station_name);
    %DEBUG:    disp(['Using barometer height ' num2str(pz)]);
  end;

  [pix,aix] = intersect_dates(stn.(pfld).date,stn.(afld).date);

  p = stn.(pfld).data(pix);
  a = stn.(afld).data(aix);
  sp = barom_to_surf(p,a,pz);

  stn.(spfld).date = stn.(pfld).date(pix);
  stn.(spfld).data = sp(:);

return;
