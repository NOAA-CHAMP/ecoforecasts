function stn = station_insol_to_par(stn,insfld,parfld,parWfld,goodix)
%function stn = station_insol_to_par(stn,insfld,parfld,parWfld,goodix)
%
% Calculate new in situ (or other) PAR light time series field STN.(PARFLD)
% in quantum units (micromol quanta/m^2/s), and optional PAR in energy units
% [W/m^2] STN.(PARWFLD), from insolation (downward shortwave radiative flux)
% STN.(INSFLD), using published conversion constants. If GOODIDX is given,
% use only those indices from time series STN.(INSFLD) to estimate PAR. If
% arg PARWFLD is absent or empty, W/m^2 time series is not added to STN.
% CALLS: INSOL_TO_PAR.
%
% Last Saved Time-stamp: <Sat 2011-10-22 18:10:43  lew.gramer>

  if ( ~exist('goodix','var') || isempty(goodix) )
    goodix = 1:length(stn.(insfld).date);
  end;

  stn.(parfld).date = stn.(insfld).date(goodix);
  [stn.(parfld).data,w] = insol_to_par(stn.(insfld).data(goodix));

  if ( exist('parWfld','var') && ~isempty(parWfld) )
    stn.(parWfld).date = stn.(insfld).date(goodix);
    stn.(parWfld).data = w;
  else
    clear w;
  end;

return;
