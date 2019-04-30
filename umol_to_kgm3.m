function kgm3 = umol_to_kgm3(umol,sp)
%function kgm3 = umol_to_kgm3(umol,sp)
%
% Convert micromolar concentrations of species SP ('N'=DEFAULT, 'P', or 'Si')
% to [kg/m^3], e.g., for upwelling nutrient flux calculations.
%
% From Charles.Fischer@noaa.gov, 2017 Aug 02:
%  14.00 kg-N/m3
%  30.97 kg-P/m3
%  60.09 kg-Si/m3
% E.g.,
%  1 u-mole P/L x 30.97 ug/1u-mole-P x 10e-9 Kg/1ug x 1L/0.001M3 = 3.097 e-5 Kg-P/M3 or 0.03097g-P/M3
%
% Last Saved Time-stamp: <Tue 2018-03-27 17:48:01 Eastern Daylight Time gramer>

  if ( ~exist('sp','var') || isempty(sp) )
    sp = 'N';
  end;
  switch (upper(sp)),
   case 'N',	atwgt = 14.00;
   case 'P',	atwgt = 30.97;
   case 'SI',	atwgt = 60.09;
   otherwise,	error('Unrecognized chemical species %s',sp);
  end;

  kgm3 = umol.*atwgt.*1e-9./0.001;

return;
