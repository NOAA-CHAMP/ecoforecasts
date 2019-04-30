function [kgm3,umol] = seatemp_to_nutes(T,sp)
%function [kgm3,umol] = seatemp_to_nutes(T[,sp])
%
% Use the numeric array of sea temperatures T to calculate mass of nutrient
% species SP (DEFAULT: 'N') within a cubic meter of water, using an empirical
% formula of Gramer et al. (2017).
%
% Last Saved Time-stamp: <Sun 2018-08-26 15:42:52 Eastern Daylight Time gramer>

  nutes = [];

  if ( ~exist('sp','var') || isempty(sp) )
    sp = 'N';
  end;

  switch (upper(sp)),
   case 'N',	a = 1.620;	b = 0.0;	eta = 2.420;
   case 'P',	a = 0.073;	b = 0.0;	eta = 0.132;
   case 'Si',	a = 0.370;	b = 0.0;	eta = 0.790;
   otherwise,	error('Unrecognized chemical species "%s"',sp);
  end;

  Tub = 23;
  Tlb = 12;

  T(T<Tlb) = Tlb;
  umol = ( (a.*(Tub-T)) + b );
  umol(umol < 0) = 0;

  kgm3 = umol_to_kgm3(umol,sp);

return;
