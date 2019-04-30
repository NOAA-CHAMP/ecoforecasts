function [p,w] = insol_to_par(g)
%function [p,w] = insol_to_par(g)
%
% G is a vector or scalar of insolation in Watts/m^2. Returned vector or
% scalar P is PAR in micromol quanta / m^2 / s, W is PAR in Watts/m^2.
% Conversion factor from Insolation to PAR is derived from Jacovides et al.,
% 2003. Conversion factor from W/m^2 to micromol quanta/m^2.s (micro-E/J) is
% derived from Dye, 2004.
%
% SEE ALSO: PAR_TO_INSOL
%
% Last Saved Time-stamp: <Fri 2011-12-16 13:11:37  lew.gramer>

  % Conversion factor from Insolation to PAR

  % % Per Papaioannou, Papanikolaou and Retalis, 1993
  % PAR_PER_INSOL = 0.473;

  % % Mean modeled J/micro-E for 4cm toal precipitable water, per Gonzalez and Calbo, 2002
  % PAR_PER_INSOL = 0.485;

  % % Mean measured J/micro-E, per Gonzalez and Calbo, 2002
  % PAR_PER_INSOL = 0.503;

  % % Hourly mean for 41% median cloud cover, per Jacovides et al., 2003
  % PAR_PER_INSOL = 0.477;

  % Hourly mean for overcast skies, per Jacovides et al., 2003
  PAR_PER_INSOL = 0.501;


  % Conversion factor from W/m^2 to micromole quanta/m^2.s

  % % Per Morel and Smith, 1974
  % UMOL_PER_WATT = 4.1513469579;

  % Dye (JGR-A, 2004) recommends 4.56 mumol/J (quantum-to-energy conversion)
  % for PAR "for a wide range of cloud conditions with little or no  error"
  UMOL_PER_WATT = 4.56;

  % Consistent with PPFD/R_b of 2.2846 Einstein/J, see Gonzalez and Calbo (2002)
  UMOL_PAR_PER_INSOL = PAR_PER_INSOL * UMOL_PAR_WATT;

  w = g .* PAR_PER_INSOL;
  p = w .* UMOL_PER_WATT;

  if ( nargout < 2 )
    clear w;
  end;

return;
