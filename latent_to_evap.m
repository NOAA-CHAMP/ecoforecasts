function [ev,L_v,rho] = latent_to_evap(lh,t,s,p)
%function [ev,L_v,rho] = latent_to_evap(lh,t,s,p)
%
% Calculate evaporation rate EV [mm/hr] by dividing latent heat flux LH
% [W/m^2] by latent heat of evaporation of pure water L_V (J/kg; DEFAULT:
% 2.44e6) and sea water density RHO [DEFAULT: 1023.3). If sea temperature T
% [oC] is also given and scalar or size(LH)==size(T), calculate L_v from time
% series T. If EITHER salinity S [psu, DEFAULT: 35] OR pressure P [dbar,
% DEFAULT: 5] are given and non-empty, calculate RHO from S,T,P.
%
% Last Saved Time-stamp: <Wed 2010-12-01 18:32:18 Eastern Standard Time gramer>

  L_v = 2.44e6;		%[J kg-1]
  rho = 1023.3;		%[kg m-3]

  if ( exist('t','var') || ~isempty(t) )
    L_v = vapor(t);
    if ( ~exist('s','var') || isempty(s) )
      s = 35;
    end;
    if ( ~exist('p','var') || isempty(p) )
      p = 5;
    end;
    rho = sw_dens(s,t,p);
  end;

  ev = (lh .* 3600) ./ (L_v .* rho);

return;
