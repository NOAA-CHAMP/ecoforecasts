function pers = natural_periods(lat)
%function pers = natural_periods([lat])
% Convenience function: returns local tidal, diurnal, ocean mesoscale,
% "weather band", biannual and annual periods in days (e.g., for DATENUM).
% If optional LAT is numeric scalar, also return local INERTIAL_PERIOD (v.)
% NOTE: Ocean "mesoscale" 4d period based on barotropic & baroclinic Rossby
% length and phase speed estimates from Chelton et al (1998) "Geographical
% Variability of the First Baroclinic Rossby Radius of Deformation", JPO. 

  %pers = [0.5 1 4 7 14 42 183 365];
  
  if ( exist('tidefreq') > 1 )
    pers = sort([(1./(tidefreq)./24),1,4,7,42,365.25636/2,365.25636]);
  else
    % K2,S2,M2,N2,K1,diurnal,P1,O1,oc.mesoscale,min.weather,MF,max.weather,biannual,annual
    pers = [0.4986,0.5,0.5175,0.5274,0.9973,1.0,1.0027,1.0758,4,7,13.6604,42,183,365];
  end;
  
  if  ( exist('lat','var') && isscalar(lat) && isnumeric(lat) && abs(lat)<=90 )
    pers = sort([pers,inertial_period(lat)]);
  end;

return;
