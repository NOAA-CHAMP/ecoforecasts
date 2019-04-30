function d = relhumid_to_dewp(t,q)
%function d = relhumid_to_dewp(t,q)
%
% August-Roche-Magnus approximation (wiki/Dewpoint, or v. Barenbrug,
% A. W. T. (1974), Psychrometry and Psychrometric Charts, 3rd ed., Cape and
% Transvaal, Cape Town, South Africa.)
%
% NOTE: This differs from the simpler algorithm currently recommended by NDBC!
%
% 't' is air temperature (scalar, vector, matrix), 'q' is the Relative Humidity
% (e.g., 100.0 for 100%) taken at same elevation. Size(t) must equal size(q).
% Returned value 'd' is dew point temperature - scalar, vector or matrix of
% size identical to inputs.
%
% Last Saved Time-stamp: <Wed 2010-01-13 09:36:25 Eastern Standard Time Lew.Gramer>

  c1 = 17.271; c2 = 237.7;
  gamma = ((c1.*t) ./ (c2+t)) + log(q./100);
  d = (c2.*gamma) ./ (c1 - gamma);

return;
