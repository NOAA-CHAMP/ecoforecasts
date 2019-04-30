function q = dewp_to_relhumid(t,d)
%function q = dewp_to_relhumid(t,d)
%
% INVERTED form of August-Roche-Magnus approximation (wiki/Dewpoint, or
% v. Barenbrug, A. W. T. (1974), Psychrometry and Psychrometric Charts, 3rd
% ed., Cape and Transvaal, Cape Town, South Africa.)
%
% NOTE: This differs from the simpler algorithm currently recommended by NDBC!
%
% 't' is air temperature (scalar, vector, matrix), 'd' is the dew point
% temperature taken at same elevation. Size(t) must equal size(d). Returned
% value 'q' is the Relative Humidity (e.g., 100.0 for 100%) - scalar, vector
% or matrix of size identical to inputs.
%
% Last Saved Time-stamp: <Wed 2010-01-13 09:36:19 Eastern Standard Time Lew.Gramer>

  c1 = 17.271; c2 = 237.7;
  gamma = (c1.*d) ./ (d + c2);
  q = 100 .* exp(gamma - ((c1.*t) ./ (c2 + t)));

  % Do basic QA on result
  q(q > 100) = nan;

return;
