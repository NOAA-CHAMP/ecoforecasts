function [wvper,wvhgt] = wind_to_wave(windspd)
%function [wvper,wvhgt] = wind_to_wave(windspd)
%
% Estimate wind-wave height WVHGT, period WVPER from wind speed vector WINDSPD.
% Formula from Fairall MATLAB script COR3_0AH.M. (Tp from Pierson-Moskowitz
% relation, significant wave height from one of several relations reviewed in
% Huang et al., 1990; both referenced in e.g., Bourassa et al, 2001. Taylor
% and Yelland, 2001, alternate formulation for Hs was also considered.)
%
% NOTE: Assumes WINDSPD is in [kts], and converts to [m/s].
%
% Last Saved Time-stamp: <Tue 2010-05-04 14:48:37 Eastern Daylight Time Lew.Gramer>

  % Convert winds from [kts] to [m/s]
  windspd = windspd .* 0.5144444444;

  Ca = 0.018;
  % Cb ~ 7.1 / g
  Cb = 0.729;
  wvper = Cb .* windspd;
  wvhgt = Ca .* (windspd.^2) .* (1 + (0.015 .* windspd));

return;
