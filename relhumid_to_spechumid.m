function spechumid = relhumid_to_spechumid(airt,relhumid)
%function spechumid = relhumid_to_spechumid(airt,relhumid)
%
% From Lenntech "Relative humidity calculator"
%  http://www.lenntech.com/calculators/humidity/relative-humidity.htm
%  http://www.lenntech.com/calculators/relative-humidity.js
% /*made by Vijay Bhagwandas, Lenntech BV*/
%
% Last Saved Time-stamp: <Thu 2010-02-18 20:33:53 Eastern Standard Time gramer>

  % absMoisture1 [g/kg]
  absMoisture1 = (0.42 .* relhumid .* exp(airt * 10 * 0.006235398) ./ 10);

  spechumid = absMoisture1 * 1e-3;


%%%% MORE PRECISE FORMULAE - to be implemented later?

% Another approximation (Magnus' formula) would be
%
%         log10(es) = -2937.4/T - 4.9283*log10(T) + 23.5470
%
% Read more: http://www.faqs.org/faqs/meteorology/temp-dewpoint/#ixzz0fwPAQsk6

%  2.  Saturated vapor pressure (es):
%
%      es = 6.1078 * 10 ** ((T * A)/(T + B))  [hPa]
%                             \
%                              temperature in C
%                    A =   7.5 } for use in vapor pressure
%                    B = 237.3 } with respect to WATER
%
% Read more: http://www.faqs.org/faqs/meteorology/temp-dewpoint/#ixzz0fwPMnJQ3


% 1)  compute e as [es(T)*rH/100]
%     where es(T) = 0.611*EXP(17.27*T/(T+237.3)) in kPa
%     T is drybulb temp in C
%
%     e = (rH/100)* 0.611*EXP(17.27*T/(T+237.3))
%     where e is ambient vapor pressure in kPa
%
%
% Read more: http://www.faqs.org/faqs/meteorology/temp-dewpoint/#ixzz0fwOEjm9r


% Calculate the saturation pressure (es) from one of the formulas given
% above.
%
% Then multiply by the relative humidity (rh). This gives you the ambient
% water vapour pressure, (e).
%
% Then the specific humidity is given by the following formula:
%
%         R
%          L        e
% rho =  --- -----------------
%         R   p + e(R / R - 1)
%          W         L   W
%
% WHERE:
%         R / R  = 0.62197   (see the example for the mixing ratio)
%          L   W
%
% Read more: http://www.faqs.org/faqs/meteorology/temp-dewpoint/#ixzz0fwQL02AW
%
% Read more: http://www.faqs.org/faqs/meteorology/temp-dewpoint/#ixzz0fwNlm8AL



return;
