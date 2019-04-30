function sd = stokes_drift(w,hs,tp,estMethod,H)
%function sd = stokes_drift(w,hs,tp,estMethod,H)
%
% Estimate Stokes drift (residual wave-driven surface circulation [m/s]) from
% wind speed [kts] W, and significant wave height [m] HS. If provided, uses
% wave period [s] TP to estimate cutoff frequency for incident waves at site.
%
% ASSUMES wind vector with given magnitude W, is parallel to peak wave dir.
% (not given). If not, project W onto wave direction before calling this.
%
% If ESTMETHOD is 'shay', use simple wind-drift estimate of Shay et al 1998;
% 'monismith' == Monismith and Fong 2004 eqn 2; 'mellor' == use Mellor 2008
% eqn 10c; 'ardhuin' == Ardhuin et al (JPO) 2009 eqn 7. (DEFAULT: 'ardhuin'.)
%
% If ESTMETHOD is 'monismith' or 'mellor', optional H specifies an assumed
% water depth (DEFAULT 3 m) used in calculating wavenumber from period TP.
%
% Last Saved Time-stamp: <Wed 2016-11-16 17:17:20 Eastern Standard Time gramer>

  if ( ~exist('tp','var') || isempty(tp) )
    tp = 2;
  end;
  if ( ~exist('estMethod','var') || isempty(estMethod) )
    estMethod = 'ardhuin';
  end;
  if ( ~exist('H','var') || isempty(H) )
    % Assumed water depth [m]
    H = 3;
  end;

  % Convert winds from [kts] to [m/s]
  w = kts2mps(w);


  switch ( estMethod ),
   case 'shay',
    %%%
    %% Simple relationship from Shay et al 1998
    sd = 0.02.*w;


   case 'ardhuin',
    %%%
    %% Empirical relationship from Ardhuin et al. 2009

    % In Ardhuin, "f_c" is Bragg wave frequency for their HF radar?? Assume it
    % is instead a "reasonable upper bound" for wave frequency at this site.
    %%%% ??? NOTE the factor 0.5 below is COMPLETELY ARBITRARY! We multiple by
    % a factor as a convenience, assuming user will specify a "dominant" or
    % "average" wave period in the field TP (if they specify anything at all).
    % But this factor should probably be determined empirically in future...
    %fc = (2.*pi)./(tp.*0.5);   % How are WW3 peak and ERAI avg "period" calculated??
    fc = 1./(tp.*0.5);

    basew = abs(w);
    basew(basew > 14.5) = 14.5;

    % % APPLY A CUTOFF FOR LOCAL WIND EFFECT??
    % w(abs(w) < 1.8) = 0;

    % Equation from Ardhuin et al. 2009
    % $$U_1_0<14.5 ms^-^1: 5\times10^-^4(1.25 - (0.25\times((0.5/f_c)^1^.^3)))U_1_0^2 + 0.025(H_s-0.4)$$
    sd = (5e-4.*(1.25 - (0.25.*((0.5./fc).^(1.3)))).*w.*basew) + (0.025.*(hs-0.4));


   case 'monismith',
    %%%
    %% Theoretical Stokes formula from Monismith and Fong 2004
    a = hs./2;
    g = sw_g(25);
    %H = 3;  % Assumed water depth [m]
    sigma = 2.*pi./tp;
    % General back-solver for k given sigma:
    test_k = 0.1:0.1:10;
    test_sigma = sqrt( test_k .* g .* tanh(test_k.*H) );
    fm = fit(test_sigma',test_k','smoothingspline');
    k = fm(sigma);
    sd = (k.*sigma.*(a.^2)) .* cosh(2.*k.*H) ./ (2.*(sinh(k.*H).^2));


   case 'mellor',
    %%%
    %% Theoretical formula for monochromatic waves, e.g., Mellor 2008, eqn 10c
    % (2k_alphaE/c) * (cosh(2k(z+h))/sinh(2kD))
    % = (2k_alpha*ga^2/2c) * (cosh(2k(z+h))/sinh(2kD))
    % = (k_alpha*ga^2/(sigma/k)) * (cosh(2k(z+h))/sinh(2kD))
    % = (k^2*ga^2/sigma) * (cosh(2k(z+h))/sinh(2kD))
    a = hs./2; % Wave amplitude estimated from H_s
    %H = 3;  % Assumed water depth [m]
    sigma = (2.*pi) ./ tp;  % Wave frequency
    g = sw_g(25);
  
    % Wavenumber: k = (sigma.^2) ./ (g.*tanh(k.*H))
    % % Shallow water:
    % k = sqrt(g.*H) .* tp;
    % Intermediate depth back-solver:
    test_k = 0.1:0.1:10;
    test_sigma = sqrt( test_k .* g .* tanh(test_k.*H) );
    fm = fit(test_sigma',test_k','smoothingspline');
    k = fm(sigma);
    %% From Stokes_Drift.m (the original)
    %sd = (k.*g.*(a.^2)) .* (cosh(2.*k.*a) ./ sinh(2.*k.*H));
    sd = ((k^2).*g.*(a.^2)) .* (cosh(2.*k.*a) ./ sinh(2.*k.*H)) ./ sigma;


   otherwise,
    error('Unknown Stokes drift estimation method "%s"',estMethod);
  end;

  % Ensure we have no unphysical speeds
  sd(sd<0) = 0;

return;
