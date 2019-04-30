1;
% SCRIPT par_model
%
% Parameters for composited "two sine-curve" model of daily total insolation,
% at various Caribbean and Gulf of Mexico sites [van Woesik et al., 2006]
%
% Last Saved Time-stamp: <Tue 2008-08-19 16:46:38 Eastern Daylight Time gramer>
%

fgbl1 = [ 6.611     0.2084  0.2675  -0.2992 0.8668  0.2826  ];
cmrc3 = [ -2.209    0.356   2.741   7.802   0.01513 6.847   ];
mlrf1 = [ 6.443     0.1815  0.4659  0.221   1.519   1.737   ];
srvi2 = [ -0.3057   1.372  -0.2698  -6.495  0.1299  4.0     ];
lppr1 = [ 6.456     0.1298  0.8587  0.3093  1.378	-3.453  ];
lcci1 = [ 0.3586    1.374   -3.546  6.404   0.1371  0.8115  ];
dbjm1 = [ 6.476     0.1292  0.8537  0.2422  1.386   -3.475  ];
% Belize numbers are actually from Veracruz, Mexico!
grbz1 = [ 0.3362    0.8195  1.404   6.181   0.178   0.4923  ];

% van Woesik curves match monthly averages - smooth 'em out a bit
t = [1:365] * (12.0/365.0);

muMol_PER_W = 4.1513469579;

W_PER_kWh_PER_D = 41.666667;

muMol_PER_kWh_PER_D = (muMol_PER_W * W_PER_kWh_PER_D);

% fI = insolation_model(t, fgbl1, muMol_PER_kWh_PER_D);
% mI = insolation_model(t, mlrf1, muMol_PER_kWh_PER_D);
% bI = insolation_model(t, grbz1, muMol_PER_kWh_PER_D);
fI = insolation_model(t, fgbl1, muMol_PER_kWh_PER_D) * 1.5;
mI = insolation_model(t, mlrf1, muMol_PER_kWh_PER_D) * 1.5;
bI = insolation_model(t, grbz1, muMol_PER_kWh_PER_D) * 1.5;

figure;
plot(t, [ fI ; mI ; bI ]);
xlim([t(1) t(end)]);
legend('FGBL1','MLRF1','GRBZ1', 'Location', 'East');
xlabel('Month');
%ylabel('kW m^-^2 day^-^1');
ylabel('Avg. \mu-mole quanta m^-^2 s^-^1');
title('van Woesik et al. Insolation Model');
print('-djpeg', 'van-Woesik-curves.jpg');


fFit = fit(t', fI', 'pchipinterp');
mFit = fit(t', mI', 'pchipinterp');
bFit = fit(t', bI', 'pchipinterp');

fDer1 = differentiate(fFit, t)';
mDer1 = differentiate(mFit, t)';
bDer1 = differentiate(bFit, t)';

figure;
plot(t, [ fDer1 ; mDer1 ; bDer1 ]);
xlim([t(1) t(end)]);
legend('dFGBL1','dMLRF1','dGRBZ1', 'Location', 'East');
xlabel('Month');
ylabel('\Delta \mu-mole quanta m^-^2 s^-^1 month^-^1');
title('van Woesik et al. Insolation Curve Derivatives');
print('-djpeg', 'van-Woesik-derivatives.jpg');


fInt1 = integrate(fFit, t, 0)';
mInt1 = integrate(mFit, t, 0)';
bInt1 = integrate(bFit, t, 0)';

figure;
plot(t, [ fInt1 ; mInt1 ; bInt1 ]);
xlim([t(1) t(end)]);
legend('\Sigma FGBL1','\Sigma MLRF1','\Sigma GRBZ1', 'Location', 'East');
xlabel('Month');
ylabel('\mu-mole quanta m^-^2');
title('van Woesik et al. Insolation Curve Integrals');
print('-djpeg', 'van-Woesik-integrals.jpg');


return;
