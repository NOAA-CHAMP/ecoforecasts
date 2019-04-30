function stn = model_surf_par(stanam, stn, begyr, endyr)
%function stn = model_surf_par(stanam, stn, begyr, endyr)
%
% Reproduce the composited "two sine-curve" model of daily total insolation,
% at various Caribbean and Gulf of Mexico sites [van Woesik et al., 2006]
%
% Last Saved Time-stamp: <Tue 2008-08-19 18:52:59 Eastern Daylight Time gramer>
%

  stn.model_surf_par.date = [];
  stn.model_surf_par.data = [];

  modelparms = [];
  switch (stanam)
   case 'fgbl1',
    modelparms = [ 6.611     0.2084  0.2675  -0.2992 0.8668  0.2826  ];
   case 'cmrc3',
    modelparms = [ -2.209    0.356   2.741   7.802   0.01513 6.847   ];
   case {'fwyf1','mlrf1','lonf1','smkf1','sanf1','plsf1','dryf1'},
    modelparms = [ 6.443     0.1815  0.4659  0.221   1.519   1.737   ];
   case 'srvi2',
    modelparms = [ -0.3057   1.372  -0.2698  -6.495  0.1299  4.0     ];
   case 'lppr1',
    modelparms = [ 6.456     0.1298  0.8587  0.3093  1.378	-3.453  ];
   case 'lcci1',
    modelparms = [ 0.3586    1.374   -3.546  6.404   0.1371  0.8115  ];
   case 'dbjm1',
    modelparms = [ 6.476     0.1292  0.8537  0.2422  1.386   -3.475  ];
   % Belize numbers are actually from Veracruz, Mexico!
   case 'grbz1',
    modelparms = [ 0.3362    0.8195  1.404   6.181   0.178   0.4923  ];
   otherwise,
    error('Unrecognized station name "%s"!', stanam);
  end;

  % van Woesik curves are for monthly averages - smooth them out a bit
  t = [0:(1.0/24.0):364.99] * (12.0/365.0);

  muMol_PER_W = 4.1513469579;
  W_PER_kWh_PER_D = 41.666667;
  muMol_PER_kWh_PER_D = (muMol_PER_W * W_PER_kWh_PER_D);

  % Holy "empirical fudge-factors", Batman!
  model = insolation_model(t, modelparms, muMol_PER_kWh_PER_D) * 1.54;

  for yr = begyr:endyr

    begdt = datenum(yr,1,1);
    enddt = datenum(yr,12,31,23,59,59);
    dts = begdt:(1.0/24.0):enddt;
    ndays = floor(enddt) - floor(begdt) +  1;

    clear modelyr;
    modelyr = model;
    if ( ndays == 366 )
      modelyr(end+1:end+24) = model(end-23:end);
    end;

    % We have to calculate an HOURLY time series from a MONTHLY model X^*
    % Figure out local mid-day
    noon = 12; noondeg = (noon/24.0)*360;
    daily_model = max(0, cosd(linspace(-noondeg, (360-noondeg), 25)));
    daily_model = daily_model(1:24);
    daily_model = repmat(daily_model, [1 round(length(dts)/24)]);
    val = modelyr .* daily_model;

    stn.model_surf_par.date = [stn.model_surf_par.date dts];
    stn.model_surf_par.data = [stn.model_surf_par.data val];

  end;

return;
