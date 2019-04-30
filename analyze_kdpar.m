function analyze_kdpar(station)

  global kddts kd udts u vdts v spd3dts spd3;
  global rawkd rawkddts spddts spd wvardts wvar;

  if ( ~exist('station', 'var') )
    % station = 'CMRC3';
    station = 'SRVI2';
  end

%   [kddts, kd] = load_g2_data([station '-KdPAR-sfc.csv']);
%   % [kddts, kd] = load_g2_data([station '-KdPAR-deep.csv']);

%   [udts, u] = load_g2_data([station '-WIND1-U.csv']);

%   [vdts, v] = load_g2_data([station '-WIND1-V.csv']);

%   [spd3dts, spd3] = load_g2_data([station '-WIND1-SPEED-3-day.csv']);

%   [wvardts,wvar] = load_g2_data([station '-wind-3-day-var.csv']);

%   spddts = intersect(udts,vdts);
%   spd = sqrt((u(ismember(udts,spddts)).*u(ismember(udts,spddts))) ...
%              + (v(ismember(vdts,spddts)).*v(ismember(vdts,spddts))));


%   rawkddts = kddts;
%   rawkd = kd;

  % Ignore all unreasonable Kd values
%   idx = find(rawkd >= 0.0);
%   kddts = rawkddts(idx);
%   kd = rawkd(idx);

  % Ignore all uninteresting Kd values
  idx = find(0.02 < rawkd & rawkd < 1.0);
  kddts = rawkddts(idx);
  kd = rawkd(idx);


%   % Consider mid-day values only!
%   kddvec = datevec(rawkddts);
%   kddtsdaily = rawkddts(kddvec(:,4) == 18);
%   kddaily = rawkd(ismember(rawkddts,kddtsdaily));
%   kddts = kddtsdaily;
%   kd = kddaily;


  regress_kd(wvardts, wvar, 'Var_3_d_a_y^w^i^n^d');
  regress_kd(udts, u, 'U^w^i^n^d');
  regress_kd(vdts, v, 'V^w^i^n^d');
  regress_kd(spddts, spd, 'Speed^w^i^n^d');
  regress_kd(spd3dts, spd3, 'Speed_3_d_a_y^w^i^n^d');

return;


function regress_kd(ydts, ydat, yttl)

  global kddts kd;

  dts = intersect(kddts,ydts);
  d1 = kd(ismember(kddts,dts)); d2 = ydat(ismember(ydts,dts));
  X = d1; X(:,2) = ones(size(X)); Y = d2;
  [B, BINT, R, RINT, RSTATS] = regress(Y, X);

  figure; hold on;
  plot(d1, d2, '.'); plot(d1, ((d1*B(1)) + B(2)), 'r');
  legend( ['K_d^P^A^R vs. ' yttl], ...
          sprintf('%g * K_d + %g', B(1), B(2)) );
  ttlclr = 'black';
  if ( RSTATS(3) < 0.08 ); ttlclr = [0.2 0.5 0.2]; end;
  title( sprintf('Regression: R^2=%g p=%g err=%g resid=[%g %g]', ...
                 RSTATS(1), RSTATS(3), RSTATS(4), RINT(1,1), RINT(1,2)), ...
         'Color', ttlclr);
  hold off;

return;
