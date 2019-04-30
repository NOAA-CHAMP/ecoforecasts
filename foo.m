function foo(arg1)
  disp(['1: ',inputname(1)]);
  disp(['2: ',inputname(2)]);
return;

function foo_was_earlier_foo(varargin)
  disp(numel(varargin));
return;

function stn = foo_was_also_foo(stn)
      z = 3.55;
      %% Ohlmann [2003] parameterization of Kd from chlorophyll a [mg/m^3]
      %% See Sweeney et al. [2005] for field comparison results
      stn = verify_variable(stn,'amodis_chlor_a_1_hour_spline');
      chl = stn.amodis_chlor_a_1_hour_spline.data/10;
      zeta1=0.015+(0.176*((0.462*chl).^0.5));
      zeta2=0.688+(0.060.*log(0.125*chl));
      V1=0.571+(0.025.*log(0.149*chl));
      V2=0.22374007 + (0.01609314.*log(2.33384554*chl));
      stn.amodis_chlor_a_Kd.date=stn.amodis_chlor_a_1_hour_spline.date;
      stn.amodis_chlor_a_Kd.data=V1.*exp(-z./zeta1) + V2.*exp(-z./zeta2);
return;



function foo_formerly_known_as_foo

cd ..

return;



w = [-30:2:30]; hs=[0.2 0.4 0.5:0.5:3]; tps = [2 4 8];

for tpix=1:length(tps)

  for hsix = 1:length(hs)
    sh(hsix,1:length(w)) = stokes_drift(w,hs(hsix),tps(tpix));
  end;

  maxigraph(figure); plot(w,sh); title(sprintf('Stokes vs. W, Hs=%g-%g [m], Tp=%g [s]',hs(1),hs(end),tps(tpix))); legend(num2str(hs')); xlabel('U_1_0 [kts]'); ylabel('Stokes [m/s]');

end;


return;

% k = (sigma.^2) ./ (g.*tanh(k.*D));

D = 3;
g = grv(25);

k = 0.1:0.1:10;
sigma = sqrt( k .* g .* tanh(k.*D) );

fm = fit(sigma',k','smoothingspline');
fitted_k = fm(sigma);

figure;
hold on;
plot(sigma, k, 'b-');
plot(sigma, fitted_k', 'r:');
xlabel('\sigma'); ylabel('k');
legend('original', 'fitted');
