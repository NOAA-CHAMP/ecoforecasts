1;

 if ( ~exist('doPrint','var') )
   doPrint = false;
 end;

 s = 35;
 t = 25;
 %h = 10.96;
 dt = 3600;

 dTPerDay = 24*dt;
 dTPerWeek = 7*24*dt;

 R = (1-0.08);

 nq0s = 60;
 nbets = 40;

 % %q0s = linspace(-500,0,nq0s);				% dissertation figure
 % %q0s = linspace(-200,0,nq0s);
 % %q0s = linspace(-250,0,nq0s);				% proposal figure
 % %q0s = linspace(-250,+250,nq0s);			% smart proposal figure
 % %q0s = linspace(-300,+300,nq0s);			% smarter proposal figure
 % q0s = linspace(-300,+300,nq0s);			% smartest proposal figure

 q0s = linspace(-500,+500,nq0s);			% Paper figure??


 % %bets = [0.002,0.02,0.20];
 % %bets = 2 .* logspace(-3,-1,nbets);			% dissertation figure
 % %bets = 5 .* logspace(-3,-2,nbets);			% proposal figure
 % %bets = 5 .* logspace(-2.3,-1.3,nbets);		% smart proposal figure: 0.5 to 5%
 % %bets = logspace(-3,0,nbets);				% smartish proposal figure
 % %bets = 5 .* logspace(-3,-1.3,nbets);			% smarter proposal figure
 % bets = logspace(-3,-1.3,nbets);			% smartest proposal figure

 bets = logspace(log10(0.0075),log10(0.15),nbets);	% Paper figure??

 %for h=[3,5,10,15,20];
 %for h=[3,10,20];
 %for h=[3,15];
 for h=[5];
 %for h=[3];

  rho = sw_dens(s,t,h);					%[kg/m^3]
  Cp = sw_cp(s,t,h);					%[J/kg*K]
  rhoCp = rho.*Cp;					%[J/K*m^3]
  g = 9.79;						%[m/s^2]
  alph = sw_alpha(s,t,h);				%[1/K]

  clear dTdthc_SS dTdthc_US dTdthc_SU dTdthc_UU
  for qix=1:nq0s;

    % Net sea-surface heat flux, Q_0
    q0 = q0s(qix);

    for bix=1:nbets;
      % Seafloor slope, \beta
      bet = bets(bix);
      % Net sea-surface buoyancy flux B_0
      res.B0 = (g.*alph.*abs(q0))./(rhoCp);
      % "Characteristic" convective velocity (Sturmann)
      res.uf = (res.B0.*h).^(1/3);

      % Volumetric flow [m^2/s]: Panel letter refers to Fig. 10, Monismith et al (2006)
      % Advective inertial (steady) momentum balance, steady thermal balance: Panel (c)
      res.Qv_SS = res.uf.*h./(bet.^(1/3));
      % Viscous (unsteady) inertial balance, balanced thermal forcing: Panel (a)
      res.Qv_US = sqrt((res.uf.^3) .* (dTPerWeek) .* h);
      % Advective inertial balance, unbalanced thermal forcing: Panel (f)
      res.Qv_SU = (bet.^(2/3)).*res.uf.*((res.uf.*(dTPerWeek)./h).^(3/2));
      % Viscous (unsteady) inertia, unbalanced thermal forcing: Panel (d)
      res.Qv_UU = bet.*(res.uf.^3).*((dTPerWeek)^2)./h;

      % Convective flow rate [m/s]
      res.u_SS  = (5.0.*res.Qv_SS./h) - 0.05;   res.u_SS(res.u_SS<0) = 0;
      res.u_US  = (3.0.*res.Qv_US./h) - 0.026;  res.u_US(res.u_US<0) = 0;
      res.u_SU  = (2.7.*res.Qv_SU./h) - 0.0224; res.u_SU(res.u_SU<0) = 0;
      res.u_UU  = (0.1.*res.Qv_UU./h);          res.u_UU(res.u_UU<0) = 0;

      res.dTdtq0 = dt.*q0./(rhoCp.*h);

      % Furthest distance traveled by convective flow in an hour [m]
      res.dx_SS = res.u_SS .* dt;
      res.dx_US = res.u_US .* dt;
      res.dx_SU = res.u_SU .* dt;
      res.dx_UU = res.u_UU .* dt;

      % Temperature change due to Q0 at furthest extent of convection
      res.dTdtx_SS = dt.*q0./(rhoCp.*(h+(bet.*res.dx_SS)));
      res.dTdtx_US = dt.*q0./(rhoCp.*(h+(bet.*res.dx_US)));
      res.dTdtx_SU = dt.*q0./(rhoCp.*(h+(bet.*res.dx_SU)));
      res.dTdtx_UU = dt.*q0./(rhoCp.*(h+(bet.*res.dx_UU)));

      % Static temperature gradient due to Q0 over depth difference between
      % observation point, and point of further extent of convection
      res.dTdx_SS=(res.dTdtq0-res.dTdtx_SS)./res.dx_SS;
      res.dTdx_SS(~isfinite(res.dTdx_SS)) = 0;
      res.dTdx_US=(res.dTdtq0-res.dTdtx_US)./res.dx_US;
      res.dTdx_US(~isfinite(res.dTdx_US)) = 0;
      res.dTdx_SU=(res.dTdtq0-res.dTdtx_SU)./res.dx_SU;
      res.dTdx_SU(~isfinite(res.dTdx_SU)) = 0;
      res.dTdx_UU=(res.dTdtq0-res.dTdtx_UU)./res.dx_UU;
      res.dTdx_UU(~isfinite(res.dTdx_UU)) = 0;

      % Rayleigh Benard instability may dampen HC during warming (Mao Lei Patterson)
      if ( res.dTdtq0 > 0 )
        RB = 0.66;
      else
        RB = 1.00;
      end;

      % Temperature change due to horizontal convection at observation point
      dTdthc_SS(bix,qix) = -RB.*R.*dTPerWeek.*res.u_SS.*res.dTdx_SS;
      dTdthc_US(bix,qix) = -RB.*R.*dTPerWeek.*res.u_US.*res.dTdx_US;
      dTdthc_SU(bix,qix) = -RB.*R.*dTPerWeek.*res.u_SU.*res.dTdx_SU;
      dTdthc_UU(bix,qix) = -RB.*R.*dTPerWeek.*res.u_UU.*res.dTdx_UU;

      % Net temperature change (HC+Q0) at observation point
      dTdtq0(bix,qix) = res.dTdtq0.*24;
      dTdt_SS(bix,qix) = dTdthc_SS(bix,qix) + dTdtq0(bix,qix);
      dTdt_US(bix,qix) = dTdthc_US(bix,qix) + dTdtq0(bix,qix);
      dTdt_SU(bix,qix) = dTdthc_SU(bix,qix) + dTdtq0(bix,qix);
      dTdt_UU(bix,qix) = dTdthc_UU(bix,qix) + dTdtq0(bix,qix);
    end;
  end;

  fh=fmg;
  surf(q0s,bets,dTdt_SS);
  surf(q0s,bets,dTdt_US);
  surf(q0s,bets,dTdt_SU);
  surf(q0s,bets,dTdt_UU);
  sfh=mesh(q0s,bets,dTdtq0);
  %lh=plot3(q0s,bets([1,end]),dTdtq0([1,end],:)');
  set(gca,'yscale','log');
  % xlim([-500,0]);
  % ylim([0.002,0.2]);
  % zlim([0,1]);
  mindT = min([min(dTdt_SS(:)),min(dTdt_SU(:)),min(dTdt_US(:)),min(dTdt_UU(:)),min(dTdtq0(:))]);
  maxdT = max([max(dTdt_SS(:)),max(dTdt_SU(:)),max(dTdt_US(:)),max(dTdt_UU(:)),max(dTdtq0(:))]);
  xlim([min(q0s),max(q0s)]);
  ylim([min(bets),max(bets)]);
  zlim([mindT,maxdT]);
  % %%%view(-70,30);
  % %%%view(95,30);
  % %%view(175,22.5);
  % %view(-40,15);
  % view(-70,30);
  % view(-70,45);
  view(-35,45); % Paper figure??
  % %shading interp;
  % shading faceted;
  shading interp; % Paper figure??
  set(sfh,'FaceColor','none');

  % q0midix = round(nq0s/2);
  % bemidix = round(nbets/2);

  % text(q0s(1),bets(midix),dTdt_SS(midix,1),'SS');
  % text(q0s(1),bets(midix),dTdt_US(midix,1),'US');
  % text(q0s(1),bets(midix),dTdt_SU(midix,1),'SU');
  % text(q0s(1),bets(midix),dTdt_UU(midix,1),'UU');
  % %text(q0s(1),bets(midix),dTdtq0(midix,1),'Q0');

  for qix=[1,nq0s]
    for bix=[1,nbets]
      text(q0s(qix),bets(bix),dTdt_SS(bix,qix),'SS');
      text(q0s(qix),bets(bix),dTdt_US(bix,qix),'US');
      text(q0s(qix),bets(bix),dTdt_SU(bix,qix),'SU');
      text(q0s(qix),bets(bix),dTdt_UU(bix,qix),'UU');
      %text(q0s(qix),bets(bix),dTdtq0(bix,qix),'Q0');
    end;
  end;

  xlabel('Sea-surface flux Q_0 [Wm^-^2]');
  ylabel('Seafloor slope \beta');
  zlabel('Total temperature change \partial_tT [K\bulletweek^-^1]');

  zlim([mindT,maxdT]);
  colorbar;

  titlename(['\partial_tT sensitivity vs. \beta, Q_0 with horizontal convection, h=',num2str(h),' m']);

  if ( doPrint )
    print('-dpng',sprintf('%s-%gm.png',mfilename,h));
    %print('-dtiff',fullfile(get_thesis_path('../figs'),'phd-hc-vs-beta-vs-q0.tif'));
    %saveas(fh,fullfile(get_thesis_path('../figs'),'phd-hc-vs-beta-vs-q0.fig'));
  end;

 end;

 clear s t h dt R rho Cp rhoCp g alph qix bix q0 bet res ig zx ix jx fh ans
 %clear q0s bets
