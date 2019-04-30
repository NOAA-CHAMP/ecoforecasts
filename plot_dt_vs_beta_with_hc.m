1;

 if ( ~exist('doPrint','var') )
   % doPrint = true;
   doPrint = false;
 end;
 
 figspath = get_ecoforecasts_path('figs');
 
 fontSize = 20;
 
 s = 35;
 t = 25;
 %h = 10.96;
 dt = 3600;
 
 R = (1-0.08);
 
 %nq0s = 60;
 nq0s = 4;
 nbets = 40;
 
 % %q0s = linspace(-500,0,nq0s);				% dissertation figure
 % %q0s = linspace(-200,0,nq0s);
 % %q0s = linspace(-250,0,nq0s);				% proposal figure
 % %q0s = linspace(-250,+250,nq0s);			% smart proposal figure
 % %q0s = linspace(-300,+300,nq0s);			% smarter proposal figure
 % q0s = linspace(-300,+300,nq0s);			% smartest proposal figure
 
 % % %q0s = linspace(-300,+300,nq0s);
 % % q0s = [-800,400];
 % q0s = [-1400,800];
 %q0s = [-1400,800];
 q0s = [-800,400];
 nq0s = numel(q0s);
 
 % %bets = [0.002,0.02,0.20];
 % %bets = 2 .* logspace(-3,-1,nbets);			% dissertation figure
 % %bets = 5 .* logspace(-3,-2,nbets);			% proposal figure
 % %bets = 5 .* logspace(-2.3,-1.3,nbets);		% smart proposal figure: 0.5 to 5%
 % %bets = logspace(-3,0,nbets);				% smartish proposal figure
 % %bets = 5 .* logspace(-3,-1.3,nbets);			% smarter proposal figure
 % bets = logspace(-3,-1.3,nbets);			% smartest proposal figure
 
 % %bets = logspace(log10(0.0010),log10(0.0405),nbets);
 % bets = logspace(log10(0.01),log10(0.10),nbets);
 bets = logspace(log10(0.0010),log10(0.0500),nbets);
 nbets = numel(bets);
 
 hoursHC = 24;
 daysAccum = 1;
 
 % %for h=[3,5,10,15,20];
 % %for h=[3];
 % % for h=[3,10,20];
 % for h=[5,15];
 for h=[3];
   
   rho = sw_dens(s,t,h);				%[kg/m^3]
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
       res.Qv_US = sqrt((res.uf.^3) .* (hoursHC.*dt) .* h);
       % Advective inertial balance, unbalanced thermal forcing: Panel (f)
       res.Qv_SU = (bet.^(2/3)).*res.uf.*((res.uf.*(hoursHC.*dt)./h).^(3/2));
       % Viscous (unsteady) inertia, unbalanced thermal forcing: Panel (d)
       res.Qv_UU = bet.*(res.uf.^3).*((hoursHC.*dt)^2)./h;
       
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
       dTdthc_SS(bix,qix) = -RB.*R.*hoursHC.*dt.*res.u_SS.*res.dTdx_SS;
       dTdthc_US(bix,qix) = -RB.*R.*hoursHC.*dt.*res.u_US.*res.dTdx_US;
       dTdthc_SU(bix,qix) = -RB.*R.*hoursHC.*dt.*res.u_SU.*res.dTdx_SU;
       dTdthc_UU(bix,qix) = -RB.*R.*hoursHC.*dt.*res.u_UU.*res.dTdx_UU;
       
       % Net temperature change (HC+Q0) at observation point
       dTdtq0(bix,qix) = res.dTdtq0.*hoursHC;
       dTdt_SS(bix,qix) = dTdthc_SS(bix,qix) + dTdtq0(bix,qix);
       dTdt_US(bix,qix) = dTdthc_US(bix,qix) + dTdtq0(bix,qix);
       dTdt_SU(bix,qix) = dTdthc_SU(bix,qix) + dTdtq0(bix,qix);
       dTdt_UU(bix,qix) = dTdthc_UU(bix,qix) + dTdtq0(bix,qix);
     end; %for bix=1:nbets;
   end; %for qix=1:nq0s;
   
   % mindT = 1.1*min([min(dTdt_SS(:)),min(dTdt_SU(:)),min(dTdt_US(:)),min(dTdt_UU(:)),min(dTdtq0(:))]);
   % maxdT = 1.1*max([max(dTdt_SS(:)),max(dTdt_SU(:)),max(dTdt_US(:)),max(dTdt_UU(:)),max(dTdtq0(:))]);
   mindT = 1.1*min([min(dTdt_SS(:)),min(dTdt_SU(:))]);
   maxdT = 1.1*max([max(dTdt_SS(:)),max(dTdt_SU(:))]);
   
   fh=fmg;
   ax(1) = spt(2,1,1);
   ax(2) = spt(2,1,2);
   for qix=1:nq0s;
     % Net sea-surface heat flux, Q_0
     q0 = q0s(qix);
     if ( dTdtq0(1,qix) > 0 )
       axes(ax(1));
       plot(bets,daysAccum.*dTdt_SS(:,qix),'LineWidth',2,'Color','r');
       %plot(bets,daysAccum.*dTdt_SU(:,qix),'LineWidth',2);
     else
       axes(ax(2));
       plot(bets,daysAccum.*dTdt_SU(:,qix),'LineWidth',2,'Color','b');
     end;
     % plot(bets,dTdt_US);
     % plot(bets,dTdt_UU);
     % set(gca,'XScale','log');
     xlim([min(bets),max(bets)]);
     %ylim([mindT,maxdT]);
     
     % for qix=[1,nq0s]
     for qix=1:nq0s
       for bix=[1,nbets]
         if ( dTdtq0(bix,qix) > 0 )
           axes(ax(1));
           %text(bets(bix),dTdt_SS(bix,qix),[num2str(dTdt_SS(bix,qix),2),' Kd^-^1']);
           th=text(bets(bix),daysAccum.*dTdt_SS(bix,qix),[num2str(dTdt_SS(bix,qix),2),'\circC']);
         else
           axes(ax(2));
           th=text(bets(bix),daysAccum.*dTdt_SU(bix,qix),[num2str(dTdt_SU(bix,qix),2),'\circC']);
         end;
         set(th,'FontSize',fontSize);
       end;
     end;
     
     for axix=2:2;
       axes(ax(axix));
       ylim([daysAccum.*mindT,daysAccum.*maxdT]);
       xlabel('Seafloor slope (dz/dx) \beta');
       ylh=ylabel('Total temperature change \partial_tT [Kd^-^1]');
       %ylabel('Diurnal range (IQR) \circC');
       set(ylh,'HorizontalAl','left');
       % set(gca,'FontSize',fontSize);
       
       text(prctile(bets,25),0,'Reef flats','Rotation',90,'HorizontalAlign','center','FontSize',fontSize);
       text(prctile(bets,83),0,'Platform slope','Rotation',90,'HorizontalAlign','center','FontSize',fontSize);
       text(prctile(bets,93),0,'Offshore slope','Rotation',90,'HorizontalAlign','center','FontSize',fontSize);
     end;
     
     axes(ax(1));
     text(prctile(bets,70),daysAccum.*maxdT*0.70,['Diurnal warming (',num2str(q0s(2)),' W\bulletm^2)'],'Rotation',0,'HorizontalAlign','center','FontSize',fontSize);
     axes(ax(2));
     text(prctile(bets,70),daysAccum.*mindT*0.70,['Cold-front passage (',num2str(q0s(1)),' W\bulletm^2)'],'Rotation',0,'HorizontalAlign','center','FontSize',fontSize);
     
     axes(ax(1));
     xlabel([]);
     xticklabels([]);
     %titlename(['\partial_tT sensitivity to \beta for Q_0 + horizontal convection, h=',num2str(h),' m']);
     th=titlename(['Total temperature change vs. \beta, depth=',num2str(h),' m']);
     set(th,'FontSize',fontSize);
     
     axes(ax(1));
     xlim([0.0000,0.055]);
     % % %ylim([-3.50,+3.50]);
     % % %ylim([-3.00,+3.75]);
     % % %ylim([0.00,+4.00]);
     % ylim([1.0,2.5]); set(gca,'yscale','log');
     ylim([0.0,2.5]);
     
     axes(ax(2));
     xlim([0.0000,0.055]);
     % % %ylim([-3.50,+3.50]);
     % % %ylim([-3.00,+3.75]);
     % % ylim([-4.00,0.00]);
     % ylim([-4.0,-0.4]); set(gca,'yscale','log');
     ylim([-4.0,-0.0]);
     
     for axix=1:2;
       axes(ax(axix));
       set(gca,'FontSize',fontSize);
       grid on;
     end;
     axes(ax(2));
     if ( doPrint )
       print('-dpng',fullfile(figspath,sprintf('%s-%gm.png',mfilename,h)));
       saveas(fh,fullfile(figspath,sprintf('%s-%gm.fig',mfilename,h)));
     end;
     
   end; %for qix=1:nq0s;
 end; %for h=[5,15];
 
 clear s t h dt R rho Cp rhoCp g alph qix bix q0 bet res ig zx ix jx fh ans
 %clear q0s bets
