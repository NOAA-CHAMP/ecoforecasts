function y=cor30a(x)
%version with shortened iteration; modified Rt and Rq
%uses wave information wave period in s and wave ht in m
%no wave, standard coare 2.6 charnock:  jwave=0 
%Oost et al.  zo=50/2/pi L (u*/c)^4.5 if jwave=1
%taylor and yelland  zo=1200 h*(L/h)^4.5 jwave=2
%x=[5.5 0 28.7 27.2 24.2 18.5 141 419 0 600 1010 15 15 15 0 1 1 5 1 ];%sample data stream
u=x(1);%wind speed (m/s)  at height zu (m)
us=x(2);%surface current speed in the wind direction (m/s)
ts=x(3);%bulk water temperature (C) if jcool=1, interface water T if jcool=0  
t=x(4);%bulk air temperature (C), height zt
Qs=x(5)/1000;%bulk water spec hum (g/kg) if jcool=1, ...
Q=x(6)/1000;%bulk air spec hum (g/kg), height zq
Rs=x(7);%downward solar flux (W/m^2)
Rl=x(8);%downard IR flux (W/m^2)
rain=x(9);%rain rate (mm/hr)
zi=x(10);%PBL depth (m)
P=x(11);%Atmos surface pressure (mb)
zu=x(12);%wind speed measurement height (m)
zt=x(13);%air T measurement height (m)
zq=x(14);%air q measurement height (m)
lat=x(15);%latitude (deg, N=+)
jcool=x(16);%implement cool calculation skin switch, 0=no, 1=yes
jwave=x(17);%implement wave dependent roughness model
twave=x(18);%wave period (s)
hwave=x(19);%wave height (m)

     %***********   set constants *************
     Beta=1.2;
     von=0.4;
     fdg=1.00;
     tdk=273.16;
     grav=grv(lat);%9.82;
     %*************  air constants ************
     Rgas=287.1;
     Le=(2.501-.00237*ts)*1e6;
     cpa=1004.67;
     cpv=cpa*(1+0.84*Q);
     rhoa=P*100/(Rgas*(t+tdk)*(1+0.61*Q));
     visa=1.326e-5*(1+6.542e-3*t+8.301e-6*t*t-4.84e-9*t*t*t);
     %************  cool skin constants  *******
     Al=2.1e-5*(ts+3.2)^0.79;
     be=0.026;
     cpw=4000;
     rhow=1022;
     visw=1e-6;
     tcw=0.6;
     bigc=16*grav*cpw*(rhow*visw)^3/(tcw*tcw*rhoa*rhoa);
     wetc=0.622*Le*Qs/(Rgas*(ts+tdk)^2);
     
     %***************   wave parameters  *********
     lwave=grav/2/pi*twave^2;
     cwave=grav/2/pi*twave;
     
     %**************  compute aux stuff *******
     Rns=Rs*.945;
     Rnl=0.97*(5.67e-8*(ts-0.3*jcool+tdk)^4-Rl);
     
     %***************   Begin bulk loop *******
     
     %***************  first guess ************
     du=u-us;
     dt=ts-t-.0098*zt;
     dq=Qs-Q;
     ta=t+tdk;
     ug=.5;
     dter=0.3; 
     dqer=wetc*dter;
     ut=sqrt(du*du+ug*ug);
	  u10=ut*log(10/1e-4)/log(zu/1e-4);
     usr=.035*u10;
	zo10=0.011*usr*usr/grav+0.11*visa/usr;
	Cd10=(von/log(10/zo10))^2;
	Ch10=0.00115;
	Ct10=Ch10/sqrt(Cd10);
	zot10=10/exp(von/Ct10);
	Cd=(von/log(zu/zo10))^2;
	Ct=von/log(zt/zot10);
	CC=von*Ct/Cd;
	Ribcu=-zu/zi/.004/Beta^3;
	Ribu=-grav*zu/ta*((dt-dter*jcool)+.61*ta*dq)/ut^2;
	nits=3;
	if Ribu<0;
		zetu=CC*Ribu/(1+Ribu/Ribcu);
	else;
		zetu=CC*Ribu*(1+27/9*Ribu/CC);
		end;		
	L10=zu/zetu;
	if zetu>50;
		nits=1;
	end;
     usr=ut*von/(log(zu/zo10)-psiu_30(zu/L10));
     tsr=-(dt-dter*jcool)*von*fdg/(log(zt/zot10)-psit_30(zt/L10));
     qsr=-(dq-wetc*dter*jcool)*von*fdg/(log(zq/zot10)-psit_30(zq/L10));

     tkt=.001;
	
   charn=0.011;
   if ut>10
      charn=0.011+(ut-10)/(18-10)*(0.018-0.011);
   end;
   if ut>18
      charn=0.018;
   end;
   
     %disp(usr)
     
     %***************  bulk loop ************
  for i=1:nits;
     
     zet=von*grav*zu/ta*(tsr*(1+0.61*Q)+.61*ta*qsr)/(usr*usr)/(1+0.61*Q);
      %disp(usr)
      %disp(zet);
      if jwave==0;zo=charn*usr*usr/grav+0.11*visa/usr;end;
      if jwave==1;zo=50/2/pi*lwave*(usr/cwave)^4.5+0.11*visa/usr;end;%Oost et al
      if jwave==2;zo=1200*hwave*(hwave/lwave)^4.5+0.11*visa/usr;end;%Taylor and Yelland
      rr=zo*usr/visa;
     L=zu/zet;
     zoq=min(1.15e-4,5.5e-5/rr^.6);
     zot=zoq;
     usr=ut*von/(log(zu/zo)-psiu_30(zu/L));
     tsr=-(dt-dter*jcool)*von*fdg/(log(zt/zot)-psit_30(zt/L));
     qsr=-(dq-wetc*dter*jcool)*von*fdg/(log(zq/zoq)-psit_30(zq/L));
     Bf=-grav/ta*usr*(tsr+.61*ta*qsr);
     if Bf>0
     ug=Beta*(Bf*zi)^.333;
     else
     ug=.2;
     end;
     ut=sqrt(du*du+ug*ug);
     Rnl=0.97*(5.67e-8*(ts-dter*jcool+tdk)^4-Rl);
     hsb=-rhoa*cpa*usr*tsr;
     hlb=-rhoa*Le*usr*qsr;
     qout=Rnl+hsb+hlb;
     dels=Rns*(.065+11*tkt-6.6e-5/tkt*(1-exp(-tkt/8.0e-4))); 	% Eq.16 Shortwave
     qcol=qout-dels;
     alq=Al*qcol+be*hlb*cpw/Le;					% Eq. 7 Buoy flux water

     if alq>0;
     		xlamx=6/(1+(bigc*alq/usr^4)^.75)^.333;				% Eq 13 Saunders
            tkt=xlamx*visw/(sqrt(rhoa/rhow)*usr);			%Eq.11 Sub. thk

     else
            xlamx=6.0;
            tkt=min(.01,xlamx*visw/(sqrt(rhoa/rhow)*usr));			%Eq.11 Sub. thk
     end;
     
      dter=qcol*tkt/tcw;%  Eq.12 Cool skin
      dqer=wetc*dter;
     
  end;%bulk iter loop
     tau=rhoa*usr*usr*du/ut;                %stress
     hsb=-rhoa*cpa*usr*tsr;
     hlb=-rhoa*Le*usr*qsr;

     
     %****************   rain heat flux ********
     
      dwat=2.11e-5*((t+tdk)/tdk)^1.94;%! water vapour diffusivity
      dtmp=(1.+3.309e-3*t-1.44e-6*t*t)*0.02411/(rhoa*cpa); 	%!heat diffusivity
      alfac= 1/(1+(wetc*Le*dwat)/(cpa*dtmp));      	%! wet bulb factor
      RF= rain*alfac*cpw*((ts-t-dter*jcool)+(Qs-Q-dqer*jcool)*Le/cpa)/3600;
     %****************   Webb et al. correection  ************
     wbar=1.61*hlb/Le/(1+1.61*Q)/rhoa+hsb/rhoa/cpa/ta;%formulation in hlb already includes webb
     %wbar=1.61*hlb/Le/rhoa+(1+1.61*Q)*hsb/rhoa/cpa/ta;
     hl_webb=rhoa*wbar*Q*Le;
     %**************   compute transfer coeffs relative to ut @meas. ht **********
     Cd=tau/rhoa/ut/max(.1,du);
     Ch=-usr*tsr/ut/(dt-dter*jcool);
     Ce=-usr*qsr/(dq-dqer*jcool)/ut;
     %************  10-m neutral coeff realtive to ut ********
     Cdn_10=von*von/log(10/zo)/log(10/zo);
     Chn_10=von*von*fdg/log(10/zo)/log(10/zot);
     Cen_10=von*von*fdg/log(10/zo)/log(10/zoq);
   
    
   y=[hsb hlb tau zo zot zoq L usr tsr qsr dter dqer tkt RF wbar Cd Ch Ce Cdn_10 Chn_10 Cen_10 ug ];
   %   1   2   3   4  5   6  7  8   9  10   11   12  13  14  15  16 17 18    19      20    21  22
