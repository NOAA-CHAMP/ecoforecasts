function out=cor30a_warm(dts,U,vus,msp,ta,vqs,qa,rs,rl,org,zi,vP,lon,lat,zu,zt,zq,ts_depth,jwarm,jcool,jwave,vtwave,vhwave,max_pwp)
%function out=cor30a_warm(dts,U,vus,msp,ta,vqs,qa,rs,rl,org,zi,vP,lon,lat,zu,zt,zq,ts_depth,jwarm,jcool,jwave,vtwave,vhwave,max_pwp)
%
% Lew.Gramer@noaa.gov, 12 Mar 2012: reformated w/tab-stop=3, untabified
% Lew.Gramer@noaa.gov, 12 Mar 2011: allow caller to limit PWP warm-layer depth
%
% Aside from a few tweaks to embed it within a proper MATLAB function that
% takes arguments and returns a value, the stringy bulk of the following code
% is from Fairall et al. script COR3_0AH.M. [Lew.Gramer@noaa.gov, 2010-04-09]
%
% Last Saved Time-stamp: <Mon 2012-03-12 13:45:53  Lew.Gramer>

  if ( ~exist('vtwave','var') || ~exist('vhwave','var') || ...
       isempty(vtwave) || isempty(vhwave) )
    warning('Estimating wave period and height from wind speed');
    [vtwave,vhwave] = wind_to_wave(U);
  end;

  if ( ~exist('max_pwp','var') || isempty(max_pwp) )
    % Default in original COARE 3.0a code was 19m
    max_pwp = 19;
  end;

% Additional comment by Gramer - these are args expected by m-function COR30A:
%
% %version with shortened iteration; modified Rt and Rq
% %uses wave information wave period in s and wave ht in m
% %no wave, standard coare 2.6 charnock:  jwave=0 
% %Oost et al.  zo=50/2/pi L (u*/c)^4.5 if jwave=1
% %taylor and yelland  zo=1200 h*(L/h)^4.5 jwave=2
% %x=[5.5 0 28.7 27.2 24.2 18.5 141 419 0 600 1010 15 15 15 0 1 1 5 1 ];%sample data stream
% u=x(1);%wind speed (m/s)  at height zu (m)
% us=x(2);%surface current speed in the wind direction (m/s)
% ts=x(3);%bulk water temperature (C) if jcool=1, interface water T if jcool=0  
% t=x(4);%bulk air temperature (C), height zt
% Qs=x(5)/1000;%bulk water spec hum (g/kg) if jcool=1, ...
% Q=x(6)/1000;%bulk air spec hum (g/kg), height zq
% Rs=x(7);%downward solar flux (W/m^2)
% Rl=x(8);%downard IR flux (W/m^2)
% rain=x(9);%rain rate (mm/hr)
% zi=x(10);%PBL depth (m)
% P=x(11);%Atmos surface pressure (mb)
% zu=x(12);%wind speed measurement height (m)
% zt=x(13);%air T measurement height (m)
% zq=x(14);%air q measurement height (m)
% lat=x(15);%latitude (deg, N=+)
% jcool=x(16);%implement cool calculation skin switch, 0=no, 1=yes
% jwave=x(17);%implement wave dependent roughness model
% twave=x(18);%wave period (s)
% hwave=x(19);%wave height (m)

%toga coare bulk flux model version 3.0a
%***************************************
%uses following matlab subroutines:
%  cor30a.m
%  psiu_30.m
%  psit_30.m
%  qsee.m
%  grv.m
%***************************************

%*********** basic specifications  *****
%  zu=         height of wind measurement
%  zt=         height of air temperature measurement
%  zq=         height of air humidity measurement
%  ts_depth    depth of water temperature measurement
%  jwarm=         0=no warm layer calc, 1 =do warm layer calc
%  jcool=         0=no cool skin calc, 1=do cool skin calc
%   jwave=      0= Charnock, 1=Oost et al, 2=Taylor and Yelland

%***********   input data **************
%  YYYYMMHHMMSS=     date in toga coare format, Y2K version
%  u=       wind speed (m/s), height zu
%  us=         surface current (m/s)
%  ts=         bulk surface sea temp (cent)
%  t=       air temp (cent), height zt
%  qs=         sea surface sat specific humidity (g/kg)
%  q=       air specific humidity (g/kg), height zq
%  Rs=         downward solar flux (w/m^2)
%  Rl=         downward IR flux (w/m^2)
%  zi=         inversion height (m)
%  P=       air pressure (mb)
%  rain=    rain rate (mm/hr)
%  lon=     longitude (deg E=+)
%  lat=     latitude (deg N=+)


%********** output data  ***************
%  hsb=        sensible heat flux (w/m^2)
%  hlb=        latent heat flux (w/m^2)
%  RF=         rain heat flux(w/m^2)
%  wbar=       webb mean w (m/s)
%  tau=        stress (nt/m^2)
%  zo=         velocity roughness length (m)
%  zot         temperature roughness length (m)
%  zoq=        moisture roughness length (m)
%  L=       Monin_Obukhov stability length
%  usr=        turbulent friction velocity (m/s), including gustiness
%  tsr         temperature scaling parameter (K)
%  qsr         humidity scaling parameter (g/g)
%  dter=       cool skin temperature depression (K)
%  dqer=       cool skin humidity depression (g/g)
%  tkt=        cool skin thickness (m)
%  Cd=         velocity drag coefficient at zu, referenced to u
%  Ch=         heat transfer coefficient at zt
%  Ce=         moisture transfer coefficient at zq
%  Cdn_10=        10-m velocity drag coeeficient, including gustiness
%  Chn_10=        10-m heat transfer coeeficient, including gustiness
%  Cen_10=        10-m humidity transfer coeeficient, including gustiness
%


%    %x=load('c:\My Documents\matlabstf\bulkalg\test2_5b.txt');%read file with half-hour average data; set your local directory 
%     x=load('c:\data\cwf\matlabstf\cwf\bulkalg\cor3_0\matlab\test3_0.txt');%read file with 1hr-average data; set your local directory 

%    jdy=x(:,1);%time in the form YYYYMMDDHHSS.SS
%    U=x(:,2); %true wind speed, m/s; etl sonic anemometer
%    tsnk=x(:,3);%sea snake temperature, C (0.05 m depth)
%    ta=x(:,4);%air temperature, C (z=14.5 m)
%    qa=x(:,5);%air specific humidity, g/kg (z=14.5  m)
%    rs=x(:,6);%downward solar flux, W/m^2 (ETL units)
%    rl=x(:,7);%downward IR flux, W/m^2 (ETL units)
%    org=x(:,8);%rainrate, mm/hr (ETL STI optical rain gauge, uncorrected)
%    lat=x(:,9);%latitude, deg  (SCS pcode)
%    lon=x(:,10);%longitude, deg (SCS pcode)
%  msp=x(:,11);%6-m deotg T from MSP, C    

% zu=15;%anemometer ht
% zt=15;%air T height
% zq=15;%humidity height
% ts_depth=6;%bulk water temperature sensor depth, U/APL MSP&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

% jcool=1;
% jwarm=1;
% jwave=0;

icount=1;
%*********************  housekeep variables  ******** 
qcol_ac=0;
tau_ac=0;
jtime=0;
jamset=0;
tau_old=.06;
hs_old=10;
hl_old=100;
RF_old=0;
dsea=0;
dt_wrm=0;
tk_pwp=max_pwp;
fxp=.5;
q_pwp=0;
jump=1;
%*******************  set constants  ****************
   tdk=273.16;
   grav=grv(lat);
   Rgas=287.1;
   cpa=1004.67;

   be=0.026;
   cpw=4000;
   rhow=1022;
% Not used?
   % visw=1e-6;
   tcw=0.6;
    dter=0.3;
    ts=msp(1);%**************************&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

%******************  setup read data loop  **********
% [nn,nx]=size(x');  %# of lines of data  
nx = numel(U);
nby10 = floor(nx/10);

disp([num2str(nx) ' data points']);

jdy = str2num(datestr(dts,'yyyymmddHHMMSS.FFF'));

%%%% Try preallocating all arrays to speed this loop up (Lew.Gramer@noaa.gov)
%DEBUG:tic,
jdx = repmat(nan,[nx 1]);
locx = repmat(nan,[nx 1]);
rnl = repmat(nan,[nx 1]);
qsx = repmat(nan,[nx 1]);
tsx = repmat(nan,[nx 1]);
Hrain = repmat(nan,[nx 1]);
%**********  new values from this code
hnet = repmat(nan,[nx 1]);
hs = repmat(nan,[nx 1]);
hl = repmat(nan,[nx 1]);
tau = repmat(nan,[nx 1]);
hl_webb = repmat(nan,[nx 1]);
%**********  outputs of the function
dt = repmat(nan,[nx 35]);
out = repmat(nan,[nx 25]);


%d1=datenum(1998,12,31);
iyr0=fix(jdy(1)/1e10);%get year for first line of data
for ibg = 1:nx          %major read loop

   % Lew.Gramer@noaa.gov: Let user know how far along we are in this looong loop
   if (mod(ibg,nby10)==0); disp(ibg); end;

   %***********   set variables not in data base  ********
   %P=1008;                      %air pressure
   P=vP(ibg);                   %air pressure
   %us=0;                        %surface current
   us=vus(ibg);                 %surface current
   %zi=600;                      %inversion ht
   % Lew.Gramer@noaa.gov: Accept either scalar, or actual PBL height (e.g., from ERA-Interim)
   if (isscalar(zi)); pblz=zi; else; pblz=zi(ibg); end;
      %*******   decode date  ************


   %[iyr mon iday ihr imin isec]=datevec(jdy(ibg)+d1);      
   st=jdy(ibg);
   iyr=fix(st/1e10);
   mon=fix(st/1e8)-iyr*100;
   iday=fix(st/1e6)-iyr*1e4-mon*100;
   ihr=fix(st/1e4)-iyr*1e6-mon*1e4-iday*100;
   imin=fix(st/100)-fix(st/1e4)*100;
   isec=0;
   %timex=[iyr mon iday ihr imin isec];
   jd=datenum(iyr,mon,iday,ihr,imin,isec)-datenum(iyr0-1,12,31);%year day number, Jan 1=1
   jdx(ibg)=jd;
      %********   decode bulk met data ****
   u=U(ibg);%wind speed
   tsea=msp(ibg);%bulk sea surface temp********************&&&&&&&&&&&&&&&&&&&&&&&&&&&
   t=ta(ibg);%air temp
%  qs=qsee([tsea P]);%bulk sea surface humidity
   qs=vqs(ibg);%bulk sea surface humidity
   q=qa(ibg);%air humidity
   Rs=rs(ibg);%downward solar flux
   Rl=rl(ibg);%doward IR flux
   rain=org(ibg);%rain rate
   % Lew.Gramer@noaa.gov: C-MAN/SEAKEYS and CREWS/ICON stations don't move!
   %grav=grv(lat(ibg));%9.72;
   %lonx=lon(ibg);%longitude
   grav=grv(lat);
   lonx=lon;%longitude

   %*****  variables for warm layer  ***
   %ntime=1e6*mon+1e4*iday+100*ihr+imin;
   time=(ihr*3600+imin*60)/24/3600;
   intime=time;
   loc=(lonx+7.5)/15;
   locx(ibg)=loc;
   %Rnl=.97*(5.67e-8*(tsea-dter*jcool+273.16)^4-Rl);
   Rnl=.97*(5.67e-8*(ts-dter*jcool+273.16)^4-Rl);%oceanic broadband emissivity=0.97
   rnl(ibg)=Rnl;
   Rns=.945*Rs;%oceanic albedo=0.055 daily average
   %*********   set condition dependent stuff ******
   Le=(2.501-.00237*tsea)*1e6;
   cpv=cpa*(1+0.84*q/1000);
   rhoa=P*100/(Rgas*(t+tdk)*(1+0.61*q/1000));
   visa=1.326e-5*(1+6.542e-3*t+8.301e-6*t*t-4.84e-9*t*t*t);
   Al=2.1e-5*(tsea+3.2)^0.79;

   %**************   apply warm layer  ***********;
   if jwarm==1;                                            %do warm layer
      chktime=loc+intime*24;
      newtime=(chktime-24*fix(chktime/24))*3600;
      if icount>1                                     %not first time thru
         if newtime>21600 & jump==1
%              e=[num2str(icount) '  ' num2str(newtime) '  ' num2str(jtime) '   ' num2str(jump) '  ' num2str(q_pwp)  ];
%              disp(e)
            %goto 16
         else
            jump=0;
            if newtime <jtime    %re-zero at midnight
               jamset=0;
               fxp=.5;
               tk_pwp=max_pwp;
               tau_ac=0;
               qcol_ac=0;
               dt_wrm=0;
               jump=0;                   %goto 16
            else
                  %****   set warm layer constants  ***
               rich=.65;                     %crit rich  
               ctd1=sqrt(2*rich*cpw/(Al*grav*rhow));
               ctd2=sqrt(2*Al*grav/(rich*rhow))/(cpw^1.5);
                  %*********************************
               dtime=newtime-jtime;       %delta time for integrals
               qr_out=Rnl+hs_old+hl_old+RF_old; %total cooling at surface
               q_pwp=fxp*Rns-qr_out;         %tot heat abs in warm layer
               if q_pwp<50 & jamset==0       %check for threshold
                  %goto 16    
               else
                  jamset=1;         %indicates threshold crossed
                  tau_ac=tau_ac+max(.002,tau_old)*dtime; %momentum integral
                  if qcol_ac+q_pwp*dtime>0   %check threshold for warm layer existence
                     for i=1:5      %loop 5 times for fxp
                        fxp=1-(0.28*0.014*(1-exp(-tk_pwp/0.014))+0.27*0.357*(1-exp(-tk_pwp/0.357))+0.45*12.82*(1-exp(-tk_pwp/12.82)))/tk_pwp;
                        %fg=fpaul(tk_pwp);fxp=fg(1);
                        qjoule=(fxp*Rns-qr_out)*dtime;
                        if qcol_ac+qjoule>0
                           tk_pwp=min(max_pwp,ctd1*tau_ac/sqrt(qcol_ac+qjoule));
                        end;
                     end;%  end i loop
                  else           %warm layer wiped out
                     fxp=0.75;
                     tk_pwp=max_pwp;
                     qjoule=(fxp*Rns-qr_out)*dtime;
                  end;%   end sign check on qcol_ac
                  qcol_ac=qcol_ac+qjoule;    %heat integral
                     %*******  compute dt_warm  ******
                  if qcol_ac>0
                     dt_wrm=ctd2*(qcol_ac)^1.5/tau_ac;
                  else 
                     dt_wrm=0;
                  end;
               end;%                    end threshold check
            end;%                            end midnight reset
            if tk_pwp<ts_depth
               dsea=dt_wrm;
            else
               dsea=dt_wrm*ts_depth/tk_pwp;
            end;
         end;%                                    end 6am start first time thru
      end;%                                            end first time thru check
      jtime=newtime;
   end;%  end jwarm,  end warm layer model appl check
      
   ts=tsea+dsea;
   qs=qsee([ts P]);
   qsx(ibg)=qs;
   tsx(ibg)=ts;
   %%Lew.Gramer@noaa.gov: Already calculated above - if not provided by user
   %a=.018;
   %b=.729;
   %twave=b*u;
   %hwave=a*u.^2.*(1+.015*u);
   twave = vtwave(ibg);
   hwave = vhwave(ibg);

   %% Lew.Gramer@noaa.gov: Adjust surface current with a quasi-Eulerian term which includes
   % Stokes drift (Arduin et al., 2009)
   %us=us+f(twave,hwave);                       %surface current

   x=[u us ts t qs q Rs Rl rain pblz  P zu zt zq lat jcool jwave twave hwave ] ;    %set data for basic flux alogithm
      %********    call modified LKB routine *******
   y=cor30a(x);
      %************* output from routine  *****************************
        hsb=y(1);                   %sensible heat flux W/m/m
        hlb=y(2);                   %latent
        taub=y(3);                   %stress
        zo=y(4);                    %vel roughness
        zot=y(5);                   %temp "
        zoq=y(6);                   %hum  "
        L=y(7);                     %Ob Length
        usr=y(8);                   %ustar
        tsr=y(9);                   %tstar
        qsr=y(10);                  %qstar  [g/g]
        dter=y(11);                 %cool skin delta T
        dqer=y(12);                 %  "   "     "   q
        tkt=y(13);                  %thickness of cool skin
        RF=y(14);                   %rain heat flux
        wbar=y(15);                 %webb mean w     
        Cd=y(16);                   %drag @ zu
        Ch=y(17);                   %
        Ce=y(18);                   %Dalton
        Cdn_10=y(19);               %neutral drag @ 10 [includes gustiness]
        Chn_10=y(20);               %
        Cen_10=y(21);               %
        Wg=y(22);
        zax(1)=jd;                  %julian day
        zax(2:10)=x(1:9);           %
        zax(4)=tsea;                %Tsea [no cool skin]
        zax(11:32)=y;               %
        zax(33:35)=[dt_wrm tk_pwp ts];  %warm layer deltaT, thickness, corrected Tsea
   %*******   previous values from cwf hp basic code *****
   
   Hrain(ibg)=RF;
   %**********  new values from this code
   hnet(ibg)=Rns-Rnl-hsb-hlb-Hrain(ibg);%total heat input to ocean
   hs(ibg)=hsb;
   hl(ibg)=hlb;
   tau(ibg)=taub;
   hl_webb(ibg)=rhoa*Le*wbar*qa(ibg)/1000;
   %********************  save various parts of data **********************************
   dt(ibg,:)=zax;

%    %*************  create Bradley type out file
%    indx(ibg)=ibg;
%    out(ibg,1)=ibg;
%    out(ibg,2)=jdy(ibg);%output old results
%    out(ibg,3)=hsb;
%    out(ibg,4)=hlb;
%    out(ibg,5)=ts-dter*jcool;
%    out(ibg,6)=taub;
%    out(ibg,7)=wbar;
%    out(ibg,8)=RF;
%    out(ibg,9)=dter;
%    out(ibg,10)=dt_wrm;
%    out(ibg,11)=tk_pwp;
%    out(ibg,12)=tkt*1e3;
%    out(ibg,13)=Wg; 
   %*************  create Bradley type out file
   out(ibg,1)=hsb;
   out(ibg,2)=hlb;
   out(ibg,3)=taub;
   out(ibg,4)=zo;
   out(ibg,5)=zot;
   out(ibg,6)=zoq;
   out(ibg,7)=L;
   out(ibg,8)=usr;
   out(ibg,9)=tsr;
   out(ibg,10)=qsr;
   out(ibg,11)=dter;
   out(ibg,12)=dqer;
   out(ibg,13)=tkt*1e3;
   out(ibg,14)=RF;
   out(ibg,15)=wbar;
   out(ibg,16)=Cd;
   out(ibg,17)=Ch;
   out(ibg,18)=Ce;
   out(ibg,19)=Cdn_10;
   out(ibg,20)=Chn_10;
   out(ibg,21)=Cen_10;
   out(ibg,22)=Wg;
   out(ibg,23)=dt_wrm;
   out(ibg,24)=tk_pwp;
   out(ibg,25)=ts-dter*jcool;

%  e=[num2str(mon) '  ' num2str(iday) '  ' num2str(ihr) '   ' num2str(dter) '  ' num2str(dt_wrm) '  ' num2str(tk_pwp) ];
%  disp(e)
   hs_old=hsb;
   hl_old=hlb;
   RF_old=RF;
   tau_old=taub;
   icount=icount+1;
end; %  data line loop

%%%% After preallocating all arrays to speed this loop up (Lew.Gramer@noaa.gov)
%DEBUG:toc,

% figure;plot(jdx,out(:,3),'-o');xlabel('Julian Day');ylabel('Sensible Heat Flux (W/m^2)');
% figure;plot(jdx,out(:,4)','-o');xlabel('Julian Day');ylabel('Latent Heat Flux (W/m^2)');
% figure;plot(jdx,out(:,6),'-o');xlabel('Julian Day');ylabel('Stress (N/m^2)');
% figure;plot(jdx,out(:,9),jdx,out(:,10));xlabel('Julian Day');ylabel('DT cool and warm (C)');

% %*****************   write output file  ******
% prtit=1;
% if prtit==1;
 
%     save 'c:\data\cwf\matlabstf\cwf\bulkalg\outs\tst3_0ah_out' out
%  e= 'c:\data\cwf\matlabstf\cwf\bulkalg\outs\tst3_0ah_out.txt' ;%&&&&&&&&&&&&&&

%    flist=fopen(e,'w');
%  fprintf(flist,'%6i %18.0f %8.2f %8.2f %8.2f %9.5f %9.5f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f \r\n',out'); 
%  %           1     2  3     4   5    6    7     8       9  10   11  12    13  
%  fclose('all');
% end;

return;
