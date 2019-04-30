1;

if ( ~exist('swtemp','var') )
  error('Run PLOT_UPWELLING_HOVMOELLER first!');
end;

% No, thank you it is very helpful! I am doing upwelling calculations - converting from a volume of water/day x a uM concentration into a kg/m3/day estimate.
%
% Lew
%
%
% On Wed, Aug 2, 2017 at 5:29 PM, Charles Fischer - NOAA Federal <charles.fischer@noaa.gov> wrote:
% Hi Lew,
% Your request is a little complicated and might require a one on one explanation.    First, the Kg/M3 are units that aren't normally use, that I'm aware of, except maybe in sewage or water treatment calculations.  We use moles because its very easy to use in calculations and for that we we the formula weight of the salt being used: for example: KH2PO4 = 136.09; KNO3 = 101.11; NaNO2 = 69.0; NaSiF6 = 188.06.   To calculate what you are asking for requires one to use the atomic weights of the  elements involved.   For example; PO4 = 94.97.  Biologist like to use just P and in that case the atomic weight would be 30.97.  And listed as Kg-P/M3.    For NO3+ NO2, what you are looking at is total NO2, because you are reducing the NO3  to NO2 by the Cd column.  (This is because NO3 can not be measured directly.)  Then, we subtract the NO3+NO2 from NO2 that was calculated on a separate channel to calculate NO3.  So, the atomic weights used are 14 (N) and 32 (O x2) = 46.  Now, if you're just looking at N, then one uses 14 and show the result as Kg-N/M3.  Finally with Si, since you are using weight and not moles, then the molecule you are using is SiO2 (Silicate).   It's atomic weight is 60.09.  When we use moles we usually just list it as Si (Silicon).  However, we call it Silicate.    
% FYI: 1 ug = 10e-9 Kg. and 1 L = .001 M3.
% An example of a calculation for PO4:  1u-mole/L  x 94.97ug/1u-mole x 10e-9 Kg/1ug x 1L/0.001M3 = 9.497 e-5 Kg/M3.
% An example of a calculation with P:  1 u-mole P/L x 30.97 ug/1u-mole-P x 10e-9 Kg/1ug x 1L/0.001M3 = 3.097 e-5 Kg-P/M3   or 0.03097g-P/M3.
%
% Hope this helps and not too confusing.
%
% On Wed, Aug 2, 2017 at 2:35 PM, Lew Gramer - NOAA Affiliate <lew.gramer@noaa.gov> wrote:
%
% Thanks, Charlie! If I wanted to convert from uM/liter to kg/m3, is there an "average" molecular weight for each of NO3+NO2, P, and Si that I could use?
%
% Lew
%
%
% On Wed, Aug 2, 2017 at 12:18 PM, Charles Fischer - NOAA Federal <charles.fischer@noaa.gov> wrote:
% Lew,
% References witch should be of help are attached.  Also, I should point out that we look at uM values rather than m-mol/L.  The papers may have been modified a little.
%
% On Wed, Aug 2, 2017 at 11:40 AM, Lew Gramer - NOAA Affiliate <lew.gramer@noaa.gov> wrote:
%
% Hi, Charlie, I had another (hopefully quick) question for you: Is there a good reference or Web page I can cite for the methods used in analyzing mmol of NO2+NO3, P, and Si for the various coastal cruises (FACE, SFP) done out of our Lab?
%
% Thanks,
% Lew


% 14.00 kg-N/m3
% 30.97 kg-P/m3
% 60.09 kg-Si/m3
%
% 1 u-mole P/L x 30.97 ug/1u-mole-P x 10e-9 Kg/1ug x 1L/0.001M3 = 3.097 e-5 Kg-P/M3   or 0.03097g-P/M3

% Nkgm3 = umol_to_kgm3(N,'N');

% [X,T] = intersect_tses(sfomc.ne_buoy.adcp_x_btm,sfomc.ne_buoy.seatemp);
% dep = sfomc.ne_buoy.depth;
[X,T] = intersect_tses(sfomc.sw_buoy.adcp_x_btm,sfomc.sw_buoy.seatemp);
dep = sfomc.ne_buoy.depth;
% % [X,T] = intersect_tses(face.deep.adcp_btm_x,face.deep.seatemp);
% % dep = face.deep.depth;

dt = median(diff(T.date))*3600*24; %[s]

nute_flux = 0;
dy = 130; %[km/reef-system]
if ( ~isfield(T,'prof') )
  dz = dep; %[m/water-column]
  Tprof = T.data;
  dN = seatemp_to_nutes(Tprof); %[K] -> [kgN/m^3]
  dN_per_s = dN.*X.data.*dz.*dy.*1e3; %[kgN/m^3]*[m/s]*m*m = [kgN/s/reef-system]
  nute_flux_per_m2 = dt.*dN.*X.data; %[s * m/s * kgN*s/m^3] = [kgN/m^2]
  nute_flux = nute_flux + (dz*1e3*nute_flux_per_m2); %[m*m/km*kgN/m^2] = [kgN/km]
else
  for binix=2:size(T.prof,2)
    Tprof=T.prof(:,binix);
    dz = abs(T.depths(binix) - T.depths(binix-1)); %[m] bin height
    dN = seatemp_to_nutes(Tprof); %[K] -> [kgN/m^3]
    dN_per_s = dN.*X.data.*dz.*dy.*1e3; %[kgN/m^3]*[m/s]*m*m = [kgN/s/reef-system]
    nute_flux_per_m2 = dt.*dN.*X.data; %[s * m/s * kgN*s/m^3] = [kgN/m^2]
    nute_flux = nute_flux + (dz*1e3*nute_flux_per_m2); %[m*m/km*kgN/m^2] = [kgN/km]
  end;
end;
fmg; plot(X.date,cumsum(dy.*nute_flux)); datetick3;
% Lee et al. 1991, SAB, p. 22,200: "about 25 kg N s-1 over the total 278-km alongshore distance"
fmg; plot(X.date,dN_per_s); datetick3;


% % IMPLIED HEAT FLUX ESTIMATES: Year:
% % 3600 s x rho_w x c_p x h x nHours x (Ta-Ts)
%
% % 1999
% (24*(1025*3999*11)/3600) * 1.8 * 4		   % 2.1643     
% (24*(1025*3999*11)/3600) * 3.6 * 3                 % 3.2464     
% % 2000                                             
% (24*(1025*3999*11)/3600) * 3.8 * 7                 % 7.9957     
% (24*(1025*3999*11)/3600) * 0.5 * 3                 % 0.45089    
% (24*(1025*3999*11)/3600) * 2.4 * 9                 % 6.4928     
% (24*(1025*3999*11)/3600) * 1.0 * 2                 % 0.601183   
% % 2001
% 0                                                  %   0        
% % 2002                                             
% 0                                                  %   0        
% % 2003                                             
% (24*(1025*3999*11)/3600) * 1.5 * 3                 % 1.3527     
% (24*(1025*3999*11)/3600) * 2.7 * 5                 % 4.0580     
% (24*(1025*3999*11)/3600) * 1.3 * 6                 % 2.3446     
% % 2004                                             
% (24*(1025*3999*11)/3600) * 1.4 * 3                 % 1.2625     
% (24*(1025*3999*11)/3600) * 2.3 * 6                 % 4.1482     
% (24*(1025*3999*11)/3600) * 7.5 * 30                % 67.633     
% % 2005                                             
% (24*(1025*3999*11)/3600) * 1.9 * 6                 % 3.4267     
% (24*(1025*3999*11)/3600) * 3.7 * 10                % 11.122     
% % 2006                                             
% (24*(1025*3999*11)/3600) * 1.3 * 7                 % 2.7354     
% (24*(1025*3999*11)/3600) * 1.7 * 9                 % 4.5990     
% (24*(1025*3999*11)/3600) * 2.1 * 7                 % 4.4187     
% (24*(1025*3999*11)/3600) * 2.2 * 5                 % 3.3065     
% % 2007                                             
% (24*(1025*3999*11)/3600) * 1.4 * 4                 % 1.6833     
% (24*(1025*3999*11)/3600) * 3.2 * 12                % 11.543     
% (24*(1025*3999*11)/3600) * 1.1 * 3                 % 0.99195    
% (24*(1025*3999*11)/3600) * 1.4 * 6                 % 2.5250     
% % 2008
% 0                                                  %   0         
% % 2009                                             
% (24*(1025*3999*11)/3600) * 3.9 * 8                 % 9.3785     
% (24*(1025*3999*11)/3600) * 4.2 * 6                 % 7.5749     
% (24*(1025*3999*11)/3600) * 1.3 * 3                 % 1.1723     
% % 2010                                             
% (24*(1025*3999*11)/3600) * 4.2 * 4                 % 5.0499     
% (24*(1025*3999*11)/3600) * 3.2 * 3                 % 2.8857     
% (24*(1025*3999*11)/3600) * 3.6 * 4                 % 4.3285     
% % 2011                                             
% (24*(1025*3999*11)/3600) * 5.4 * 4                 % 6.4928     
% (24*(1025*3999*11)/3600) * 1.5 * 4                 % 1.8035     
% (24*(1025*3999*11)/3600) * 4.3 * 4                 % 5.1702     
% (24*(1025*3999*11)/3600) * 2.5 * 6                 % 4.5089	


%% ACROSS-SHORE MIXING RATE

% Recalculate mixing rates using method from PLOT_UPWELLING_HOVMOELLER

mixrate=[]; clear mixrate

mixrate.date = swtemp.date;

mixrate.data = (1-(nwtemp.data-swtemp.prof(:,end))./(Tsfcpt-swtemp.prof(:,end)));
mixrate.data(0>mixrate.data) = 0;
mixrate.data(mixrate.data>1) = nan;

mixrate.prof = (1-(swtemp.prof-repmat(Tendpt,[1,size(swtemp.prof,2)])) ...
           ./ repmat((Tsfcpt-Tendpt),[1,size(swtemp.prof,2)]));
mixrate.prof(0>mixrate.prof) = 0;
mixrate.prof(mixrate.prof>1) = nan;

mixrate.depths = swtemp.depths;


fmg;
boxplot(mixrate.prof(:,4),get_month(mixrate.date),'positions',get_month(mixrate.date));
xticklabels({'Jan','Feb','Mar','Jul','Aug','Sep','Oct','Nov','Dec'});
axis([6.5,10.5,0,1]);
if doPrint; print('-dpng',fullfile(figspath,'upwelling_mixrate_cross_boxplot.png')); end;


%% ALONG-SHORE MIXING RATE

4*(17-14)/station_dist(sefcri.updb,sefcri.pb2)
% 0.20
4*(22-17)/station_dist(sefcri.pb2,sefcri.bc3)
% 0.34
4*(25-22)/station_dist(sefcri.bc3,sefcri.dc3)
% 0.34

fmg; plot(1:3,[0.2,0.34,0.34],'-p'); axis([0,4,0,1]); xticks([1:3]); xticklabels({'Far North-North','North-Central','Central-South'});
if doPrint; print('-dpng',fullfile(figspath,'upwelling_mixrate_along_boxplot.png')); end;



N=25; %Gmol N/yr: SAB=10^11 mol N/yr, NFRT=SAB/4 = 25e9 mol N/yr

Nfld = lkwf1.ngdc_hires_bathy.field;
Nfld(-50>=Nfld | Nfld>=0) = nan;
Nfld = (1 - ((-50-Nfld)/(-50)) ).*N;

%{
fmg;
contourf(lkwf1.ngdc_hires_bathy.lon,lkwf1.ngdc_hires_bathy.lat,Nfld);
colorbar;
plot_hires_coastline(lkwf1.ngdc_hires_bathy);
axis([-80.2,-79.9,25.54,27.25]);
daspect([1,5,1]);
titlename('NO LATITUDE DILUTION');
%}

%attenuate with latitude from a peak at 27.23
%(27.23-26.67)*111 * 0.2
[Nlon,Nlat] = meshgrid(lkwf1.ngdc_hires_bathy.lon,lkwf1.ngdc_hires_bathy.lat);
Nlat = (27.69 - Nlat)./(27.69 - 25.0);
Nfld = Nfld.*(1-Nlat);

fmg;
contourf(lkwf1.ngdc_hires_bathy.lon,lkwf1.ngdc_hires_bathy.lat,Nfld);
cbh=colorbar;
plot_hires_coastline(lkwf1.ngdc_hires_bathy);
axis([-80.2,-79.9,25.54,27.25]);
xlabel(cbh,'Gmol NO_x/yr');
if doPrint; print('-dpng',fullfile(figspath,'upwelling_NOx_flux_map.png')); end;
pause;
daspect([1,5,1]);
if doPrint; print('-dpng',fullfile(figspath,'upwelling_NOx_flux_map_exaggerated.png')); end;

