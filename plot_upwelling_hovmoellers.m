1;
%% SCRIPT plot_upwelling_hovmoellers.m:
%
% Display and optionally (DEFAULT: doPrint=False) print Hovmoeller diagrams
% (vs. time and depth) of interpolated sea temperature, ocean currents along-
% and cross-shore, and nutrient fluxes from SFOMC moorings, especially during
% periods of expected upwelling (summer) onto the southeast Florida shelf.
%
% PARAMETERS [DEFAULTS]:
%
%  basenm - root filename for printed figures ['upwelling_hovmoellers_SFOMC']
%  do2000 - plot the upwelling season for 2000? [FALSE - plot for 1999 instead]
%  doCums - plot cumulative sums of nutrient inputs? [TRUE]
%  doHovs - plot Hovmoellers (using SURF)? [TRUE]
%  doKgs  - plot onshore fluxes in [kg/m] (using SURF)? [TRUE]
%  doPrint - PRINT plotted figures? [FALSE]
%  figspath - path in which to PRINT figures [.../CRCP/Upwelling/CoRIS]
%
% Last Saved Time-stamp: <Fri 2018-04-13 14:40:44 Eastern Daylight Time gramer>


if ( ~exist('figspath','var') )
  figspath = fullfile(get_coral_path,'CRCP','Upwelling','CoRIS');
end;
if ( ~exist('doHovs','var') || isempty(doHovs) )
  doHovs = true;
end;
if ( ~exist('doCums','var') || isempty(doCums) )
  doCums = true;
end;
if ( ~exist('doKgs','var') || isempty(doKgs) )
  doKgs = true;
end;
if ( ~exist('doPrint','var') || isempty(doPrint) )
  doPrint = false;
end;

if ( ~exist('do2000','var') || isempty(do2000) )
  % Shall we plot the upwelling season for 1999 (DEFAULT) or 2000?
  do2000 = false;
end;

if ( ~exist('fontsz','var') || isempty(fontsz) )
  fontsz = 16;
end;

if ( ~exist('xl','var') || isempty(xl) )
  if ( do2000 )
    xl = datenum(2000,[6,7],20); %Upwelling period
  else
    %xl = sfomc.ne_buoy.adcp_u.date([1,end]); %Whole record
    %xl = datenum(1999,7,[25,31]); %Transition to high-frequency
    %xl = datenum(1999,7,[15,22]); %Start of records
    %xl = datenum(1999,[7,8],20); %Upwelling period
    xl = datenum(1999,[7,8],[15,29]); %Broader upwelling period
  end;
end;

%  minspd - speed below which SURF should not contour data (0 == no min. speed)
if ( ~exist('minspd','var') || isempty(minspd) )
  %minspd = 0.02; %Roughly the measurement accuracy of the ADCP
  %minspd = 0.01; %Less than the measurement accuracy of the ADCP
  minspd = 0; %Turn off this control
end;

if ( ~exist('basenm','var') || isempty(basenm) )
  basenm = ['upwelling_hovmoellers_SFOMC'];
end;
fignm = [basenm,datestr(xl(1),'_yyyymmmdd'),datestr(xl(end),'_yyyymmmdd')];
dtstr = [datestr(xl(1),'mmm-dd'),' to ',datestr(xl(end),'mmm-dd-yyyy')];


if ( doPrint )
  disp(['Will print figures to ',fullfile(figspath,[fignm,'*'])]);
else
  disp(['Will not print figures']);
end;


% Subset all SFOMC mooring data to the same period for speed of plotting
if ( ~exist('nw','var') || any(nw.xl~=xl) || (nw.minspd~=minspd) )
  if ( do2000 )
    subfn = @(x)(ts_date_range(x,sfomc.e_buoy.adcp_l.date([1,end])));
  else
    subfn = @(x)(ts_date_range(x,sfomc.ne_buoy.adcp_l.date([1,end])));
  end;

  nw=[]; clear nw
  nw.xl = xl; nw.minspd = minspd;
  nw.adcp_depths = sfomc.nw_w_btm.adcp_depths;
  nw.adcp_speed = subset_ts(sfomc.nw_w_btm.adcp_speed,subfn); %nw.adcp_speed.prof(abs(nw.adcp_speed.prof)<minspd) = nan;
  nw.adcp_u = subset_ts(sfomc.nw_w_btm.adcp_u,subfn); %nw.adcp_u.prof(isnan(nw.adcp_speed.prof)) = nan;
  nw.adcp_v = subset_ts(sfomc.nw_w_btm.adcp_v,subfn); %nw.adcp_v.prof(isnan(nw.adcp_speed.prof)) = nan;
  nw.adcp_x = subset_ts(sfomc.nw_w_btm.adcp_x,subfn); %nw.adcp_x.prof(isnan(nw.adcp_speed.prof)) = nan;
  nw.adcp_l = subset_ts(sfomc.nw_w_btm.adcp_l,subfn); %nw.adcp_l.prof(isnan(nw.adcp_speed.prof)) = nan;
  nw.adcp_seatemp = subset_ts(sfomc.nw_w_btm.adcp_seatemp,subfn);
  nw.sbe_seatemp = subset_ts(sfomc.nw_w_btm.sbe_seatemp,subfn);
  [nw.seatemp,s1] = intersect_tses(nw.adcp_seatemp,nw.sbe_seatemp);
  nw.seatemp.depths = [-5,-11];
  nw.seatemp.prof(:,2) = nw.seatemp.data;
  nw.seatemp.prof(:,1) = s1.data;
end;
if ( ~exist('sw','var') || any(sw.xl~=xl) || (sw.minspd~=minspd) )
  if ( do2000 )
    sw=[]; clear sw
    sw.xl = xl; sw.minspd = minspd;
    sw.adcp_depths = sfomc.c_buoy.adcp_depths;
    sw.adcp_speed = sfomc.c_buoy.adcp_speed; %sw.adcp_speed.prof(abs(sw.adcp_speed.prof)<minspd) = nan;
    sw.adcp_u = sfomc.c_buoy.adcp_u; %sw.adcp_u.prof(isnan(sw.adcp_speed.prof)) = nan;
    sw.adcp_v = sfomc.c_buoy.adcp_v; %sw.adcp_v.prof(isnan(sw.adcp_speed.prof)) = nan;
    sw.adcp_x = sfomc.c_buoy.adcp_x; %sw.adcp_x.prof(isnan(sw.adcp_speed.prof)) = nan;
    sw.adcp_l = sfomc.c_buoy.adcp_l; %sw.adcp_l.prof(isnan(sw.adcp_speed.prof)) = nan;
    sw.sfctemp = sfomc.c_buoy.at_0_6m.sbe_seatemp;
    sw.seatemp = sfomc.c_buoy.seatemp;
  else
    sw=[]; clear sw
    sw.xl = xl; sw.minspd = minspd;
    sw.adcp_depths = sfomc.sw_buoy.adcp_depths;
    sw.adcp_speed = sfomc.sw_buoy.adcp_speed; %sw.adcp_speed.prof(abs(sw.adcp_speed.prof)<minspd) = nan;
    sw.adcp_u = sfomc.sw_buoy.adcp_u; %sw.adcp_u.prof(isnan(sw.adcp_speed.prof)) = nan;
    sw.adcp_v = sfomc.sw_buoy.adcp_v; %sw.adcp_v.prof(isnan(sw.adcp_speed.prof)) = nan;
    sw.adcp_x = sfomc.sw_buoy.adcp_x; %sw.adcp_x.prof(isnan(sw.adcp_speed.prof)) = nan;
    sw.adcp_l = sfomc.sw_buoy.adcp_l; %sw.adcp_l.prof(isnan(sw.adcp_speed.prof)) = nan;
    sw.sfctemp = sfomc.sw_buoy.at_0_6m.sbe_seatemp;
    sw.seatemp = sfomc.sw_buoy.seatemp;
  end;
end;
if ( ~exist('ne','var') || any(ne.xl~=xl) || (ne.minspd~=minspd) )
  if ( do2000 )
    ne=[]; clear ne
    ne.xl = xl; ne.minspd = minspd;
    ne.adcp_depths = sfomc.e_buoy.adcp_depths;
    ne.adcp_speed = sfomc.e_buoy.adcp_speed; %ne.adcp_speed.prof(abs(ne.adcp_speed.prof)<minspd) = nan;
    ne.adcp_dir = sfomc.e_buoy.adcp_dir; %ne.adcp_speed.prof(abs(ne.adcp_speed.prof)<minspd) = nan;
    ne.adcp_u = sfomc.e_buoy.adcp_u; %ne.adcp_u.prof(isnan(ne.adcp_speed.prof)) = nan;
    ne.adcp_v = sfomc.e_buoy.adcp_v; %ne.adcp_v.prof(isnan(ne.adcp_speed.prof)) = nan;
    ne.adcp_x = sfomc.e_buoy.adcp_x; %ne.adcp_x.prof(isnan(ne.adcp_speed.prof)) = nan;
    ne.adcp_l = sfomc.e_buoy.adcp_l; %ne.adcp_l.prof(isnan(ne.adcp_speed.prof)) = nan;
    ne.sfctemp = sfomc.e_buoy.at_5m.sbe_seatemp;
    ne.seatemp = sfomc.e_buoy.seatemp;
  else
    ne=[]; clear ne
    ne.xl = xl; ne.minspd = minspd;
    ne.adcp_depths = sfomc.ne_buoy.adcp_depths;
    ne.adcp_speed = sfomc.ne_buoy.adcp_speed; %ne.adcp_speed.prof(abs(ne.adcp_speed.prof)<minspd) = nan;
    ne.adcp_dir = sfomc.ne_buoy.adcp_dir; %ne.adcp_speed.prof(abs(ne.adcp_speed.prof)<minspd) = nan;
    ne.adcp_u = sfomc.ne_buoy.adcp_u; %ne.adcp_u.prof(isnan(ne.adcp_speed.prof)) = nan;
    ne.adcp_v = sfomc.ne_buoy.adcp_v; %ne.adcp_v.prof(isnan(ne.adcp_speed.prof)) = nan;
    ne.adcp_x = sfomc.ne_buoy.adcp_x; %ne.adcp_x.prof(isnan(ne.adcp_speed.prof)) = nan;
    ne.adcp_l = sfomc.ne_buoy.adcp_l; %ne.adcp_l.prof(isnan(ne.adcp_speed.prof)) = nan;
    ne.sfctemp = sfomc.ne_buoy.at_0_6m.sbe_seatemp;
    ne.seatemp = sfomc.ne_buoy.seatemp;
  end;
end;


if ( doHovs )

  fmg;
  tpax(1) = subplot(8,1,1);
  plot_ts(fwyf1.ndbc_air_t,'k-',nw.adcp_seatemp);
  ylim([26.5,31.5]); %colorbar; %Adding useless COLORBAR just lines up the X axes nicely
  xlim(xl); %datetick3;
  grid on;
  ylabel('NW'); set(gca,'FontSize',fontsz);
  titlename(['Temperature Hovmuellers ',dtstr]);
  
  tpax(2) = subplot(8,1,2:3);
  surf(sw.seatemp.date,sw.seatemp.depths,sw.seatemp.prof');
  view(0,270); shading interp; set(gca,'YDir','reverse');
  hold on;
  contour(sw.seatemp.date,sw.seatemp.depths,sw.seatemp.prof',[28.8,28.8],'LineWidth',0.5,'Color','k');
  caxis([25.8,31.6]); colorbar;
  ylim([-sfomc.sw_buoy.depth,0]); 
  for dp=sw.seatemp.depths(:)'; annotline([],dp); end;
  xlim(xl); %datetick3;
  ylabel('Mid (SW/C)'); set(gca,'FontSize',fontsz);

  tpax(3) = subplot(8,1,4:8);
  surf(ne.seatemp.date,ne.seatemp.depths,ne.seatemp.prof'); 
  view(0,270); shading interp; set(gca,'YDir','reverse');
  hold on;
  [c,h]=contour(ne.seatemp.date,ne.seatemp.depths,ne.seatemp.prof',[23,23],'LineWidth',0.5,'Color','k');
  caxis([12,32]); colorbar;
  ylim([-sfomc.ne_buoy.depth,0]); 
  for dp=ne.seatemp.depths(:)'; annotline([],dp); end;
  xlim(xl); datetick('x',2,'keeplimits'); %datetick3;
  ylabel('Outer (NE/E) Sea Temp. [\circC]'); set(gca,'FontSize',fontsz);

  set(tpax(1:2),'XTickLabel','');
  align_subplot_to_subplot(tpax(1:2));
  if doPrint; print('-dpng',fullfile(figspath,[fignm,'_temps.png'])); end;


  fmg;
  xsax(1) = subplot(8,1,1);
  surf(nw.adcp_x.date,nw.adcp_depths,nw.adcp_x.prof'); 
  view(0,270); shading interp; set(gca,'YDir','reverse');
  hold on;
  contour(nw.adcp_x.date,nw.adcp_depths,nw.adcp_x.prof',[0,0],'LineWidth',0.5,'Color','k');
  colorbar;
  ylim([-sfomc.nw_w_btm.depth,0]); 
  xlim(xl); %datetick3;
  ylabel('NW'); set(gca,'FontSize',fontsz);
  titlename(['Cross-shore current Hovmuellers ',dtstr]);
  
  xsax(2) = subplot(8,1,2:3);
  surf(sw.adcp_x.date,sw.adcp_depths,sw.adcp_x.prof'); 
  view(0,270); shading interp; set(gca,'YDir','reverse');
  hold on;
  contour(sw.adcp_x.date,sw.adcp_depths,sw.adcp_x.prof',[0,0],'LineWidth',0.5,'Color','k');
  %caxis([-0.1,+0.1]);
  colorbar;
  ylim([-sfomc.sw_buoy.depth,0]); 
  xlim(xl); %datetick3;
  ylabel('SW'); set(gca,'FontSize',fontsz);
  
  xsax(3) = subplot(8,1,4:8);
  surf(ne.adcp_x.date,ne.adcp_depths,ne.adcp_x.prof'); 
  view(0,270); shading interp; set(gca,'YDir','reverse');
  hold on;
  contour(ne.adcp_x.date,ne.adcp_depths,ne.adcp_x.prof',[0,0],'LineWidth',0.5,'Color','k');
  %caxis([-0.1,+0.1]);
  caxis([-0.9,+0.9]);
  colorbar;
  ylim([-sfomc.ne_buoy.depth,0]); 
  xlim(xl); datetick('x',2,'keeplimits'); %datetick3;
  ylabel('NE (Outer) Cross'); set(gca,'FontSize',fontsz);
  
  set(xsax(1:2),'XTickLabel','');
  if doPrint; print('-dpng',fullfile(figspath,[fignm,'_across.png'])); end;


  fmg;
  lsax(1) = subplot(8,1,1);
  surf(nw.adcp_l.date,nw.adcp_depths,nw.adcp_l.prof'); 
  view(0,270); shading interp; set(gca,'YDir','reverse');
  hold on;
  contour(nw.adcp_l.date,nw.adcp_depths,nw.adcp_l.prof',[0,0],'LineWidth',0.5,'Color','k');
  %caxis([-0.1,+0.1]);
  colorbar;
  ylim([-sfomc.nw_w_btm.depth,0]); 
  xlim(xl); %datetick3;
  ylabel('NW'); set(gca,'FontSize',fontsz);
  titlename(['Alongshore current Hovmuellers ',dtstr]);
  
  lsax(2) = subplot(8,1,2:3);
  surf(sw.adcp_l.date,sw.adcp_depths,sw.adcp_l.prof'); 
  view(0,270); shading interp; set(gca,'YDir','reverse');
  hold on;
  contour(sw.adcp_l.date,sw.adcp_depths,sw.adcp_l.prof',[0,0],'LineWidth',0.5,'Color','k');
  %caxis([-0.1,+0.1]);
  colorbar;
  ylim([-sfomc.sw_buoy.depth,0]); 
  xlim(xl); %datetick3;
  ylabel('SW'); set(gca,'FontSize',fontsz);
  
  lsax(3) = subplot(8,1,4:8);
  surf(ne.adcp_l.date,ne.adcp_depths,ne.adcp_l.prof');
  view(0,270); shading interp; set(gca,'YDir','reverse');
  hold on;
  contour(ne.adcp_l.date,ne.adcp_depths,ne.adcp_l.prof',[0,0],'LineWidth',0.5,'Color','k');
  %caxis([-0.1,+0.1]);
  colorbar;
  ylim([-sfomc.ne_buoy.depth,0]); 
  xlim(xl); datetick('x',2,'keeplimits'); %datetick3;
  ylabel('NE (Outer) Along'); set(gca,'FontSize',fontsz);

  set(lsax(1:2),'XTickLabel','');
  if doPrint; print('-dpng',fullfile(figspath,[fignm,'_along.png'])); end;

end; %if ( doHovs )



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CALCULATE NUTRIENT CONCENTRATIONS and MIXING RATES:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Nitrate + Nitrite
ne.N = ne.seatemp;
ne.N.prof = (23-ne.seatemp.prof).*1.62; ne.N.prof(ne.N.prof<0) = 0;
ne.N.data = nansum(ne.N.prof,2);

% Phosphates
ne.P = ne.seatemp;
ne.P.prof = (23-ne.seatemp.prof).*0.073; ne.P.prof(ne.P.prof<0) = 0;
ne.P.data = nansum(ne.P.prof,2);

% Silicates
ne.Si = ne.seatemp;
ne.Si.prof = (23-ne.seatemp.prof).*0.37; ne.Si.prof(ne.Si.prof<0) = 0;
ne.Si.data = nansum(ne.Si.prof,2);


% Concentrations at other sites are approximated from assumed heat mixing rates

% Which is the right temperature end-member?
% >> fmg; plot_ts(sfomc.ne_buoy.at_0_6m.sbe_seatemp,sfomc.sw_buoy.at_0_6m.sbe_seatemp,sfomc.nw_w_btm.sbe_seatemp,fwyf1.ndbc_air_t,fwyf1.ndbc_sea_t); xlim(sfomc.ne_buoy.at_0_6m.sbe_seatemp.date([1,end])); datetick3; legend('NE','SW','NW','FWY','T_s');

if 0;
  plot_hires_bathymetry(sfomc.c_buoy,-[0:120]);
  plot_all_se_florida_upwelling_sites;
  th=text(sfomc.lons([1,3,4,6])+0.001,sfomc.lats([1,3,4,6]),textize(upper(sfomc.station_names([1,3,4,6])))); set(th,'Color','w');
  axis([sfomc.c_buoy.ngdc_hires_bathy.lon(1),sfomc.c_buoy.ngdc_hires_bathy.lon(end),26.015,26.095]);
end;
if 0;
  fmg;
  spt(3,1,1); plot_bathy_transect(sfomc.c_buoy,[],[45,90,135]); ylim([-120,5]);
  spt(3,1,2); plot_bathy_transect({sfomc.c_buoy,sfomc.sw_buoy},[],[45,90,135]); ylim([-120,5]);
  spt(3,1,3); plot_bathy_transect({sfomc.c_buoy,face.deep},[],[45,90,135]); ylim([-120,5]);
end;

if ( do2000 )
  % % [nenitr,nephos,nesili, netemp,swtemp,nwtemp, sfctemp, airtemp] = ...
  % %     intersect_tses(ne.N,ne.P,ne.Si, ne.seatemp,sw.seatemp,nw.adcp_seatemp, fwyf1.ndbc_air_t, fwyf1.ndbc_air_t);
  % [nenitr,nephos,nesili, netemp,swtemp,nwtemp, sfctemp, airtemp] = ...
  %     intersect_tses(ne.N,ne.P,ne.Si, ne.seatemp,sw.seatemp,nw.adcp_seatemp, sw.sfctemp, fwyf1.ndbc_air_t);
  [nenitr,nephos,nesili, netemp,swtemp,nwtemp, sfctemp, airtemp] = ...
      intersect_tses(ne.N,ne.P,ne.Si, ne.seatemp,sw.seatemp,nw.seatemp, sw.sfctemp, fwyf1.ndbc_air_t);
else
  % % [nenitr,nephos,nesili,netemp,swtemp,nwtemp,sfctemp, airtemp] = ...
  % %     intersect_tses(ne.N,ne.P,ne.Si,ne.seatemp,sw.seatemp,nw.adcp_seatemp,sfomc.ne_buoy.at_0_6m.sbe_seatemp, fwyf1.ndbc_air_t);
  % % [nenitr,nephos,nesili, netemp,swtemp,nwtemp, sfctemp, airtemp] = ...
  % %     intersect_tses(ne.N,ne.P,ne.Si, ne.seatemp,sw.seatemp,nw.adcp_seatemp, sw.sfctemp, fwyf1.ndbc_air_t);
  % [nenitr,nephos,nesili,netemp,swtemp,nwtemp,sfctemp, airtemp] = ...
  %     intersect_tses(ne.N,ne.P,ne.Si, ne.seatemp,sw.seatemp,nw.adcp_seatemp, nw.sbe_seatemp, fwyf1.ndbc_air_t);
  [nenitr,nephos,nesili,netemp,swtemp,nwtemp,sfctemp, airtemp] = ...
      intersect_tses(ne.N,ne.P,ne.Si, ne.seatemp,sw.seatemp,nw.seatemp, nw.sbe_seatemp, fwyf1.ndbc_air_t);
end;

Tendpt = squeeze(netemp.prof(:,end));
%Tsfcpt = sfctemp.data - 0.8; %0.8=mean summer difference between ADCP and SBE Ts
Tsfcpt = sfctemp.data;
Nendpt = squeeze(nenitr.prof(:,end));
Pendpt = squeeze(nephos.prof(:,end));
Sendpt = squeeze(nesili.prof(:,end));


%%%
% SW/C buoy (20 m isobath) nutrient concentrations

% What is the right mixing ratio for cool water moving onshore over rugose bottom?

% Munk 1966: vertical eddy diffusivity K=1.3e-4 [m^2/s].
% Garrett & Gilbert 2008, in Small-Scale Turbulence...: "diapycnal mixing rate in the deep ocean of 3 to 4e-4 [m^2/s]"
% Lentz et al. 2017: "Synthesis of results from the new field observations with estimates from previous field and laboratory studies indicate that coral reef drag coefficients range from 0.2 to 0.005 and hydrodynamic roughnesses generally range from 2 to 8 cm. While coral reef drag coefficients depend on factors such as physical roughness and surface waves, a substantial fraction of the scatter in estimates of coral reef drag coefficients is due to variations in water depth."

sw.mixprof = (1-(swtemp.prof-repmat(Tendpt,[1,size(swtemp.prof,2)])) ...
              ./ repmat((Tsfcpt-Tendpt),[1,size(swtemp.prof,2)]));

% % fmg; lh=plot(sfctemp.date,sfctemp.data,nwtemp.date,nwtemp.data, swtemp.date,swtemp.prof(:,[1:end]),netemp.date,netemp.prof(:,[1:end])); xlim(datenum(1999,[7,10],15)); datetick3; legend( num2str([-5 ; -11 ; swtemp.depths([1:end])' ; netemp.depths([1:end])']) , 'Location','SouthWest'); set(lh(1:2),'Color','k','LineWidth',2.5); set(lh(2),'LineStyle',':'); set(lh(3:6),'LineWidth',1.5);
% % axis([datenum('19-Jul-1999 21:19:07'),datenum('03-Aug-1999 04:15:02'),14.1,31.2]); datetick3;
% % >> 29.64-26.76
% % ans =
% %     2.8800
% % >> 28.2-18.0
% % ans =
% %    10.2000
% % >> 2.88/10.2
% % ans =
% %     0.2824
% % axis([datenum('19-Jul-1999 21:19:07'),datenum('03-Aug-1999 04:15:02'),14.1,31.2]); datetick3;

% sw.mixprof = repmat(0.28,size(swtemp.prof));


% Ignore non-physical mixing rate profiles
sw.mixprof(0>sw.mixprof | sw.mixprof>1) = 0;

swNendpt = repmat(Nendpt,[1,size(swtemp.prof,2)]);
sw.N = swtemp;
sw.N.prof = swNendpt .* sw.mixprof;
sw.N.data = nansum(sw.N.prof,2);

swPendpt = repmat(Pendpt,[1,size(swtemp.prof,2)]);
sw.P = swtemp;
sw.P.prof = swPendpt .* sw.mixprof;
sw.P.data = nansum(sw.P.prof,2);

swSendpt = repmat(Sendpt,[1,size(swtemp.prof,2)]);
sw.Si = swtemp;
sw.Si.prof = swSendpt .* sw.mixprof;
sw.Si.data = nansum(sw.Si.prof,2);


%%%
% NW/W mooring (11 m isobath) nutrient concentrations

nw.mixprof = (1-(nwtemp.prof-repmat(Tendpt,[1,size(nwtemp.prof,2)])) ...
              ./ repmat((Tsfcpt-Tendpt),[1,size(nwtemp.prof,2)]));
nw.mixprof(0>nw.mixprof | nw.mixprof>1) = 0;

nwNendpt = repmat(Nendpt,[1,size(nwtemp.prof,2)]);
nw.N = nwtemp;
nw.N.prof = nwNendpt .* nw.mixprof;
nw.N.data = nansum(nw.N.prof,2);

nwPendpt = repmat(Pendpt,[1,size(nwtemp.prof,2)]);
nw.P = nwtemp;
nw.P.prof = nwPendpt .* nw.mixprof;
nw.P.data = nansum(nw.P.prof,2);

nwSendpt = repmat(Sendpt,[1,size(nwtemp.prof,2)]);
nw.Si = nwtemp;
nw.Si.prof = nwSendpt .* nw.mixprof;
nw.Si.data = nansum(nw.Si.prof,2);



if ( doHovs )
  fmg;
  ncax(1) = subplot(8,1,1);
  plot_ts(nw.N);
  ylim([0.0,1.0]); %colorbar; %Adding useless COLORBAR just lines up the X axes nicely
  xlim(xl); %datetick3;
  ylabel('NW'); set(gca,'FontSize',fontsz);
  titlename(['Nitrate+Nitrite Concentration Hovmoellers ',dtstr]);

  ncax(2) = subplot(8,1,2:3);
  surf(sw.N.date,sw.N.depths,sw.N.prof');
  view(0,270); shading interp; set(gca,'YDir','reverse');
  hold on;
  caxis([0,1.5]); colorbar;
  ylim([-sfomc.sw_buoy.depth,0]); 
  for dp=sw.seatemp.depths(:)'; annotline([],dp); end;
  xlim(xl); %datetick3;
  ylabel('SW'); set(gca,'FontSize',fontsz);

  ncax(3) = subplot(8,1,4:8);
  surf(ne.N.date,ne.N.depths,ne.N.prof'); 
  view(0,270); shading interp; set(gca,'YDir','reverse');
  hold on;
  caxis([0,15]); colorbar;
  ylim([-sfomc.ne_buoy.depth,0]); 
  for dp=ne.seatemp.depths(:)'; annotline([],dp); end;
  xlim(xl); datetick('x',2,'keeplimits'); %datetick3;
  ylabel('NE (Outer) [NO2+NO3] [\muM]'); set(gca,'FontSize',fontsz);

  set(ncax(1:2),'XTickLabel','');
  align_subplot_to_subplot(ncax(1:2));
  if doPrint; print('-dpng',fullfile(figspath,[fignm,'_nitr.png'])); end;


  fmg;
  pcax(1) = subplot(8,1,1);
  plot_ts(nw.P);
  ylim([0,0.06]); %colorbar; %Adding useless COLORBAR just lines up the X axes nicely
  xlim(xl); %datetick3;
  ylabel('NW'); set(gca,'FontSize',fontsz);
  titlename(['Phosphate Concentration Hovmoellers ',dtstr]);

  pcax(2) = subplot(8,1,2:3);
  surf(sw.P.date,sw.P.depths,sw.P.prof');
  view(0,270); shading interp; set(gca,'YDir','reverse');
  hold on;
  caxis([0,0.1]); colorbar;
  ylim([-sfomc.sw_buoy.depth,0]); 
  for dp=sw.seatemp.depths(:)'; annotline([],dp); end;
  xlim(xl); %datetick3;
  ylabel('SW'); set(gca,'FontSize',fontsz);

  pcax(3) = subplot(8,1,4:8);
  surf(ne.P.date,ne.P.depths,ne.P.prof'); 
  view(0,270); shading interp; set(gca,'YDir','reverse');
  hold on;
  caxis([0,0.6]); colorbar;
  ylim([-sfomc.ne_buoy.depth,0]); 
  for dp=ne.seatemp.depths(:)'; annotline([],dp); end;
  xlim(xl); datetick('x',2,'keeplimits'); %datetick3;
  ylabel('NE (Outer) [PO_4] [\muM]'); set(gca,'FontSize',fontsz);

  set(pcax(1:2),'XTickLabel','');
  align_subplot_to_subplot(pcax(1:2));
  if doPrint; print('-dpng',fullfile(figspath,[fignm,'_phos.png'])); end;



  fmg;
  scax(1) = subplot(8,1,1);
  plot_ts(nw.Si);
  ylim([0,0.25]); %colorbar; %Adding useless COLORBAR just lines up the X axes nicely
  xlim(xl); %datetick3;
  ylabel('NW'); set(gca,'FontSize',fontsz);
  titlename(['Silicate Concentration Hovmoellers ',dtstr]);

  scax(2) = subplot(8,1,2:3);
  surf(sw.Si.date,sw.Si.depths,sw.Si.prof');
  view(0,270); shading interp; set(gca,'YDir','reverse');
  hold on;
  caxis([0,0.3]); colorbar;
  ylim([-sfomc.sw_buoy.depth,0]); 
  for dp=sw.seatemp.depths(:)'; annotline([],dp); end;
  xlim(xl); %datetick3;
  ylabel('SW'); set(gca,'FontSize',fontsz);

  scax(3) = subplot(8,1,4:8);
  surf(ne.Si.date,ne.Si.depths,ne.Si.prof'); 
  view(0,270); shading interp; set(gca,'YDir','reverse');
  hold on;
  caxis([0,2]); colorbar;
  ylim([-sfomc.ne_buoy.depth,0]); 
  for dp=ne.seatemp.depths(:)'; annotline([],dp); end;
  xlim(xl); datetick('x',2,'keeplimits'); %datetick3;
  ylabel('NE (Outer) [Si] [\muM]'); set(gca,'FontSize',fontsz);

  set(scax(1:2),'XTickLabel','');
  align_subplot_to_subplot(scax(1:2));
  if doPrint; print('-dpng',fullfile(figspath,[fignm,'_sili.png'])); end;

end; %if ( doHovs )

if ( doCums )
  fmg;
  naax(1) = subplot(8,1,1);
  plot(nw.N.date,cumsum(nw.N.data)); datetick('x',2);
  %ylim([0,0.5]);
  %colorbar; %Adding useless COLORBAR just lines up the X axes nicely
  xlim(xl); %datetick3;
  ylabel('NW'); set(gca,'FontSize',fontsz);
  titlename(['Nitrate+Nitrite Accumulation Hovmoellers ',dtstr]);

  naax(2) = subplot(8,1,2:3);
  surf(sw.N.date,sw.N.depths,cumsum(sw.N.prof',2,'omitnan'));
  view(0,270); shading interp; set(gca,'YDir','reverse');
  hold on;
  %caxis([0,1.5]);
  colorbar;
  ylim([-sfomc.sw_buoy.depth,0]); 
  for dp=sw.seatemp.depths(:)'; annotline([],dp); end;
  xlim(xl); %datetick3;
  ylabel('SW'); set(gca,'FontSize',fontsz);

  naax(3) = subplot(8,1,4:8);
  surf(ne.N.date,ne.N.depths,cumsum(ne.N.prof',2,'omitnan')); 
  view(0,270); shading interp; set(gca,'YDir','reverse');
  hold on;
  %caxis([0,15]);
  colorbar;
  ylim([-sfomc.ne_buoy.depth,0]); 
  for dp=ne.seatemp.depths(:)'; annotline([],dp); end;
  xlim(xl); datetick('x',2,'keeplimits'); %datetick3;
  ylabel('NE (Outer) [NO2+NO3] [\muM]'); set(gca,'FontSize',fontsz);

  set(naax(1:2),'XTickLabel','');
  align_subplot_to_subplot(naax(1:2));
  if doPrint; print('-dpng',fullfile(figspath,[fignm,'_nitr_cumsum.png'])); end;
end; %if ( doCums )


% Calculate {nw,sw,ne}.{N,P,Si}x.cprof: Profile of nutrient mass flux, per
% profile bin, per meter of reef-line, per ADCP sample period:
%  Cross-shore-velocity x nutrient-conc. x height-per-bin x time-per-sample
%   [m/s] * [kg/m^3] * [m] * [s] = [kg/m]
%
% Then also calculate {nw,sw,ne}.int{N,P,Si}x: Time-cumulative, water-column
% integrated nutrient flux per meter of reef-line over each record [kg/m]

% For comparison, see:
%  2008 Longdill Healy Black - Transient wind-driven coastal upwelling on a shelf with varying width and orientation.pdf
% Upwelled water in New Zealand, T ~ 15oC, NO2+NO3 ~ 80 micrograms / liter
%   80 mu-g/L = 80e-6 kg / m^3 = 8.0e-5 kg / m^3


% Track only movement of water "onshore" (isobath-oriented 'u' < 0)
nw.adcp_onshore = nw.adcp_x;
nw.adcp_onshore.prof(nw.adcp_onshore.prof>0)=0;
nw.adcp_onshore.prof(nw.adcp_onshore.prof<0)=-nw.adcp_onshore.prof(nw.adcp_onshore.prof<0);

sw.adcp_onshore = sw.adcp_x;
sw.adcp_onshore.prof(sw.adcp_onshore.prof>0)=0;
sw.adcp_onshore.prof(sw.adcp_onshore.prof<0)=-sw.adcp_onshore.prof(sw.adcp_onshore.prof<0);

ne.adcp_onshore = ne.adcp_x;
ne.adcp_onshore.prof(ne.adcp_onshore.prof>0)=0;
ne.adcp_onshore.prof(ne.adcp_onshore.prof<0)=-ne.adcp_onshore.prof(ne.adcp_onshore.prof<0);


%%%
% NW/W buoy (11 m isobath) nutrient fluxes

nw.Nx.date = nw.adcp_onshore.date;
nw.Nx.depths = nw.adcp_depths(end:-1:1);

nw.Nx.prof = interp2(nw.N.date,nw.N.depths,nw.N.prof',nw.Nx.date,nw.Nx.depths,'linear',NaN)';

nw.Nx.xprof = nw.adcp_onshore.prof .* umol_to_kgm3(nw.Nx.prof,'N') ...
    .* abs(median(diff(nw.Nx.depths))) .* median(diff(nw.Nx.date)*60*60*24);

nw.Nx.cprof = cumsum(nw.Nx.xprof,'omitnan');

nw.intNx.date = nw.Nx.date; nw.intNx.data = sum(nw.Nx.cprof,2,'omitnan');


%%%
% SW/C buoy (20 m isobath) nutrient fluxes

sw.Nx.date = sw.adcp_onshore.date;
sw.Nx.depths = sw.adcp_depths(end:-1:1);
btmix = find(sw.Nx.depths<sw.N.depths(end));

sw.Nx.prof = interp2(sw.N.date,sw.N.depths,sw.N.prof',sw.Nx.date,sw.Nx.depths,'linear',NaN)';
sw.Nx.prof(:,btmix) = repmat(sw.Nx.prof(:,btmix(1)-1),[1,numel(btmix)]);

sw.Nx.xprof = sw.adcp_onshore.prof .* umol_to_kgm3(sw.Nx.prof,'N') ...
    .* abs(median(diff(sw.Nx.depths))) .* median(diff(sw.Nx.date)*60*60*24);
sw.Nx.xprof(:,btmix) = repmat(sw.Nx.xprof(:,btmix(1)-1),[1,numel(btmix)]);

sw.Nx.cprof = cumsum(sw.Nx.xprof,'omitnan');
sw.Nx.cprof(:,btmix) = repmat(sw.Nx.cprof(:,btmix(1)-1),[1,numel(btmix)]);

sw.intNx.date = sw.Nx.date; sw.intNx.data = sum(sw.Nx.cprof,2,'omitnan');



%%%
% NE/E buoy (50 m isobath) nutrient fluxes

ne.Nx.date = ne.adcp_onshore.date;
ne.Nx.depths = ne.adcp_depths(end:-1:1);
btmix = find(ne.Nx.depths<ne.N.depths(end));

ne.Nx.prof = interp2(ne.N.date,ne.N.depths,ne.N.prof',ne.Nx.date,ne.Nx.depths,'linear',NaN)';
ne.Nx.prof(:,btmix) = repmat(ne.Nx.prof(:,btmix(1)-1),[1,numel(btmix)]);

ne.Nx.xprof = ne.adcp_onshore.prof .* umol_to_kgm3(ne.Nx.prof,'N') ...
    .* abs(median(diff(ne.Nx.depths))) .* median(diff(ne.Nx.date)*60*60*24);
ne.Nx.xprof(:,btmix) = repmat(ne.Nx.xprof(:,btmix(1)-1),[1,numel(btmix)]);

ne.Nx.cprof = cumsum(ne.Nx.xprof,'omitnan');
ne.Nx.cprof(:,btmix) = repmat(ne.Nx.cprof(:,btmix(1)-1),[1,numel(btmix)]);

ne.intNx.date = ne.Nx.date; ne.intNx.data = sum(ne.Nx.cprof,2,'omitnan');



if ( doKgs )
  fmg;
  nkax(1) = subplot(8,1,1);
  surf(nw.Nx.date,nw.Nx.depths,nw.Nx.cprof');
  view(0,270); shading interp; set(gca,'YDir','reverse');
  set(gca,'CLim',[0,0.0675]);
  hold on;
  colorbar;
  ylim([-sfomc.nw_w_btm.depth,0]); 
  %for dp=nw.Nx.depths(:)'; annotline([],dp); end;
  xlim(xl); datetick('x',2,'keeplimits'); %datetick3;
  ylabel('Inshore'); set(gca,'FontSize',fontsz);
  titlename(['Cumulative Nitrate+Nitrite [kg/m/bin] Hovmuellers ',dtstr]);

  nkax(2) = subplot(8,1,2:3);
  %surf(sw.Nx.date,sw.Nx.depths,sw.Nx.prof');
  %surf(sw.Nx.date,sw.Nx.depths,sw.Nx.xprof');
  surf(sw.Nx.date,sw.Nx.depths,sw.Nx.cprof');
  view(0,270); shading interp; set(gca,'YDir','reverse');
  set(gca,'CLim',[0,0.0675]);
  hold on;
  colorbar;
  ylim([-sfomc.sw_buoy.depth,0]); 
  %for dp=sw.Nx.depths(:)'; annotline([],dp); end;
  xlim(xl); datetick('x',2,'keeplimits'); %datetick3;
  ylabel('Mid (SW/C)'); set(gca,'FontSize',fontsz);

  nkax(3) = subplot(8,1,4:8);
  % surf(ne.Nx.date,ne.Nx.depths,ne.Nx.prof');
  % surf(ne.Nx.date,ne.Nx.depths,ne.Nx.xprof');
  surf(ne.Nx.date,ne.Nx.depths,ne.Nx.cprof');
  view(0,270); shading interp; set(gca,'YDir','reverse');
  hold on;
  colorbar;
  ylim([-sfomc.ne_buoy.depth,0]); 
  %for dp=ne.Nx.depths(:)'; annotline([],dp); end;
  xlim(xl); datetick('x',2,'keeplimits'); %datetick3;
  ylabel('Outer (NE/E) u \bullet [NO2+NO3] [kg\bulletm^-^1]'); set(gca,'FontSize',fontsz);

  set(nkax(1:2),'XTickLabel','');
  if doPrint; print('-dpng',fullfile(figspath,[fignm,'_N_cumkg_prof.png'])); end;


  fmg; plot_ts(ne.intNx,sw.intNx,nw.intNx); legend('Outer','Mid','Inshore','Location','East');
  xlim(xl); datetick('x',2,'keeplimits'); %datetick3;
  ylabel('\int^0_-_h (u_{onshore}\bullet[NO2+NO3])dz   [kg\bulletm^-^1]'); set(gca,'FontSize',fontsz);
  titlename('Accumulated Nitrate+Nitrite per meter of reef');
  if doPrint; print('-dpng',fullfile(figspath,[fignm,'_N_cumkg.png'])); end;
end; %if ( doKgs )


[ig,ix11m] = min(abs((-sfomc.ne_buoy.adcp_depths)-11));
[ig,ix30m] = min(abs((-sfomc.ne_buoy.adcp_depths)-30));
[ig,ix45m] = min(abs((-sfomc.ne_buoy.adcp_depths)-45));

sfomc.ne_buoy.adcp_x_11m.date = sfomc.ne_buoy.adcp_x.date;
sfomc.ne_buoy.adcp_x_11m.data = sfomc.ne_buoy.adcp_x.prof(:,ix11m);
sfomc.ne_buoy.adcp_x_30m.date = sfomc.ne_buoy.adcp_x.date;
sfomc.ne_buoy.adcp_x_30m.data = sfomc.ne_buoy.adcp_x.prof(:,ix30m);
sfomc.ne_buoy.adcp_x_45m.date = sfomc.ne_buoy.adcp_x.date;
sfomc.ne_buoy.adcp_x_45m.data = sfomc.ne_buoy.adcp_x.prof(:,ix45m);


% % OOPS! Much of the time, upwelled water seems to be moving OFFshore?!
% fmg; boxplot(sw.adcp_x.prof(sw.Nx.xprof(:,end)>0,:),sw.adcp_depths,'notch','on');
%  titlename('Mid-Station upwelling u\bullet\nablah profile');
% fmg; boxplot(ne.adcp_x.prof(ne.Nx.xprof(:,end)>0,:),ne.adcp_depths,'notch','on');
%  titlename('East-Station upwelling u\bullet\nablah profile');

% Let's focus on the moments right at the start (and end) of an upwelling event...
%upwix = find(abs(sfomc.ne_buoy.at_45m.sbe_seatemp.data-23)<2e-4);
%upwix = find((0 < ne.Nx.prof(:,end)) & (ne.Nx.prof(:,end)<1e-4));
upwix = find(ne.Nx.prof(1:end-1,end)==0 & ne.Nx.prof(2:end,end)>0);
dwnix = find(ne.Nx.prof(1:end-1,end)>0 & ne.Nx.prof(2:end,end)==0);


% % And sure enough, at the START of upwelling, water is moving ONshore...
% fmg; plot(ne.Nx.date(upwix),ne.adcp_x.prof(upwix,[ix45m,ix30m,ix11m]));
%  xlim(datenum(1999,[7,8],[21,1])); datetick3; legend('45m','35m','11m');
%  titlename('Currents at upwelling ONSET');
% % ... while it is actually flowing OFFSHORE only as upwelling relaxes.
% fmg; plot(ne.Nx.date(dwnix),ne.adcp_x.prof(dwnix,[ix45m,ix30m,ix11m]));
%  xlim(datenum(1999,[7,8],[21,1])); datetick3; legend('45m','35m','11m'); ylim([-0.5,0.3]);
%  titlename('Currents at upwelling RELAXATION');
% scatter_fit(ne.Nx.prof(ne.Nx.prof(:,end)>0,end),ne.adcp_x.prof(ne.Nx.prof(:,end)>0,ix45m))

 
fmg;
ax(1)=subplot_tight(4,1,[1:3],'sharex');
 plot_ts(sfomc.ne_buoy.at_45m.sbe_seatemp,sfomc.ne_buoy.at_30m.sbe_seatemp,sfomc.nw_w_btm.adcp_seatemp,...
         ts_nanify_gaps(airtemp,3/24),'k');
 ylim([15,31]); grid on; hold on;
 plot(ne.Nx.date(upwix),repmat(23,size(upwix)),'r^');
 plot(ne.Nx.date(dwnix),repmat(23,size(dwnix)),'bv');
 xlim(xl); datetick('x',2,'keeplimits'); %datetick3;
 ylabel('Temperature [\circC]'); set(gca,'FontSize',fontsz);
titlename(['Upwelling Mechanisms - ',dtstr]);
ax(2)=subplot_tight(4,1,4,'sharex');
 plot_ts(sfomc.ne_buoy.adcp_x_45m,sfomc.ne_buoy.adcp_x_30m,sfomc.nw_w_btm.adcp_x_btm,sfomc.ne_buoy.adcp_x_11m);
 ylim([-0.5,+0.5]); grid on; hold on;
 plot(ne.Nx.date(upwix),repmat(-0.5,size(upwix)),'r^');
 plot(ne.Nx.date(dwnix),repmat(-0.5,size(dwnix)),'bv');
 ylabel('u_x [m/s]'); set(gca,'FontSize',fontsz);
xlim(xl); datetick('x',2,'keeplimits'); %datetick3;
set(ax(1),'XTickLabel','');
if doPrint; print('-dpng',fullfile(figspath,[fignm,'_dynamics_1999.png'])); end;

%{
disp('Hit enter to continue...'); pause;
xlim([datenum(1999,7,15),datenum(1999,7,20)]); datetick('x',2,'keeplimits'); %datetick3;
if doPrint; print('-dpng',fullfile(figspath,[fignm,'_dynamics_1_1999.png'])); end;

disp('Hit enter to continue...'); pause;
xlim([datenum(1999,7,20),datenum(1999,7,25)]); datetick('x',2,'keeplimits'); %datetick3;
if doPrint; print('-dpng',fullfile(figspath,[fignm,'_dynamics_2_1999.png'])); end;

disp('Hit enter to continue...'); pause;
xlim([datenum(1999,7,25),datenum(1999,8,1)]); datetick('x',2,'keeplimits'); %datetick3;
if doPrint; print('-dpng',fullfile(figspath,[fignm,'_dynamics_3_1999.png'])); end;

disp('Hit enter to continue...'); pause;
xlim([datenum(1999,8,1,17,0,0),datenum(1999,8,9,19,0,0)]); datetick('x',2,'keeplimits'); %datetick3;
if doPrint; print('-dpng',fullfile(figspath,[fignm,'_dynamics_4_1999.png'])); end;

disp('Hit enter to continue...'); pause;
xlim([datenum(1999,8,13,0,0,0),datenum(1999,8,21,0,0,0)]); datetick('x',2,'keeplimits'); %datetick3;
if doPrint; print('-dpng',fullfile(figspath,[fignm,'_dynamics_5_1999.png'])); end;

disp('Hit enter to continue...'); pause;
xlim([datenum(1999,8,21,0,0,0),datenum(1999,8,28,0,0,0)]); datetick('x',2,'keeplimits'); %datetick3;
if doPrint; print('-dpng',fullfile(figspath,[fignm,'_dynamics_6_1999.png'])); end;
%}


fmg;
 subplot(2,1,1); plot(get_yearday(sfomc.ne_buoy.at_20m.sbe_seatemp.date),sfomc.ne_buoy.at_20m.sbe_seatemp.data,get_yearday(sfomc.e_buoy.at_20m.sbe_seatemp.date),sfomc.e_buoy.at_20m.sbe_seatemp.data); legend('1999','2000','Location','SouthEast'); ylabel('20 m'); axis([datenum(0,[7,9],15),19,31]); datetick3('x','dd-mmm','keeplimits'); grid on; titlename('Offshore (NE/E) sea temperature'); 
set(gca,'FontSize',fontsz);
 subplot(2,1,2); plot(get_yearday(sfomc.ne_buoy.at_30m.sbe_seatemp.date),sfomc.ne_buoy.at_30m.sbe_seatemp.data,get_yearday(sfomc.e_buoy.at_30m.sbe_seatemp.date),sfomc.e_buoy.at_30m.sbe_seatemp.data); legend('1999','2000','Location','SouthEast'); ylabel('30 m'); axis([datenum(0,[7,9],15),19,31]); datetick3('x','dd-mmm','keeplimits'); grid on;
set(gca,'FontSize',fontsz);
if doPrint; print('-dpng',fullfile(figspath,[fignm,'_offshore_1999_vs_2000.png'])); end;
