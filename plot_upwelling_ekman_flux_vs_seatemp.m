1;
%% SCRIPT plot_upwelling_ekman_flux_vs_seatemp.m
%
% (FORMERLY "like_with_jack_with_sfomc.m")
%
% Display and optionally (DEFAULT: doPrint=False) print a series of plots
% showing sea temperature during upwelling events and their corresponding
% ocean current profiles, coastal Ekman flux, and nutrient flux estimates.
%
% PARAMETERS (defaults):
%   doSFOMC = true; % Plot 1999-2000 temp./current/forcing for SFOMC stations
%   doSFOMCevents = false; % Plot 2000-2011 temp./current/fcg. for SFOMC
%   doSFOMCextras = false; % Plot 1999 closeup, 1999 vs. 2000 intercomparison
%   doSEFCRI = false; % Plot SEFCRI ("Far North" etc.) sea temperature
%   doFACE = false; % Plot FACE ("Deep", "Shallow", etc.) sea temperature & currents
%   doNF09 = false; % Plot Nancy Foster '09 mooring temperatures, etc.
%   doHanS = false; % Plot Hazen & Sawyer (and "FACE AOML") ADCP data
%   doNutPlots = false; % Plot nutrient flux accumulation time series
%   doNutMaps = false; % Plot nutrient flux accumulation regional maps (SLOW!!!)
%   doPause = false; % Pause between sets of plots (useful for user review)
%   doPrint = false; % Print plots
%   figspath = fullfile(get_coral_path,'CRCP','Upwelling','CoRIS');
%   fontsz = 24; % Default font size for plots
%   linewid = 2.5; % Default line width for plots
%   maxgp = 3; % Maximum gap (days) for TS_NANIFY_GAPS
%   upwmT = 24.0; % Near-bottom sea temperature below which upwelling may be occurring
%   upwMr = @(x)(find(ismember(get_month(x.date),[5:10]))):
%           Months during which (detectable) upwelling may be occurring
%
% Last Saved Time-stamp: <Wed 2018-06-20 22:11:16 Eastern Daylight Time gramer>

disp('LOOK AT: 1999 low- then high-, and 2007-7-29 - 8-8'); pause(1);

set_more off;

if ( ~exist('doSFOMC','var') || isempty(doSFOMC) )
  doSFOMC = true;
end;
if ( ~exist('doSFOMCevents','var') || isempty(doSFOMCevents) )
  doSFOMCevents = false;
end;
if ( ~exist('doSFOMCextras','var') || isempty(doSFOMCextras) )
  doSFOMCextras = false;
end;
if ( ~exist('doSEFCRI','var') || isempty(doSEFCRI) )
  doSEFCRI = false;
end;
if ( ~exist('doFACE','var') || isempty(doFACE) )
  doFACE = false;
end;
if ( ~exist('doNF09','var') || isempty(doNF09) )
  doNF09 = false;
end;
if ( ~exist('doHanS','var') || isempty(doHanS) )
  doHanS = false;
end;
if ( ~exist('doNutPlots','var') || isempty(doNutPlots) )
  doNutPlots = false;
end;
if ( ~exist('doNutMaps','var') || isempty(doNutMaps) )
  doNutMaps = false;
end;

if ( ~exist('doPause','var') || isempty(doPause) )
  doPause = false;
end;

if ( ~exist('doPrint','var') || isempty(doPrint) )
  doPrint = false;
end;
if ( ~exist('figspath','var') )
  figspath = fullfile(get_coral_path,'CRCP','Upwelling','CoRIS');
end;
if ( doPrint )
  disp(['Will print figures to ',figspath]);
else
  disp('Figures will not be printed');
end;

if ( ~exist('fontsz','var') || isempty(fontsz) )
  %% For narrow graphs
  fontsz = 24;
  %% For medium graphs
  %fontsz = 20;
  %% For wide graphs
  %fontsz = 16;
end;

if ( ~exist('linewid','var') || isempty(linewid) )
  %linewid = 0.5;
  linewid = 2.5;
end;
% Width of gridlines and boundary for each AXES
if ( ~exist('gridwid','var') || isempty(gridwid) )
  gridwid = 1.5;
end;

maxgap = 3; % Maximum gap (days) for TS_NANIFY_GAPS
GAPFUN = @(x)(ts_nanify_gaps(x,maxgap));

upwmT = 24.0; % Near-bottom sea temperature below which upwelling may be occurring

%upwMr = @ts_jas; % Months during which (detectable) upwelling may be occurring
upwMr = @(x)(find(ismember(get_month(x.date),[5:10]))); % Months during which (detectable) upwelling may be occurring

% Alternate indicator of upwelling - high-frequency near-bottom sea temperature variability

% upwdts = union( ...
% floor(sfomc.nw_w_btm.adcp_seatemp_3_h_hp_9_h_lp.date(abs(sfomc.nw_w_btm.adcp_seatemp_3_h_hp_9_h_lp.data)>0.2)), ...
% floor(sfomc.nw_w_btm.adcp_seatemp_11_h_hp_14_h_lp.date(abs(sfomc.nw_w_btm.adcp_seatemp_11_h_hp_14_h_lp.data)>0.2)) ...
%     );
% upwdts = union( ...
%     floor(sfomc.nw_w_btm.adcp_seatemp_3_h_hp_9_h_lp.date(sfomc.nw_w_btm.adcp_seatemp_3_h_hp_9_h_lp.data<-0.2)), ...
%     floor(sfomc.nw_w_btm.adcp_seatemp_11_h_hp_14_h_lp.date(sfomc.nw_w_btm.adcp_seatemp_11_h_hp_14_h_lp.data<-0.2)) ...
%     );
upwdts = union( ...
    sfomc.nw_w_btm.adcp_seatemp_3_h_hp_9_h_lp.date(sfomc.nw_w_btm.adcp_seatemp_3_h_hp_9_h_lp.data<-0.2), ...
    sfomc.nw_w_btm.adcp_seatemp_11_h_hp_14_h_lp.date(sfomc.nw_w_btm.adcp_seatemp_11_h_hp_14_h_lp.data<-0.2) ...
    );
upwdts = upwdts(ismember(get_month(upwdts),[5:10]));


dT = GAPFUN(ts_d_dt(interp_ts(sfomc.nw_w_btm.adcp_seatemp_hourly),1/24));
dsp = GAPFUN(ts_d_dt(interp_ts(sfomc.nw_w_btm.adcp_speed_btm_hourly),1/24));
dxs = GAPFUN(ts_d_dt(interp_ts(sfomc.nw_w_btm.adcp_x_btm_hourly),1/24));
dls = GAPFUN(ts_d_dt(interp_ts(sfomc.nw_w_btm.adcp_l_btm_hourly),1/24));


if ( doSFOMC )
  legargs = { 'Orientation','horiz','Location','SouthEast' };

  % %xl = sfomc.ne_buoy.seatemp.date([1,end]);
  % %xl = [datenum(1999,6,26),sfomc.ne_buoy.seatemp.date(end)];
  % xl = [datenum(1999,6,15),datenum(2000,4,1)];
  %xl = [datenum(1999,6,25),datenum(2000,4,1)];
  xl = [datenum(1999,6,25),datenum(2000,3,30)];


  % Some day, may be worth analyzing with Lanczos 3- or 6-hlp? Or Butterworth
  % 14-hlp? But that would mask these cool, high-frequency variability events!


  % SHALLOW site currents

  shfh = fmg;
  shax(1) = subplot_tight(7,1,1:3);
  titlename('SFOMC SHALLOW site (NW/W) currents vs. temperatures');
  title('SFOMC temperatures');
  %sh_lhs=plot_ts(fwyf1.ndbc_air_t,'k-',sfomc.ne_buoy.seatemp,sfomc.sw_buoy.seatemp,sfomc.nw_w_btm.adcp_seatemp);
  sh_lhs = plot_ts(GAPFUN(fwyf1.ndbc_air_t),'k-',GAPFUN(sfomc.ne_buoy.at_30m.sbe_seatemp),GAPFUN(sfomc.sw_buoy.at_15m.sbe_seatemp),GAPFUN(sfomc.nw_w_btm.adcp_seatemp));
  plot_ts(GAPFUN(sfomc.e_buoy.at_30m.sbe_seatemp),'specix',2,GAPFUN(sfomc.c_buoy.at_15m.sbe_seatemp),'specix',3,GAPFUN(sfomc.nw_w_btm.adcp_seatemp),'specix',4);
  %ylim([16.5,32.5]); xlim_datetick(xl);
  ylim([16.5,34.5]); xlim_datetick(xl); 
  grid on; grid minor;
  %annotate_ts_event(upwMr,sfomc.ne_buoy.at_45m.sbe_seatemp,upwmT,17.5,'k.');
  sh_lhs(end+1) = annotate_ts_event(upwMr,sfomc.ne_buoy.at_45m.sbe_seatemp,upwmT,-inf);
  %plot(upwdts,repmat(32.4,size(upwdts)),'r.');
  plot(upwdts,repmat(34.4,size(upwdts)),'r.');
  legend(sh_lhs,'Air (FWYF1)','Deep (NE/E 30 m)','Mid (SW/C 15 m)','Shallow (NW 11 m)','(Upwelling)', legargs{:});
  panel_label_subplot(gca,'a','ur','FontSize',fontsz,'BackgroundColor','w');
  set(gca,'FontSize',fontsz); set(gca,'LineWidth',gridwid); set(get(gca,'Child'),'LineWidth',linewid);

  shax(2) = subplot_tight(7,1,4:5); title('SHALLOW site cross-shore current'); %ylabel('Cross-shore');
  lh=plot_ts(GAPFUN(sfomc.nw_w_btm.adcp_x_uwc),GAPFUN(sfomc.nw_w_btm.adcp_x_btm));
  ylim([-0.425,+0.425]); xlim_datetick(xl);
  center_axes; grid on; grid minor;
  %xlim_datetick('Y',xl); %datetick('x',2,'keeplimits'); %datetick3;
  %annotate_ts_event(upwMr,sfomc.ne_buoy.at_45m.sbe_seatemp,upwmT,-0.35,'k.');
  lh(end+1)=annotate_ts_event(upwMr,sfomc.ne_buoy.at_45m.sbe_seatemp,upwmT,-inf);
  legend(lh,'Upper Water Column (NW)','Near Seafloor (NW)','(Upwelling)', legargs{:});
  panel_label_subplot(gca,'b','ur','FontSize',fontsz,'BackgroundColor','w');
  set(gca,'FontSize',fontsz); set(gca,'LineWidth',gridwid); set(get(gca,'Child'),'LineWidth',linewid);

  shax(3) = subplot_tight(7,1,6:7); title('SHALLOW site along-shore current'); %ylabel('Alongshore');
  lh=plot_ts(GAPFUN(sfomc.nw_w_btm.adcp_l_uwc),GAPFUN(sfomc.nw_w_btm.adcp_l_btm));
  ylim([-1.5,+1.5]); xlim_datetick(xl);
  center_axes; grid on; grid minor;
  xlim_datetick('Y',xl); %datetick('x',2,'keeplimits'); %datetick3;
  %annotate_ts_event(upwMr,sfomc.ne_buoy.at_45m.sbe_seatemp,upwmT,-1.25,'k.');
  lh(end+1)=annotate_ts_event(upwMr,sfomc.ne_buoy.at_45m.sbe_seatemp,upwmT,-inf);
  legend(lh,'Upper Water Column (NW)','Near Seafloor (NW)','(Upwelling)', legargs{:});
  panel_label_subplot(gca,'c','ur','FontSize',fontsz,'BackgroundColor','w');
  set(gca,'FontSize',fontsz); set(gca,'LineWidth',gridwid); set(get(gca,'Child'),'LineWidth',linewid);

  %set(shax(1:2),'XTickLabel','');


  % DEEP site currents

  dpfh = fmg;
  dpax(1) = subplot_tight(7,1,1:3);
  titlename('SFOMC DEEP site (NE/E) currents vs. temperatures');
  title('SFOMC temperatures');
  %lh=plot_ts(fwyf1.ndbc_air_t,'k-',sfomc.ne_buoy.seatemp,sfomc.sw_buoy.seatemp,sfomc.nw_w_btm.adcp_seatemp);
  lh=plot_ts(GAPFUN(fwyf1.ndbc_air_t),'k-',GAPFUN(sfomc.ne_buoy.at_30m.sbe_seatemp),GAPFUN(sfomc.sw_buoy.at_15m.sbe_seatemp),GAPFUN(sfomc.nw_w_btm.adcp_seatemp));
  plot_ts(GAPFUN(sfomc.e_buoy.at_30m.sbe_seatemp),'specix',2,GAPFUN(sfomc.c_buoy.at_15m.sbe_seatemp),'specix',3,GAPFUN(sfomc.nw_w_btm.adcp_seatemp),'specix',4);
  %ylim([16.5,32.5]); xlim_datetick(xl); 
  ylim([16.5,34.5]); xlim_datetick(xl); 
  grid on; grid minor;
  %annotate_ts_event(upwMr,sfomc.ne_buoy.at_45m.sbe_seatemp,upwmT,17.5,'k.');
  lh(end+1)=annotate_ts_event(upwMr,sfomc.ne_buoy.at_45m.sbe_seatemp,upwmT,-inf);
  %plot(upwdts,repmat(17.5,size(upwdts)),'r.');
  %plot(upwdts,repmat(32.4,size(upwdts)),'r.');
  plot(upwdts,repmat(34.4,size(upwdts)),'r.');
  legend(lh,'Air (FWYF1)','Deep (NE/E 30 m)','Mid (SW/C 15 m)','Shallow (NW 11 m)','(Upwelling)', legargs{:});
  panel_label_subplot(gca,'a','ur','FontSize',fontsz,'BackgroundColor','w');
  set(gca,'FontSize',fontsz); set(gca,'LineWidth',gridwid); set(get(gca,'Child'),'LineWidth',linewid);

  dpax(2) = subplot_tight(7,1,4:5); title('DEEP site cross-shore current'); %ylabel('Cross-shore');
  lh=plot_ts(GAPFUN(sfomc.ne_buoy.adcp_x_uwc),GAPFUN(sfomc.ne_buoy.adcp_x_btm));
  plot_ts(GAPFUN(sfomc.e_buoy.adcp_x_uwc),GAPFUN(sfomc.e_buoy.adcp_x_btm));
  ylim([-0.425,+0.425]); xlim_datetick(xl);
  center_axes; grid on; grid minor;
  %xlim_datetick('Y',xl); %datetick('x',2,'keeplimits'); %datetick3;
  %annotate_ts_event(upwMr,sfomc.ne_buoy.at_45m.sbe_seatemp,upwmT,-0.35,'k.');
  lh(end+1)=annotate_ts_event(upwMr,sfomc.ne_buoy.at_45m.sbe_seatemp,upwmT,-inf);
  legend(lh,'Upper Water Column (NE/E)','Near Seafloor (NE/E)','(Upwelling)', legargs{:});
  panel_label_subplot(gca,'d','ur','FontSize',fontsz,'BackgroundColor','w');
  set(gca,'FontSize',fontsz); set(gca,'LineWidth',gridwid); set(get(gca,'Child'),'LineWidth',linewid);

  dpax(3) = subplot_tight(7,1,6:7); title('DEEP site along-shore current'); %ylabel('Alongshore');
  lh=plot_ts(GAPFUN(sfomc.ne_buoy.adcp_l_uwc),GAPFUN(sfomc.ne_buoy.adcp_l_btm));
  plot_ts(GAPFUN(sfomc.e_buoy.adcp_l_uwc),GAPFUN(sfomc.e_buoy.adcp_l_btm));
  ylim([-1.5,+1.5]); %xlim_datetick(xl); 
  center_axes; grid on; grid minor;
  xlim_datetick('Y',xl); %datetick('x',2,'keeplimits'); %datetick3;
  %annotate_ts_event(upwMr,sfomc.ne_buoy.at_45m.sbe_seatemp,upwmT,-1.25,'k.');
  lh(end+1)=annotate_ts_event(upwMr,sfomc.ne_buoy.at_45m.sbe_seatemp,upwmT,-inf);
  legend(lh,'Upper Water Column (NE/E)','Near Seafloor (NE/E)','(Upwelling)', legargs{:});
  panel_label_subplot(gca,'e','ur','FontSize',fontsz,'BackgroundColor','w');
  set(gca,'FontSize',fontsz); set(gca,'LineWidth',gridwid); set(get(gca,'Child'),'LineWidth',linewid);

  %set(dpax(1:2),'XTickLabel','');


  % Ekman flux

  ekfh = fmg;
  ekax(1) = subplot_tight(7,1,1:4); %ylabel('Cross-shore');
  titlename('Ekman volumetric flux [m^2s^-^1] near SFOMC section - cross-shore');
  title('Cross-shore volume flux');
  lh=plot_ts(GAPFUN(canf1.ndbc_ekman_flux_x_volume),GAPFUN(lkwf1.ndbc_ekman_flux_x_volume),GAPFUN(pvgf1.ndbc_ekman_flux_x_volume),GAPFUN(rsm.rsmas_ekman_flux_x_volume),GAPFUN(fwyf1.ndbc_ekman_flux_x_volume));
  %xlim(xl);
  xlim_datetick(xl);
  %%ylim([-4.9,+4.9]);
  %ylim([-25.9,+25.9]); % Accommodate tropical weather systems!
  ylim([-9.9,+9.9]); % Partially accommodate tropical weather systems
  center_axes; grid on; grid minor;
  %annotate_ts_event(upwMr,sfomc.ne_buoy.at_45m.sbe_seatemp,upwmT,-4.5,'k.');
  annotate_ts_event(upwMr,sfomc.ne_buoy.at_45m.sbe_seatemp,upwmT,-inf);
  legh=legend('Cape Can. (''88-''16)','Lake Worth (''01-''16)','Port Ev. (''09-''17)','RSMAS Roof (''03-''11)','Fowey (''91-''16)','(Upw.)', legargs{:}); set(legh,'FontSize',19);
  panel_label_subplot(gca,'f','ur','FontSize',fontsz,'BackgroundColor','w');
  set(gca,'FontSize',fontsz); set(gca,'LineWidth',gridwid); set(get(gca,'Child'),'LineWidth',linewid);

  ekax(2) = subplot_tight(7,1,5:6); %ylabel('Alongshore');
  title('Along-shore volume flux');
  lh=plot_ts(GAPFUN(canf1.ndbc_ekman_flux_l_volume),GAPFUN(lkwf1.ndbc_ekman_flux_l_volume),GAPFUN(pvgf1.ndbc_ekman_flux_l_volume),GAPFUN(rsm.rsmas_ekman_flux_l_volume),GAPFUN(fwyf1.ndbc_ekman_flux_l_volume));
  %xlim(xl);
  xlim_datetick(xl);
  %%ylim([-4.9,+4.9]);
  %ylim([-25.9,+25.9]); % Accommodate tropical weather systems!
  ylim([-9.9,+9.9]); % Partially accommodate tropical weather systems
  center_axes; grid on; grid minor;
  %annotate_ts_event(upwMr,sfomc.ne_buoy.at_45m.sbe_seatemp,upwmT,-4.5,'k.');
  annotate_ts_event(upwMr,sfomc.ne_buoy.at_45m.sbe_seatemp,upwmT,-inf);
  panel_label_subplot(gca,'g','ur','FontSize',fontsz,'BackgroundColor','w');
  set(gca,'FontSize',fontsz); set(gca,'LineWidth',gridwid); set(get(gca,'Child'),'LineWidth',linewid);

  ekax(3) = subplot_tight(7,1,7); %ylabel('Barom.');
  title('Barometric pressure');
  lh=plot_ts(GAPFUN(canf1.ndbc_barom),GAPFUN(lkwf1.ndbc_barom),GAPFUN(pvgf1.ndbc_barom),GAPFUN(rsm.rsmas_barom),GAPFUN(fwyf1.ndbc_barom));
  %xlim(xl);
  ylim([980,1036]);
  grid on; grid minor;
  xlim_datetick('Y',xl); %datetick('x',2,'keeplimits'); %datetick3;
  %annotate_ts_event(upwMr,sfomc.ne_buoy.at_45m.sbe_seatemp,upwmT,985,'k.');
  annotate_ts_event(upwMr,sfomc.ne_buoy.at_45m.sbe_seatemp,upwmT,-inf);
  panel_label_subplot(gca,'h','ur','FontSize',fontsz,'BackgroundColor','w');
  set(gca,'FontSize',fontsz); set(gca,'LineWidth',gridwid); set(get(gca,'Child'),'LineWidth',linewid);
  %set(ekax(1:2),'XTickLabel','');

  % ekax(3) = subplot_tight(7,1,7); ylabel('Barom.');
  % axen = plotyy_ts(GAPFUN(canf1.ndbc_barom),GAPFUN(canf1.ndbc_sigwavehgt));
  % lh=plot_ts(axen(1),GAPFUN(lkwf1.ndbc_barom),GAPFUN(pvgf1.ndbc_barom),GAPFUN(rsm.rsmas_barom),GAPFUN(fwyf1.ndbc_barom));
  % xlim(axen,xl);
  % ylim(axen(1),[980,1040]);
  % grid on; grid minor;
  % xlim_datetick('Y',xl); %datetick('x',2,'keeplimits'); %datetick3;
  % set(axen,'FontSize',fontsz); set(get(axen(1),'Child'),'LineWidth',linewid); set(get(axen(2),'Child'),'LineWidth',linewid);
  % axes(axen(2));



  if 0;
    % Wind stress

    wsfh = fmg;
    wsax(1) = subplot_tight(7,1,1:4); ylabel('Alongshore');
    titlename('Wind stress [Nm^-^2] near SFOMC section');
    lh=plot_ts(GAPFUN(canf1.ndbc_bulk_windstress_ls),GAPFUN(lkwf1.ndbc_bulk_windstress_ls),GAPFUN(pvgf1.ndbc_bulk_windstress_ls),GAPFUN(rsm.rsmas_bulk_windstress_ls),GAPFUN(fwyf1.ndbc_bulk_windstress_ls));
    xlim_datetick(xl);
    ylim([-1.5,+1.5]);
    center_axes; grid on; grid minor;
    %annotate_ts_event(upwMr,sfomc.ne_buoy.at_45m.sbe_seatemp,upwmT,-1.25,'k.');
    annotate_ts_event(upwMr,sfomc.ne_buoy.at_45m.sbe_seatemp,upwmT,-inf);
    legh=legend('Cape Can. (''88-''16)','Lake Worth (''01-''16)','Port Ev. (''09-''17)','RSMAS Roof (''03-''11)','Fowey (''91-''16)','(Upw.)', legargs{:}); set(legh,'FontSize',19);
    set(gca,'FontSize',fontsz); set(gca,'LineWidth',gridwid); set(get(gca,'Child'),'LineWidth',linewid);

    wsax(2) = subplot_tight(7,1,5:6); ylabel('Cross-shore');
    lh=plot_ts(GAPFUN(canf1.ndbc_bulk_windstress_xs),GAPFUN(lkwf1.ndbc_bulk_windstress_xs),GAPFUN(pvgf1.ndbc_bulk_windstress_xs),GAPFUN(rsm.rsmas_bulk_windstress_xs),GAPFUN(fwyf1.ndbc_bulk_windstress_xs));
    xlim_datetick(xl);
    ylim([-1.5,+1.5]);
    center_axes; grid on; grid minor;
    %annotate_ts_event(upwMr,sfomc.ne_buoy.at_45m.sbe_seatemp,upwmT,-1.25,'k.');
    annotate_ts_event(upwMr,sfomc.ne_buoy.at_45m.sbe_seatemp,upwmT,-inf);
    set(gca,'FontSize',fontsz); set(gca,'LineWidth',gridwid); set(get(gca,'Child'),'LineWidth',linewid);

    % wsax(3) = subplot_tight(7,1,7); ylabel('H_s,Pcp.');
    % lh=plot_ts(GAPFUN(canf1.ndbc_sigwavehgt),GAPFUN(cnnf1.ndbc_sigwavehgt),GAPFUN(ftpf1.ndbc_sigwavehgt),GAPFUN(rsm.rsmas_precip),GAPFUN(fdepk.fdep_rainfall));
    wsax(3) = subplot_tight(7,1,7); ylabel('Waves');
    lh=plot_ts(GAPFUN(canf1.ndbc_sigwavehgt),GAPFUN(cnnf1.ndbc_sigwavehgt),GAPFUN(ftpf1.ndbc_sigwavehgt));
    %xlim(xl);
    ylim([0,10]);
    grid on; grid minor;
    xlim_datetick('Y',xl); %datetick('x',2,'keeplimits'); %datetick3;
    %annotate_ts_event(upwMr,sfomc.ne_buoy.at_45m.sbe_seatemp,upwmT,0.5,'k.');
    annotate_ts_event(upwMr,sfomc.ne_buoy.at_45m.sbe_seatemp,upwmT,-inf);
    set(gca,'FontSize',fontsz); set(gca,'LineWidth',gridwid); set(get(gca,'Child'),'LineWidth',linewid);
    %set(wsax(1:2),'XTickLabel','');


    % Shallow-site variability metrics

    vmfh = fmg;
    vmax(1) = subplot_tight(7,1,1:4); ylabel('\Delta_t Temperature [\circC]');
    titlename('Sea temperature variability [T\bullethr^-^1] at NW/W Bottom');
    lh=plot_ts(GAPFUN(dT));
    xlim_datetick(xl);
    ylim([-2.0,+2.0]);
    center_axes; grid on; grid minor;
    set(gca,'FontSize',fontsz); set(gca,'LineWidth',gridwid); set(get(gca,'Child'),'LineWidth',linewid);
    %annotate_ts_event(upwMr,sfomc.ne_buoy.at_45m.sbe_seatemp,upwmT,-1.75,'k.');
    annotate_ts_event(upwMr,sfomc.ne_buoy.at_45m.sbe_seatemp,upwmT,-inf);
    
    vmax(2) = subplot_tight(7,1,5:6); ylabel('\Delta_tu near btm.');
    lh=plot_ts(GAPFUN(dls),GAPFUN(dxs));
    xlim_datetick(xl);
    ylim([-1.5,+1.5]);
    center_axes; grid on; grid minor;
    %annotate_ts_event(upwMr,sfomc.ne_buoy.at_45m.sbe_seatemp,upwmT,-1.25,'k.');
    annotate_ts_event(upwMr,sfomc.ne_buoy.at_45m.sbe_seatemp,upwmT,-inf);
    legend('Along','Cross','(Upwelling)');
    set(gca,'FontSize',fontsz); set(gca,'LineWidth',gridwid); set(get(gca,'Child'),'LineWidth',linewid);
    
    vmax(3) = subplot_tight(7,1,7); ylabel('\Delta_t Spd.');
    lh=plot_ts(GAPFUN(dsp));
    %xlim(xl);
    ylim([0,+2.5]);
    grid on; grid minor;
    xlim_datetick('Y',xl); %datetick('x',2,'keeplimits'); %datetick3;
    %annotate_ts_event(upwMr,sfomc.ne_buoy.at_45m.sbe_seatemp,upwmT,0.25,'k.');
    annotate_ts_event(upwMr,sfomc.ne_buoy.at_45m.sbe_seatemp,upwmT,-inf);
    set(gca,'FontSize',fontsz); set(gca,'LineWidth',gridwid); set(get(gca,'Child'),'LineWidth',linewid);
    %set(wsax(1:2),'XTickLabel','');
  end;


  prfh = fmg;
  prax(1) = subplot_tight(7,1,1:3); %ylabel('Sea temp.');
  titlename('Sea temperature variability by period at NW/W Bottom');
  lh=plot_ts(GAPFUN(sfomc.nw_w_btm.adcp_seatemp_3_h_hp_9_h_lp),GAPFUN(sfomc.nw_w_btm.adcp_seatemp_11_h_hp_14_h_lp),GAPFUN(sfomc.nw_w_btm.adcp_seatemp_17_h_hp_28_h_lp),GAPFUN(sfomc.nw_w_btm.adcp_seatemp_48_h_hp_240_h_lp));
  ylim([-1.0,+1.0]);
  xlim_datetick(xl);
  center_axes; grid on; grid minor;
  %annotate_ts_event(upwMr,sfomc.ne_buoy.at_45m.sbe_seatemp,upwmT,-0.75,'k.');
  annotate_ts_event(upwMr,sfomc.ne_buoy.at_45m.sbe_seatemp,upwmT,-inf);
  legend('3-9h','11-14h','17-28h','2-10d','(Upw.)', legargs{:});
  panel_label_subplot(gca,'a','ur','FontSize',fontsz,'BackgroundColor','w');
  set(gca,'FontSize',fontsz); set(gca,'LineWidth',gridwid); set(get(gca,'Child'),'LineWidth',linewid);

  prax(2) = subplot_tight(7,1,4:5); %ylabel('Btm. u_x');
  title('Cross-shore current variability by period');
  lh=plot_ts(GAPFUN(sfomc.nw_w_btm.adcp_x_btm_3_h_hp_9_h_lp),GAPFUN(sfomc.nw_w_btm.adcp_x_btm_11_h_hp_14_h_lp),GAPFUN(sfomc.nw_w_btm.adcp_x_btm_17_h_hp_28_h_lp),GAPFUN(sfomc.nw_w_btm.adcp_x_btm_48_h_hp_240_h_lp));
  ylim([-0.1,+0.1]);
  xlim_datetick(xl);
  center_axes; grid on; grid minor;
  %annotate_ts_event(upwMr,sfomc.ne_buoy.at_45m.sbe_seatemp,upwmT,-0.075,'k.');
  annotate_ts_event(upwMr,sfomc.ne_buoy.at_45m.sbe_seatemp,upwmT,-inf);
  %legend('3-9h','11-14h','17-28h','2-10d','(Upw.)', 'Orientation','horiz');
  panel_label_subplot(gca,'b','ur','FontSize',fontsz,'BackgroundColor','w');
  set(gca,'FontSize',fontsz); set(gca,'LineWidth',gridwid); set(get(gca,'Child'),'LineWidth',linewid);

  prax(3) = subplot_tight(7,1,6:7); %ylabel('Btm. u_l');
  title('Along-shore current variability by period');
  lh=plot_ts(GAPFUN(sfomc.nw_w_btm.adcp_l_btm_3_h_hp_9_h_lp),GAPFUN(sfomc.nw_w_btm.adcp_l_btm_11_h_hp_14_h_lp),GAPFUN(sfomc.nw_w_btm.adcp_l_btm_17_h_hp_28_h_lp),GAPFUN(sfomc.nw_w_btm.adcp_l_btm_48_h_hp_240_h_lp));
  %xlim(datenum([2013,2015],[11,11],[12,12])); datetick3;
  %xlim(xl); datetick3;
  xlim_datetick('Y',xl);
  ylim([-0.425,+0.425]);
  center_axes; grid on; grid minor;
  %annotate_ts_event(upwMr,sfomc.ne_buoy.at_45m.sbe_seatemp,upwmT,-0.35,'k.');
  annotate_ts_event(upwMr,sfomc.ne_buoy.at_45m.sbe_seatemp,upwmT,-inf);
  %legend('3-9h','11-14h','17-28h','2-10d','(Upw.)', 'Orientation','horiz');
  panel_label_subplot(gca,'c','ur','FontSize',fontsz,'BackgroundColor','w');
  set(gca,'FontSize',fontsz); set(gca,'LineWidth',gridwid); set(get(gca,'Child'),'LineWidth',linewid);

  if 0;
    prfh = fmg;
    prax(1) = subplot_tight(7,1,1:3); %ylabel('Sea temp.');
    titlename('Sea temperature variability by period at NE Buoy');
    lh=plot_ts(GAPFUN(sfomc.ne_buoy.seatemp_3_h_hp_9_h_lp),GAPFUN(sfomc.ne_buoy.seatemp_11_h_hp_14_h_lp),GAPFUN(sfomc.ne_buoy.seatemp_17_h_hp_28_h_lp),GAPFUN(sfomc.ne_buoy.seatemp_48_h_hp_240_h_lp));
    xlim_datetick(xl);
    %annotate_ts_event(upwMr,sfomc.ne_buoy.at_45m.sbe_seatemp,upwmT,-0.75,'k.');
    annotate_ts_event(upwMr,sfomc.ne_buoy.at_45m.sbe_seatemp,upwmT,-inf);
    legend('3-9h','11-14h','17-28h','2-10d','(Upw.)', legargs{:});
    set(gca,'FontSize',fontsz); set(gca,'LineWidth',gridwid); set(get(gca,'Child'),'LineWidth',linewid);

    prax(2) = subplot_tight(7,1,4:5); %ylabel('Btm. u_x');
    title('Cross-shore current variability by period');
    lh=plot_ts(GAPFUN(sfomc.ne_buoy.adcp_x_btm_3_h_hp_9_h_lp),GAPFUN(sfomc.ne_buoy.adcp_x_btm_11_h_hp_14_h_lp),GAPFUN(sfomc.ne_buoy.adcp_x_btm_17_h_hp_28_h_lp),GAPFUN(sfomc.ne_buoy.adcp_x_btm_48_h_hp_240_h_lp));
    center_axes; grid on; grid minor;
    xlim_datetick(xl);
    %annotate_ts_event(upwMr,sfomc.ne_buoy.at_45m.sbe_seatemp,upwmT,-0.075,'k.');
    annotate_ts_event(upwMr,sfomc.ne_buoy.at_45m.sbe_seatemp,upwmT,-inf);
    %legend('3-9h','11-14h','17-28h','2-10d','(Upw.)', 'Orientation','horiz');
    set(gca,'FontSize',fontsz); set(gca,'LineWidth',gridwid); set(get(gca,'Child'),'LineWidth',linewid);

    prax(3) = subplot_tight(7,1,6:7); %ylabel('Btm. u_l');
    title('Along-shore current variability by period');
    lh=plot_ts(GAPFUN(sfomc.ne_buoy.adcp_l_btm_3_h_hp_9_h_lp),GAPFUN(sfomc.ne_buoy.adcp_l_btm_11_h_hp_14_h_lp),GAPFUN(sfomc.ne_buoy.adcp_l_btm_17_h_hp_28_h_lp),GAPFUN(sfomc.ne_buoy.adcp_l_btm_48_h_hp_240_h_lp));
    %xlim(datenum([2013,2015],[11,11],[12,12])); datetick3;
    %xlim(xl); datetick3;
    center_axes; grid on; grid minor;
    xlim_datetick('Y',xl);
    %annotate_ts_event(upwMr,sfomc.ne_buoy.at_45m.sbe_seatemp,upwmT,-0.35,'k.');
    annotate_ts_event(upwMr,sfomc.ne_buoy.at_45m.sbe_seatemp,upwmT,-inf);
    %legend('3-9h','11-14h','17-28h','2-10d','(Upw.)', 'Orientation','horiz');
    set(gca,'FontSize',fontsz); set(gca,'LineWidth',gridwid); set(get(gca,'Child'),'LineWidth',linewid);
  end;

  annot_dts = [datenum(1999,[7,8],[8,19]) ; datenum(1999,[8,10],[20,8]) ; datenum([1999,2000],[10,2],[11,29])];
  annot_h = textbox_range([shax,dpax,ekax,prax],annot_dts,{'Upwelling','Mixed','Atmospheric Cooling'},'uc','EdgeColor','k','LineWidth',2);

  if ( doPrint )
    axes(shax(end));
    print('-dpng',fullfile(figspath,'upwelling-sfomc-1999-2000-shallow-whole.png'));
    axes(dpax(end));
    print('-dpng',fullfile(figspath,'upwelling-sfomc-1999-2000-deep-whole.png'));
    axes(ekax(end));
    print('-dpng',fullfile(figspath,'upwelling-sfomc-1999-2000-ekman-whole.png'));
    % axes(wsax(end));
    % print('-dpng',fullfile(figspath,'upwelling-sfomc-1999-2000-stress-whole.png'));
    % axes(vmax(end));
    % print('-dpng',fullfile(figspath,'upwelling-sfomc-1999-2000-variab-whole.png'));
    axes(prax(end));
    print('-dpng',fullfile(figspath,'upwelling-sfomc-1999-2000-partn-whole.png'));
  end;
  disp('Whole records');
  %keyboard;

end; %if ( doSFOMC )


if ( ~doSFOMCevents )
  if ( doSFOMC )
    disp('Not plotting individual events (doSFOMCevents FALSE)');
  end;

elseif ( ~doSFOMC )
  disp('OOPS! doSFOMCevents=True requires doSFOMCevents=True! IGNORED...');

else
  if doPause; pause; else; pause(0.5); end;

  if 0;
    disp('Broader focus on periods of high-frequency variability');
    axes(shax(end)); xlim_datetick(shax,datenum(1999,[7,8],[15,29])); datetick(shax(3),'x',2,'keeplimits');
    if doPrint; print('-dpng',fullfile(figspath,'upwelling-sfomc-1999-shallow-Jul-Aug.png')); end;
    axes(dpax(end)); xlim_datetick(dpax,datenum(1999,[7,8],[15,29])); datetick(dpax(3),'x',2,'keeplimits');
    if doPrint; print('-dpng',fullfile(figspath,'upwelling-sfomc-1999-deep-Jul-Aug.png')); end;
    axes(ekax(end)); xlim_datetick(ekax,datenum(1999,[7,8],[15,29])); datetick(ekax(3),'x',2,'keeplimits');
    if doPrint; print('-dpng',fullfile(figspath,'upwelling-sfomc-1999-ekman-Jul-Aug.png')); end;
    % axes(wsax(end)); xlim_datetick(wsax,datenum(1999,[7,8],[15,29])); datetick(wsax(3),'x',2,'keeplimits');
    % if doPrint; print('-dpng',fullfile(figspath,'upwelling-sfomc-1999-stress-Jul-Aug.png')); end;
    % axes(vmax(end)); xlim_datetick(vmax,datenum(1999,[7,8],[15,29])); datetick(vmax(3),'x',2,'keeplimits');
    % if doPrint; print('-dpng',fullfile(figspath,'upwelling-sfomc-1999-variab-Jul-Aug.png')); end;
    axes(prax(end)); xlim_datetick(prax,datenum(1999,[7,8],[15,29])); datetick(prax(3),'x',2,'keeplimits');
    if doPrint; print('-dpng',fullfile(figspath,'upwelling-sfomc-1999-partn-Jul-Aug.png')); end;
    if doPause; pause; else; pause(0.5); end;

    disp('Close focus on period of low-frequency (initial) variability');
    axes(shax(end)); xlim_datetick(shax,datenum(1999,7,[22,28])); datetick(shax(3),'x',2,'keeplimits');
    if doPrint; print('-dpng',fullfile(figspath,'upwelling-sfomc-1999-shallow-Jul-lowfreq.png')); end;
    axes(dpax(end)); xlim_datetick(dpax,datenum(1999,7,[22,28])); datetick(dpax(3),'x',2,'keeplimits');
    if doPrint; print('-dpng',fullfile(figspath,'upwelling-sfomc-1999-deep-Jul-lowfreq.png')); end;
    axes(ekax(end)); xlim_datetick(ekax,datenum(1999,7,[22,28])); datetick(ekax(3),'x',2,'keeplimits');
    if doPrint; print('-dpng',fullfile(figspath,'upwelling-sfomc-1999-ekman-Jul-lowfreq.png')); end;
    % axes(wsax(end)); xlim_datetick(wsax,datenum(1999,7,[22,28])); datetick(wsax(3),'x',2,'keeplimits');
    % if doPrint; print('-dpng',fullfile(figspath,'upwelling-sfomc-1999-stress-Jul-lowfreq.png')); end;
    % axes(vmax(end)); xlim_datetick(vmax,datenum(1999,7,[22,28])); datetick(vmax(3),'x',2,'keeplimits');
    % if doPrint; print('-dpng',fullfile(figspath,'upwelling-sfomc-1999-variab-Jul-lowfreq.png')); end;
    axes(prax(end)); xlim_datetick(prax,datenum(1999,7,[22,28])); datetick(prax(3),'x',2,'keeplimits');
    if doPrint; print('-dpng',fullfile(figspath,'upwelling-sfomc-1999-partn-Jul-lowfreq.png')); end;
    if doPause; pause; else; pause(0.5); end;

    disp('CLOSE focus on period of high-frequency variability');
    axes(shax(end)); xlim_datetick(shax,datenum(1999,7,[28,32])); datetick(shax(3),'x',2,'keeplimits');
    if doPrint; print('-dpng',fullfile(figspath,'upwelling-sfomc-1999-shallow-Jul-highfreq.png')); end;
    axes(dpax(end)); xlim_datetick(dpax,datenum(1999,7,[28,32])); datetick(dpax(3),'x',2,'keeplimits');
    if doPrint; print('-dpng',fullfile(figspath,'upwelling-sfomc-1999-deep-Jul-highfreq.png')); end;
    axes(ekax(end)); xlim_datetick(ekax,datenum(1999,7,[28,32])); datetick(ekax(3),'x',2,'keeplimits');
    if doPrint; print('-dpng',fullfile(figspath,'upwelling-sfomc-1999-ekman-Jul-highfreq.png')); end;
    % axes(wsax(end)); xlim_datetick(wsax,datenum(1999,7,[28,32])); datetick(wsax(3),'x',2,'keeplimits');
    % if doPrint; print('-dpng',fullfile(figspath,'upwelling-sfomc-1999-stress-Jul-highfreq.png')); end;
    % axes(vmax(end)); xlim_datetick(vmax,datenum(1999,7,[28,32])); datetick(vmax(3),'x',2,'keeplimits');
    % if doPrint; print('-dpng',fullfile(figspath,'upwelling-sfomc-1999-variab-Jul-highfreq.png')); end;
    axes(prax(end)); xlim_datetick(prax,datenum(1999,7,[28,32])); datetick(prax(3),'x',2,'keeplimits');
    if doPrint; print('-dpng',fullfile(figspath,'upwelling-sfomc-1999-partn-Jul-highfreq.png')); end;
    if doPause; pause; else; pause(0.5); end;

    % disp('Whole period when all current meters were operating');
    % axes(shax(end)); xlim_datetick(shax,datenum(1999,[9,11],[18,12])); datetick(shax(3),'x',2,'keeplimits');
    % axes(dpax(end)); xlim_datetick(dpax,datenum(1999,[9,11],[18,12])); datetick(dpax(3),'x',2,'keeplimits');
    % axes(ekax(end)); xlim_datetick(ekax,datenum(1999,[9,11],[18,12])); datetick(ekax(3),'x',2,'keeplimits');
    % % axes(wsax(end)); xlim_datetick(wsax,datenum(1999,[9,11],[18,12])); datetick(wsax(3),'x',2,'keeplimits');
    % % axes(vmax(end)); xlim_datetick(vmax,datenum(1999,[9,11],[18,12])); datetick(vmax(3),'x',2,'keeplimits');
    % axes(prax(end)); xlim_datetick(prax,datenum(1999,[9,11],[18,12])); datetick(prax(3),'x',2,'keeplimits');
    % if doPause; pause; else; pause(0.5); end;
    
    disp('Upwelling event when three current meters were operating');
    axes(shax(end)); xlim_datetick(shax,datenum(1999,7,[22,28])); datetick(shax(3),'x',2,'keeplimits');
    axes(dpax(end)); xlim_datetick(dpax,datenum(1999,7,[22,28])); datetick(dpax(3),'x',2,'keeplimits');
    axes(ekax(end)); xlim_datetick(ekax,datenum(1999,7,[22,28])); datetick(ekax(3),'x',2,'keeplimits');
    % axes(wsax(end)); xlim_datetick(wsax,datenum(1999,7,[22,28])); datetick(wsax(3),'x',2,'keeplimits');
    % axes(vmax(end)); xlim_datetick(vmax,datenum(1999,7,[22,28])); datetick(vmax(3),'x',2,'keeplimits');
    axes(prax(end)); xlim_datetick(prax,datenum(1999,7,[22,28])); datetick(prax(3),'x',2,'keeplimits');
    if doPause; pause; else; pause(0.5); end;
  end;

  % %% 1999, 2000, 2003, 2004, 2005?, 2006, 2007, 2009?, 2010?, 2011
  % modys = datenum(0,[7,9],[5,5]);
  % for yr=[1999, 2000, 2003, 2004, 2005, 2006, 2007, 2011];
  modys = datenum(0,[6,9],[5,5]);
  for yr=[1999, 2000, 2003, 2004, 2005, 2006, 2007, 2009, 2010, 2011];
    yrstr = num2str(yr);
    disp(['Focus on high-frequency variability in 5-min sea temperature in *',yrstr,'*']);

    xl = datenum(yr,1,0)+modys;

    ylim(shax(1),[21.5,33.5]);
    axes(shax(end)); xlim_datetick(shax,xl); datetick(shax(3),'x',2,'keeplimits');
    if (yr>2000); legend(sh_lhs([1,end-1,end]),{'Air (FWYF1)','Shallow (NW 11 m)','(Upwelling)'}, legargs{:}); end;
    if ( yr == 1999 || yr == 2000 )
      ylim(dpax(1),[21.5,33.5]);
      axes(dpax(end)); xlim_datetick(dpax,xl); datetick(dpax(3),'x',2,'keeplimits');
    end;
    axes(ekax(end)); xlim_datetick(ekax,xl); datetick(ekax(3),'x',2,'keeplimits');
    if ( yr > 1999 )
      panel_label_subplot(ekax(1),'d','ur','FontSize',fontsz,'BackgroundColor','w');
      panel_label_subplot(ekax(2),'e','ur','FontSize',fontsz,'BackgroundColor','w');
      panel_label_subplot(ekax(3),'f','ur','FontSize',fontsz,'BackgroundColor','w');
    end;
    % axes(wsax(end)); xlim_datetick(wsax,xl); datetick(wsax(3),'x',2,'keeplimits');
    % axes(vmax(end)); xlim_datetick(vmax,xl); datetick(vmax(3),'x',2,'keeplimits');
    axes(prax(end)); xlim_datetick(prax,xl); datetick(prax(3),'x',2,'keeplimits');
    if ( yr > 1999 )
      panel_label_subplot(prax(1),'g','ur','FontSize',fontsz,'BackgroundColor','w');
      panel_label_subplot(prax(2),'h','ur','FontSize',fontsz,'BackgroundColor','w');
      panel_label_subplot(prax(3),'i','ur','FontSize',fontsz,'BackgroundColor','w');
    end;
    % if ( yr == 1999 )
    %   axes(prax(end)); xlim_datetick(prax,xl); datetick(prax(3),'x',2,'keeplimits');
    % end;

    try, delete(annot_h); catch, end;
    if ( yr == 1999 || yr == 2000 )
      annot_h = textbox_range([shax,dpax,ekax,prax],sfomc_upwelling_dates(xl),'ROMAN','um','EdgeColor','k','LineWidth',2);
    else
      annot_h = textbox_range([shax,ekax,prax],sfomc_upwelling_dates(xl),'ROMAN','um','EdgeColor','k','LineWidth',2);
    end;

    if ( doPrint )
      axes(shax(end)); print('-dpng',fullfile(figspath,['upwelling-sfomc-',yrstr,'-shallow-wide.png']));
      if ( yr == 1999 || yr == 2000 )
        axes(dpax(end)); print('-dpng',fullfile(figspath,['upwelling-sfomc-',yrstr,'-deep-wide.png']));
      end;
      axes(ekax(end)); print('-dpng',fullfile(figspath,['upwelling-sfomc-',yrstr,'-ekman-wide.png']));
      %axes(wsax(end)); print('-dpng',fullfile(figspath,['upwelling-sfomc-',yrstr,'-stress-wide.png']));
      %axes(vmax(end)); print('-dpng',fullfile(figspath,['upwelling-sfomc-',yrstr,'-variab-wide.png']));
      axes(prax(end)); print('-dpng',fullfile(figspath,['upwelling-sfomc-',yrstr,'-partn-wide.png']));
    end;

    if doPause; pause; else; pause(0.5); end;

  end;
  %xlim([shax,dpax,ekax,wsax],datenum(1999,1,0)+modys); jumpax([shax,dpax,ekax,wsax],nan,365);
  %keyboard;

  disp('Done with SFOMC');
  if doPause; pause; else; pause(0.5); end;



  % 1999 TRANSIENTS
  % Narrow focus on periods of high-frequency variability
  xl = datenum(1999,7,[9,29]);

  fmg;
  subplot_tight(7,1,1:3);
  titlename('Section of sea temperatures and currents (1999 TRANSIENTS)');
  title('Temperatures (1999 TRANSIENTS)');
  lh=plot_ts(GAPFUN(fwyf1.ndbc_air_t),'k-',GAPFUN(sfomc.se_buoy.seatemp),GAPFUN(sfomc.ne_buoy.seatemp),GAPFUN(sfomc.sw_buoy.seatemp),GAPFUN(sfomc.nw_w_btm.adcp_seatemp));
  ylim([16.5,32.5]);
  grid on; grid minor;
  xlim_datetick(xl);
  legend(lh,'Air (FWYF1)','Offshore (SE)','Deep (NE)','Mid (SW)','Shallow (NW)', legargs{:});
  set(gca,'FontSize',fontsz); set(gca,'LineWidth',gridwid); set(get(gca,'Child'),'LineWidth',linewid);

  subplot_tight(7,1,4:5); %ylabel('Cross-shore');
  title('Cross-shore currents');
  lh=plot_ts(GAPFUN(sfomc.ne_buoy.adcp_u),GAPFUN(sfomc.sw_buoy.adcp_u),GAPFUN(sfomc.nw_w_btm.adcp_u));
  ylim([-0.425,+0.425]);
  grid on; grid minor;
  xlim_datetick(xl);
  set(gca,'FontSize',fontsz); set(gca,'LineWidth',gridwid); set(get(gca,'Child'),'LineWidth',linewid);

  subplot_tight(7,1,6:7); %ylabel('Alongshore');
  title('Along-shore currents');
  lh=plot_ts(GAPFUN(sfomc.ne_buoy.adcp_v),GAPFUN(sfomc.sw_buoy.adcp_v),GAPFUN(sfomc.nw_w_btm.adcp_v));
  ylim([-1.5,+1.5]);
  grid on; grid minor;
  xlim_datetick('Y',xl);
  legend(lh,'Deep (NE)','Mid (SW)','Shallow (NW)', legargs{:});
  set(gca,'FontSize',fontsz); set(gca,'LineWidth',gridwid); set(get(gca,'Child'),'LineWidth',linewid);

  annot_h = textbox_range(gcf,sfomc_upwelling_dates(xl),'ROMAN');

  if doPrint; print('-dpng',fullfile(figspath,'upwelling-sfomc-1999-cross-vs-alongshore-Jul-transients.png')); end;


  % For analysis of counter-current and possible wastewater onshore transport (winters)
  sfomc.ne_buoy = verify_variable(sfomc.ne_buoy,{'adcp_u_36_h_sum','adcp_v_36_h_sum'});
  sfomc.ne_buoy = filter_gaps(sfomc.ne_buoy,'adcp_u','adcp_u_36_h_sum');
  sfomc.ne_buoy = filter_gaps(sfomc.ne_buoy,'adcp_v','adcp_v_36_h_sum');

  sfomc.sw_buoy = verify_variable(sfomc.sw_buoy,{'adcp_u_36_h_sum','adcp_v_36_h_sum'});
  sfomc.sw_buoy = filter_gaps(sfomc.sw_buoy,'adcp_u','adcp_u_36_h_sum');
  sfomc.sw_buoy = filter_gaps(sfomc.sw_buoy,'adcp_v','adcp_v_36_h_sum');

  sfomc.nw_w_btm = verify_variable(sfomc.nw_w_btm,{'adcp_u_36_h_sum','adcp_v_36_h_sum'});
  sfomc.nw_w_btm = filter_gaps(sfomc.nw_w_btm,'adcp_u','adcp_u_36_h_sum');
  sfomc.nw_w_btm = filter_gaps(sfomc.nw_w_btm,'adcp_v','adcp_v_36_h_sum');


  % 1999 COLD FRONTS
  % Focus on period of cold-front passage and counter-current
  xl = datenum(1999,10,[13,32]);

  fmg;
  subplot_tight(7,1,1:3);
  titlename('Section of sea temperatures and currents (1999 COLD FRONTS)');
  title('Temperatures (1999 COLD FRONTS)');
  lh=plot_ts(GAPFUN(fwyf1.ndbc_air_t),'k-',GAPFUN(sfomc.se_buoy.seatemp),GAPFUN(sfomc.ne_buoy.seatemp),GAPFUN(sfomc.sw_buoy.seatemp),GAPFUN(sfomc.nw_w_btm.adcp_seatemp));
  ylim([16.5,32.5]);
  grid on; grid minor;
  xlim_datetick(xl);
  legend(lh,'Air (FWYF1)','Offshore (SE)','Deep (NE)','Mid (SW)','Shallow (NW)', legargs{:});
  set(gca,'FontSize',fontsz); set(gca,'LineWidth',gridwid); set(get(gca,'Child'),'LineWidth',linewid);

  subplot_tight(7,1,4:5); %ylabel('Alongshore');
  title('Along-shore currents');
  lh=plot_ts(GAPFUN(sfomc.ne_buoy.adcp_v_36_h_sum),GAPFUN(sfomc.sw_buoy.adcp_v_36_h_sum),GAPFUN(sfomc.nw_w_btm.adcp_v_36_h_sum));
  grid on; grid minor;
  xlim_datetick(xl);
  legend(lh,'\Sigma_3_6_h Deep (NE)','\Sigma_3_6_h Mid (SW)','\Sigma_3_6_h Shallow (NW)', legargs{:});
  set(gca,'FontSize',fontsz); set(gca,'LineWidth',gridwid); set(get(gca,'Child'),'LineWidth',linewid);

  subplot_tight(7,1,6:7); %ylabel('Cross-shore');
  title('Cross-shore currents');
  lh=plot_ts(GAPFUN(sfomc.ne_buoy.adcp_u_36_h_sum),GAPFUN(sfomc.sw_buoy.adcp_u_36_h_sum),GAPFUN(sfomc.nw_w_btm.adcp_u_36_h_sum));
  grid on; grid minor;
  xlim_datetick('Y',xl);
  set(gca,'FontSize',fontsz); set(gca,'LineWidth',gridwid); set(get(gca,'Child'),'LineWidth',linewid);

  if doPrint; print('-dpng',fullfile(figspath,'upwelling-sfomc-1999-cross-vs-alongshore-Oct-36_h_sum-cold-fronts.png')); end;

  % clear ans dpax dpfh ekfh legargs shax shfh

  disp('Done with SFOMC events');
  if doPause; pause; else; pause(0.5); end;

end; %if ( ~doSFOMCevents ) else


if ( ~doSFOMCextras )
  if ( doSFOMC )
    disp('Not plotting extras for SFOMC (doSFOMCextras FALSE)');
  end;

else
  if doPause; pause; else; pause(0.5); end;

  xl = datenum(1999,[7,8],[15,28]); %Broader upwelling period
  basenm = ['upwelling-sfomc-mechanisms'];
  fignm = [basenm,datestr(xl(1),'_yyyymmmdd'),datestr(xl(end),'_yyyymmmdd')];
  dtstr = [datestr(xl(1),'mmm-dd'),' to ',datestr(xl(end),'mmm-dd-yyyy')];

  upwix = find(sfomc.ne_buoy.Nx.prof(1:end-1,end)==0 & sfomc.ne_buoy.Nx.prof(2:end,end)>0);
  dwnix = find(sfomc.ne_buoy.Nx.prof(1:end-1,end)>0 & sfomc.ne_buoy.Nx.prof(2:end,end)==0);


  fmg;

  ax(1)=subplot_tight(4,1,[1:3]);
  titlename(['Upwelling Mechanisms - ',dtstr]);
  title('Temperatures by depth [\circC]');
  plot_ts(GAPFUN(sfomc.ne_buoy.at_45m.sbe_seatemp),...
          GAPFUN(sfomc.ne_buoy.at_30m.sbe_seatemp),...
          GAPFUN(sfomc.nw_w_btm.adcp_seatemp),...
          GAPFUN(fwyf1.ndbc_air_t),'k');
  xlim_datetick(xl);
  ylim([15,31]); grid on; grid minor; hold on;
  %plot(sfomc.ne_buoy.Nx.date(upwix),repmat(23,size(upwix)),'r^');
  %plot(sfomc.ne_buoy.Nx.date(dwnix),repmat(23,size(dwnix)),'bv');
  %ylabel('Temperature [\circC]');
  annotate_ts_event(upwMr,sfomc.ne_buoy.at_45m.sbe_seatemp,upwmT,-inf);
  legend('45 m','30 m','11 m','Air','(Upw.)', 'Location','SouthWest');
  panel_label_subplot(gca,'a','ur','FontSize',fontsz,'BackgroundColor','w');
  set(gca,'FontSize',fontsz); set(gca,'LineWidth',gridwid); set(get(gca,'Child'),'LineWidth',linewid);
  
  ax(2)=subplot_tight(4,1,4);
  %title('Cross-shore currents by bin');
  title('Cross-shore curr.');
  plot_ts(GAPFUN(sfomc.ne_buoy.adcp_x_45m),...
          GAPFUN(sfomc.ne_buoy.adcp_x_30m),...
          GAPFUN(sfomc.ne_buoy.adcp_x_11m));
  ylim([-0.5,+0.5]); grid on; grid minor; hold on;
  %plot(sfomc.ne_buoy.Nx.date(upwix),repmat(-0.5,size(upwix)),'r^');
  %plot(sfomc.ne_buoy.Nx.date(dwnix),repmat(-0.5,size(dwnix)),'bv');
  %ylabel('u_x [m/s]');
  xlim_datetick('Y',xl);
  annotate_ts_event(upwMr,sfomc.ne_buoy.at_45m.sbe_seatemp,upwmT,-inf);
  %legend('45 m','30 m','11 m','(Upw.)', 'Location','SouthWest');
  panel_label_subplot(gca,'b','ur','FontSize',fontsz,'BackgroundColor','w');
  set(gca,'FontSize',fontsz); set(gca,'LineWidth',gridwid); set(get(gca,'Child'),'LineWidth',linewid);
  
  % annotation('doublearrow',[0.1809 0.2530],[0.2483 0.2483],'Color','r','LineWidth',2,'Head1Width',14,'Head2Width',14);
  % annotation('doublearrow',[0.2552 0.3178],[0.2483 0.2483],'Color','b','LineWidth',2,'Head1Width',14,'Head2Width',14);
  % annotation('doublearrow',[0.3193 0.4024],[0.2483 0.2483],'Color','m','LineWidth',2,'Head1Width',14,'Head2Width',14);
  phase_dts = datenum(1999,7,[21,24; 24,27; 28,31],[0,12; 12,23; 0,23],0,0);
  annot_h = textbox_range(ax,phase_dts,'ROMAN','ul','EdgeColor','k','LineWidth',2);

  if doPrint; print('-dpng',fullfile(figspath,[fignm,'_dynamics_1999.png'])); end;


  fmg;

  subplot(2,1,1);
  titlename('Offshore (NE/E) sea temperatures');
  plot(get_yearday(sfomc.ne_buoy.at_20m.sbe_seatemp.date),sfomc.ne_buoy.at_20m.sbe_seatemp.data,get_yearday(sfomc.e_buoy.at_20m.sbe_seatemp.date),sfomc.e_buoy.at_20m.sbe_seatemp.data);
  legend('1999','2000','Location','SouthEast');
  ylabel('20 m');
  axis([datenum(0,[7,9],15),19,31]);
  datetick3('x','dd-mmm','keeplimits');
  grid on; grid minor;
  panel_label_subplot(gca,'a','ll','FontSize',fontsz,'BackgroundColor','w');
  set(gca,'FontSize',fontsz); set(gca,'LineWidth',gridwid); set(get(gca,'Child'),'LineWidth',linewid);

  subplot(2,1,2);
  plot(get_yearday(sfomc.ne_buoy.at_30m.sbe_seatemp.date),sfomc.ne_buoy.at_30m.sbe_seatemp.data,get_yearday(sfomc.e_buoy.at_30m.sbe_seatemp.date),sfomc.e_buoy.at_30m.sbe_seatemp.data);
  legend('1999','2000','Location','SouthEast');
  ylabel('30 m');
  axis([datenum(0,[7,9],15),19,31]);
  datetick3('x','dd-mmm','keeplimits');
  grid on; grid minor;
  panel_label_subplot(gca,'b','ll','FontSize',fontsz,'BackgroundColor','w');
  set(gca,'FontSize',fontsz); set(gca,'LineWidth',gridwid); set(get(gca,'Child'),'LineWidth',linewid);
  if doPrint; print('-dpng',fullfile(figspath,[fignm,'_offshore_1999_vs_2000.png'])); end;

end; %if ( ~doSFOMCextras ) else



if ( doSEFCRI )
  legargs = { 'Orientation','horiz','Location','NorthWest' };

  dtlm = [datenum(2010,5,1),datenum(2013,10,1)];
  subfn = @(x)(ts_date_range(x,dtlm(1),dtlm(2)));

  if ( ~exist('fdepk_smooth_airtemp') )
    disp('Smoothing air temperature time series');
    % fdepk_smooth_airtemp = GAPFUN(interp_ts(subset_ts(fdepk.fdep_airtemp,subfn)));
    % lkwf1_smooth_airtemp = GAPFUN(interp_ts(subset_ts(lkwf1.ndbc_air_t,subfn)));
    % pvgf1_smooth_airtemp = GAPFUN(interp_ts(subset_ts(pvgf1.ndbc_air_t,subfn)));
    fdepk = verify_variable(fdepk,'fdep_airtemp_1_d_avg');
    lkwf1 = verify_variable(lkwf1,'ndbc_air_t_1_d_avg');
    pvgf1 = verify_variable(pvgf1,'ndbc_air_t_1_d_avg');
    fdepk_smooth_airtemp = GAPFUN(subset_ts(fdepk.fdep_airtemp_1_d_avg,subfn));
    lkwf1_smooth_airtemp = GAPFUN(subset_ts(lkwf1.ndbc_air_t_1_d_avg,subfn));
    pvgf1_smooth_airtemp = GAPFUN(subset_ts(pvgf1.ndbc_air_t_1_d_avg,subfn));
  end;

  fmg;
  jax(1) = subplot_tight(3,1,1); title('Far North (Martin County)');
   %lh=plot_ts(subset_ts(fwyf1.erai_air_t,subfn),'k-',sefcri.updb.hourly_t,sefcri.mc2.hourly_t);
   lh=plot_ts(fdepk_smooth_airtemp,'k-','LineWidth',1,sefcri.updb.hourly_t,sefcri.mc2.hourly_t);
   xlim_datetick(dtlm); %datetick3;
   ylim([12,32]);
   grid on; grid minor;
   %annotate_ts_event(upwMr,sefcri.updb.hourly_t,upwmT,12.5);
   annotate_ts_event(upwMr,sefcri.updb.hourly_t,upwmT,-inf);
   legend(['Air temp. (FDEPK)'],...
          ['Offshore (',num2str(sefcri.updb.depth),' m)'],...
          ['Near-shore (',num2str(sefcri.mc2.depth),' m)'],...
          '(Upwelling)', ...
          'Location','SouthWest', 'Orientation','horiz');
   set(gca,'FontSize',fontsz); set(gca,'LineWidth',gridwid); set(get(gca,'Child'),'LineWidth',linewid);
  jax(2) = subplot_tight(3,1,2); title('Northern (Palm Beach)');
   lh=plot_ts(lkwf1_smooth_airtemp,'k-','LineWidth',1,sefcri.pb2.hourly_t,sefcri.pb1.hourly_t);
   xlim_datetick(dtlm); %datetick3;
   ylim([12,32]);
   grid on; grid minor;
   %annotate_ts_event(upwMr,sefcri.pb2.hourly_t,upwmT,12.5);
   annotate_ts_event(upwMr,sefcri.pb2.hourly_t,upwmT,-inf);
   legend(['Air temp. (LKWF1)'],...
          ['Offshore (',num2str(sefcri.pb2.depth),' m)'],...
          ['Near-shore (',num2str(sefcri.pb1.depth),' m)'],...
          '(Upwelling)', ...
          'Location','SouthWest', 'Orientation','horiz');
   set(gca,'FontSize',fontsz); set(gca,'LineWidth',gridwid); set(get(gca,'Child'),'LineWidth',linewid);
  jax(3) = subplot_tight(3,1,3); title('Central (Ft. Lauderdale)');
   lh=plot_ts(pvgf1_smooth_airtemp,'k-','LineWidth',1,sefcri.bc3.hourly_t,sefcri.bc1.hourly_t);
   xlim_datetick('Y',dtlm); %datetick3;
   ylim([12,32]);
   grid on; grid minor;
   %annotate_ts_event(upwMr,sefcri.bc3.hourly_t,upwmT,12.5);
   annotate_ts_event(upwMr,sefcri.bc3.hourly_t,upwmT,-inf);
   legend(['Air temp. (PVGF1)'],...
          ['Offshore (',num2str(sefcri.bc3.depth),' m)'],...
          ['Near-shore (',num2str(sefcri.bc1.depth),' m)'],...
          '(Upwelling)', ...
          'Location','SouthWest', 'Orientation','horiz');
   set(gca,'FontSize',fontsz); set(gca,'LineWidth',gridwid); set(get(gca,'Child'),'LineWidth',linewid);
  if doPrint; print('-dpng',fullfile(figspath,'upwelling-sefcri-2010-2013-air.png')); end;


  prfh = fmg;

  prax(1) = subplot_tight(3,1,1); title('Far North (Martin County)');
   lh=plot_ts(sefcri.updb.hourly_t_3_h_hp_9_h_lp,sefcri.updb.hourly_t_11_h_hp_14_h_lp,sefcri.updb.hourly_t_17_h_hp_28_h_lp,sefcri.updb.hourly_t_48_h_hp_240_h_lp);
   xlim_datetick(dtlm); %datetick3;
   ylim([-0.075,+0.075]);
   %annotate_ts_event(upwMr,sefcri.updb.hourly_t,upwmT,-0.075);
   annotate_ts_event(upwMr,sefcri.updb.hourly_t,upwmT,-inf);
   %legend('3-9h','11-14h','17-28h','2-10d','(Upw.)', legargs{:});
   grid on; grid minor;
   set(gca,'FontSize',fontsz); set(gca,'LineWidth',gridwid); set(get(gca,'Child'),'LineWidth',linewid);
  prax(2) = subplot_tight(3,1,2); title('Northern (Palm Beach)');
   lh=plot_ts(sefcri.pb2.hourly_t_3_h_hp_9_h_lp,sefcri.pb2.hourly_t_11_h_hp_14_h_lp,sefcri.pb2.hourly_t_17_h_hp_28_h_lp,sefcri.pb2.hourly_t_48_h_hp_240_h_lp);
   xlim_datetick(dtlm); %datetick3;
   ylim([-0.075,+0.075]);
   %annotate_ts_event(upwMr,sefcri.pb2.hourly_t,upwmT,-0.075);
   annotate_ts_event(upwMr,sefcri.pb2.hourly_t,upwmT,-inf);
   %legend('3-9h','11-14h','17-28h','2-10d','(Upw.)', legargs{:});
   grid on; grid minor;
   set(gca,'FontSize',fontsz); set(gca,'LineWidth',gridwid); set(get(gca,'Child'),'LineWidth',linewid);
  prax(3) = subplot_tight(3,1,3); title('Central (Ft. Lauderdale)');
   lh=plot_ts(sefcri.bc3.hourly_t_3_h_hp_9_h_lp,sefcri.bc3.hourly_t_11_h_hp_14_h_lp,sefcri.bc3.hourly_t_17_h_hp_28_h_lp,sefcri.bc3.hourly_t_48_h_hp_240_h_lp);
   xlim_datetick('Y',dtlm); %datetick3;
   ylim([-0.075,+0.075]);
   %annotate_ts_event(upwMr,sefcri.bc3.hourly_t,upwmT,-0.075);
   annotate_ts_event(upwMr,sefcri.bc3.hourly_t,upwmT,-inf);
   legend('3-9h','11-14h','17-28h','2-10d','(Upw.)', legargs{:});
   grid on; grid minor;
   set(gca,'FontSize',fontsz); set(gca,'LineWidth',gridwid); set(get(gca,'Child'),'LineWidth',linewid);
  if doPrint; print('-dpng',fullfile(figspath,'upwelling-sefcri-2010-2013-partn.png')); end;


  for yr=2010:2013
    yrstr = num2str(yr);
    disp(['Hit Enter to display ',yrstr]);
    if doPause; pause; else; pause(0.5); end;
    xlim_datetick(jax(1:2),datenum(yr,[5,10],1)); xlim_datetick('Y',jax(3),datenum(yr,[5,10],1)); %datetick3;
    if doPrint; print('-dpng',fullfile(figspath,['upwelling-sefcri-',yrstr,'-air.png'])); end;
    xlim_datetick(prax(1:2),datenum(yr,[5,10],1)); xlim_datetick('Y',prax(3),datenum(yr,[5,10],1)); %datetick3;
    if doPrint; print('-dpng',fullfile(figspath,['upwelling-sefcri-',yrstr,'-partn.png'])); end;
  end;

end; %if ( doSEFCRI )


if ( doFACE )

  % FACE temperature/currents (2014)

  if 1;
    %%xl = datenum(2014,[7,10],[30,1]);
    %xl = datenum(2014,[7,10],[18,1]);
    xl = datenum(2014,[7,10],[25,1]);

    % Ocean deep response
    dpfh = fmg;

    dpax(1) = subplot_tight(7,1,1:3);
    titlename('FACE HWD moorings - air/sea temperature [\circC] (2014)');
    title('Air/sea temperatures [\circC]');
    % Match Deep and Shallow colors with the 2015 plot (below)
    lh=plot_ts(GAPFUN(pvgf1.ndbc_air_t),'k-',...
               GAPFUN(face.deep.seatemp),'specix',1,...
               GAPFUN(face.shallow.seatemp),'specix',2);
    xlim_datetick(xl);
    ylim([20.0,33.0]);
    %annotate_ts_event(upwMr,face.deep.seatemp,upwmT,20.0);
    annotate_ts_event(upwMr,face.deep.seatemp,upwmT,-inf);
    legend('Air (PVGF1)','Deep','Shallow','(Upwelling)',...
           'Location','SouthEast', 'Orientation','horiz');
    grid on; grid minor;
    panel_label_subplot(gca,'a','ur','FontSize',fontsz,'BackgroundColor','w');
    set(gca,'FontSize',fontsz); set(gca,'LineWidth',gridwid); set(get(gca,'Child'),'LineWidth',linewid);

    dpax(2) = subplot_tight(7,1,4:5);
    title('Cross-shore near-bottom currents [ms^-^1]');
    lh=plot_ts(GAPFUN(face.deep.adcp_btm_x),'specix',1,...
               GAPFUN(face.shallow.adcp_btm_x),'specix',2);
    xlim_datetick(xl);
    ylim([-0.425,+0.425]);
    %annotate_ts_event(upwMr,face.deep.seatemp,upwmT,-0.425);
    annotate_ts_event(upwMr,face.deep.seatemp,upwmT,-inf);
    %legend('Deep','Shallow','(Upwelling)', 'Location','SouthEast','Orientation','horiz');
    grid on; grid minor;
    panel_label_subplot(gca,'b','ur','FontSize',fontsz,'BackgroundColor','w');
    set(gca,'FontSize',fontsz); set(gca,'LineWidth',gridwid); set(get(gca,'Child'),'LineWidth',linewid);

    dpax(3) = subplot_tight(7,1,6:7);
    title('Along-shore near-bottom currents [ms^-^1]');
    lh=plot_ts(GAPFUN(face.deep.adcp_btm_l),'specix',1,...
               GAPFUN(face.shallow.adcp_btm_l),'specix',2);
    xlim_datetick('Y',xl);
    ylim([-1.25,+1.25]);
    %annotate_ts_event(upwMr,face.deep.seatemp,upwmT,-1.25);
    annotate_ts_event(upwMr,face.deep.seatemp,upwmT,-inf);
    %legend('Deep','Shallow','(Upwelling)', 'Location','SouthEast','Orientation','horiz');
    grid on; grid minor;
    panel_label_subplot(gca,'c','ur','FontSize',fontsz,'BackgroundColor','w');
    set(gca,'FontSize',fontsz); set(gca,'LineWidth',gridwid); set(get(gca,'Child'),'LineWidth',linewid);

    annot_h(1) = dsuannotation('textarrow',...
                               [double(dpax(2)),datenum(2014,8,8,3,20,0),double(dpax(1)),datenum(2014,8,10,3,0,0)],...
                               [double(dpax(2)),0.35,double(dpax(1)),25.5],'String','Onset',...
                               'FontSize',fontsz,'VerticalAlignment','middle',...
                               'LineWidth',5,'LineStyle',':','HeadWidth',fontsz,'HeadLength',fontsz);
    annot_h(2) = dsuannotation('textarrow',...
                               [double(dpax(3)),datenum(2014,8,22,18,20,0),double(dpax(1)),datenum(2014,8,22,13,0,0)],...
                               [double(dpax(3)),0.85,double(dpax(1)),25.0],'String','Relaxation',...
                               'FontSize',fontsz,'VerticalAlignment','middle',...
                               'LineWidth',5,'LineStyle',':','HeadWidth',fontsz,'HeadLength',fontsz);

    if doPrint; print('-dpng',fullfile(figspath,'upwelling-face-2014.png')); end;


    % Ekman response
    ekfh = fmg;

    ekax(1) = subplot_tight(7,1,1:3);
    titlename('Cross-shore Ekman volumetric flux [m^2s^-^1] near FACE section (2014)');
    title('Cross-shore Ekman volumetric flux [m^2s^-^1]');
    lh=plot_ts(GAPFUN(canf1.ndbc_ekman_flux_x_volume),GAPFUN(lkwf1.ndbc_ekman_flux_x_volume),GAPFUN(pvgf1.ndbc_ekman_flux_x_volume),GAPFUN(rsm.rsmas_ekman_flux_x_volume),GAPFUN(fwyf1.ndbc_ekman_flux_x_volume));
    xlim_datetick(xl);
    ylim([-4.9,+4.9]);
    %annotate_ts_event(upwMr,face.deep.seatemp,upwmT,-4.9);
    annotate_ts_event(upwMr,face.deep.seatemp,upwmT,-inf);
    center_axes; grid on; grid minor;
    legh=legend('Cape Can. (''88-''16)','Lake Worth (''01-''16)','Port Ev. (''09-''17)','RSMAS Roof (''03-''11)','Fowey (''91-''16)','(Upw.)', 'Orientation','horiz','Location','SouthEast'); set(legh,'FontSize',19)
    panel_label_subplot(gca,'d','ur','FontSize',fontsz,'BackgroundColor','w');
    set(gca,'FontSize',fontsz); set(gca,'LineWidth',gridwid); set(get(gca,'Child'),'LineWidth',linewid);

    ekax(2) = subplot_tight(7,1,4:6);
    title('Along-shore Ekman volumetric flux [m^2s^-^1]');
    lh=plot_ts(GAPFUN(canf1.ndbc_ekman_flux_l_volume),GAPFUN(lkwf1.ndbc_ekman_flux_l_volume),GAPFUN(pvgf1.ndbc_ekman_flux_l_volume),GAPFUN(rsm.rsmas_ekman_flux_l_volume),GAPFUN(fwyf1.ndbc_ekman_flux_l_volume));
    xlim_datetick(xl);
    ylim([-4.9,+4.9]);
    %annotate_ts_event(upwMr,face.deep.seatemp,upwmT,-4.9);
    annotate_ts_event(upwMr,face.deep.seatemp,upwmT,-inf);
    center_axes; grid on; grid minor;
    panel_label_subplot(gca,'e','ur','FontSize',fontsz,'BackgroundColor','w');
    set(gca,'FontSize',fontsz); set(gca,'LineWidth',gridwid); set(get(gca,'Child'),'LineWidth',linewid);

    ekax(3) = subplot_tight(7,1,7);
    title('Barometric pressure [mbar]');
    lh=plot_ts(GAPFUN(canf1.ndbc_barom),GAPFUN(lkwf1.ndbc_barom),GAPFUN(pvgf1.ndbc_barom),GAPFUN(rsm.rsmas_barom),GAPFUN(fwyf1.ndbc_barom));
    xlim_datetick('Y',xl);
    %ylim([980,1040]);
    ylim([1006,1023])
    %annotate_ts_event(upwMr,face.deep.seatemp,upwmT,1006);
    annotate_ts_event(upwMr,face.deep.seatemp,upwmT,-inf);
    grid on; grid minor;
    panel_label_subplot(gca,'f','ur','FontSize',fontsz,'BackgroundColor','w');
    set(gca,'FontSize',fontsz); set(gca,'LineWidth',gridwid); set(get(gca,'Child'),'LineWidth',linewid);

    if doPrint; print('-dpng',fullfile(figspath,'upwelling-face-2014-ekman.png')); end;


    % Variability partition
    prfh = fmg;

    prax(1) = subplot_tight(7,1,1:3);
    titlename('Sea temperature variability [\circC]: FACE deep mooring (2014)');
    title('Sea temperature variability [\circC]');
    lh=plot_ts(GAPFUN(face.deep.seatemp_3_h_hp_9_h_lp),...
               GAPFUN(face.deep.seatemp_11_h_hp_14_h_lp),...
               GAPFUN(face.deep.seatemp_17_h_hp_28_h_lp),...
               GAPFUN(face.deep.seatemp_48_h_hp_240_h_lp));
    xlim_datetick(xl);
    center_axes;
    annotate_ts_event(upwMr,face.deep.seatemp,upwmT,-inf);
    legend('3-9h','11-14h','17-28h','2-10d','(Upw.)', 'Orientation','horiz','Location','SouthEast');
    grid on; grid minor;
    panel_label_subplot(gca,'g','ur','FontSize',fontsz,'BackgroundColor','w');
    set(gca,'FontSize',fontsz); set(gca,'LineWidth',gridwid); set(get(gca,'Child'),'LineWidth',linewid);

    prax(2) = subplot_tight(7,1,4:5);
    title('Cross-shore near-bottom current var. [ms^-^1]');
    lh=plot_ts(GAPFUN(face.deep.adcp_btm_x_3_h_hp_9_h_lp),...
               GAPFUN(face.deep.adcp_btm_x_11_h_hp_14_h_lp),...
               GAPFUN(face.deep.adcp_btm_x_17_h_hp_28_h_lp),...
               GAPFUN(face.deep.adcp_btm_x_48_h_hp_240_h_lp));
    xlim_datetick(xl); %datetick3;
    %center_axes;
    ylim([-0.1,+0.1]);
    annotate_ts_event(upwMr,face.deep.seatemp,upwmT,-inf);
    grid on; grid minor;
    %legend('3-9h','11-14h','17-28h','2-10d','(Upw.)', 'Orientation','horiz');
    panel_label_subplot(gca,'h','ur','FontSize',fontsz,'BackgroundColor','w');
    set(gca,'FontSize',fontsz); set(gca,'LineWidth',gridwid); set(get(gca,'Child'),'LineWidth',linewid);

    prax(3) = subplot_tight(7,1,6:7);
    title('Along-shore near-bottom current var. [ms^-^1]');
    lh=plot_ts(GAPFUN(face.deep.adcp_btm_l_3_h_hp_9_h_lp),...
               GAPFUN(face.deep.adcp_btm_l_11_h_hp_14_h_lp),...
               GAPFUN(face.deep.adcp_btm_l_17_h_hp_28_h_lp),...
               GAPFUN(face.deep.adcp_btm_l_48_h_hp_240_h_lp));
    %center_axes;
    ylim([-0.3,+0.3]);
    annotate_ts_event(upwMr,face.deep.seatemp,upwmT,-inf);
    grid on; grid minor;
    %legend('3-9h','11-14h','17-28h','2-10d','(Upw.)', 'Orientation','horiz');
    panel_label_subplot(gca,'i','ur','FontSize',fontsz,'BackgroundColor','w');
    set(gca,'FontSize',fontsz); set(gca,'LineWidth',gridwid); set(get(gca,'Child'),'LineWidth',linewid);
    xlim_datetick('Y',xl);

    if doPrint; print('-dpng',fullfile(figspath,'upwelling-face-2014-partn.png')); end;
  end;


  % FACE temperature/currents (2015)

  if 1;
    % %xl = datenum(2015,[7,10],[30,1]);
    % xl = datenum(2015,[7,10],[18,1]);
    % %xl = datenum(2015,[7,10],[25,1]);
    xl = datenum(2015,[7,9],[22,28]);

    % Ocean deep response
    dpfh = fmg;

    dpax(1) = subplot_tight(7,1,1:3);
    titlename('FACE HWD moorings - air/sea temperature [\circC] (2015)');
    title('Air/sea temperatures [\circC]');
    lh=plot_ts(GAPFUN(pvgf1.ndbc_air_t),'k-',...
               GAPFUN(face.deep.seatemp),'specix',1,...
               GAPFUN(face.tcm2.seatemp),'specix',3,...
               GAPFUN(face.tcm1.seatemp),'specix',4,...
               GAPFUN(face.shallow.seatemp),'specix',2);
    xlim_datetick(xl);
    ylim([20.0,33.0]);
    %annotate_ts_event(upwMr,face.deep.seatemp,upwmT,20.0);
    annotate_ts_event(upwMr,face.deep.seatemp,upwmT,-inf);
    legend('Air (PVGF1)','Deep','TCM2','TCM1','Shallow','(Upwelling)',...
           'Location','SouthEast', 'Orientation','horiz');
    grid on; grid minor;
    panel_label_subplot(gca,'d','ur','FontSize',fontsz,'BackgroundColor','w');
    set(gca,'FontSize',fontsz); set(gca,'LineWidth',gridwid); set(get(gca,'Child'),'LineWidth',linewid);

    dpax(2) = subplot_tight(7,1,4:5);
    title('Cross-shore near-bottom currents [ms^-^1]');
    lh=plot_ts(GAPFUN(face.deep.adcp_btm_x),'specix',1,...
               GAPFUN(face.tcm2.alt_x),'specix',3,...
               GAPFUN(face.tcm1.alt_x),'specix',4,...
               GAPFUN(face.shallow.adcp_btm_x),'specix',2);
    xlim_datetick(xl);
    ylim([-0.425,+0.425]);
    %annotate_ts_event(upwMr,face.deep.seatemp,upwmT,-0.425);
    annotate_ts_event(upwMr,face.deep.seatemp,upwmT,-inf);
    %legend('Deep','TCM2','TCM1','Shallow','(Upwelling)', 'Location','SouthEast','Orientation','horiz');
    grid on; grid minor;
    panel_label_subplot(gca,'e','ur','FontSize',fontsz,'BackgroundColor','w');
    set(gca,'FontSize',fontsz); set(gca,'LineWidth',gridwid); set(get(gca,'Child'),'LineWidth',linewid);

    dpax(3) = subplot_tight(7,1,6:7);
    title('Along-shore near-bottom currents [ms^-^1]');
    lh=plot_ts(GAPFUN(face.deep.adcp_btm_l),'specix',1,...
               GAPFUN(face.tcm2.alt_l),'specix',3,...
               GAPFUN(face.tcm1.alt_l),'specix',4,...
               GAPFUN(face.shallow.adcp_btm_l),'specix',2);
    xlim_datetick('Y',xl);
    ylim([-1.25,+1.25]);
    %annotate_ts_event(upwMr,face.deep.seatemp,upwmT,-1.25);
    annotate_ts_event(upwMr,face.deep.seatemp,upwmT,-inf);
    %legend('Deep','TCM2','TCM1','Shallow','(Upwelling)', 'Location','SouthEast','Orientation','horiz');
    grid on; grid minor;
    panel_label_subplot(gca,'f','ur','FontSize',fontsz,'BackgroundColor','w');
    set(gca,'FontSize',fontsz); set(gca,'LineWidth',gridwid); set(get(gca,'Child'),'LineWidth',linewid);

    annot_h=[];
    annot_h(end+1) = ...
        dsuannotation('textarrow',...
                      [double(dpax(2)),datenum(2015,7,25,04,20,0),double(dpax(1)),datenum(2015,7,28,0,0,0)],...
                      [double(dpax(2)),0.35,double(dpax(1)),25.5],'String','Onset',...
                      'FontSize',fontsz,'VerticalAlignment','middle',...
                      'LineWidth',5,'LineStyle',':','HeadWidth',fontsz,'HeadLength',fontsz);
    annot_h(end+1) = ...
        dsuannotation('textarrow',...
                      [double(dpax(3)),datenum(2015,8,2,18,0,0),double(dpax(1)),datenum(2015,8,3,12,0,0)],...
                      [double(dpax(3)),0.85,double(dpax(1)),25.5],'String','Relaxation',...
                      'FontSize',fontsz,'VerticalAlignment','middle',...
                      'LineWidth',5,'LineStyle',':','HeadWidth',fontsz,'HeadLength',fontsz);
    annot_h(end+1) = ...
        dsuannotation('textarrow',...
                      [double(dpax(2)),datenum(2015,8,10,3,0,0),double(dpax(1)),datenum(2015,8,13,0,20,0)],...
                      [double(dpax(2)),0.35,double(dpax(1)),25.5],'String','Onset',...
                      'FontSize',fontsz,'VerticalAlignment','middle',...
                      'LineWidth',5,'LineStyle',':','HeadWidth',fontsz,'HeadLength',fontsz);
    annot_h(end+1) = ...
        dsuannotation('textarrow',...
                      [double(dpax(3)),datenum(2015,8,21,08,20,0),double(dpax(1)),datenum(2015,8,22,17,0,0)],...
                      [double(dpax(3)),0.85,double(dpax(1)),25.0],'String','Relaxation',...
                      'FontSize',fontsz,'VerticalAlignment','middle',...
                      'LineWidth',5,'LineStyle',':','HeadWidth',fontsz,'HeadLength',fontsz);
    annot_h(end+1) = ...
        dsuannotation('textarrow',...
                      [double(dpax(2)),datenum(2015,9,5,7,20,0),double(dpax(1)),datenum(2015,9,8,11,0,0)],...
                      [double(dpax(2)),0.35,double(dpax(1)),25.5],'String','Onset',...
                      'FontSize',fontsz,'VerticalAlignment','middle',...
                      'LineWidth',5,'LineStyle',':','HeadWidth',fontsz,'HeadLength',fontsz);
    annot_h(end+1) = ...
        dsuannotation('textarrow',...
                      [double(dpax(3)),datenum(2015,9,15,16,40,0),double(dpax(1)),datenum(2015,9,16,4,42,0)],...
                      [double(dpax(3)),0.85,double(dpax(1)),25.0],'String','Relaxation',...
                      'FontSize',fontsz,'VerticalAlignment','middle',...
                      'LineWidth',5,'LineStyle',':','HeadWidth',fontsz,'HeadLength',fontsz);

    if doPrint; print('-dpng',fullfile(figspath,'upwelling-face-2015.png')); end;


    % Ekman response
    ekfh = fmg;

    ekax(1) = subplot_tight(7,1,1:3);
    titlename('Cross-shore Ekman volumetric flux [m^2s^-^1] near FACE section (2015)');
    title('Cross-shore Ekman volumetric flux [m^2s^-^1]');
    lh=plot_ts(GAPFUN(canf1.ndbc_ekman_flux_x_volume),GAPFUN(lkwf1.ndbc_ekman_flux_x_volume),GAPFUN(pvgf1.ndbc_ekman_flux_x_volume),GAPFUN(rsm.rsmas_ekman_flux_x_volume),GAPFUN(fwyf1.ndbc_ekman_flux_x_volume));
    xlim_datetick(xl);
    ylim([-4.9,+4.9]);
    %annotate_ts_event(upwMr,face.deep.seatemp,upwmT,-4.9);
    annotate_ts_event(upwMr,face.deep.seatemp,upwmT,-inf);
    center_axes; grid on; grid minor;
    legh=legend('Cape Can. (''88-''16)','Lake Worth (''01-''16)','Port Ev. (''09-''17)','RSMAS Roof (''03-''11)','Fowey (''91-''16)','(Upw.)', 'Orientation','horiz','Location','SouthEast'); set(legh,'FontSize',19);
    panel_label_subplot(gca,'g','ur','FontSize',fontsz,'BackgroundColor','w');
    set(gca,'FontSize',fontsz); set(gca,'LineWidth',gridwid); set(get(gca,'Child'),'LineWidth',linewid);

    ekax(2) = subplot_tight(7,1,4:6);
    title('Along-shore Ekman volumetric flux [m^2s^-^1]');
    lh=plot_ts(GAPFUN(canf1.ndbc_ekman_flux_l_volume),GAPFUN(lkwf1.ndbc_ekman_flux_l_volume),GAPFUN(pvgf1.ndbc_ekman_flux_l_volume),GAPFUN(rsm.rsmas_ekman_flux_l_volume),GAPFUN(fwyf1.ndbc_ekman_flux_l_volume));
    xlim_datetick(xl);
    ylim([-4.9,+4.9]);
    %annotate_ts_event(upwMr,face.deep.seatemp,upwmT,-4.9);
    annotate_ts_event(upwMr,face.deep.seatemp,upwmT,-inf);
    center_axes; grid on; grid minor;
    panel_label_subplot(gca,'h','ur','FontSize',fontsz,'BackgroundColor','w');
    set(gca,'FontSize',fontsz); set(gca,'LineWidth',gridwid); set(get(gca,'Child'),'LineWidth',linewid);

    ekax(3) = subplot_tight(7,1,7);
    title('Barometric pressure [mbar]');
    lh=plot_ts(GAPFUN(canf1.ndbc_barom),GAPFUN(lkwf1.ndbc_barom),GAPFUN(pvgf1.ndbc_barom),GAPFUN(rsm.rsmas_barom),GAPFUN(fwyf1.ndbc_barom));
    xlim_datetick('Y',xl);
    %ylim([980,1040]);
    ylim([1006,1023])
    %annotate_ts_event(upwMr,face.deep.seatemp,upwmT,1006);
    annotate_ts_event(upwMr,face.deep.seatemp,upwmT,-inf);
    grid on; grid minor;
    panel_label_subplot(gca,'i','ur','FontSize',fontsz,'BackgroundColor','w');
    set(gca,'FontSize',fontsz); set(gca,'LineWidth',gridwid); set(get(gca,'Child'),'LineWidth',linewid);

    if doPrint; print('-dpng',fullfile(figspath,'upwelling-face-2015-ekman.png')); end;


    % Variability partition
    prfh = fmg;

    prax(1) = subplot_tight(7,1,1:3);
    titlename('Sea temperature variability [\circC]: FACE deep mooring (2015)');
    title('Sea temperature variability [\circC]');
    lh=plot_ts(GAPFUN(face.deep.seatemp_3_h_hp_9_h_lp),...
               GAPFUN(face.deep.seatemp_11_h_hp_14_h_lp),...
               GAPFUN(face.deep.seatemp_17_h_hp_28_h_lp),...
               GAPFUN(face.deep.seatemp_48_h_hp_240_h_lp));
    xlim_datetick(xl);
    center_axes;
    annotate_ts_event(upwMr,face.deep.seatemp,upwmT,-inf);
    legend('3-9h','11-14h','17-28h','2-10d','(Upw.)', 'Orientation','horiz','Location','SouthEast');
    grid on; grid minor;
    panel_label_subplot(gca,'j','ur','FontSize',fontsz,'BackgroundColor','w');
    set(gca,'FontSize',fontsz); set(gca,'LineWidth',gridwid); set(get(gca,'Child'),'LineWidth',linewid);

    prax(2) = subplot_tight(7,1,4:5);
    title('Cross-shore near-bottom current var. [ms^-^1]');
    lh=plot_ts(GAPFUN(face.deep.adcp_btm_x_3_h_hp_9_h_lp),...
               GAPFUN(face.deep.adcp_btm_x_11_h_hp_14_h_lp),...
               GAPFUN(face.deep.adcp_btm_x_17_h_hp_28_h_lp),...
               GAPFUN(face.deep.adcp_btm_x_48_h_hp_240_h_lp));
    xlim_datetick(xl); %datetick3;
    %center_axes;
    ylim([-0.1,+0.1]);
    annotate_ts_event(upwMr,face.deep.seatemp,upwmT,-inf);
    %legend('3-9h','11-14h','17-28h','2-10d','(Upw.)', 'Orientation','horiz');
    grid on; grid minor;
    panel_label_subplot(gca,'k','ur','FontSize',fontsz,'BackgroundColor','w');
    set(gca,'FontSize',fontsz); set(gca,'LineWidth',gridwid); set(get(gca,'Child'),'LineWidth',linewid);

    prax(3) = subplot_tight(7,1,6:7);
    title('Along-shore near-bottom current var. [ms^-^1]');
    lh=plot_ts(GAPFUN(face.deep.adcp_btm_l_3_h_hp_9_h_lp),...
               GAPFUN(face.deep.adcp_btm_l_11_h_hp_14_h_lp),...
               GAPFUN(face.deep.adcp_btm_l_17_h_hp_28_h_lp),...
               GAPFUN(face.deep.adcp_btm_l_48_h_hp_240_h_lp));
    %center_axes;
    ylim([-0.3,+0.3]);
    annotate_ts_event(upwMr,face.deep.seatemp,upwmT,-inf);
    %legend('3-9h','11-14h','17-28h','2-10d','(Upw.)', 'Orientation','horiz');
    grid on; grid minor;
    panel_label_subplot(gca,'l','ur','FontSize',fontsz,'BackgroundColor','w');
    set(gca,'FontSize',fontsz); set(gca,'LineWidth',gridwid); set(get(gca,'Child'),'LineWidth',linewid);
    xlim_datetick('Y',xl);

    if doPrint; print('-dpng',fullfile(figspath,'upwelling-face-2015-partn.png')); end;
  end;


  % FACE temperature/currents (Sep 2015)

  if 1;
    xl = datenum(2015,[8,10],[31,1]);

    % Ocean deep response
    dpfh = fmg;

    dpax(1) = subplot_tight(7,1,1:3);
    titlename('FACE HWD moorings - air/sea temp. & near-bottom curr. [\circC] (Sep 2015)');
    title('Air/sea temperatures [\circC]');
    lh=plot_ts(GAPFUN(pvgf1.ndbc_air_t),'k-',...
               GAPFUN(face.deep.seatemp),'specix',1,...
               GAPFUN(face.tcm2.seatemp),'specix',3,...
               GAPFUN(face.tcm1.seatemp),'specix',4,...
               GAPFUN(face.shallow.seatemp),'specix',2);
    xlim_datetick(xl);
    ylim([20.0,33.0]);
    %annotate_ts_event(upwMr,face.deep.seatemp,upwmT,20.0);
    annotate_ts_event(upwMr,face.deep.seatemp,upwmT,-inf);
    legend('Air (PVGF1)','Deep','TCM2','TCM1','Shallow','(Upwelling)',...
           'Location','SouthEast', 'Orientation','horiz');
    grid on; grid minor;
    panel_label_subplot(gca,'a','ur','FontSize',fontsz,'BackgroundColor','w');
    set(gca,'FontSize',fontsz); set(gca,'LineWidth',gridwid); set(get(gca,'Child'),'LineWidth',linewid);

    dpax(2) = subplot_tight(7,1,4:5);
    title('Cross-shore near-bottom currents [ms^-^1]');
    lh=plot_ts(GAPFUN(face.deep.adcp_btm_x),'specix',1,...
               GAPFUN(face.tcm2.alt_x),'specix',3,...
               GAPFUN(face.tcm1.alt_x),'specix',4,...
               GAPFUN(face.shallow.adcp_btm_x),'specix',2);
    xlim_datetick(xl);
    ylim([-0.425,+0.425]);
    %annotate_ts_event(upwMr,face.deep.seatemp,upwmT,-0.425);
    annotate_ts_event(upwMr,face.deep.seatemp,upwmT,-inf);
    %legend('Deep','TCM2','TCM1','Shallow','(Upwelling)', 'Location','SouthEast','Orientation','horiz');
    grid on; grid minor;
    panel_label_subplot(gca,'b','ur','FontSize',fontsz,'BackgroundColor','w');
    set(gca,'FontSize',fontsz); set(gca,'LineWidth',gridwid); set(get(gca,'Child'),'LineWidth',linewid);

    dpax(3) = subplot_tight(7,1,6:7);
    title('Along-shore near-bottom currents [ms^-^1]');
    lh=plot_ts(GAPFUN(face.deep.adcp_btm_l),'specix',1,...
               GAPFUN(face.tcm2.alt_l),'specix',3,...
               GAPFUN(face.tcm1.alt_l),'specix',4,...
               GAPFUN(face.shallow.adcp_btm_l),'specix',2);
    xlim_datetick('Y',xl);
    ylim([-1.25,+1.25]);
    %annotate_ts_event(upwMr,face.deep.seatemp,upwmT,-1.25);
    annotate_ts_event(upwMr,face.deep.seatemp,upwmT,-inf);
    %legend('Deep','TCM2','TCM1','Shallow','(Upwelling)', 'Location','SouthEast','Orientation','horiz');
    grid on; grid minor;
    panel_label_subplot(gca,'c','ur','FontSize',fontsz,'BackgroundColor','w');
    set(gca,'FontSize',fontsz); set(gca,'LineWidth',gridwid); set(get(gca,'Child'),'LineWidth',linewid);

    annot_h=[];
    annot_h(end+1) = ...
        dsuannotation('textarrow',...
                      [double(dpax(2)),datenum(2015,9,5,7,20,0),double(dpax(1)),datenum(2015,9,8,11,0,0)],...
                      [double(dpax(2)),0.35,double(dpax(1)),25.5],'String','Onset',...
                      'FontSize',fontsz,'VerticalAlignment','middle',...
                      'LineWidth',5,'LineStyle',':','HeadWidth',fontsz,'HeadLength',fontsz);
    annot_h(end+1) = ...
        dsuannotation('textarrow',...
                      [double(dpax(3)),datenum(2015,9,15,16,40,0),double(dpax(1)),datenum(2015,9,16,4,42,0)],...
                      [double(dpax(3)),0.85,double(dpax(1)),25.0],'String','Relaxation',...
                      'FontSize',fontsz,'VerticalAlignment','middle',...
                      'LineWidth',5,'LineStyle',':','HeadWidth',fontsz,'HeadLength',fontsz);


    % Ocean near-surface response
    shfh = fmg;

    shax(1) = subplot_tight(7,1,1:3);
    titlename('FACE HWD moorings - air/sea temp. & near-surface curr. [\circC] (Sep 2015)');
    title('Air/sea temperatures [\circC]');
    lh=plot_ts(GAPFUN(pvgf1.ndbc_air_t),'k-',...
               GAPFUN(face.deep.seatemp),'specix',1,...
               GAPFUN(face.tcm2.seatemp),'specix',3,...
               GAPFUN(face.tcm1.seatemp),'specix',4,...
               GAPFUN(face.shallow.seatemp),'specix',2);
    xlim_datetick(xl);
    ylim([20.0,33.0]);
    %annotate_ts_event(upwMr,face.deep.seatemp,upwmT,20.0);
    annotate_ts_event(upwMr,face.deep.seatemp,upwmT,-inf);
    legend('Air (PVGF1)','Deep','TCM2','TCM1','Shallow','(Upwelling)', 'Location','SouthEast','Orientation','horiz');
    grid on; grid minor;
    panel_label_subplot(gca,'d','ur','FontSize',fontsz,'BackgroundColor','w');
    set(gca,'FontSize',fontsz); set(gca,'LineWidth',gridwid); set(get(gca,'Child'),'LineWidth',linewid);

    shax(2) = subplot_tight(7,1,4:5);
    title('Cross-shore near-surface currents [ms^-^1]');
    lh=plot_ts(GAPFUN(face.deep.adcp_sfc_x),'specix',1,...
               GAPFUN(face.shallow.adcp_sfc_x),'specix',2);
    xlim_datetick(xl);
    %ylim([-0.25,+0.25]);
    %ylim([-0.42,+0.42]);
    ylim([-0.425,+0.425]);
    %annotate_ts_event(upwMr,face.deep.seatemp,upwmT,-0.425);
    annotate_ts_event(upwMr,face.deep.seatemp,upwmT,-inf);
    %legend('Deep','Shallow','(Upwelling)', 'Location','SouthEast','Orientation','horiz');
    grid on; grid minor;
    panel_label_subplot(gca,'e','ur','FontSize',fontsz,'BackgroundColor','w');
    set(gca,'FontSize',fontsz); set(gca,'LineWidth',gridwid); set(get(gca,'Child'),'LineWidth',linewid);

    shax(3) = subplot_tight(7,1,6:7);
    title('Along-shore near-surface currents [ms^-^1]');
    lh=plot_ts(GAPFUN(face.deep.adcp_sfc_l),'specix',1,...
               GAPFUN(face.shallow.adcp_sfc_l),'specix',2);
    xlim_datetick('Y',xl);
    ylim([-1.25,+1.25]);
    %annotate_ts_event(upwMr,face.deep.seatemp,upwmT,-1.25);
    annotate_ts_event(upwMr,face.deep.seatemp,upwmT,-inf);
    %legend('Deep','Shallow','(Upwelling)', 'Location','SouthEast','Orientation','horiz');
    grid on; grid minor;
    panel_label_subplot(gca,'f','ur','FontSize',fontsz,'BackgroundColor','w');
    set(gca,'FontSize',fontsz); set(gca,'LineWidth',gridwid); set(get(gca,'Child'),'LineWidth',linewid);


    % Ekman response
    ekfh = fmg;

    ekax(1) = subplot_tight(7,1,1:3);
    titlename('Ekman volumetric flux [m^2s^-^1] & barom. near FACE section (Sep 2015)');
    title('Cross-shore Ekman volumetric flux [m^2s^-^1]');
    lh=plot_ts(GAPFUN(canf1.ndbc_ekman_flux_x_volume),GAPFUN(lkwf1.ndbc_ekman_flux_x_volume),GAPFUN(pvgf1.ndbc_ekman_flux_x_volume),GAPFUN(rsm.rsmas_ekman_flux_x_volume),GAPFUN(fwyf1.ndbc_ekman_flux_x_volume));
    xlim_datetick(xl);
    ylim([-4.9,+4.9]);
    %annotate_ts_event(upwMr,face.deep.seatemp,upwmT,-4.9);
    annotate_ts_event(upwMr,face.deep.seatemp,upwmT,-inf);
    center_axes; grid on; grid minor;
    legh=legend('Cape Can. (''88-''16)','Lake Worth (''01-''16)','Port Ev. (''09-''17)','RSMAS Roof (''03-''11)','Fowey (''91-''16)','(Upw.)', 'Orientation','horiz','Location','SouthEast'); set(legh,'FontSize',19)
    panel_label_subplot(gca,'d','ur','FontSize',fontsz,'BackgroundColor','w');
    set(gca,'FontSize',fontsz); set(gca,'LineWidth',gridwid); set(get(gca,'Child'),'LineWidth',linewid);

    ekax(2) = subplot_tight(7,1,4:6);
    title('Along-shore Ekman volumetric flux [m^2s^-^1]');
    lh=plot_ts(GAPFUN(canf1.ndbc_ekman_flux_l_volume),GAPFUN(lkwf1.ndbc_ekman_flux_l_volume),GAPFUN(pvgf1.ndbc_ekman_flux_l_volume),GAPFUN(rsm.rsmas_ekman_flux_l_volume),GAPFUN(fwyf1.ndbc_ekman_flux_l_volume));
    xlim_datetick(xl);
    ylim([-4.9,+4.9]);
    %annotate_ts_event(upwMr,face.deep.seatemp,upwmT,-4.9);
    annotate_ts_event(upwMr,face.deep.seatemp,upwmT,-inf);
    center_axes; grid on; grid minor;
    panel_label_subplot(gca,'e','ur','FontSize',fontsz,'BackgroundColor','w');
    set(gca,'FontSize',fontsz); set(gca,'LineWidth',gridwid); set(get(gca,'Child'),'LineWidth',linewid);

    ekax(3) = subplot_tight(7,1,7);
    title('Barometric pressure [mbar]');
    lh=plot_ts(GAPFUN(canf1.ndbc_barom),GAPFUN(lkwf1.ndbc_barom),GAPFUN(pvgf1.ndbc_barom),GAPFUN(rsm.rsmas_barom),GAPFUN(fwyf1.ndbc_barom));
    xlim_datetick('Y',xl);
    %ylim([980,1040]);
    ylim([1006,1023])
    %annotate_ts_event(upwMr,face.deep.seatemp,upwmT,1006);
    annotate_ts_event(upwMr,face.deep.seatemp,upwmT,-inf);
    grid on; grid minor;
    panel_label_subplot(gca,'f','ur','FontSize',fontsz,'BackgroundColor','w');
    set(gca,'FontSize',fontsz); set(gca,'LineWidth',gridwid); set(get(gca,'Child'),'LineWidth',linewid);


    % Variability partition
    prfh = fmg;

    prax(1) = subplot_tight(7,1,1:3);
    titlename('Sea temperature variability [\circC]: FACE deep mooring (Sep 2015)');
    title('Sea temperature variability [\circC]');
    lh=plot_ts(GAPFUN(face.deep.seatemp_3_h_hp_9_h_lp),...
               GAPFUN(face.deep.seatemp_11_h_hp_14_h_lp),...
               GAPFUN(face.deep.seatemp_17_h_hp_28_h_lp),...
               GAPFUN(face.deep.seatemp_48_h_hp_240_h_lp));
    xlim_datetick(xl);
    center_axes;
    annotate_ts_event(upwMr,face.deep.seatemp,upwmT,-inf);
    legend('3-9h','11-14h','17-28h','2-10d','(Upw.)', 'Orientation','horiz','Location','SouthEast');
    grid on; grid minor;
    panel_label_subplot(gca,'g','ur','FontSize',fontsz,'BackgroundColor','w');
    set(gca,'FontSize',fontsz); set(gca,'LineWidth',gridwid); set(get(gca,'Child'),'LineWidth',linewid);

    prax(2) = subplot_tight(7,1,4:5);
    title('Cross-shore near-bottom current var. [ms^-^1]');
    lh=plot_ts(GAPFUN(face.deep.adcp_btm_x_3_h_hp_9_h_lp),...
               GAPFUN(face.deep.adcp_btm_x_11_h_hp_14_h_lp),...
               GAPFUN(face.deep.adcp_btm_x_17_h_hp_28_h_lp),...
               GAPFUN(face.deep.adcp_btm_x_48_h_hp_240_h_lp));
    xlim_datetick(xl); %datetick3;
    %center_axes;
    ylim([-0.1,+0.1]);
    annotate_ts_event(upwMr,face.deep.seatemp,upwmT,-inf);
    %legend('3-9h','11-14h','17-28h','2-10d','(Upw.)', 'Orientation','horiz');
    grid on; grid minor;
    panel_label_subplot(gca,'h','ur','FontSize',fontsz,'BackgroundColor','w');
    set(gca,'FontSize',fontsz); set(gca,'LineWidth',gridwid); set(get(gca,'Child'),'LineWidth',linewid);

    prax(3) = subplot_tight(7,1,6:7);
    title('Along-shore near-bottom current var. [ms^-^1]');
    lh=plot_ts(GAPFUN(face.deep.adcp_btm_l_3_h_hp_9_h_lp),...
               GAPFUN(face.deep.adcp_btm_l_11_h_hp_14_h_lp),...
               GAPFUN(face.deep.adcp_btm_l_17_h_hp_28_h_lp),...
               GAPFUN(face.deep.adcp_btm_l_48_h_hp_240_h_lp));
    %center_axes;
    ylim([-0.3,+0.3]);
    annotate_ts_event(upwMr,face.deep.seatemp,upwmT,-inf);
    %legend('3-9h','11-14h','17-28h','2-10d','(Upw.)', 'Orientation','horiz');
    grid on; grid minor;
    panel_label_subplot(gca,'i','ur','FontSize',fontsz,'BackgroundColor','w');
    set(gca,'FontSize',fontsz); set(gca,'LineWidth',gridwid); set(get(gca,'Child'),'LineWidth',linewid);
    xlim_datetick('Y',xl);

    annot_h = textbox_range([dpax,shax,ekax,prax],datenum(2015,9,[7,18],12,0,0),'I');
    annot_h = textbox_range([dpax,shax,ekax,prax],datenum(2015,9,[22,25],12,0,0),'II');

    if doPrint; axes(dpax(1)); print('-dpng',fullfile(figspath,'upwelling-face-2015-Sep.png')); end;
    if doPrint; axes(shax(1)); print('-dpng',fullfile(figspath,'upwelling-face-2015-Sep-upper.png')); end;
    if doPrint; axes(ekax(1)); print('-dpng',fullfile(figspath,'upwelling-face-2015-Sep-ekman.png')); end;
    if doPrint; axes(prax(1)); print('-dpng',fullfile(figspath,'upwelling-face-2015-Sep-partn.png')); end;
  end;

end; %if ( doFACE )


if ( doNF09 )

  fmg;
  titlename('FACE 2009: High frequency variability at Broward');
  lh=plot_ts(face.brwd100.brwd100_upp.seatemp,face.brwd40.brwd40_btm.seatemp,face.brwdbe.brwdbenth5.seatemp,face.brwdbe.brwdbenth4.seatemp,face.brwd20.brwd20_btm.seatemp,fwyf1.ndbc_air_t,'k-');
  legend(lh,'100m @38 m','40m @38 m','Be @33 m','Be @32 m','20m @19 m','Air', 'Location','SouthWest');
  %xlim(face.brwd100.brwd100_upp.seatemp.date([1,end]));
  axis([datenum(2009,[10,11],[27,8]),22.0,29.0]);
  datetick3('x','mmm-dd','keeplimits');
  set(gca,'FontSize',fontsz); set(gca,'LineWidth',gridwid); set(get(gca,'Child'),'LineWidth',linewid);
  if doPrint; print('-dpng',fullfile(figspath,'upwelling-nf09-brwd-variability.png')); end;

  fmg;
  titlename('FACE 2009: Onshore propagation at Broward');
  lh=plot_ts(face.brwd100.brwd100_upp.seatemp,face.brwd40.brwd40_btm.seatemp,face.brwdbe.brwdbenth5.seatemp,face.brwdbe.brwdbenth4.seatemp,face.brwd20.brwd20_btm.seatemp,fwyf1.ndbc_air_t,'k-');
  legend(lh,'100m @38 m','40m @38 m','Be @33 m','Be @32 m','20m @19 m','Air', 'Location','SouthWest');
  %xlim(face.brwd100.brwd100_upp.seatemp.date([1,end]));
  axis([datenum(2009,10,30,13,0,0),datenum(2009,10,30,23,0,0),24.0,29.0]);
  datetick3('x','mmm-dd HH:MM','keeplimits');
  set(gca,'FontSize',fontsz); set(gca,'LineWidth',gridwid); set(get(gca,'Child'),'LineWidth',linewid);
  annotation('textarrow','String',sprintf('940 m  \n 25 cm/s  '),'Position',[0.1466 0.6044 0.0884 0.1882],'TextRotation',45,'LineWidth',1,'FontSize',fontsz);
  annotation('textarrow','String',sprintf('320 m  \n 30 cm/s  '),'Position',[0.2473 0.7815 0.0268 0.0610],'TextRotation',45,'LineWidth',1,'FontSize',fontsz);
  annotation('textarrow','String',sprintf('940 m  \n 13 cm/s  '),'Position',[0.7014 0.2488 0.1786 0.5123],'TextRotation',45,'LineWidth',1,'FontSize',fontsz);
  if doPrint; print('-dpng',fullfile(figspath,'upwelling-nf09-brwd-onshore-Oct-30.png')); end;

  fmg;
  titlename('FACE 2009: Alongshore propagation from Broward to Boca');
  lh=plot_ts(face.brwd100.brwd100_low.seatemp,face.brwd100.brwd100_upp.seatemp,face.boca40.boca40_btm.seatemp,fwyf1.ndbc_air_t,'k-');
  legend(lh,'Broward @66 m','Broward @38 m','Boca @41 m','Air', 'Location','SouthWest');
  %xlim(face.brwd100.brwd100_upp.seatemp.date([1,end]));
  axis([datenum(2009,[10,11],[20,8]),20.0,29.0]);
  datetick3('x','mmm-dd','keeplimits');
  set(gca,'FontSize',fontsz); set(gca,'LineWidth',gridwid); set(get(gca,'Child'),'LineWidth',linewid);
  if doPrint; print('-dpng',fullfile(figspath,'upwelling-nf09-alongshore.png')); end;

end; %if ( doNF09 )


if ( doHanS )
  fmg; plot_ts(face.AOML_HW.adcp_seatemp,face.AOML_BR.adcp_seatemp,face.HanS_HW.adcp_seatemp,face.HanS_BR.adcp_seatemp); legend('AOML HW 7.3 m','AOML BR 8.2 m','H&S HW 26.2 m','H&S BR 32.4 m'); set(get(gca,'Child'),'LineWidth',2);
end; %if ( doHanS )


if ( doNutPlots )

  fmg; plot_ts(sfomc.ne_buoy.intNx,sfomc.sw_buoy.intNx,sfomc.nw_w_btm.intNx); legend('Outer','Mid','Inshore','Location','East');
  axis([datenum(1999,[7,8],[15,29]),0,8]); datetick('x',2,'keeplimits'); %datetick3;
  ylabel('\int^0_-_h (u_{onshore}\bullet[NO2+NO3])dz   [kg\bulletm^-^1]'); set(gca,'FontSize',fontsz);
  titlename('Accumulated nitrogen per meter of reef');
  if doPrint; print('-dpng',fullfile(figspath,['upwelling-sfomc-1999-N_cumkg.png'])); end;

  fmg; plot_ts(sfomc.ne_buoy.intPx,sfomc.sw_buoy.intPx,sfomc.nw_w_btm.intPx); legend('Outer','Mid','Inshore','Location','East');
  axis([datenum(1999,[7,8],[15,29]),0,0.8]); datetick('x',2,'keeplimits'); %datetick3;
  ylabel('\int^0_-_h (u_{onshore}\bullet[P])dz   [kg\bulletm^-^1]'); set(gca,'FontSize',fontsz);
  titlename('Accumulated phosphorous per meter of reef');
  if doPrint; print('-dpng',fullfile(figspath,['upwelling-sfomc-1999-P_cumkg.png'])); end;

  fmg; plot_ts(sfomc.ne_buoy.intSix,sfomc.sw_buoy.intSix,sfomc.nw_w_btm.intSix); legend('Outer','Mid','Inshore','Location','East');
  axis([datenum(1999,[7,8],[15,29]),0,8]); datetick('x',2,'keeplimits'); %datetick3;
  ylabel('\int^0_-_h (u_{onshore}\bullet[Si])dz   [kg\bulletm^-^1]'); set(gca,'FontSize',fontsz);
  titlename('Accumulated silicon per meter of reef');
  if doPrint; print('-dpng',fullfile(figspath,['upwelling-sfomc-1999-Si_cumkg.png'])); end;


  fmg;
  spt(1,3,1);
  plot_ts(sfomc.ne_buoy.intNx,sfomc.sw_buoy.intNx,sfomc.nw_w_btm.intNx); legend('Outer','Mid','Inshore','Location','East');
  axis([datenum(1999,[7,8],[15,29]),0,8]); xticks('auto'); datetick('x',6,'keeplimits','keepticks'); %datetick3;
  ylabel('\int^0_-_h (u_{onshore}\bullet[NO2+NO3])dz   [kg\bulletm^-^1]'); set(gca,'FontSize',16);
  title('Accumulated N m^-^1 (1999)');

  spt(1,3,2);
  plot_ts(sfomc.ne_buoy.intPx,sfomc.sw_buoy.intPx,sfomc.nw_w_btm.intPx); legend('Outer','Mid','Inshore','Location','East');
  axis([datenum(1999,[7,8],[15,29]),0,8]); xticks('auto'); datetick('x',6,'keeplimits','keepticks'); %datetick3;
  ylabel('\int^0_-_h (u_{onshore}\bullet[P])dz   [kg\bulletm^-^1]'); set(gca,'FontSize',16);
  title('Accumulated P m^-^1');

  spt(1,3,3);
  plot_ts(sfomc.ne_buoy.intSix,sfomc.sw_buoy.intSix,sfomc.nw_w_btm.intSix); legend('Outer','Mid','Inshore','Location','East');
  axis([datenum(1999,[7,8],[15,29]),0,8]); xticks('auto'); datetick('x',6,'keeplimits','keepticks'); %datetick3;
  ylabel('\int^0_-_h (u_{onshore}\bullet[Si])dz   [kg\bulletm^-^1]'); set(gca,'FontSize',16);
  title('Accumulated Si m^-^1');

  %th=suptitlename('Accumulated nutrients per meter of reef (1999)');
  %set(th,'FontSize',fontsz);
  if doPrint; print('-dpng',fullfile(figspath,['upwelling-sfomc-1999-nuts_cumkg.png'])); end;

end; %if ( doNutPlots )


if ( doNutMaps )

  % Max accumulation of N [kg/m] at 50 m isobath from summer 1999
  N = max(sfomc.ne_buoy.intNx.data(sfomc.ne_buoy.intNx.date<datenum(1999,11,1)));

  % Create field of nitrogen accumulation mirroring hires bathymetry field
  Nfld = lkwf1.ngdc_hires_bathy.field;
  Nfld(-50>=Nfld | Nfld>=0) = nan;
  % % (OLD) Assume constant cross-shore mixing rate: peak at 50 m isobath, zero at shore
  % Nfld = (1 - ((-50-Nfld)/(-50)) ).*N;
  % Use cross-shore mixing rate estimated from sea temperature troughs
  Nfld = (1 - ((-50-Nfld)/(-50)) ).*N;

  [Nlon,Nlat] = meshgrid(lkwf1.ngdc_hires_bathy.lon,lkwf1.ngdc_hires_bathy.lat);
  % % Attenuate with latitude from a peak at 27.23N
  % %(27.23-26.67)*111 * 0.2
  % Nlat = (27.69 - Nlat)./(27.69 - 25.0);
  % Nfld = Nfld.*(1-Nlat);
  % Attenuate with latitude from estimated value (above) for NE buoy latitude, 26.0695N
  Nlat = (Nlat - 25.00)/(26.0695 - 25.000);
  Nfld = Nfld.*Nlat;

  fmg;
  contourf(lkwf1.ngdc_hires_bathy.lon,lkwf1.ngdc_hires_bathy.lat,Nfld);
  cbh=colorbar;
  plot_hires_coastline(lkwf1.ngdc_hires_bathy);
  %axis([-80.20,-79.9,25.54,27.25]);
  axis([-80.20,-79.95,25.80,27.25]);
  xlabel(cbh,'kg\bulletN\bulletm^-^1yr^-^1');
  if doPrint; print('-dpng',fullfile(figspath,'upwelling-N-flux-map.png')); end;
  pause;
  daspect([1,5,1]);
  if doPrint; print('-dpng',fullfile(figspath,'upwelling-N-flux-map-exaggerated.png')); end;

  % Bathymetry is 10 m resolution, max along dim 2 is peak at 50 m isobath
  min(diff(lkwf1.ngdc_hires_bathy.lat))*111e3,

  %[NOx]
  nansum(nanmax(Nfld,[],2))*10,
  %[P]
  nansum(nanmax(Nfld,[],2))*10 * max(sfomc.ne_buoy.intPx.data(sfomc.ne_buoy.intPx.date<datenum(1999,11,1))) / max(sfomc.ne_buoy.intNx.data(sfomc.ne_buoy.intNx.date<datenum(1999,11,1))),
  %[P]
  nansum(nanmax(Nfld,[],2))*10 * max(sfomc.ne_buoy.intSix.data(sfomc.ne_buoy.intSix.date<datenum(1999,11,1))) / max(sfomc.ne_buoy.intNx.data(sfomc.ne_buoy.intNx.date<datenum(1999,11,1))),

end; %if ( doNutMaps )


if 0;
  fmg; spt(2,1,1); plot(get_monthdayhour(sfomc.nw_w_btm.adcp_seatemp_3_h_hp_9_h_lp.date),sfomc.nw_w_btm.adcp_seatemp_3_h_hp_9_h_lp.data.^3,'r.',get_monthdayhour(sfomc.nw_w_btm.adcp_seatemp_11_h_hp_14_h_lp.date),sfomc.nw_w_btm.adcp_seatemp_11_h_hp_14_h_lp.data.^3,'b.'); grid on; grid minor; legend('3-9h','11-14h', 'Location','southwest'); spt(2,1,2); plot(get_monthdayhour(sfomc.nw_w_btm.adcp_x_btm_3_h_hp_9_h_lp.date),sfomc.nw_w_btm.adcp_x_btm_3_h_hp_9_h_lp.data.^3,'r.',get_monthdayhour(sfomc.nw_w_btm.adcp_x_btm_11_h_hp_14_h_lp.date),sfomc.nw_w_btm.adcp_x_btm_11_h_hp_14_h_lp.data.^3,'b.'); datetick3('x',6); xlim(datenum(0,[5,11],1)); grid on; grid minor;
end;

set_more;

