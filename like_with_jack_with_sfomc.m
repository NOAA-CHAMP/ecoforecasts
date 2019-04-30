1;

error('THIS SCRIPT HAS NOW BEEN RENAMED ekman_flux_vs_seatemp.m');

% fmg; plot_ts(fwyf1.ndbc_ekman_flux_volume,rsm.rsmas_ekman_flux_volume,lkwf1.ndbc_ekman_flux_volume,fdepk.fdep_ekman_flux_volume,canf1.ndbc_ekman_flux_volume); xlim(datenum(2007,[7,9],[5,5])); datetick3;
% pause;


if ( ~exist('doPrint','var') || isempty(doPrint) )
  doPrint = false;
end;
if ( ~exist('figspath','var') || isempty(figspath) )
  figspath = fullfile(get_coral_path,'CRCP','Upwelling','CoRIS');
end;
if ( doPrint )
  disp(['Will print figures to ',figspath]);
else
  disp('Figures will not be printed');
end;

legargs = { 'Orientation','horiz','Location','southeast' };

% Some day, may be worth analyzing with Lanczos 3- or 6-hlp? Or Butterworth
% 14-hlp? But that would mask these weird high-frequency variability events!

dpfh = fmg;
dpax(1) = subplot(7,1,1:3); title('Sea temperatures vs. deep site currents');
plot_ts(fwyf1.ndbc_air_t,'k-',sfomc.se_buoy.seatemp,sfomc.ne_buoy.seatemp,sfomc.sw_buoy.seatemp,sfomc.nw_w_btm.adcp_seatemp);
ylim([19.5,31.5]);
grid on; set(gca,'XTickLabel','');
legend('Air','Offshore (SE)','Deep (NE)','Mid (SW)','Shallow (NW)', legargs{:});

dpax(2) = subplot(7,1,4:5); ylabel('Alongshore'); %title('Deep alongshore current');
plot_ts(sfomc.ne_buoy.adcp_l_uwc,sfomc.ne_buoy.adcp_l_btm);
ylim([-1.5,+1.5]); xlim(sfomc.ne_buoy.seatemp.date([1,end])); 
grid on; center_axes; set(gca,'XTickLabel','');
legend('Along Upper Water Column (NE)','Along Near Seafloor (NE)', legargs{:});

dpax(3) = subplot(7,1,6:7); ylabel('Cross-shore'); %title('Deep cross-shore current');
plot_ts(sfomc.ne_buoy.adcp_x_uwc,sfomc.ne_buoy.adcp_x_btm);
ylim([-0.4,+0.4]); xlim(sfomc.ne_buoy.seatemp.date([1,end]));
grid on; center_axes;
legend('Cross Upper Water Column (NE)','Cross Near Seafloor (NE)', legargs{:});
datetick3;


shfh = fmg;
shax(1) = subplot(7,1,1:3); title('Sea temperatures vs. shallow site currents');
plot_ts(fwyf1.ndbc_air_t,'k-',sfomc.se_buoy.seatemp,sfomc.ne_buoy.seatemp,sfomc.sw_buoy.seatemp,sfomc.nw_w_btm.adcp_seatemp);
ylim([19.5,31.5]); xlim(sfomc.ne_buoy.seatemp.date([1,end]));
grid on; set(gca,'XTickLabel','');
legend('Air','Offshore (SE)','Deep (NE)','Mid (SW)','Shallow (NW)', legargs{:});

shax(2) = subplot(7,1,4:5); ylabel('Alongshore'); %title('Shallow alongshore current');
plot_ts(sfomc.nw_w_btm.adcp_l_uwc,sfomc.nw_w_btm.adcp_l_btm);
ylim([-1.5,+1.5]); xlim(sfomc.ne_buoy.seatemp.date([1,end]));
grid on; center_axes; set(gca,'XTickLabel','');
legend('Along Upper Water Column (NW)','Along Near Seafloor (NW)', legargs{:});

shax(3) = subplot(7,1,6:7); ylabel('Cross-shore'); %title('Shallow cross-shore current');
plot_ts(sfomc.nw_w_btm.adcp_x_uwc,sfomc.nw_w_btm.adcp_x_btm);
ylim([-0.4,+0.4]); xlim(sfomc.ne_buoy.seatemp.date([1,end]));
grid on; center_axes;
legend('Cross Upper Water Column (NW)','Cross Near Seafloor (NW)', legargs{:});
datetick3;


ekfh = fmg;
plot_ts(lkwf1.ndbc_ekman_flux_volume,pvgf1.ndbc_ekman_flux_volume,rsm.rsmas_ekman_flux_volume,fwyf1.ndbc_ekman_flux_volume);
ylim([-3.0,+3.0]); xlim(sfomc.ne_buoy.seatemp.date([1,end]));
titlename('Ekman volumetric flux near SFOMC section');
legend('Lake Worth (2001-2016)','Port Everglades (2009-2017)','RSMAS Rooftop (2003-2011)','Fowey Rocks (1991-2016)', legargs{:});
datetick3;


if ( doPrint )
  figure(dpfh);
  print('-dpng',fullfile(figspath,'upwelling-sfomc-1999-2000-deep-whole.png'));
  figure(shfh);
  print('-dpng',fullfile(figspath,'upwelling-sfomc-1999-2000-shallow-whole.png'));
  figure(ekfh);
  print('-dpng',fullfile(figspath,'upwelling-sfomc-1999-2000-ekman-whole.png'));
end;
disp('Whole records');
pause;


disp('Focus on high-frequency variability in 5-min sea temperature in *2007*');
figure(dpfh); xlim(datenum(2007,[7,9],[5,5])); datetick3;
if doPrint; print('-dpng',fullfile(figspath,'upwelling-sfomc-1999-2000-deep-wide.png')); end;
figure(shfh); xlim(datenum(2007,[7,9],[5,5])); datetick3;
if doPrint; print('-dpng',fullfile(figspath,'upwelling-sfomc-1999-2000-shallow-wide.png')); end;
figure(ekfh); xlim(datenum(2007,[7,9],[5,5])); datetick3;
if doPrint; print('-dpng',fullfile(figspath,'upwelling-sfomc-1999-2000-ekman-wide.png')); end;
keyboard;
pause;


disp('Broader focus on periods of high-frequency variability');
figure(dpfh); xlim(datenum(1999,[7,8],[15,29])); datetick3;
if doPrint; print('-dpng',fullfile(figspath,'upwelling-sfomc-1999-2000-deep-Jul-Aug.png')); end;
figure(shfh); xlim(datenum(1999,[7,8],[15,29])); datetick3;
if doPrint; print('-dpng',fullfile(figspath,'upwelling-sfomc-1999-2000-shallow-Jul-Aug.png')); end;
figure(ekfh); xlim(datenum(1999,[7,8],[15,29])); datetick3;
if doPrint; print('-dpng',fullfile(figspath,'upwelling-sfomc-1999-2000-ekman-Jul-Aug.png')); end;
pause;

disp('CLOSE focus on period of high-frequency variability');
figure(dpfh); xlim(datenum(1999,7,[28,32])); datetick3;
if doPrint; print('-dpng',fullfile(figspath,'upwelling-sfomc-1999-2000-deep-Jul.png')); end;
figure(shfh); xlim(datenum(1999,7,[28,32])); datetick3;
if doPrint; print('-dpng',fullfile(figspath,'upwelling-sfomc-1999-2000-shallow-Jul.png')); end;
figure(ekfh); xlim(datenum(1999,7,[28,32])); datetick3;
if doPrint; print('-dpng',fullfile(figspath,'upwelling-sfomc-1999-2000-ekman-Jul.png')); end;
pause;

% disp('Whole period when all current meters were operating');
% figure(dpfh); xlim(datenum(1999,[9,11],[18,12])); datetick3;
% figure(shfh); xlim(datenum(1999,[9,11],[18,12])); datetick3;
% figure(ekfh); xlim(datenum(1999,[9,11],[18,12])); datetick3;
% pause;

disp('Upwelling event when three current meters were operating');
figure(dpfh); xlim(datenum(1999,7,[22,28])); datetick3;
figure(shfh); xlim(datenum(1999,7,[22,28])); datetick3;
figure(ekfh); xlim(datenum(1999,7,[22,28])); datetick3;
pause;

disp('Done');
pause;


% Cross-shore section: cross- vs. along-shore currents

fmg;
subplot(7,1,1:3); title('Section of sea temperatures and currents (TRANSIENTS)');
plot_ts(fwyf1.ndbc_air_t,'k-',sfomc.se_buoy.seatemp,sfomc.ne_buoy.seatemp,sfomc.sw_buoy.seatemp,sfomc.nw_w_btm.adcp_seatemp);
ylim([19.5,31.5]);
grid on;
legend('Air (FWY)','Offshore (SE)','Deep (NE)','Mid (SW)','Shallow (NW)', legargs{:});

subplot(7,1,4:5); ylabel('Alongshore'); %title('Alongshore');
plot_ts(sfomc.ne_buoy.adcp_v,sfomc.sw_buoy.adcp_v,sfomc.nw_w_btm.adcp_v);
ylim([-1.5,+1.5]);
legend('Deep (NE)','Mid (SW)','Shallow (NW)', legargs{:});

subplot(7,1,6:7); ylabel('Cross-shore'); %title('Cross-shore');
plot_ts(sfomc.ne_buoy.adcp_u,sfomc.sw_buoy.adcp_u,sfomc.nw_w_btm.adcp_u);
ylim([-0.4,+0.4]);

% Broader focus on periods of high-frequency variability
xlim(datenum(1999,7,[9,29])); datetick3;

if doPrint; print('-dpng',fullfile(figspath,'upwelling-sfomc-1999-2000-cross-vs-alongshore.png')); end;



% For analysis of counter-current and possible wastewater onshore transport (winters)
sfomc.ne_buoy = verify_variable(sfomc.ne_buoy,{'adcp_u_36_h_sum','adcp_v_36_h_sum'});
sfomc.sw_buoy = verify_variable(sfomc.sw_buoy,{'adcp_u_36_h_sum','adcp_v_36_h_sum'});
sfomc.nw_w_btm = verify_variable(sfomc.nw_w_btm,{'adcp_u_36_h_sum','adcp_v_36_h_sum'});

fmg;
subplot(7,1,1:3); title('Section of sea temperatures and currents (COLD FRONT)');
plot_ts(fwyf1.ndbc_air_t,'k-',sfomc.se_buoy.seatemp,sfomc.ne_buoy.seatemp,sfomc.sw_buoy.seatemp,sfomc.nw_w_btm.adcp_seatemp);
ylim([19.5,31.5]);
grid on;
legend('Air (FWY)','Offshore (SE)','Deep (NE)','Mid (SW)','Shallow (NW)', legargs{:});

subplot(7,1,4:5); ylabel('Alongshore'); %title('Alongshore');
plot_ts(sfomc.ne_buoy.adcp_v_36_h_sum,sfomc.sw_buoy.adcp_v_36_h_sum,sfomc.nw_w_btm.adcp_v_36_h_sum);
grid on;
legend('\Sigma_3_6_h Deep (NE)','\Sigma_3_6_h Mid (SW)','\Sigma_3_6_h Shallow (NW)', legargs{:});

subplot(7,1,6:7); ylabel('Cross-shore'); %title('Cross-shore');
plot_ts(sfomc.ne_buoy.adcp_u_36_h_sum,sfomc.sw_buoy.adcp_u_36_h_sum,sfomc.nw_w_btm.adcp_u_36_h_sum);
grid on;

% Focus on period of cold-front passage and counter-current
xlim(datenum(1999,10,[13,32])); datetick3;

if doPrint; print('-dpng',fullfile(figspath,'upwelling-sfomc-1999-2000-cross-vs-alongshore-36_h_sum.png')); end;

% clear ans dpax dpfh ekfh legargs shax shfh
