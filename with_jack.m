1;

legargs = { 'Orientation','horiz','Location','southeast' };

% Could maybe try Lanczos 3- or 6-hlp? Or Butterworth 14-hlp?
% But that would mask these weird events of high-frequency variability!

dpfh = fmg;
dpax(1) = subplot(7,1,1:3); title('Sea temperatures');
plot_ts(jack.deep.seatemp,jack.tcm2.seatemp,jack.tcm1.seatemp,jack.shallow.seatemp);
ylim([19.5,31.5]);
grid on; set(gca,'XTickLabel','');
legend('Deep','TCM2','TCM1','Shallow', legargs{:});

dpax(2) = subplot(7,1,4:5); title('Deep alongshore current');
plot_ts(jack.deep.adcp_v_sfc,jack.deep.adcp_v_mid,jack.deep.adcp_v_btm);
ylim([-1.5,+1.5]); xlim(jack.deep.seatemp.date([1,end])); 
grid on; center_axes; set(gca,'XTickLabel','');
legend('Surface','Mid-column','Bottom', legargs{:});

dpax(3) = subplot(7,1,6:7); title('Deep cross-shore current');
plot_ts(jack.deep.adcp_u_sfc,jack.deep.adcp_u_mid,jack.deep.adcp_u_btm);
ylim([-0.4,+0.4]); xlim(jack.deep.seatemp.date([1,end]));
grid on; center_axes;
legend('Surface','Mid-column','Bottom', legargs{:});


shfh = fmg;
shax(1) = subplot(7,1,1:3); title('Sea temperatures');
plot_ts(pvgf1.ndbc_air_t,'k-',jack.deep.seatemp,jack.tcm2.seatemp,jack.tcm1.seatemp,jack.shallow.seatemp);
ylim([19.5,31.5]); xlim(jack.deep.seatemp.date([1,end]));
grid on; set(gca,'XTickLabel','');
legend('Air','Deep','TCM2','TCM1','Shallow', legargs{:});

shax(2) = subplot(7,1,4:5); title('Shallow alongshore current');
plot_ts(jack.shallow.adcp_v_sfc,jack.shallow.adcp_v_mid,jack.shallow.adcp_v_btm);
ylim([-1.5,+1.5]); xlim(jack.deep.seatemp.date([1,end]));
grid on; center_axes; set(gca,'XTickLabel','');
legend('Surface','Mid-column','Bottom', legargs{:});

shax(3) = subplot(7,1,6:7); title('Shallow cross-shore current');
plot_ts(jack.shallow.adcp_u_sfc,jack.shallow.adcp_u_mid,jack.shallow.adcp_u_btm);
ylim([-0.5,+0.5]); xlim(jack.deep.seatemp.date([1,end]));
grid on; center_axes;
legend('Surface','Mid-column','Bottom', legargs{:});


ekfh = fmg;
plot_ts(lkwf1.ndbc_ekman_flux_volume,pvgf1.ndbc_ekman_flux_volume);
ylim([-3.0,+3.0]); xlim(jack.deep.seatemp.date([1,end]));
titlename('Ekman volumetric flux near FACE HWD section');
legend('Lake Worth','Port Everglades','ERA-Interim', legargs{:});


disp('Whole records');
pause;


disp('Broader focus on periods of high-frequency variability');
figure(dpfh); xlim(datenum(2015,[6,8],[10,29])); datetick3;
figure(shfh); xlim(datenum(2015,[6,8],[10,29])); datetick3;
figure(ekfh); xlim(datenum(2015,[6,8],[10,29])); datetick3;
pause;

disp('CLOSE focus on period of high-frequency variability');
figure(dpfh); xlim(datenum(2015,7,[10,19])); datetick3;
figure(shfh); xlim(datenum(2015,7,[10,19])); datetick3;
figure(ekfh); xlim(datenum(2015,7,[10,19])); datetick3;
pause;

disp('Whole period when all four current meters were operating');
figure(dpfh); xlim(datenum(2015,[9,11],[18,12])); datetick3;
figure(shfh); xlim(datenum(2015,[9,11],[18,12])); datetick3;
figure(ekfh); xlim(datenum(2015,[9,11],[18,12])); datetick3;
pause;

disp('Upwelling event when all four current meters were operating');
figure(dpfh); xlim(datenum(2015,9,[18,29])); datetick3;
figure(shfh); xlim(datenum(2015,9,[18,29])); datetick3;
figure(ekfh); xlim(datenum(2015,9,[18,29])); datetick3;
pause;

disp('Done');
pause;


% Cross- vs. along-shore currents

fmg;
subplot(7,1,1:3); title('Sea temperatures');
plot_ts(pvgf1.ndbc_air_t,'k-',jack.deep.seatemp,jack.tcm2.seatemp,jack.tcm1.seatemp,jack.shallow.seatemp);
legend('Air','Deep','TCM2','TCM1','Shallow', legargs{:});

subplot(7,1,4:5); title('Alongshore');
plot_ts(jack.deep.adcp_v,jack.tcm2.v,jack.tcm1.v,jack.shallow.adcp_v);
legend('Deep','TCM2','TCM1','Shallow', legargs{:});

subplot(7,1,6:7); title('Cross-shore');
plot_ts(jack.deep.adcp_u,jack.tcm2.u,jack.tcm1.u,jack.shallow.adcp_u);

% Broader focus on periods of high-frequency variability
xlim(datenum(2015,7,[9,29])); datetick3;




% For analysis of counter-current and possible wastewater onshore transport (winters)
jack.deep = verify_variable(jack.deep,{'adcp_u_36_h_sum','adcp_v_36_h_sum'});
jack.tcm2 = verify_variable(jack.tcm2,{'u_36_h_sum','v_36_h_sum'});
jack.tcm1 = verify_variable(jack.tcm1,{'u_36_h_sum','v_36_h_sum'});
jack.shallow = verify_variable(jack.shallow,{'adcp_u_36_h_sum','adcp_v_36_h_sum'});

fmg;
subplot(7,1,1:3); title('Sea temperature');
plot_ts(pvgf1.ndbc_air_t,'k-',jack.deep.seatemp,jack.tcm2.seatemp,jack.tcm1.seatemp,jack.shallow.seatemp);
legend('Air','Deep','TCM2','TCM1','Shallow', legargs{:});

subplot(7,1,4:5); title('Alongshore');
plot_ts(jack.deep.adcp_v_36_h_sum,jack.tcm2.v_36_h_sum,jack.tcm1.v_36_h_sum,jack.shallow.adcp_v_36_h_sum);
legend('Deep','TCM2','TCM1','Shallow', legargs{:});

subplot(7,1,6:7); title('Cross-shore');
plot_ts(jack.deep.adcp_u_36_h_sum,jack.tcm2.u_36_h_sum,jack.tcm1.u_36_h_sum,jack.shallow.adcp_u_36_h_sum);

% Focus on period of cold-front passage and counter-current
xlim(datenum(2015,10,[13,32])); datetick3;
