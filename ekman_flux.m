1;
error('THIS SCRIPT HAS BEEN RENAMED upwelling_ekman_flux.m');
% SCRIPT ekman_flux.m
% Analyze in situ sea temperature extremes in summer vs. alongshore wind
% stress and nominal coastal Ekman divergence (i.e., examines feasibility of
% explaining upwelling in terms of wind) on the SE Florida shelf. Relies on
% data from NOAA-AOML's CHAMP and FACE Programs, RSMAS rooftop station, NOAA
% NDBC, Florida DEP, SEFCRI/FWC/FDEP, and NSUOC/SFOMC. Produces figures in
% the FIGSPATH directory (DEFAULT: GET_CORAL_PATH('CRCP/Upwelling/CoRIS')).
%
% Last Saved Time-stamp: <Fri 2017-06-09 10:17:09 Eastern Daylight Time gramer>

if ( ~exist('doPlots','var') || isempty(doPlots) )
  doPlots = false;
  %doPlots = true;
end;
if ( ~exist('doPrint','var') || isempty(doPrint) )
  doPrint = false;
  %doPrint = true;
end;
if ( doPrint && ~exist('figspath','var') )
  figspath = get_coral_path('CRCP/Upwelling/CoRIS');
end;

% Coast-line data
if ( ~exist('fwyf1','var') )
  fwyf1 = get_station_from_station_name('fwyf1'); fwyf1 = load_all_ndbc_data(fwyf1);
end;
if ( ~exist('rsm','var') )
  rsm = read_rsmas_weatherpack_data;
end;
if ( ~exist('lkwf1','var') )
  lkwf1 = get_station_from_station_name('lkwf1'); lkwf1 = load_all_ndbc_data(lkwf1);
end;
if ( ~exist('pvgf1','var') )
  pvgf1 = get_station_from_station_name('pvgf1'); pvgf1 = load_all_ndbc_data(pvgf1);
end;
if ( ~exist('ftpf1','var') )
  ftpf1 = get_station_from_station_name('41114'); ftpf1 = load_all_ndbc_data(ftpf1);
end;
if ( ~exist('fdepk','var') )
  fdepk = read_fdep_stevens_data('k');
end;
if ( ~exist('cnnf1','var') )
  cnnf1 = get_station_from_station_name('41113'); cnnf1 = load_all_ndbc_data(cnnf1);
end;
if ( ~exist('canf1','var') )
  canf1 = get_station_from_station_name('41009'); canf1 = load_all_ndbc_data(canf1);
end;

% Reef-line data - just temperature
if ( ~exist('sefcri','var') )
  sefcri = load(get_coral_path('CRCP/Upwelling/SEFCRI.mat'));
end;

% Reef-line data - currents
if ( ~exist('jack','var') )
  jack = get_face_jack_currents;
end;
if ( ~exist('sfomc','var') )
  [ig,h]= system('hostname');
  if ( strncmpi(h,'manannan',length('manannan')) && ~exist('fullsfomc','var') )
    fullsfomc = true;
  end;
  if ( exist('fullsfomc','var') && fullsfomc )
    sfomc = get_sfomc_data;
  else
    sfomc = get_sfomc_trim_data;
  end;
end;


% Bathymetry
if ( ~isfield(lkwf1,'ngdc_hires_bathy') )
  lkwf1 = read_hires_bathymetry(lkwf1,[30e3,120e3]);
end;
if ( ~isfield(sefcri.dc3,'ngdc_hires_bathy') )
  sefcri.dc3 = read_hires_bathymetry(sefcri.dc3,[30e3,50e3],[],false);
end;
if ( exist('sfomc','var') && ~isfield(sfomc.c_buoy,'ngdc_hires_bathy') )
  % NOTE: Reads HIGH-resolution (10 m) bathymetry - keep the bbox tight!
  sfomc.c_buoy = read_hires_bathymetry(sfomc.c_buoy,[3e3,11e3],[],true);
end;
if ( ~isfield(fdepk,'ngdc_hires_bathy') )
  fdepk = read_hires_bathymetry(fdepk,[30e3,50e3]);
end;
if ( ~isfield(canf1,'ngdc_hires_bathy') )
  canf1 = read_hires_bathymetry(canf1,[30e3,50e3]);
end;

if ( ~isfield(lkwf1,'slope') )
  [lkwf1.slope,lkwf1.slope_orientation,lkwf1.isobath_orientation,lkwf1.ngdc_hires_bathy] = ...
      find_ngdc_slope(lkwf1.ngdc_hires_bathy,lkwf1.lon,lkwf1.lat,7);
end;
if ( ~isfield(fwyf1,'slope') )
  [fwyf1.slope,fwyf1.slope_orientation,fwyf1.isobath_orientation,lkwf1.ngdc_hires_bathy] = ...
      find_ngdc_slope(lkwf1.ngdc_hires_bathy,fwyf1.lon,fwyf1.lat,7);
end;
if ( ~isfield(rsm,'slope') )
  [rsm.slope,rsm.slope_orientation,rsm.isobath_orientation,lkwf1.ngdc_hires_bathy] = ...
      find_ngdc_slope(lkwf1.ngdc_hires_bathy,fwyf1.lon,fwyf1.lat,7);
end;
if ( ~isfield(pvgf1,'slope') )
  [pvgf1.slope,pvgf1.slope_orientation,pvgf1.isobath_orientation,lkwf1.ngdc_hires_bathy] = ...
      find_ngdc_slope(lkwf1.ngdc_hires_bathy,pvgf1.lon,pvgf1.lat,7);
end;
if ( ~isfield(ftpf1,'slope') )
  [ftpf1.slope,ftpf1.slope_orientation,ftpf1.isobath_orientation,lkwf1.ngdc_hires_bathy] = ...
      find_ngdc_slope(lkwf1.ngdc_hires_bathy,ftpf1.lon,ftpf1.lat,7);
end;
if ( ~isfield(fdepk,'slope') )
  [fdepk.slope,fdepk.slope_orientation,fdepk.isobath_orientation,lkwf1.ngdc_hires_bathy] = ...
      find_ngdc_slope(lkwf1.ngdc_hires_bathy,fdepk.lon,fdepk.lat,7);
end;
if ( ~isfield(canf1,'slope') )
  [canf1.slope,canf1.slope_orientation,canf1.isobath_orientation,canf1.ngdc_hires_bathy] = ...
      find_ngdc_slope(canf1.ngdc_hires_bathy,canf1.lon,canf1.lat,7);
end;
if ( ~isfield(cnnf1,'slope') )
  [cnnf1.slope,cnnf1.slope_orientation,cnnf1.isobath_orientation,canf1.ngdc_hires_bathy] = ...
      find_ngdc_slope(canf1.ngdc_hires_bathy,cnnf1.lon,cnnf1.lat,7);
end;

if ( ~isfield(jack.shallow,'slope') )
  [jack.shallow.slope,jack.shallow.slope_orientation,jack.shallow.isobath_orientation,lkwf1.ngdc_hires_bathy] = ...
      find_ngdc_slope(lkwf1.ngdc_hires_bathy,jack.shallow.lon,jack.shallow.lat,7);
  [jack.tcm1.slope,jack.tcm1.slope_orientation,jack.tcm1.isobath_orientation,lkwf1.ngdc_hires_bathy] = ...
      find_ngdc_slope(lkwf1.ngdc_hires_bathy,jack.tcm1.lon,jack.tcm1.lat,7);
  [jack.tcm2.slope,jack.tcm2.slope_orientation,jack.tcm2.isobath_orientation,lkwf1.ngdc_hires_bathy] = ...
      find_ngdc_slope(lkwf1.ngdc_hires_bathy,jack.tcm2.lon,jack.tcm2.lat,7);
  [jack.deep.slope,jack.deep.slope_orientation,jack.deep.isobath_orientation,lkwf1.ngdc_hires_bathy] = ...
      find_ngdc_slope(lkwf1.ngdc_hires_bathy,jack.deep.lon,jack.deep.lat,7);
end;

if ( exist('sfomc','var') && ~isfield(sfomc.nw_w_btm,'slope') )
  [sfomc.nw_w_btm.slope,sfomc.nw_w_btm.slope_orientation,sfomc.nw_w_btm.isobath_orientation,lkwf1.ngdc_hires_bathy] = ...
      find_ngdc_slope(lkwf1.ngdc_hires_bathy,sfomc.nw_w_btm.lon,sfomc.nw_w_btm.lat,7);
  [sfomc.ne_buoy.slope,sfomc.ne_buoy.slope_orientation,sfomc.ne_buoy.isobath_orientation,lkwf1.ngdc_hires_bathy] = ...
      find_ngdc_slope(lkwf1.ngdc_hires_bathy,sfomc.ne_buoy.lon,sfomc.ne_buoy.lat,7);
  [sfomc.e_buoy.slope,sfomc.e_buoy.slope_orientation,sfomc.e_buoy.isobath_orientation,lkwf1.ngdc_hires_bathy] = ...
      find_ngdc_slope(lkwf1.ngdc_hires_bathy,sfomc.e_buoy.lon,sfomc.e_buoy.lat,7);
  [sfomc.c_buoy.slope,sfomc.c_buoy.slope_orientation,sfomc.c_buoy.isobath_orientation,lkwf1.ngdc_hires_bathy] = ...
      find_ngdc_slope(lkwf1.ngdc_hires_bathy,sfomc.c_buoy.lon,sfomc.c_buoy.lat,7);
  [sfomc.sw_buoy.slope,sfomc.sw_buoy.slope_orientation,sfomc.sw_buoy.isobath_orientation,lkwf1.ngdc_hires_bathy] = ...
      find_ngdc_slope(lkwf1.ngdc_hires_bathy,sfomc.sw_buoy.lon,sfomc.sw_buoy.lat,7);
end;


if ( ~isfield(fwyf1,'ndbc_wind1_x') )
  fwyf1 = station_spddir_to_uv(fwyf1,'ndbc_wind1_speed','ndbc_wind1_dir','ndbc_wind1_u','ndbc_wind1_v');
  rsm = station_spddir_to_uv(rsm,  'rsmas_wind1_speed','rsmas_wind1_dir','rsmas_wind1_u','rsmas_wind1_v');
  pvgf1 = station_spddir_to_uv(pvgf1,'ndbc_wind1_speed','ndbc_wind1_dir','ndbc_wind1_u','ndbc_wind1_v');
  lkwf1 = station_spddir_to_uv(lkwf1,'ndbc_wind1_speed','ndbc_wind1_dir','ndbc_wind1_u','ndbc_wind1_v');
  fdepk = station_spddir_to_uv(fdepk,'fdep_wind_speed','fdep_wind_dir','fdep_wind_u','fdep_wind_v');
  ftpf1 = station_spddir_to_uv(ftpf1,'ndbc_sigwavehgt','ndbc_avgwavedir','ndbc_wave_u','ndbc_wave_v');
  cnnf1 = station_spddir_to_uv(cnnf1,'ndbc_sigwavehgt','ndbc_avgwavedir','ndbc_wave_u','ndbc_wave_v');
  canf1 = station_spddir_to_uv(canf1,'ndbc_wind1_speed','ndbc_wind1_dir','ndbc_wind1_u','ndbc_wind1_v');

  fwyf1 = station_reorient_vectors(fwyf1,'isobath_orientation','ndbc_wind1_u','ndbc_wind1_v','ndbc_wind1_x','ndbc_wind1_l');
  rsm = station_reorient_vectors(rsm,'isobath_orientation','rsmas_wind1_u','rsmas_wind1_v','rsmas_wind1_x','rsmas_wind1_l');
  pvgf1 = station_reorient_vectors(pvgf1,'isobath_orientation','ndbc_wind1_u','ndbc_wind1_v','ndbc_wind1_x','ndbc_wind1_l');
  lkwf1 = station_reorient_vectors(lkwf1,'isobath_orientation','ndbc_wind1_u','ndbc_wind1_v','ndbc_wind1_x','ndbc_wind1_l');
  fdepk = station_reorient_vectors(fdepk,'isobath_orientation','fdep_wind_u','fdep_wind_v','fdep_wind_x','fdep_wind_l');
  ftpf1 = station_reorient_vectors(ftpf1,'isobath_orientation','ndbc_wave_u','ndbc_wave_v','ndbc_wave_x','ndbc_wave_l');
  cnnf1 = station_reorient_vectors(cnnf1,'isobath_orientation','ndbc_wave_u','ndbc_wave_v','ndbc_wave_x','ndbc_wave_l');
  canf1 = station_reorient_vectors(canf1,'isobath_orientation','ndbc_wind1_u','ndbc_wind1_v','ndbc_wind1_x','ndbc_wind1_l');
end;

if ( ~isfield(jack.shallow,'adcp_x') )
  jack.shallow = station_reorient_vectors(jack.shallow,'isobath_orientation','adcp_u','adcp_v','adcp_x','adcp_l');
  jack.tcm1 = station_reorient_vectors(jack.tcm1,'isobath_orientation','u','v','x','l');
  jack.tcm2 = station_reorient_vectors(jack.tcm2,'isobath_orientation','u','v','x','l');
  jack.deep = station_reorient_vectors(jack.deep,'isobath_orientation','adcp_u','adcp_v','adcp_x','adcp_l');
end;


if ( ~isfield(fwyf1,'ndbc_ekman_flux') )
  disp('Calculating Ekman fluxes...');
  fwyf1 = station_ekman_flux(fwyf1,[],'ndbc','erai');
  pvgf1 = station_ekman_flux(pvgf1,[],'ndbc');
  lkwf1 = station_ekman_flux(lkwf1,[],'ndbc');
  fdepk = station_ekman_flux(fdepk,[],'fdep');
  canf1 = station_ekman_flux(canf1,[],'ndbc');

  rsm.lon = fwyf1.lon; rsm.lat = fwyf1.lat;
  rsm.ndbc_sea_t = fwyf1.ndbc_sea_t;
  rsm = station_ekman_flux(rsm,[],'rsmas');
end;

for hrs = 12:12:36
  for cderf={'lp','sum','avg'};
    fldstub = sprintf('_%d_h_%s',hrs,cderf{:});
    fld = ['ndbc_ekman_flux_volume',fldstub];
    disp(fld);
    fwyf1 = verify_variable(fwyf1,fld);
    rsm = verify_variable(rsm,fld);
    pvgf1 = verify_variable(pvgf1,fld);
    lkwf1 = verify_variable(lkwf1,fld);
    canf1 = verify_variable(canf1,fld);
    jack.deep = verify_variable(jack.deep,['adcp_x_',fldstub]);
    jack.deep = verify_variable(jack.deep,['adcp_l_',fldstub]);
    jack.tcm2 = verify_variable(jack.tcm2,['x_',fldstub]);
    jack.tcm2 = verify_variable(jack.tcm2,['l_',fldstub]);
    jack.tcm1 = verify_variable(jack.tcm1,['x_',fldstub]);
    jack.tcm1 = verify_variable(jack.tcm1,['l_',fldstub]);
    jack.shallow = verify_variable(jack.shallow,['adcp_x_',fldstub]);
    jack.shallow = verify_variable(jack.shallow,['adcp_l_',fldstub]);
  end;
end;


if ( doPlots )

  fhs=[];
  [lkwf1,ig,ig,ig,fhs(1)] = plot_hires_bathymetry(lkwf1,-[0:2:20,30:20:200]);
  titlename('SE Florida shelf - USGS Bathymetry - instrumentation');

  [sefcri.dc3,ig,ig,ig,fhs(end+1)] = plot_hires_bathymetry(sefcri.dc3,-[0:2:20,30:10:50]);
  titlename('SE Florida shelf (southern) - USGS Bathymetry - instrumentation');

  [fdepk,ig,ig,ig,fhs(end+1)] = plot_hires_bathymetry(fdepk,-[0:2:20,30:10:50]);
  titlename('SE Florida shelf (northern) - USGS Bathymetry - instrumentation');

  %[canf1,ig,ig,ig,fhs(end+1)] = plot_hires_bathymetry(canf1,-[0:2:20,30:10:50]);
  %titlename('Central Florida shelf - USGS Bathymetry - instrumentation');

  fhix = 1;
  figure(fhs(fhix));
  plot_all_se_florida_upwelling_sites;
  if doPrint; print('-dpng',fullfile(figspath,['upwelling-bathymetry.png'])); end;

  for fhix = 2:numel(fhs)
    figure(fhs(fhix));
    plot_all_se_florida_upwelling_sites;
    if doPrint; print('-dpng',fullfile(figspath,['upwelling-bathymetry-',num2str(fhix),'.png'])); end;
  end; %for fhix = 2:numel(fhs)

  clear fhix

  % High-resolution map of Alex's SFOMC section
  if ( exist('sfomc','var') )
    sfomc.c_buoy = plot_hires_bathymetry(sfomc.c_buoy,-[0:1:30,40:10:120],[3e3,11e3]);
    plot_all_se_florida_upwelling_sites;
    if doPrint; print('-dpng',fullfile(figspath,['upwelling-bathymetry-sfomc.png'])); end;
  end;

end; %if ( doPlots )


if ( doPlots )

  if 1;
    fmg;
    subplot(3,1,1); title('Far North (Martin County)');
     plot_ts(sefcri.updb.hourly_t,sefcri.mc2.hourly_t); 
     ylim([12,32]);
     legend(['Offshore (',num2str(sefcri.updb.depth),' m)'],...
            ['Near-shore (',num2str(sefcri.mc2.depth),' m)'],'Location','SouthWest');
     grid on;
    subplot(3,1,2); title('Northern (Palm Beach)');
     plot_ts(sefcri.pb2.hourly_t,sefcri.pb1.hourly_t);
     ylim([12,32]);
     legend(['Offshore (',num2str(sefcri.pb2.depth),' m)'],...
            ['Near-shore (',num2str(sefcri.pb1.depth),' m)'],'Location','SouthWest');
     grid on;
    subplot(3,1,3); title('Central (Ft. Lauderdale)');
     plot_ts(sefcri.bc3.hourly_t,sefcri.bc1.hourly_t);
     ylim([12,32]);
     legend(['Offshore (',num2str(sefcri.bc3.depth),' m)'],...
            ['Near-shore (',num2str(sefcri.bc1.depth),' m)'],'Location','SouthWest');
     grid on;
    if doPrint; print('-dpng',fullfile(figspath,'upwelling-sefcri-2010-2013.png')); end;
  end;

  if 1;
    fmg;
    %spt(2,1,1);
    subplot(2,1,1);
    shax = plotyy_ts(jack.shallow.seatemp,jack.shallow.adcp_x);
    legend('11m Ts','11m u','Location','SouthWest');
    grid on;
    titlename('FACE HWD moorings - upwelling and cross-shore currents');
    %spt(2,1,2);
    subplot(2,1,2);
    dpax = plotyy_ts(jack.deep.seatemp,jack.deep.adcp_x);
    legend('27m Ts','27m u','Location','SouthWest');
    grid on;
    if doPrint; print('-dpng',fullfile(figspath,'upwelling-jack.png')); end;
  end;

  if 1;
    fmg;
    [ax,lh(1),lh(2)] = plotyy_ts(jack.deep.adcp_x,jack.deep.seatemp);
    lh(3) = plot(ax(2),jack.shallow.seatemp.date,jack.shallow.seatemp.data,'-','Color',[0,0.5,0]);
    xlim(datenum(2015,[1,10],1)); datetick3;
    ylim(ax(1),[-0.35,+0.35]);
    ylim(ax(2),[20,32]);
    legend(lh,'27m u','27m T','11m T','Location','SouthWest');
    grid on;
    titlename('FACE HWD moorings - upwelling and cross-shore currents');
    if doPrint; print('-dpng',fullfile(figspath,'upwelling-jack-2015.png')); end;
  end;

  if 1;
    fmg;
    subplot(7,1,1:5);
    plot_ts(jack.deep.seatemp,jack.tcm2.seatemp,jack.tcm1.seatemp,jack.shallow.seatemp);
    ylim([20,32]);
    xlim(datenum(2015,[8,10],1));
    legend('Deep','TCM2','TCM1','Shallow','Location','SouthEast');
    grid on;
    titlename('FACE HWD moorings - upwelling and cross-shore currents');

    subplot(7,1,6:7);
    plot_ts(jack.deep.adcp_l,jack.deep.adcp_x,jack.tcm2.x,jack.tcm1.x,jack.shallow.adcp_x);
    %ylim([-0.35,+0.35]);
    %ylim([-0.60,+0.60]);
    ylim([-1.10,+1.10]);
    xlim(datenum(2015,[8,10],1)); datetick3;
    legend('Deep ALONG','Deep cross','TCM2 cross','TCM1 cross','Shallow cross','Location','SouthEast');
    grid on;

    if doPrint; print('-dpng',fullfile(figspath,'upwelling-jack-2015-onshore.png')); end;
  end;

  if 1;
    fmg;
    subplot(7,1,1:3);
    plot_ts(jack.deep.seatemp,jack.tcm2.seatemp,jack.tcm1.seatemp,jack.shallow.seatemp);
    legend('Deep','TCM2','TCM1','Shallow','Location','SouthEast');
    ylim([27,31]);
    ylabel('T [^oC]');
    xlim(datenum(2015,9,[18,31])); datetick3;
    set(get(gca,'XAxis'),'Vis','off');
    grid on;
    titlename('FACE HWD moorings - upwelling and cross-shore currents');

    subplot(7,1,4:5);
    plot_ts(jack.deep.adcp_l,jack.tcm2.l,jack.tcm1.l,jack.shallow.adcp_l);
    legend('Deep','TCM2','TCM1','Shallow','Location','SouthEast');
    ylim([-0.80,+0.80]);
    ylabel('v [ms^-^1]');
    xlim(datenum(2015,9,[18,31])); datetick3;
    set(get(gca,'XAxis'),'Vis','off');
    grid on;

    subplot(7,1,6:7);
    plot_ts(jack.deep.adcp_x,jack.tcm2.x,jack.tcm1.x,jack.shallow.adcp_x);
    legend('Deep','TCM2','TCM1','Shallow','Location','SouthEast');
    ylim([-0.20,+0.20]);
    ylabel('u [ms^-^1]');
    xlim(datenum(2015,9,[18,31])); datetick3;
    grid on;

    if doPrint; print('-dpng',fullfile(figspath,'upwelling-jack-2015-onshore-Sep.png')); end;
  end;

end;

clear ans ax cderf dpax fhs fld fldstub hrs ig lh shax
