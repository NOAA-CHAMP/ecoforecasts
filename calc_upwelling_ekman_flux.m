1;
% SCRIPT calc_upwelling_ekman_flux.m
%
% (FORMERLY just "ekman_flux.m")
%
% Analyze in situ sea temperature extremes in summer vs. alongshore wind
% stress and nominal coastal Ekman divergence (i.e., examines feasibility of
% explaining upwelling in terms of wind) on the SE Florida shelf. Relies on
% data from NOAA-AOML's CHAMP and FACE Programs, RSMAS rooftop station, NOAA
% NDBC, Florida DEP, SEFCRI/FWC/FDEP, and NSUOC/SFOMC. Produces figures in
% the FIGSPATH directory (DEFAULT: GET_CORAL_PATH('CRCP/Upwelling/CoRIS')).
%
% PARAMETERS (DEFAULTS):
%   doCSV - Create CSV file with DUMP_UPWELLING_CSV.m (false)
%   doNuts - Estimate nutrient concentrations and fluxes (true)
%   doRi - Estimate Nyqvist and shear frequencies and Richardson number (true)
%   doGrid - Created gridded dataset from all (ship or moored) station data (true)
%   doMaps - Plot bathymetric maps showing station locations (false, CAN BE SLOW!)
%   doPrint - Print bathymetric maps (false)
%   figspath - Directory in which to plot figures (default above)
%   datapath - Directory to save data files to (same as FIGSPATH)
%   fontsz - Default font point size for plotted figures (24)
%
% Last Saved Time-stamp: <Mon 2018-07-09 12:05:54 Eastern Daylight Time gramer>

if ( ~exist('doCSV','var') || isempty(doCSV) )
  doCSV = false;
end;
if ( ~exist('doNuts','var') || isempty(doNuts) )
  doNuts = true;
end;
if ( ~exist('doRi','var') || isempty(doRi) )
  doRi = true;
end;
if ( ~exist('doGrid','var') || isempty(doGrid) )
  %doGrid = false;
  doGrid = true;
end;
if ( ~exist('doMaps','var') || isempty(doMaps) )
  doMaps = false;
  %doMaps = true;
end;
if ( ~exist('doPrint','var') || isempty(doPrint) )
  doPrint = false;
  %doPrint = true;
end;
if ( ~exist('figspath','var') )
  figspath = fullfile(get_coral_path,'CRCP','Upwelling','CoRIS');
end;
if ( ~exist('datapath','var') )
  datapath = fullfile(get_coral_path,'CRCP','Upwelling','CoRIS');
end;
if ( doMaps )
  if ( doPrint )
    disp(['Will print figures to ',figspath]);
  else
    disp('Figures will not be printed');
  end;
end;

upwpath = get_coral_path('CRCP/Upwelling');

if ( ~exist('fontsz','var') || isempty(fontsz) )
  %%fontsz = 16;
  %fontsz = 20;
  fontsz = 24;
end;


% NEW DATA!!

%{
error('D:\coral\FACE\BR-HW\Nuts_Master-Data-Sheet');
error('Add to and analyze with d:/coral/FACE/analyze_face_nutrients.m');

error('Analyze D:\coral\FACE\Partners\Broward');

% Cell 1	1	7/1/2010 20:00	0.13
error('d:/coral/CRCP/Upwelling/NSUOC/NW_W-Calypso/Data/AWAC/01Jul10_26Jul10/010710_northvelocity.dat');
error('04Jun10_13Jun10');
error('13Jun10_30Jun10');

fid = fopen(fullfile(upwpath,'NSUOC','NW_W-Calypso','Data','AWAC','01Jul10_26Jul10','010710_northvelocity.dat'),'r');
C = fscanf(fid,'Cell %d\t%d\t%s %s\t%g\n');
fclose(fid);

fwc = load(fullfile(upwpath,'fwc_fdep.mat'));
fmg;
plot_ts(fwc.UPDB.sea_t,fwc.TUCA.sea_t);
legend(sprintf('UPDB %g m, %g ^oN',fwc.UPDB.depth,fwc.UPDB.lat),...
       sprintf('TUCA %g m, %g ^oN',fwc.TUCA.depth,fwc.TUCA.lat));
%}



% Coast-line data
if ( ~exist('fwyf1','var') )
  fwyf1 = get_station_from_station_name('fwyf1'); fwyf1 = load_all_ndbc_data(fwyf1);
end;
if ( ~exist('rsm','var') )
  rsm = read_rsmas_weatherpack_data;
  % RSMAS rooftop STRUCT has no coordinates, no sea temperature
  rsm.lon = fwyf1.lon;  rsm.lat = fwyf1.lat;
  %%%% COMMENTED OUT: 2018 Apr 15 - allow volume flux to use default sea temp.!
  %%%% rsm.ndbc_sea_t = fwyf1.ndbc_sea_t;

  % RSMAS rooftop also has some funky negative wind speeds!
  rsm.rsmas_wind1_speed.data(rsm.rsmas_wind1_speed.data<0) = 0;
  rsm.rsmas_wind2_speed.data(rsm.rsmas_wind2_speed.data<0) = 0;
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

% Straits of Florida data - nutrients and temperatures (sub-thermocline empirical relationships)
if ( ~exist('sfp','var') )
  sfp = load_SFP_nutrients;
end;

% Reef-line data - just temperature
if ( ~exist('sefcri','var') )
  sefcri = load(fullfile(upwpath,'SEFCRI.mat'));

  sefcri.station_names = {...
      'sefcri_dc1','sefcri_dc3','sefcri_bc1','sefcri_bc3', ...
      'sefcri_pb1','sefcri_pb2','sefcri_mc2', 'sefcri_pela','sefcri_evca', ...
      'sefcri_slib','sefcri_updb', ...
                   };
  sefcri.lats = [sefcri.dc1.lat,sefcri.dc3.lat,sefcri.bc1.lat,sefcri.bc3.lat, ...
                 sefcri.pb1.lat,sefcri.pb2.lat,sefcri.mc2.lat, sefcri.pela.lat,sefcri.evca.lat, ...
                 sefcri.slib.lat, sefcri.updb.lat, ...
              ];
  sefcri.lons = [sefcri.dc1.lon,sefcri.dc3.lon,sefcri.bc1.lon,sefcri.bc3.lon, ...
                 sefcri.pb1.lon,sefcri.pb2.lon,sefcri.mc2.lon, sefcri.pela.lon,sefcri.evca.lon, ...
                 sefcri.slib.lon, sefcri.updb.lon, ...
              ];
end;
if ( ~exist('fwc','var') )
  fwc = load(fullfile(upwpath,'fwc_fdep.mat'));

  fwc.station_names = {'fwc_UPDB','fwc_TUCA','fwc_Noula','fwc_Alpha','fwc_Rodeo'};
  fwc.lats = [fwc.UPDB.lat,fwc.TUCA.lat,fwc.Noula.lat,fwc.Alpha.lat,fwc.Rodeo.lat];
  fwc.lons = [fwc.UPDB.lon,fwc.TUCA.lon,fwc.Noula.lon,fwc.Alpha.lon,fwc.Rodeo.lon];
end;

% Reef-line data - currents
if ( ~exist('face','var') )
  face = get_face_jack_currents;
  % Eddy09 experiment moorings
  face = read_NF09(face);
  % Nancy Foster Cruises NF08 and NF09
  face = load_FACE_nutrients(face);

  % Outfall ADCPs of Hazen & Sawyer and "AOML" (FACE)
  face = get_FACE_Hazen_currents(face);

  face.station_names = {...
      'face_shallow','face_tcm1','face_tcm2','face_deep', ...
      'face_brwd20','face_brwd40','face_brwd100', 'face_boca20','face_boca40', ...
      'face_AOML_HW','face_HanS_HW','face_AOML_BR','face_HanS_BR', ...
                   };
  face.lats = [face.shallow.lat,face.tcm1.lat,face.tcm2.lat,face.deep.lat, ...
               face.brwd20.lat,face.brwd40.lat,face.brwd100.lat, face.boca20.lat,face.boca40.lat, ...
               face.AOML_HW.lat, face.HanS_HW.lat, face.AOML_BR.lat, face.HanS_BR.lat, ...
              ];
  face.lons = [face.shallow.lon,face.tcm1.lon,face.tcm2.lon,face.deep.lon, ...
               face.brwd20.lon,face.brwd40.lon,face.brwd100.lon, face.boca20.lon,face.boca40.lon, ...
               face.AOML_HW.lon, face.HanS_HW.lon, face.AOML_BR.lon, face.HanS_BR.lon, ...
              ];
end;

if ( ~exist('sfomc','var') )
  if ( ~exist('fullsfomc','var') )
    [ig,h]= system('hostname');
    if ( ~strncmpi(h,'titan',length('titan')) )
      fullsfomc = true;
    end;
  end;
  if ( fullsfomc )
    sfomc = get_sfomc_data;
  else
    sfomc = get_sfomc_trim_data;
  end;

  sfomc.station_names = cellstr(char(sfomc.mc_sites{1:4:end}));
  sfomc.metadata = [sfomc.mc_sites{2:4:end}];
  sfomc.lats = sfomc.metadata(1:4:end);
  sfomc.lons = sfomc.metadata(2:4:end);

  if ( ~isfield(sfomc.nw_w_btm,'adcp_speed_uwc') )
    sfomc.nw_w_btm.adcp_speed_uwc.date = sfomc.nw_w_btm.adcp_u_uwc.date;
    sfomc.nw_w_btm.adcp_speed_uwc.data = uv_to_spd(sfomc.nw_w_btm.adcp_u_uwc.data,sfomc.nw_w_btm.adcp_v_uwc.data);
  end;
  if ( ~isfield(sfomc.nw_w_btm,'adcp_speed_btm') )
    sfomc.nw_w_btm.adcp_speed_btm.date = sfomc.nw_w_btm.adcp_u_btm.date;
    sfomc.nw_w_btm.adcp_speed_btm.data = uv_to_spd(sfomc.nw_w_btm.adcp_u_btm.data,sfomc.nw_w_btm.adcp_v_btm.data);
  end;

  [ig,ix11m] = min(abs((-sfomc.ne_buoy.adcp_depths)-11));
  [ig,ix30m] = min(abs((-sfomc.ne_buoy.adcp_depths)-30));
  [ig,ix45m] = min(abs((-sfomc.ne_buoy.adcp_depths)-45));
  
  sfomc.ne_buoy.adcp_x_11m.date = sfomc.ne_buoy.adcp_x.date;
  sfomc.ne_buoy.adcp_x_11m.data = sfomc.ne_buoy.adcp_x.prof(:,ix11m);
  sfomc.ne_buoy.adcp_x_30m.date = sfomc.ne_buoy.adcp_x.date;
  sfomc.ne_buoy.adcp_x_30m.data = sfomc.ne_buoy.adcp_x.prof(:,ix30m);
  sfomc.ne_buoy.adcp_x_45m.date = sfomc.ne_buoy.adcp_x.date;
  sfomc.ne_buoy.adcp_x_45m.data = sfomc.ne_buoy.adcp_x.prof(:,ix45m);

  if ( ~isfield(sfomc.ne_buoy,'salin') )
    [c0,c5,c10,c15] = ...
        intersect_tses([],sfomc.c_buoy.at_0_6m.sbe_salin,sfomc.c_buoy.at_5m.sbe_salin,sfomc.c_buoy.at_10m.sbe_salin,sfomc.c_buoy.at_15m.sbe_salin);
    sfomc.c_buoy.salin.date = c0.date;
    sfomc.c_buoy.salin.prof = [c0.data,c5.data,c10.data,c15.data];
    sfomc.c_buoy.salin.data = nanmean(sfomc.c_buoy.salin.prof,2);
    sfomc.c_buoy.salin.depths = -[0.6,5,10,15];
    clear c0 c5 c10 c15
    
    [sw0,sw5,sw10,sw15] = ...
        intersect_tses([],sfomc.sw_buoy.at_0_6m.sbe_salin,sfomc.sw_buoy.at_5m.sbe_salin,sfomc.sw_buoy.at_10m.sbe_salin,sfomc.sw_buoy.at_15m.sbe_salin);
    sfomc.sw_buoy.salin.date = sw0.date;
    sfomc.sw_buoy.salin.prof = [sw0.data,sw5.data,sw10.data,sw15.data];
    sfomc.sw_buoy.salin.data = nanmean(sfomc.sw_buoy.salin.prof,2);
    sfomc.sw_buoy.salin.depths = -[0.6,5,10,15];
    clear sw0 sw5 sw10 sw15
    
    [ne0,ne5,ne10,ne15,ne20,ne30,ne40,ne45] = ...
        intersect_tses([],sfomc.ne_buoy.at_0_6m.sbe_salin,sfomc.ne_buoy.at_5m.sbe_salin,sfomc.ne_buoy.at_10m.sbe_salin,sfomc.ne_buoy.at_15m.sbe_salin,sfomc.ne_buoy.at_20m.sbe_salin,sfomc.ne_buoy.at_30m.sbe_salin,sfomc.ne_buoy.at_40m.sbe_salin,sfomc.ne_buoy.at_45m.sbe_salin);
    sfomc.ne_buoy.salin.date = ne0.date;
    sfomc.ne_buoy.salin.prof = [ne0.data,ne5.data,ne10.data,ne15.data,ne20.data,ne30.data,ne40.data,ne45.data];
    sfomc.ne_buoy.salin.data = nanmean(sfomc.ne_buoy.salin.prof,2);
    sfomc.ne_buoy.salin.depths = -[0.6,5,10,15,20,30,40,45];
    clear ne0 ne5 ne10 ne15 ne20 ne30 ne40 ne45
    
    [e5,e10,e15,e20,e30] = ...
        intersect_tses([],sfomc.e_buoy.at_5m.sbe_salin,sfomc.e_buoy.at_10m.sbe_salin,sfomc.e_buoy.at_15m.sbe_salin,sfomc.e_buoy.at_20m.sbe_salin,sfomc.e_buoy.at_30m.sbe_salin);
    sfomc.e_buoy.salin.date = e5.date;
    sfomc.e_buoy.salin.prof = [e5.data,e10.data,e15.data,e20.data,e30.data];
    sfomc.e_buoy.salin.data = nanmean(sfomc.e_buoy.salin.prof,2);
    sfomc.e_buoy.salin.depths = -[5,10,15,20,30];
    clear e0 e5 e10 e15 e20 e30
    
    [se26,se35,se98] = ...
        intersect_tses([],sfomc.se_buoy.at_26m.sbe_salin,sfomc.se_buoy.at_35m.sbe_salin,sfomc.se_buoy.at_98m.sbe_salin);
    sfomc.se_buoy.salin.date = se26.date;
    sfomc.se_buoy.salin.prof = [se26.data,se35.data,se98.data];
    sfomc.se_buoy.salin.data = nanmean(sfomc.se_buoy.salin.prof,2);
    sfomc.se_buoy.salin.depths = -[26,35,98];
    clear se26 se35 se98
  end;

  if ( ~isfield(sfomc.nw_w_btm,'adcp_seatemp_hourly') )
    ws = warning('OFF','MATLAB:interp1:NaNstrip');

    sfomc.nw_w_btm.adcp_seatemp_hourly = ts_nanify_gaps(interp_ts(sfomc.nw_w_btm.adcp_seatemp),10);

    sfomc.nw_w_btm.adcp_speed_hourly = ts_nanify_gaps(interp_ts(sfomc.nw_w_btm.adcp_speed),10);
    sfomc.nw_w_btm.adcp_x_hourly = ts_nanify_gaps(interp_ts(sfomc.nw_w_btm.adcp_x),10);
    sfomc.nw_w_btm.adcp_l_hourly = ts_nanify_gaps(interp_ts(sfomc.nw_w_btm.adcp_l),10);
    sfomc.nw_w_btm.adcp_dir_hourly.date = sfomc.nw_w_btm.adcp_x_hourly.date;
    sfomc.nw_w_btm.adcp_dir_hourly.data = uv_to_dir_curr(sfomc.nw_w_btm.adcp_x_hourly.data,sfomc.nw_w_btm.adcp_l_hourly.data);

    sfomc.nw_w_btm.adcp_speed_uwc_hourly = ts_nanify_gaps(interp_ts(sfomc.nw_w_btm.adcp_speed_uwc),10);
    sfomc.nw_w_btm.adcp_x_uwc_hourly = ts_nanify_gaps(interp_ts(sfomc.nw_w_btm.adcp_x_uwc),10);
    sfomc.nw_w_btm.adcp_l_uwc_hourly = ts_nanify_gaps(interp_ts(sfomc.nw_w_btm.adcp_l_uwc),10);
    sfomc.nw_w_btm.adcp_dir_uwc_hourly.date = sfomc.nw_w_btm.adcp_x_uwc_hourly.date;
    sfomc.nw_w_btm.adcp_dir_uwc_hourly.data = uv_to_dir_curr(sfomc.nw_w_btm.adcp_x_uwc_hourly.data,sfomc.nw_w_btm.adcp_l_uwc_hourly.data);

    sfomc.nw_w_btm.adcp_speed_btm_hourly = ts_nanify_gaps(interp_ts(sfomc.nw_w_btm.adcp_speed_btm),10);
    sfomc.nw_w_btm.adcp_x_btm_hourly = ts_nanify_gaps(interp_ts(sfomc.nw_w_btm.adcp_x_btm),10);
    sfomc.nw_w_btm.adcp_l_btm_hourly = ts_nanify_gaps(interp_ts(sfomc.nw_w_btm.adcp_l_btm),10);
    sfomc.nw_w_btm.adcp_dir_btm_hourly.date = sfomc.nw_w_btm.adcp_x_btm_hourly.date;
    sfomc.nw_w_btm.adcp_dir_btm_hourly.data = uv_to_dir_curr(sfomc.nw_w_btm.adcp_x_btm_hourly.data,sfomc.nw_w_btm.adcp_l_btm_hourly.data);

    warning(ws); clear ws
  end;

end;


% Bathymetry
if ( ~isfield(lkwf1,'ngdc_hires_bathy') )
  lkwf1 = read_hires_bathymetry(lkwf1,[30e3,120e3]);
end;
if ( exist('sefcri','var') && ~isfield(sefcri.dc3,'ngdc_hires_bathy') )
  sefcri.dc3 = read_hires_bathymetry(sefcri.dc3,[30e3,50e3],[],false);
end;
if ( exist('sfomc','var') && ~isfield(sfomc.c_buoy,'ngdc_hires_bathy') )
  % NOTE: Reads HIGH-resolution (10 m) bathymetry - keep the bbox tight!
  sfomc.c_buoy = read_hires_bathymetry(sfomc.c_buoy,[3e3,11e3],[],true);
end;
if ( ~isfield(fdepk,'ngdc_hires_bathy') )
  fdepk = read_hires_bathymetry(fdepk,[30e3,50e3]);
end;
if ( ~isfield(ftpf1,'ngdc_hires_bathy') )
  ftpf1 = read_hires_bathymetry(ftpf1,[6e3,10e3]);
end;
if ( ~isfield(canf1,'ngdc_hires_bathy') )
  canf1 = read_hires_bathymetry(canf1,[30e3,50e3]);
end;

if ( ~isfield(lkwf1,'slope') )
  [lkwf1.slope,lkwf1.slope_orientation,lkwf1.isobath_orientation,lkwf1.ngdc_hires_bathy] = ...
      find_ngdc_slope(lkwf1.ngdc_hires_bathy,lkwf1.lon,lkwf1.lat,3);
end;
if ( ~isfield(fwyf1,'slope') )
  [fwyf1.slope,fwyf1.slope_orientation,fwyf1.isobath_orientation,lkwf1.ngdc_hires_bathy] = ...
      find_ngdc_slope(lkwf1.ngdc_hires_bathy,fwyf1.lon,fwyf1.lat,3);
end;
if ( ~isfield(rsm,'slope') )
  [rsm.slope,rsm.slope_orientation,rsm.isobath_orientation,lkwf1.ngdc_hires_bathy] = ...
      find_ngdc_slope(lkwf1.ngdc_hires_bathy,rsm.lon,rsm.lat,3);
end;
if ( ~isfield(pvgf1,'slope') )
  [pvgf1.slope,pvgf1.slope_orientation,pvgf1.isobath_orientation,lkwf1.ngdc_hires_bathy] = ...
      find_ngdc_slope(lkwf1.ngdc_hires_bathy,pvgf1.lon,pvgf1.lat,3);
end;
if ( ~isfield(ftpf1,'slope') )
  % [ftpf1.slope,ftpf1.slope_orientation,ftpf1.isobath_orientation,lkwf1.ngdc_hires_bathy] = ...
  %     find_ngdc_slope(lkwf1.ngdc_hires_bathy,ftpf1.lon,ftpf1.lat,3);
  % Once we fixed a bug in READ_HIRES_BATHYMETRY, suddenly station FTPF1 was
  % no longer within "120e3" meters of LKWF1!
  [ftpf1.slope,ftpf1.slope_orientation,ftpf1.isobath_orientation,ftpf1.ngdc_hires_bathy] = ...
      find_ngdc_slope(ftpf1.ngdc_hires_bathy,ftpf1.lon,ftpf1.lat,3);
end;
if ( ~isfield(fdepk,'slope') )
  [fdepk.slope,fdepk.slope_orientation,fdepk.isobath_orientation,lkwf1.ngdc_hires_bathy] = ...
      find_ngdc_slope(lkwf1.ngdc_hires_bathy,fdepk.lon,fdepk.lat,3);
end;
if ( ~isfield(canf1,'slope') )
  [canf1.slope,canf1.slope_orientation,canf1.isobath_orientation,canf1.ngdc_hires_bathy] = ...
      find_ngdc_slope(canf1.ngdc_hires_bathy,canf1.lon,canf1.lat,3);
end;
if ( ~isfield(cnnf1,'slope') )
  [cnnf1.slope,cnnf1.slope_orientation,cnnf1.isobath_orientation,canf1.ngdc_hires_bathy] = ...
      find_ngdc_slope(canf1.ngdc_hires_bathy,cnnf1.lon,cnnf1.lat,3);
end;

if ( exist('face','var') && ~isfield(face.shallow,'slope') )
  [face.shallow.slope,face.shallow.slope_orientation,face.shallow.isobath_orientation,lkwf1.ngdc_hires_bathy] = ...
      find_ngdc_slope(lkwf1.ngdc_hires_bathy,face.shallow.lon,face.shallow.lat,3);
  [face.tcm1.slope,face.tcm1.slope_orientation,face.tcm1.isobath_orientation,lkwf1.ngdc_hires_bathy] = ...
      find_ngdc_slope(lkwf1.ngdc_hires_bathy,face.tcm1.lon,face.tcm1.lat,3);
  [face.tcm2.slope,face.tcm2.slope_orientation,face.tcm2.isobath_orientation,lkwf1.ngdc_hires_bathy] = ...
      find_ngdc_slope(lkwf1.ngdc_hires_bathy,face.tcm2.lon,face.tcm2.lat,3);
  [face.deep.slope,face.deep.slope_orientation,face.deep.isobath_orientation,lkwf1.ngdc_hires_bathy] = ...
      find_ngdc_slope(lkwf1.ngdc_hires_bathy,face.deep.lon,face.deep.lat,3);

  if ( isfield(face,'HanS_HW') )
    [face.HanS_HW.slope,face.HanS_HW.slope_orientation,face.HanS_HW.isobath_orientation,lkwf1.ngdc_hires_bathy] = ...
        find_ngdc_slope(lkwf1.ngdc_hires_bathy,face.HanS_HW.lon,face.HanS_HW.lat,3);
    [face.HanS_BR.slope,face.HanS_BR.slope_orientation,face.HanS_BR.isobath_orientation,lkwf1.ngdc_hires_bathy] = ...
        find_ngdc_slope(lkwf1.ngdc_hires_bathy,face.HanS_BR.lon,face.HanS_BR.lat,3);
    [face.AOML_HW.slope,face.AOML_HW.slope_orientation,face.AOML_HW.isobath_orientation,lkwf1.ngdc_hires_bathy] = ...
        find_ngdc_slope(lkwf1.ngdc_hires_bathy,face.AOML_HW.lon,face.AOML_HW.lat,3);
    [face.AOML_BR.slope,face.AOML_BR.slope_orientation,face.AOML_BR.isobath_orientation,lkwf1.ngdc_hires_bathy] = ...
        find_ngdc_slope(lkwf1.ngdc_hires_bathy,face.AOML_BR.lon,face.AOML_BR.lat,3);
  end;

  if ( isfield(face,'brwd20') )
    [face.brwd20.slope,face.brwd20.slope_orientation,face.brwd20.isobath_orientation,lkwf1.ngdc_hires_bathy] = ...
        find_ngdc_slope(lkwf1.ngdc_hires_bathy,face.brwd20.lon,face.brwd20.lat,3);
    [face.brwd40.slope,face.brwd40.slope_orientation,face.brwd40.isobath_orientation,lkwf1.ngdc_hires_bathy] = ...
        find_ngdc_slope(lkwf1.ngdc_hires_bathy,face.brwd40.lon,face.brwd40.lat,3);
    [face.brwd100.slope,face.brwd100.slope_orientation,face.brwd100.isobath_orientation,lkwf1.ngdc_hires_bathy] = ...
        find_ngdc_slope(lkwf1.ngdc_hires_bathy,face.brwd100.lon,face.brwd100.lat,3);
    [face.boca20.slope,face.boca20.slope_orientation,face.boca20.isobath_orientation,lkwf1.ngdc_hires_bathy] = ...
        find_ngdc_slope(lkwf1.ngdc_hires_bathy,face.boca20.lon,face.boca20.lat,3);
    [face.boca40.slope,face.boca40.slope_orientation,face.boca40.isobath_orientation,lkwf1.ngdc_hires_bathy] = ...
        find_ngdc_slope(lkwf1.ngdc_hires_bathy,face.boca40.lon,face.boca40.lat,3);
  end;
end;

if ( exist('sefcri','var') && ~isfield(sefcri.updb,'slope') )
  [sefcri.dc1.slope,sefcri.dc1.slope_orientation,sefcri.dc1.isobath_orientation,lkwf1.ngdc_hires_bathy] = ...
      find_ngdc_slope(lkwf1.ngdc_hires_bathy,sefcri.dc1.lon,sefcri.dc1.lat,3);
  [sefcri.dc3.slope,sefcri.dc3.slope_orientation,sefcri.dc3.isobath_orientation,lkwf1.ngdc_hires_bathy] = ...
      find_ngdc_slope(lkwf1.ngdc_hires_bathy,sefcri.dc3.lon,sefcri.dc3.lat,3);
  [sefcri.bc1.slope,sefcri.bc1.slope_orientation,sefcri.bc1.isobath_orientation,lkwf1.ngdc_hires_bathy] = ...
      find_ngdc_slope(lkwf1.ngdc_hires_bathy,sefcri.bc1.lon,sefcri.bc1.lat,3);
  [sefcri.bc3.slope,sefcri.bc3.slope_orientation,sefcri.bc3.isobath_orientation,lkwf1.ngdc_hires_bathy] = ...
      find_ngdc_slope(lkwf1.ngdc_hires_bathy,sefcri.bc3.lon,sefcri.bc3.lat,3);
  [sefcri.pb1.slope,sefcri.pb1.slope_orientation,sefcri.pb1.isobath_orientation,lkwf1.ngdc_hires_bathy] = ...
      find_ngdc_slope(lkwf1.ngdc_hires_bathy,sefcri.pb1.lon,sefcri.pb1.lat,3);
  [sefcri.pb2.slope,sefcri.pb2.slope_orientation,sefcri.pb2.isobath_orientation,lkwf1.ngdc_hires_bathy] = ...
      find_ngdc_slope(lkwf1.ngdc_hires_bathy,sefcri.pb2.lon,sefcri.pb2.lat,3);
  [sefcri.mc2.slope,sefcri.mc2.slope_orientation,sefcri.mc2.isobath_orientation,lkwf1.ngdc_hires_bathy] = ...
      find_ngdc_slope(lkwf1.ngdc_hires_bathy,sefcri.mc2.lon,sefcri.mc2.lat,3);
  [sefcri.pela.slope,sefcri.pela.slope_orientation,sefcri.pela.isobath_orientation,lkwf1.ngdc_hires_bathy] = ...
      find_ngdc_slope(lkwf1.ngdc_hires_bathy,sefcri.pela.lon,sefcri.pela.lat,3);
  [sefcri.evca.slope,sefcri.evca.slope_orientation,sefcri.evca.isobath_orientation,lkwf1.ngdc_hires_bathy] = ...
      find_ngdc_slope(lkwf1.ngdc_hires_bathy,sefcri.evca.lon,sefcri.evca.lat,3);
  [sefcri.slib.slope,sefcri.slib.slope_orientation,sefcri.slib.isobath_orientation,lkwf1.ngdc_hires_bathy] = ...
      find_ngdc_slope(lkwf1.ngdc_hires_bathy,sefcri.slib.lon,sefcri.slib.lat,3);
  [sefcri.updb.slope,sefcri.updb.slope_orientation,sefcri.updb.isobath_orientation,lkwf1.ngdc_hires_bathy] = ...
      find_ngdc_slope(lkwf1.ngdc_hires_bathy,sefcri.updb.lon,sefcri.updb.lat,3);
end;

if ( exist('sfomc','var') && ~isfield(sfomc.nw_w_btm,'slope') )
  [sfomc.nw_w_btm.slope,sfomc.nw_w_btm.slope_orientation,sfomc.nw_w_btm.isobath_orientation,lkwf1.ngdc_hires_bathy] = ...
      find_ngdc_slope(lkwf1.ngdc_hires_bathy,sfomc.nw_w_btm.lon,sfomc.nw_w_btm.lat,3);
  [sfomc.ne_buoy.slope,sfomc.ne_buoy.slope_orientation,sfomc.ne_buoy.isobath_orientation,lkwf1.ngdc_hires_bathy] = ...
      find_ngdc_slope(lkwf1.ngdc_hires_bathy,sfomc.ne_buoy.lon,sfomc.ne_buoy.lat,3);
  [sfomc.e_buoy.slope,sfomc.e_buoy.slope_orientation,sfomc.e_buoy.isobath_orientation,lkwf1.ngdc_hires_bathy] = ...
      find_ngdc_slope(lkwf1.ngdc_hires_bathy,sfomc.e_buoy.lon,sfomc.e_buoy.lat,3);
  [sfomc.c_buoy.slope,sfomc.c_buoy.slope_orientation,sfomc.c_buoy.isobath_orientation,lkwf1.ngdc_hires_bathy] = ...
      find_ngdc_slope(lkwf1.ngdc_hires_bathy,sfomc.c_buoy.lon,sfomc.c_buoy.lat,3);
  [sfomc.sw_buoy.slope,sfomc.sw_buoy.slope_orientation,sfomc.sw_buoy.isobath_orientation,lkwf1.ngdc_hires_bathy] = ...
      find_ngdc_slope(lkwf1.ngdc_hires_bathy,sfomc.sw_buoy.lon,sfomc.sw_buoy.lat,3);
  [sfomc.se_buoy.slope,sfomc.se_buoy.slope_orientation,sfomc.se_buoy.isobath_orientation,lkwf1.ngdc_hires_bathy] = ...
      find_ngdc_slope(lkwf1.ngdc_hires_bathy,sfomc.se_buoy.lon,sfomc.se_buoy.lat,3);
  [sfomc.pier_cc.slope,sfomc.pier_cc.slope_orientation,sfomc.pier_cc.isobath_orientation,lkwf1.ngdc_hires_bathy] = ...
      find_ngdc_slope(lkwf1.ngdc_hires_bathy,sfomc.pier_cc.lon,sfomc.pier_cc.lat,3);
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

if ( exist('face','var') && ~isfield(face.shallow,'adcp_x') )
  face.shallow = station_reorient_vectors(face.shallow,'isobath_orientation','adcp_u','adcp_v','adcp_x','adcp_l');
  face.tcm1 = station_reorient_vectors(face.tcm1,'isobath_orientation','u','v','x','l');
  face.tcm2 = station_reorient_vectors(face.tcm2,'isobath_orientation','u','v','x','l');
  face.deep = station_reorient_vectors(face.deep,'isobath_orientation','adcp_u','adcp_v','adcp_x','adcp_l');

  %%%% ??? Could FACE TCMs have had MIS-ORIENTED components??
  face.tcm1 = station_reorient_vectors(face.tcm1,-90,'x','l','alt_x','alt_l');
  face.tcm1.alt_x.data = -face.tcm1.alt_x.data;
  face.tcm2 = station_reorient_vectors(face.tcm2,-90,'x','l','alt_x','alt_l');
  face.tcm2.alt_x.data = -face.tcm2.alt_x.data;

  % % From READ_FACE_JACK_CURRENTS.m:
  %adcp_u_sfc  
  %adcp_v_sfc  
  % shallow_top_good_bin = 10;
  % shallow_sfc_bins = shallow_top_good_bin-2:shallow_top_good_bin;
  % shallow_mid_bins = 5:7;
  % shallow_btm_bins = 1:4;
  % deep_top_good_bin = 18;
  % deep_sfc_bins = deep_top_good_bin-5:deep_top_good_bin;
  % deep_mid_bins = 7:12;
  % deep_btm_bins = 1:6;
  face.shallow = station_reorient_vectors(face.shallow,'isobath_orientation',...
                                          'adcp_u_sfc','adcp_v_sfc','adcp_sfc_x','adcp_sfc_l');
  face.shallow = station_reorient_vectors(face.shallow,'isobath_orientation',...
                                          'adcp_u_mid','adcp_v_mid','adcp_mid_x','adcp_mid_l');
  face.shallow = station_reorient_vectors(face.shallow,'isobath_orientation',...
                                          'adcp_u_btm','adcp_v_btm','adcp_btm_x','adcp_btm_l');
  face.deep = station_reorient_vectors(face.deep,'isobath_orientation',...
                                       'adcp_u_sfc','adcp_v_sfc','adcp_sfc_x','adcp_sfc_l');
  face.deep = station_reorient_vectors(face.deep,'isobath_orientation',...
                                       'adcp_u_mid','adcp_v_mid','adcp_mid_x','adcp_mid_l');
  face.deep = station_reorient_vectors(face.deep,'isobath_orientation',...
                                       'adcp_u_btm','adcp_v_btm','adcp_btm_x','adcp_btm_l');

  if ( isfield(face,'HanS_HW') )
    face.HanS_HW = station_reorient_vectors(face.HanS_HW,'isobath_orientation','adcp_u','adcp_v','adcp_x','adcp_l');
    face.HanS_BR = station_reorient_vectors(face.HanS_BR,'isobath_orientation','adcp_u','adcp_v','adcp_x','adcp_l');
    face.AOML_HW = station_reorient_vectors(face.AOML_HW,'isobath_orientation','adcp_u','adcp_v','adcp_x','adcp_l');
    face.AOML_BR = station_reorient_vectors(face.AOML_BR,'isobath_orientation','adcp_u','adcp_v','adcp_x','adcp_l');
    
    if ( isfield(face.HanS_HW,'adcp_u_sfc') )
      face.HanS_HW = station_reorient_vectors(face.HanS_HW,'isobath_orientation',...
                                              'adcp_u_sfc','adcp_v_sfc','adcp_x_sfc','adcp_l_sfc');
      face.HanS_HW = station_reorient_vectors(face.HanS_HW,'isobath_orientation',...
                                              'adcp_u_mid','adcp_v_mid','adcp_x_mid','adcp_l_mid');
      face.HanS_HW = station_reorient_vectors(face.HanS_HW,'isobath_orientation',...
                                              'adcp_u_btm','adcp_v_btm','adcp_x_btm','adcp_l_btm');
      face.HanS_BR = station_reorient_vectors(face.HanS_BR,'isobath_orientation',...
                                              'adcp_u_sfc','adcp_v_sfc','adcp_x_sfc','adcp_l_sfc');
      face.HanS_BR = station_reorient_vectors(face.HanS_BR,'isobath_orientation',...
                                              'adcp_u_mid','adcp_v_mid','adcp_x_mid','adcp_l_mid');
      face.HanS_BR = station_reorient_vectors(face.HanS_BR,'isobath_orientation',...
                                              'adcp_u_btm','adcp_v_btm','adcp_x_btm','adcp_l_btm');

      face.AOML_HW = station_reorient_vectors(face.AOML_HW,'isobath_orientation',...
                                              'adcp_u_sfc','adcp_v_sfc','adcp_x_sfc','adcp_l_sfc');
      face.AOML_HW = station_reorient_vectors(face.AOML_HW,'isobath_orientation',...
                                              'adcp_u_mid','adcp_v_mid','adcp_x_mid','adcp_l_mid');
      face.AOML_HW = station_reorient_vectors(face.AOML_HW,'isobath_orientation',...
                                              'adcp_u_btm','adcp_v_btm','adcp_x_btm','adcp_l_btm');
      face.AOML_BR = station_reorient_vectors(face.AOML_BR,'isobath_orientation',...
                                              'adcp_u_sfc','adcp_v_sfc','adcp_x_sfc','adcp_l_sfc');
      face.AOML_BR = station_reorient_vectors(face.AOML_BR,'isobath_orientation',...
                                              'adcp_u_mid','adcp_v_mid','adcp_x_mid','adcp_l_mid');
      face.AOML_BR = station_reorient_vectors(face.AOML_BR,'isobath_orientation',...
                                              'adcp_u_btm','adcp_v_btm','adcp_x_btm','adcp_l_btm');
    end;
  end;

end;


if ( ~isfield(fwyf1,'ndbc_ekman_flux') )
  disp('Calculating Ekman fluxes...');
  fwyf1 = station_ekman_flux(fwyf1,[],'ndbc','erai');
  pvgf1 = station_ekman_flux(pvgf1,[],'ndbc');
  lkwf1 = station_ekman_flux(lkwf1,[],'ndbc');
  fdepk = station_ekman_flux(fdepk,[],'fdep');
  canf1 = station_ekman_flux(canf1,[],'ndbc','erai');
  rsm = station_ekman_flux(rsm,[],'rsmas');
end;

for hrs = 12:12:36
  for cderf={'lp','sum','avg'};
    fldstub = sprintf('_%d_h_%s',hrs,cderf{:});

    fld = ['ndbc_ekman_flux_x_volume',fldstub];
    disp(fld);
    fwyf1 = verify_variable(fwyf1,fld);
    fwyf1 = filter_gaps(fwyf1,'ndbc_ekman_flux_x_volume',fld);

    rsm = verify_variable(rsm,['rsmas_ekman_flux_x_volume',fldstub]);
    rsm = filter_gaps(rsm,'rsmas_ekman_flux_x_volume',['rsmas_ekman_flux_x_volume',fldstub]);

    pvgf1 = verify_variable(pvgf1,fld);
    pvgf1 = filter_gaps(pvgf1,'ndbc_ekman_flux_x_volume',fld);

    lkwf1 = verify_variable(lkwf1,fld);
    lkwf1 = filter_gaps(lkwf1,'ndbc_ekman_flux_x_volume',fld);

    canf1 = verify_variable(canf1,fld);
    canf1 = filter_gaps(canf1,'ndbc_ekman_flux_x_volume',fld);

    if ( exist('face','var') )
      face.deep = verify_variable(face.deep,['adcp_x',fldstub]);
      face.deep = verify_variable(face.deep,['adcp_l',fldstub]);
      face.deep = filter_gaps(face.deep,'adcp_x',['adcp_x',fldstub]);
      face.deep = filter_gaps(face.deep,'adcp_l',['adcp_l',fldstub]);

      face.tcm2 = verify_variable(face.tcm2,['x',fldstub]);
      face.tcm2 = verify_variable(face.tcm2,['l',fldstub]);
      face.tcm2 = filter_gaps(face.tcm2,'x',['x',fldstub]);
      face.tcm2 = filter_gaps(face.tcm2,'l',['l',fldstub]);

      face.tcm1 = verify_variable(face.tcm1,['x',fldstub]);
      face.tcm1 = verify_variable(face.tcm1,['l',fldstub]);
      face.tcm1 = filter_gaps(face.tcm1,'x',['x',fldstub]);
      face.tcm1 = filter_gaps(face.tcm1,'l',['l',fldstub]);

      face.shallow = verify_variable(face.shallow,['adcp_x',fldstub]);
      face.shallow = verify_variable(face.shallow,['adcp_l',fldstub]);
      face.shallow = filter_gaps(face.shallow,'adcp_x',['adcp_x',fldstub]);
      face.shallow = filter_gaps(face.shallow,'adcp_l',['adcp_l',fldstub]);
    end; %if ( exist('face','var') )
  end;
end;


%% Variability frequency partitions
%pers = [ 3,11 ; 11,16 ; 17,28 ; 2*24,10*24 ];
pers = [ 3,9 ; 11,14 ; 17,28 ; 2*24,10*24 ];

if ( exist('canf1','var') && ~isfield(canf1,'ndbc_ekman_flux_x_volume_3_h_hp') )
  canf1 = station_partition_periods(canf1,'ndbc_ekman_flux_x_volume',pers);
  canf1 = station_partition_periods(canf1,'ndbc_ekman_flux_l_volume',pers);
end;
if ( exist('lkwf1','var') && ~isfield(lkwf1,'ndbc_ekman_flux_x_volume_3_h_hp') )
  lkwf1 = station_partition_periods(lkwf1,'ndbc_ekman_flux_x_volume',pers);
  lkwf1 = station_partition_periods(lkwf1,'ndbc_ekman_flux_l_volume',pers);
end;
if ( exist('pvgf1','var') && ~isfield(pvgf1,'ndbc_ekman_flux_x_volume_3_h_hp') )
  pvgf1 = station_partition_periods(pvgf1,'ndbc_ekman_flux_x_volume',pers);
  pvgf1 = station_partition_periods(pvgf1,'ndbc_ekman_flux_l_volume',pers);
end;
if ( exist('rsm','var') && ~isfield(rsm,'rsmas_ekman_flux_x_volume_3_h_hp') )
  rsm = station_partition_periods(rsm,'rsmas_ekman_flux_x_volume',pers);
  rsm = station_partition_periods(rsm,'rsmas_ekman_flux_l_volume',pers);
end;
if ( exist('fwyf1','var') && ~isfield(fwyf1,'ndbc_ekman_flux_x_volume_3_h_hp') )
  fwyf1 = station_partition_periods(fwyf1,'ndbc_ekman_flux_x_volume',pers);
  fwyf1 = station_partition_periods(fwyf1,'ndbc_ekman_flux_l_volume',pers);
end;


if ( exist('sefcri','var') && ~isfield(sefcri.updb,'hourly_t_3_h_hp') )
  sefcri.updb = station_partition_periods(sefcri.updb,'hourly_t',pers);
  sefcri.pb2 = station_partition_periods(sefcri.pb2,'hourly_t',pers);
  sefcri.bc3 = station_partition_periods(sefcri.bc3,'hourly_t',pers);
  sefcri.dc3 = station_partition_periods(sefcri.dc3,'hourly_t',pers);
end;

if ( exist('face','var') && ~isfield(face.deep,'seatemp_3_h_hp') )
  face.deep = station_partition_periods(face.deep,'seatemp',pers);
  face.deep = station_partition_periods(face.deep,'adcp_btm_x',pers);
  face.deep = station_partition_periods(face.deep,'adcp_btm_l',pers);
  face.deep = station_partition_periods(face.deep,'adcp_sfc_x',pers);
  face.deep = station_partition_periods(face.deep,'adcp_sfc_l',pers);
end;

if ( exist('sfomc','var') && ~isfield(sfomc.ne_buoy,'seatemp_3_h_hp') )
  sfomc.ne_buoy = station_partition_periods(sfomc.ne_buoy,'seatemp',pers);
  sfomc.ne_buoy = station_partition_periods(sfomc.ne_buoy,'adcp_x_btm',pers);
  sfomc.ne_buoy = station_partition_periods(sfomc.ne_buoy,'adcp_l_btm',pers);
  sfomc.ne_buoy = station_partition_periods(sfomc.ne_buoy,'adcp_x_uwc',pers);
  sfomc.ne_buoy = station_partition_periods(sfomc.ne_buoy,'adcp_l_uwc',pers);
  sfomc.e_buoy = station_partition_periods(sfomc.e_buoy,'seatemp',pers);
  sfomc.e_buoy = station_partition_periods(sfomc.e_buoy,'adcp_x_btm',pers);
  sfomc.e_buoy = station_partition_periods(sfomc.e_buoy,'adcp_l_btm',pers);
  sfomc.e_buoy = station_partition_periods(sfomc.e_buoy,'adcp_x_uwc',pers);
  sfomc.e_buoy = station_partition_periods(sfomc.e_buoy,'adcp_l_uwc',pers);
  sfomc.sw_buoy = station_partition_periods(sfomc.sw_buoy,'seatemp',pers);
  sfomc.sw_buoy = station_partition_periods(sfomc.sw_buoy,'adcp_x_btm',pers);
  sfomc.sw_buoy = station_partition_periods(sfomc.sw_buoy,'adcp_l_btm',pers);
  sfomc.sw_buoy = station_partition_periods(sfomc.sw_buoy,'adcp_x_uwc',pers);
  sfomc.sw_buoy = station_partition_periods(sfomc.sw_buoy,'adcp_l_uwc',pers);
  sfomc.c_buoy = station_partition_periods(sfomc.c_buoy,'seatemp',pers);
  sfomc.c_buoy = station_partition_periods(sfomc.c_buoy,'adcp_x_btm',pers);
  sfomc.c_buoy = station_partition_periods(sfomc.c_buoy,'adcp_l_btm',pers);
  sfomc.c_buoy = station_partition_periods(sfomc.c_buoy,'adcp_x_uwc',pers);
  sfomc.c_buoy = station_partition_periods(sfomc.c_buoy,'adcp_l_uwc',pers);
  % % What was our fundamental sampling period?
  % fmg; plot(sfomc.nw_w_btm.adcp_x_btm.date(1:end-1),24*diff(sfomc.nw_w_btm.adcp_x_btm.date),'.'); datetick3; ylim([0,8]);
  sfomc.nw_w_btm = station_partition_periods(sfomc.nw_w_btm,'adcp_seatemp',pers);
  sfomc.nw_w_btm = station_partition_periods(sfomc.nw_w_btm,'adcp_x_btm',pers);
  sfomc.nw_w_btm = station_partition_periods(sfomc.nw_w_btm,'adcp_l_btm',pers);
  sfomc.nw_w_btm = station_partition_periods(sfomc.nw_w_btm,'adcp_x_uwc',pers);
  sfomc.nw_w_btm = station_partition_periods(sfomc.nw_w_btm,'adcp_l_uwc',pers);
end;


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% NUTRIENTS and NUTRIENT FLUXES

% Adapated from plot_upwelling_hovmoellers.m, 2018 Mar 06

if ( doNuts )
  
  sfomc.ne_buoy.sfctemp = sfomc.ne_buoy.at_0_6m.sbe_seatemp;
  sfomc.sw_buoy.sfctemp = sfomc.sw_buoy.at_0_6m.sbe_seatemp;
  [sfomc.nw_w_btm.seatemp,s1] = intersect_tses(sfomc.nw_w_btm.adcp_seatemp,sfomc.nw_w_btm.sbe_seatemp);
  sfomc.nw_w_btm.seatemp.depths = [-5,-11];
  sfomc.nw_w_btm.seatemp.prof(:,2) = sfomc.nw_w_btm.seatemp.data;
  sfomc.nw_w_btm.seatemp.prof(:,1) = s1.data;



  % Nitrate + Nitrite
  sfomc.ne_buoy.N = sfomc.ne_buoy.seatemp;
  sfomc.ne_buoy.N.prof = (23-sfomc.ne_buoy.seatemp.prof).*1.62;
  sfomc.ne_buoy.N.prof(sfomc.ne_buoy.N.prof<0) = 0;
  sfomc.ne_buoy.N.data = nansum(sfomc.ne_buoy.N.prof,2);
  
  % Phosphates
  sfomc.ne_buoy.P = sfomc.ne_buoy.seatemp;
  sfomc.ne_buoy.P.prof = (23-sfomc.ne_buoy.seatemp.prof).*0.073;
  sfomc.ne_buoy.P.prof(sfomc.ne_buoy.P.prof<0) = 0;
  sfomc.ne_buoy.P.data = nansum(sfomc.ne_buoy.P.prof,2);

  % Silicates
  sfomc.ne_buoy.Si = sfomc.ne_buoy.seatemp;
  sfomc.ne_buoy.Si.prof = (23-sfomc.ne_buoy.seatemp.prof).*0.37;
  sfomc.ne_buoy.Si.prof(sfomc.ne_buoy.Si.prof<0) = 0;
  sfomc.ne_buoy.Si.data = nansum(sfomc.ne_buoy.Si.prof,2);


  % Concentrations at other sites are approximated from assumed heat mixing rates
  
  % Which is the right temperature end-member?
  % >> fmg; plot_ts(sfomc.ne_buoy.at_0_6m.sbe_seatemp,sfomc.sw_buoy.at_0_6m.sbe_seatemp,sfomc.nw_w_btm.sbe_seatemp,fwyf1.ndbc_air_t,fwyf1.ndbc_sea_t); xlim(sfomc.ne_buoy.at_0_6m.sbe_seatemp.date([1,end])); datetick3; legend('NE','SW','NW','FWY','T_s');
  
  % % For estimating onshore nutrient fluxes from *2000* data
  % [nenitr,nephos,nesili, netemp,swtemp,nwtemp, sfctemp, airtemp] = ...
  %     intersect_tses(sfomc.ne_buoy.N,sfomc.ne_buoy.P,sfomc.ne_buoy.Si, sfomc.ne_buoy.seatemp,sfomc.sw_buoy.seatemp,sfomc.nw_w_btm.seatemp, sfomc.sw_buoy.sfctemp, fwyf1.ndbc_air_t);

  % For estimating onshore nutrient fluxes from *1999* data
  [nenitr,nephos,nesili,netemp,swtemp,nwtemp,sfctemp, airtemp] = ...
      intersect_tses(sfomc.ne_buoy.N,sfomc.ne_buoy.P,sfomc.ne_buoy.Si, sfomc.ne_buoy.seatemp,sfomc.sw_buoy.seatemp,sfomc.nw_w_btm.seatemp, sfomc.nw_w_btm.sbe_seatemp, fwyf1.ndbc_air_t);

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
  
  sfomc.sw_buoy.mixprof = (1-(swtemp.prof-repmat(Tendpt,[1,size(swtemp.prof,2)])) ...
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
  
  % sfomc.sw_buoy.mixprof = repmat(0.28,size(swtemp.prof));
  
  
  % Ignore non-physical mixing rate profiles
  sfomc.sw_buoy.mixprof(0>sfomc.sw_buoy.mixprof | sfomc.sw_buoy.mixprof>1) = 0;
  
  swNendpt = repmat(Nendpt,[1,size(swtemp.prof,2)]);
  sfomc.sw_buoy.N = swtemp;
  sfomc.sw_buoy.N.prof = swNendpt .* sfomc.sw_buoy.mixprof;
  sfomc.sw_buoy.N.data = nansum(sfomc.sw_buoy.N.prof,2);
  
  swPendpt = repmat(Pendpt,[1,size(swtemp.prof,2)]);
  sfomc.sw_buoy.P = swtemp;
  sfomc.sw_buoy.P.prof = swPendpt .* sfomc.sw_buoy.mixprof;
  sfomc.sw_buoy.P.data = nansum(sfomc.sw_buoy.P.prof,2);
  
  swSendpt = repmat(Sendpt,[1,size(swtemp.prof,2)]);
  sfomc.sw_buoy.Si = swtemp;
  sfomc.sw_buoy.Si.prof = swSendpt .* sfomc.sw_buoy.mixprof;
  sfomc.sw_buoy.Si.data = nansum(sfomc.sw_buoy.Si.prof,2);
  
  
  %%%
  % NW/W mooring (11 m isobath) nutrient concentrations
  
  sfomc.nw_w_btm.mixprof = (1-(nwtemp.prof-repmat(Tendpt,[1,size(nwtemp.prof,2)])) ...
                            ./ repmat((Tsfcpt-Tendpt),[1,size(nwtemp.prof,2)]));
  sfomc.nw_w_btm.mixprof(0>sfomc.nw_w_btm.mixprof | sfomc.nw_w_btm.mixprof>1) = 0;
  
  nwNendpt = repmat(Nendpt,[1,size(nwtemp.prof,2)]);
  sfomc.nw_w_btm.N = nwtemp;
  sfomc.nw_w_btm.N.prof = nwNendpt .* sfomc.nw_w_btm.mixprof;
  sfomc.nw_w_btm.N.data = nansum(sfomc.nw_w_btm.N.prof,2);
  
  nwPendpt = repmat(Pendpt,[1,size(nwtemp.prof,2)]);
  sfomc.nw_w_btm.P = nwtemp;
  sfomc.nw_w_btm.P.prof = nwPendpt .* sfomc.nw_w_btm.mixprof;
  sfomc.nw_w_btm.P.data = nansum(sfomc.nw_w_btm.P.prof,2);
  
  nwSendpt = repmat(Sendpt,[1,size(nwtemp.prof,2)]);
  sfomc.nw_w_btm.Si = nwtemp;
  sfomc.nw_w_btm.Si.prof = nwSendpt .* sfomc.nw_w_btm.mixprof;
  sfomc.nw_w_btm.Si.data = nansum(sfomc.nw_w_btm.Si.prof,2);
  
  
  
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
  
  
  % % OOPS! Much of the time, upwelled water seems to be moving OFFshore?!
  % fmg; boxplot(sfomc.sw_buoy.adcp_x.prof(sfomc.sw_buoy.Nx.xprof(:,end)>0,:),sfomc.sw_buoy.adcp_depths,'notch','on');
  %  titlename('Mid-Station upwelling u\bullet\nablah profile');
  % fmg; boxplot(sfomc.ne_buoy.adcp_x.prof(sfomc.ne_buoy.Nx.xprof(:,end)>0,:),sfomc.ne_buoy.adcp_depths,'notch','on');
  %  titlename('East-Station upwelling u\bullet\nablah profile');
  
  % Track only movement of water "onshore" (isobath-oriented 'u' < 0)
  sfomc.nw_w_btm.adcp_onshore = sfomc.nw_w_btm.adcp_x;
  sfomc.nw_w_btm.adcp_onshore.prof(sfomc.nw_w_btm.adcp_onshore.prof>0)=0;
  sfomc.nw_w_btm.adcp_onshore.prof(sfomc.nw_w_btm.adcp_onshore.prof<0)=-sfomc.nw_w_btm.adcp_onshore.prof(sfomc.nw_w_btm.adcp_onshore.prof<0);
  
  sfomc.sw_buoy.adcp_onshore = sfomc.sw_buoy.adcp_x;
  sfomc.sw_buoy.adcp_onshore.prof(sfomc.sw_buoy.adcp_onshore.prof>0)=0;
  sfomc.sw_buoy.adcp_onshore.prof(sfomc.sw_buoy.adcp_onshore.prof<0)=-sfomc.sw_buoy.adcp_onshore.prof(sfomc.sw_buoy.adcp_onshore.prof<0);
  
  sfomc.ne_buoy.adcp_onshore = sfomc.ne_buoy.adcp_x;
  sfomc.ne_buoy.adcp_onshore.prof(sfomc.ne_buoy.adcp_onshore.prof>0)=0;
  sfomc.ne_buoy.adcp_onshore.prof(sfomc.ne_buoy.adcp_onshore.prof<0)=-sfomc.ne_buoy.adcp_onshore.prof(sfomc.ne_buoy.adcp_onshore.prof<0);
  
  
  %%%
  % NW/W buoy (11 m isobath) nutrient fluxes
  
  [nwix,neix] = intersect_dates(sfomc.nw_w_btm.adcp_onshore.date,sfomc.ne_buoy.adcp_onshore.date);

  % [NO3+NO2]
  sfomc.nw_w_btm.Nx.date = sfomc.nw_w_btm.adcp_onshore.date(nwix);
  sfomc.nw_w_btm.Nx.depths = sfomc.nw_w_btm.adcp_depths(end:-1:1);
  sfomc.nw_w_btm.Nx.prof = interp2(sfomc.nw_w_btm.N.date,sfomc.nw_w_btm.N.depths,sfomc.nw_w_btm.N.prof',...
                                   sfomc.nw_w_btm.Nx.date,sfomc.nw_w_btm.Nx.depths,'linear',NaN)';
  sfomc.nw_w_btm.Nx.xprof = sfomc.nw_w_btm.adcp_onshore.prof(nwix,:) .* umol_to_kgm3(sfomc.nw_w_btm.Nx.prof,'N') ...
      .* abs(median(diff(sfomc.nw_w_btm.Nx.depths))) .* median(diff(sfomc.nw_w_btm.Nx.date)*60*60*24);
  sfomc.nw_w_btm.Nx.cprof = cumsum(sfomc.nw_w_btm.Nx.xprof,'omitnan');
  
  sfomc.nw_w_btm.intNx.date = sfomc.nw_w_btm.Nx.date;
  sfomc.nw_w_btm.intNx.data = sum(sfomc.nw_w_btm.Nx.cprof,2,'omitnan');
  
  % [P]
  sfomc.nw_w_btm.Px.date = sfomc.nw_w_btm.adcp_onshore.date(nwix);
  sfomc.nw_w_btm.Px.depths = sfomc.nw_w_btm.adcp_depths(end:-1:1);
  sfomc.nw_w_btm.Px.prof = interp2(sfomc.nw_w_btm.P.date,sfomc.nw_w_btm.P.depths,sfomc.nw_w_btm.P.prof',...
                                   sfomc.nw_w_btm.Px.date,sfomc.nw_w_btm.Px.depths,'linear',NaN)';
  sfomc.nw_w_btm.Px.xprof = sfomc.nw_w_btm.adcp_onshore.prof(nwix,:) .* umol_to_kgm3(sfomc.nw_w_btm.Px.prof,'P') ...
      .* abs(median(diff(sfomc.nw_w_btm.Px.depths))) .* median(diff(sfomc.nw_w_btm.Px.date)*60*60*24);
  sfomc.nw_w_btm.Px.cprof = cumsum(sfomc.nw_w_btm.Px.xprof,'omitnan');
  
  sfomc.nw_w_btm.intPx.date = sfomc.nw_w_btm.Px.date;
  sfomc.nw_w_btm.intPx.data = sum(sfomc.nw_w_btm.Px.cprof,2,'omitnan');
  
  % [Si]
  sfomc.nw_w_btm.Six.date = sfomc.nw_w_btm.adcp_onshore.date(nwix);
  sfomc.nw_w_btm.Six.depths = sfomc.nw_w_btm.adcp_depths(end:-1:1);
  sfomc.nw_w_btm.Six.prof = interp2(sfomc.nw_w_btm.Si.date,sfomc.nw_w_btm.Si.depths,sfomc.nw_w_btm.Si.prof',...
                                   sfomc.nw_w_btm.Six.date,sfomc.nw_w_btm.Six.depths,'linear',NaN)';
  sfomc.nw_w_btm.Six.xprof = sfomc.nw_w_btm.adcp_onshore.prof(nwix,:) .* umol_to_kgm3(sfomc.nw_w_btm.Six.prof,'Si') ...
      .* abs(median(diff(sfomc.nw_w_btm.Six.depths))) .* median(diff(sfomc.nw_w_btm.Six.date)*60*60*24);
  sfomc.nw_w_btm.Six.cprof = cumsum(sfomc.nw_w_btm.Six.xprof,'omitnan');
  
  sfomc.nw_w_btm.intSix.date = sfomc.nw_w_btm.Six.date;
  sfomc.nw_w_btm.intSix.data = sum(sfomc.nw_w_btm.Six.cprof,2,'omitnan');
  
  
  %%%
  % SW/C buoy (20 m isobath) nutrient fluxes
  
  % [NO3+NO2]
  sfomc.sw_buoy.Nx.date = sfomc.sw_buoy.adcp_onshore.date;
  sfomc.sw_buoy.Nx.depths = sfomc.sw_buoy.adcp_depths(end:-1:1);
  btmix = find(sfomc.sw_buoy.Nx.depths<sfomc.sw_buoy.N.depths(end));
  sfomc.sw_buoy.Nx.prof = interp2(sfomc.sw_buoy.N.date,sfomc.sw_buoy.N.depths,sfomc.sw_buoy.N.prof',...
                                  sfomc.sw_buoy.Nx.date,sfomc.sw_buoy.Nx.depths,'linear',NaN)';
  sfomc.sw_buoy.Nx.prof(:,btmix) = repmat(sfomc.sw_buoy.Nx.prof(:,btmix(1)-1),[1,numel(btmix)]);
  sfomc.sw_buoy.Nx.xprof = sfomc.sw_buoy.adcp_onshore.prof .* umol_to_kgm3(sfomc.sw_buoy.Nx.prof,'N') ...
      .* abs(median(diff(sfomc.sw_buoy.Nx.depths))) .* median(diff(sfomc.sw_buoy.Nx.date)*60*60*24);
  sfomc.sw_buoy.Nx.xprof(:,btmix) = repmat(sfomc.sw_buoy.Nx.xprof(:,btmix(1)-1),[1,numel(btmix)]);
  sfomc.sw_buoy.Nx.cprof = cumsum(sfomc.sw_buoy.Nx.xprof,'omitnan');
  sfomc.sw_buoy.Nx.cprof(:,btmix) = repmat(sfomc.sw_buoy.Nx.cprof(:,btmix(1)-1),[1,numel(btmix)]);
  
  sfomc.sw_buoy.intNx.date = sfomc.sw_buoy.Nx.date;
  sfomc.sw_buoy.intNx.data = sum(sfomc.sw_buoy.Nx.cprof,2,'omitnan');
  
  % [P]
  sfomc.sw_buoy.Px.date = sfomc.sw_buoy.adcp_onshore.date;
  sfomc.sw_buoy.Px.depths = sfomc.sw_buoy.adcp_depths(end:-1:1);
  btmix = find(sfomc.sw_buoy.Px.depths<sfomc.sw_buoy.P.depths(end));
  sfomc.sw_buoy.Px.prof = interp2(sfomc.sw_buoy.P.date,sfomc.sw_buoy.P.depths,sfomc.sw_buoy.P.prof',...
                                  sfomc.sw_buoy.Px.date,sfomc.sw_buoy.Px.depths,'linear',NaN)';
  sfomc.sw_buoy.Px.prof(:,btmix) = repmat(sfomc.sw_buoy.Px.prof(:,btmix(1)-1),[1,numel(btmix)]);
  sfomc.sw_buoy.Px.xprof = sfomc.sw_buoy.adcp_onshore.prof .* umol_to_kgm3(sfomc.sw_buoy.Px.prof,'P') ...
      .* abs(median(diff(sfomc.sw_buoy.Px.depths))) .* median(diff(sfomc.sw_buoy.Px.date)*60*60*24);
  sfomc.sw_buoy.Px.xprof(:,btmix) = repmat(sfomc.sw_buoy.Px.xprof(:,btmix(1)-1),[1,numel(btmix)]);
  sfomc.sw_buoy.Px.cprof = cumsum(sfomc.sw_buoy.Px.xprof,'omitnan');
  sfomc.sw_buoy.Px.cprof(:,btmix) = repmat(sfomc.sw_buoy.Px.cprof(:,btmix(1)-1),[1,numel(btmix)]);
  
  sfomc.sw_buoy.intPx.date = sfomc.sw_buoy.Px.date;
  sfomc.sw_buoy.intPx.data = sum(sfomc.sw_buoy.Px.cprof,2,'omitnan');
  
  % [Si]
  sfomc.sw_buoy.Six.date = sfomc.sw_buoy.adcp_onshore.date;
  sfomc.sw_buoy.Six.depths = sfomc.sw_buoy.adcp_depths(end:-1:1);
  btmix = find(sfomc.sw_buoy.Six.depths<sfomc.sw_buoy.Si.depths(end));
  sfomc.sw_buoy.Six.prof = interp2(sfomc.sw_buoy.Si.date,sfomc.sw_buoy.Si.depths,sfomc.sw_buoy.Si.prof',...
                                  sfomc.sw_buoy.Six.date,sfomc.sw_buoy.Six.depths,'linear',NaN)';
  sfomc.sw_buoy.Six.prof(:,btmix) = repmat(sfomc.sw_buoy.Six.prof(:,btmix(1)-1),[1,numel(btmix)]);
  sfomc.sw_buoy.Six.xprof = sfomc.sw_buoy.adcp_onshore.prof .* umol_to_kgm3(sfomc.sw_buoy.Six.prof,'Si') ...
      .* abs(median(diff(sfomc.sw_buoy.Six.depths))) .* median(diff(sfomc.sw_buoy.Six.date)*60*60*24);
  sfomc.sw_buoy.Six.xprof(:,btmix) = repmat(sfomc.sw_buoy.Six.xprof(:,btmix(1)-1),[1,numel(btmix)]);
  sfomc.sw_buoy.Six.cprof = cumsum(sfomc.sw_buoy.Six.xprof,'omitnan');
  sfomc.sw_buoy.Six.cprof(:,btmix) = repmat(sfomc.sw_buoy.Six.cprof(:,btmix(1)-1),[1,numel(btmix)]);
  
  sfomc.sw_buoy.intSix.date = sfomc.sw_buoy.Six.date;
  sfomc.sw_buoy.intSix.data = sum(sfomc.sw_buoy.Six.cprof,2,'omitnan');
  
  
  
  %%%
  % NE/E buoy (50 m isobath) nutrient fluxes
  
  % [NO3+NO2]
  sfomc.ne_buoy.Nx.date = sfomc.ne_buoy.adcp_onshore.date;
  sfomc.ne_buoy.Nx.depths = sfomc.ne_buoy.adcp_depths(end:-1:1);
  btmix = find(sfomc.ne_buoy.Nx.depths<sfomc.ne_buoy.N.depths(end));
  sfomc.ne_buoy.Nx.prof = interp2(sfomc.ne_buoy.N.date,sfomc.ne_buoy.N.depths,sfomc.ne_buoy.N.prof',...
                                  sfomc.ne_buoy.Nx.date,sfomc.ne_buoy.Nx.depths,'linear',NaN)';
  sfomc.ne_buoy.Nx.prof(:,btmix) = repmat(sfomc.ne_buoy.Nx.prof(:,btmix(1)-1),[1,numel(btmix)]);
  sfomc.ne_buoy.Nx.xprof = sfomc.ne_buoy.adcp_onshore.prof .* umol_to_kgm3(sfomc.ne_buoy.Nx.prof,'N') ...
      .* abs(median(diff(sfomc.ne_buoy.Nx.depths))) .* median(diff(sfomc.ne_buoy.Nx.date)*60*60*24);
  sfomc.ne_buoy.Nx.xprof(:,btmix) = repmat(sfomc.ne_buoy.Nx.xprof(:,btmix(1)-1),[1,numel(btmix)]);
  sfomc.ne_buoy.Nx.cprof = cumsum(sfomc.ne_buoy.Nx.xprof,'omitnan');
  sfomc.ne_buoy.Nx.cprof(:,btmix) = repmat(sfomc.ne_buoy.Nx.cprof(:,btmix(1)-1),[1,numel(btmix)]);
  
  sfomc.ne_buoy.intNx.date = sfomc.ne_buoy.Nx.date;
  sfomc.ne_buoy.intNx.data = sum(sfomc.ne_buoy.Nx.cprof,2,'omitnan');
  
  % [P]
  sfomc.ne_buoy.Px.date = sfomc.ne_buoy.adcp_onshore.date;
  sfomc.ne_buoy.Px.depths = sfomc.ne_buoy.adcp_depths(end:-1:1);
  btmix = find(sfomc.ne_buoy.Px.depths<sfomc.ne_buoy.P.depths(end));
  sfomc.ne_buoy.Px.prof = interp2(sfomc.ne_buoy.P.date,sfomc.ne_buoy.P.depths,sfomc.ne_buoy.P.prof',...
                                  sfomc.ne_buoy.Px.date,sfomc.ne_buoy.Px.depths,'linear',NaN)';
  sfomc.ne_buoy.Px.prof(:,btmix) = repmat(sfomc.ne_buoy.Px.prof(:,btmix(1)-1),[1,numel(btmix)]);
  sfomc.ne_buoy.Px.xprof = sfomc.ne_buoy.adcp_onshore.prof .* umol_to_kgm3(sfomc.ne_buoy.Px.prof,'P') ...
      .* abs(median(diff(sfomc.ne_buoy.Px.depths))) .* median(diff(sfomc.ne_buoy.Px.date)*60*60*24);
  sfomc.ne_buoy.Px.xprof(:,btmix) = repmat(sfomc.ne_buoy.Px.xprof(:,btmix(1)-1),[1,numel(btmix)]);
  sfomc.ne_buoy.Px.cprof = cumsum(sfomc.ne_buoy.Px.xprof,'omitnan');
  sfomc.ne_buoy.Px.cprof(:,btmix) = repmat(sfomc.ne_buoy.Px.cprof(:,btmix(1)-1),[1,numel(btmix)]);
  
  sfomc.ne_buoy.intPx.date = sfomc.ne_buoy.Px.date;
  sfomc.ne_buoy.intPx.data = sum(sfomc.ne_buoy.Px.cprof,2,'omitnan');
  
  % [Si]
  sfomc.ne_buoy.Six.date = sfomc.ne_buoy.adcp_onshore.date;
  sfomc.ne_buoy.Six.depths = sfomc.ne_buoy.adcp_depths(end:-1:1);
  btmix = find(sfomc.ne_buoy.Six.depths<sfomc.ne_buoy.Si.depths(end));
  sfomc.ne_buoy.Six.prof = interp2(sfomc.ne_buoy.Si.date,sfomc.ne_buoy.Si.depths,sfomc.ne_buoy.Si.prof',...
                                  sfomc.ne_buoy.Six.date,sfomc.ne_buoy.Six.depths,'linear',NaN)';
  sfomc.ne_buoy.Six.prof(:,btmix) = repmat(sfomc.ne_buoy.Six.prof(:,btmix(1)-1),[1,numel(btmix)]);
  sfomc.ne_buoy.Six.xprof = sfomc.ne_buoy.adcp_onshore.prof .* umol_to_kgm3(sfomc.ne_buoy.Six.prof,'Si') ...
      .* abs(median(diff(sfomc.ne_buoy.Six.depths))) .* median(diff(sfomc.ne_buoy.Six.date)*60*60*24);
  sfomc.ne_buoy.Six.xprof(:,btmix) = repmat(sfomc.ne_buoy.Six.xprof(:,btmix(1)-1),[1,numel(btmix)]);
  sfomc.ne_buoy.Six.cprof = cumsum(sfomc.ne_buoy.Six.xprof,'omitnan');
  sfomc.ne_buoy.Six.cprof(:,btmix) = repmat(sfomc.ne_buoy.Six.cprof(:,btmix(1)-1),[1,numel(btmix)]);
  
  sfomc.ne_buoy.intSix.date = sfomc.ne_buoy.Six.date;
  sfomc.ne_buoy.intSix.data = sum(sfomc.ne_buoy.Six.cprof,2,'omitnan');
  
  
  % % OOPS! Much of the time, upwelled water seems to be moving OFFshore?!
  % fmg; boxplot(sfomc.sw_buoy.adcp_x.prof(sfomc.sw_buoy.Nx.xprof(:,end)>0,:),sfomc.sw_buoy.adcp_depths,'notch','on');
  %  titlename('Mid-Station upwelling u\bullet\nablah profile');
  % fmg; boxplot(sfomc.ne_buoy.adcp_x.prof(sfomc.ne_buoy.Nx.xprof(:,end)>0,:),sfomc.ne_buoy.adcp_depths,'notch','on');
  %  titlename('East-Station upwelling u\bullet\nablah profile');
  
  % Let's focus on the moments right at the start (and end) of an upwelling event...
  %upwix = find(abs(sfomc.ne_buoy.at_45m.sbe_seatemp.data-23)<2e-4);
  %upwix = find((0 < sfomc.ne_buoy.Nx.prof(:,end)) & (sfomc.ne_buoy.Nx.prof(:,end)<1e-4));
  upwix = find(sfomc.ne_buoy.Nx.prof(1:end-1,end)==0 & sfomc.ne_buoy.Nx.prof(2:end,end)>0);
  dwnix = find(sfomc.ne_buoy.Nx.prof(1:end-1,end)>0 & sfomc.ne_buoy.Nx.prof(2:end,end)==0);
  
  % Sanity check (plotting code was MOVED to plot_upwelling_ekman_flux_vs_seatemp.m)
  
end; %if ( doNuts )



% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SHEAR AND STRATIFICATION: Richardson Number (Ri)

% Adapated from plot_upwelling_N2.m, 2018 Mar 06

if ( doRi )

  %for cstnm = {'sw_buoy'};
  for cstnm = {'ne_buoy','e_buoy','c_buoy','sw_buoy','nw_w_btm'};
    stnm = cstnm{:};
    
    % Calculate from hourly-interpolated values
    ws = warning('OFF','MATLAB:interp1:NaNstrip');
    spd = interp_ts(sfomc.(stnm).adcp_speed);
    warning(ws);
    tmp = interp_ts(sfomc.(stnm).seatemp);
    if ( isfield(sfomc.(stnm),'salin') )
      sal = interp_ts(sfomc.(stnm).salin);
    else
      sal = tmp;
      sal.data(:) = 35;
      sal.prof(:) = 35;
    end;
    [spd,tmp,sal] = intersect_tses(spd,tmp,sal);
    % tmp = interp_ts(tmp);
    % spd = interp_ts(spd);
    
    %[N2,Pv] = gsw_Nsquared(repmat(35,[numel(15:30),4]),repmat([15:30]',[1,4]),[5,15,25,35],26)
    t = tmp.date;
    T = tmp.prof(:,end:-1:1);
    %S = repmat(35,size(T));
    S = sal.prof(:,end:-1:1);
    z = -sfomc.(stnm).seatemp.depths(end:-1:1);
    
    
    % Calculate Brunt-Väisälä (buoyancy) frequency squared: N^2
    sfomc.(stnm).brunt_vaisala_squared.date = tmp.date;
    [n2,Pv] = gsw_Nsquared(S',T',z',sfomc.(stnm).lat);
    sfomc.(stnm).brunt_vaisala_squared.prof = n2';
    sfomc.(stnm).brunt_vaisala_squared.depths = -Pv(:,1)';

    spd.depths = -z;
    wd = warning('OFF','MATLAB:interp2:NaNstrip');
    spd.interp_prof = interp2(spd.date,sfomc.(stnm).adcp_depths,spd.prof',spd.date,spd.depths,'linear')';
    warning(wd);

    % Calculate vertical velocity shear squared: (du/dz)^2
    sfomc.(stnm).shear_squared.date = spd.date;
    sfomc.(stnm).shear_squared.prof = ((diff(spd.interp_prof')./diff(spd.depths)')').^2;
    sfomc.(stnm).shear_squared.depths = sfomc.(stnm).brunt_vaisala_squared.depths;

    % Calculate Richardson number => Ri < 0.25 suggesting Kelvin-Helmholtz instability
    sfomc.(stnm).richardson_number.date = sfomc.(stnm).brunt_vaisala_squared.date;
    sfomc.(stnm).richardson_number.prof = sfomc.(stnm).brunt_vaisala_squared.prof./sfomc.(stnm).shear_squared.prof;
    sfomc.(stnm).richardson_number.depths = sfomc.(stnm).brunt_vaisala_squared.depths;
    sfomc.(stnm).richardson_number.prof(sfomc.(stnm).richardson_number.prof<0) = 0;

    % % Sanity check
    % if ( ~isscalar(sfomc.(stnm).richardson_number.depths) )
    %   fmg; surf(sfomc.(stnm).richardson_number.date,sfomc.(stnm).seatemp.depths(2:end)',sfomc.(stnm).richardson_number.prof'); datetick3; shading interp; set(gca,'CLim',[0,0.25]); colorbar; titlename(textize(stnm));
    % end;

  end; %for cstnm = {'ne_buoy','e_buoy','c_buoy','sw_buoy','nw_w_btm'};

end; %if ( doRi )


if ( doMaps )

  disp('Plotting bathymetry...');

  clear blowupmap;

  fhs=[];
  [lkwf1,ig,ig,ig,fhs(1)] = plot_hires_bathymetry(lkwf1,-[0:2:20,30:20:200]);
  set(gca,'FontSize',fontsz);
  titlename('SE Florida shelf - USGS Bathymetry - instrumentation');

  % Draw red rectangle around close-up area (see last map below)
  if ( exist('sfomc','var') )
    rectangle('Position',bbox2rect(field_bbox(sfomc.c_buoy.ngdc_hires_bathy)),...
              'EdgeColor','r','LineWidth',2);
  end;

  warning('Skipping close-up bathymetry maps for now - too memory-intensive');
  % if ( exist('sefcri','var') )
  %   [sefcri.dc3,ig,ig,ig,fhs(end+1)] = plot_hires_bathymetry(sefcri.dc3,-[0:2:20,30:10:50]);
  %   set(gca,'FontSize',fontsz);
  %   titlename('SE Florida shelf (southern) - USGS Bathymetry - instrumentation');
  % end;
  %
  % [fdepk,ig,ig,ig,fhs(end+1)] = plot_hires_bathymetry(fdepk,-[0:2:20,30:10:50]);
  % set(gca,'FontSize',fontsz);
  % titlename('SE Florida shelf (northern) - USGS Bathymetry - instrumentation');
  %
  % %[canf1,ig,ig,ig,fhs(end+1)] = plot_hires_bathymetry(canf1,-[0:2:20,30:10:50]);
  % %titlename('Central Florida shelf - USGS Bathymetry - instrumentation');

  fhix = 1;
  figure(fhs(fhix));
  plot_all_se_florida_upwelling_sites;
  set(gca,'FontSize',fontsz);
  title([]);
  xlim([-80.37,-79.85]); daspect([1,cosd(25),1]);
  grid on; grid minor;
keyboard;
  if doPrint; print('-dpng',fullfile(figspath,['upwelling-bathymetry.png'])); end;

  % for fhix = 2:numel(fhs)
  %   figure(fhs(fhix));
  %   set(gca,'FontSize',fontsz);
  %   plot_all_se_florida_upwelling_sites;
  %   title([]);
  %   if doPrint; print('-dpng',fullfile(figspath,['upwelling-bathymetry-',num2str(fhix),'.png'])); end;
  % end; %for fhix = 2:numel(fhs)

  clear fhix

  % High-resolution map of Alex Soloviev's SFOMC section
  if ( exist('sfomc','var') )
    sfomc.c_buoy = plot_hires_bathymetry(sfomc.c_buoy,-[0:1:30,40:10:120],[3e3,11e3]);
    set(gca,'FontSize',fontsz);
    blowupmap=1; plot_all_se_florida_upwelling_sites; clear blowupmap;
    %titlename('Area of Port Everglades - USGS Bathymetry - instrumentation');
    %ylim([26.000,26.165]); daspect([1,cosd(25),1]);
    ylim([26.005,26.165]); daspect([1,cosd(25),1]);
    title([]);
    grid on; grid minor;
    if doPrint; print('-dpng',fullfile(figspath,['upwelling-bathymetry-sfomc.png'])); end;
  end;

end; %if ( doMaps )


% NOTE: GROUP OF TIME SERIES PLOTS WITH, e.g., TITLE 'Far North (Martin County)'
% OR LEGEND 'TCM1','TCM2', etc. WERE MOVED TO TEMP.m on 2017-Jul-04


if ( doCSV )
  dump_upwelling_csv;
end;


if ( doGrid )
  matfname = fullfile(datapath,'upwelling_ekman_flux_grid.mat');
  if ( exist(matfname,'file') )
    disp(['Loading ',matfname]);
    grd = load(matfname);
    
  else
    disp('Extracting grid');
    
    grd=[]; clear grd
    grd.all_lons = [canf1.lon,cnnf1.lon,face.lons,face.NF08.lon',face.NF09.lon',fdepk.lon,ftpf1.lon,fwc.lons,fwyf1.lon,lkwf1.lon,pvgf1.lon,rsm.lon,sefcri.lons,sfomc.lons,sfp.longitude.data'];
    grd.all_lats = [canf1.lat,cnnf1.lat,face.lats,face.NF08.lat',face.NF09.lat',fdepk.lat,ftpf1.lat,fwc.lats,fwyf1.lat,lkwf1.lat,pvgf1.lat,rsm.lat,sefcri.lats,sfomc.lats,sfp.latitude.data'];

    grd.all_dates = [];
    stcs = {canf1,cnnf1,face,fdepk,ftpf1,fwc,fwyf1,lkwf1,pvgf1,rsm,sefcri,sfomc};
    for stix=1:numel(stcs)
      stc = stcs{stix};
      if ( isfield(stc,'station_name') )
        disp(upper(stc.station_name));
      end;
      flds = fieldnames(stc);
      for fldix = 1:numel(flds)
        fld = flds{fldix};
        if ( isfield(stc.(fld),'date') )
          %DEBUG:      disp(fld);
          grd.all_dates = [grd.all_dates,stc.(fld).date(:)'];
        elseif ( isstruct(stc.(fld)) && numel(stc.(fld)) == 1 )
          subflds = fieldnames(stc.(fld));
          for subfldix = 1:numel(subflds)
            subfld = subflds{subfldix};
            if ( isstruct(stc.(fld).(subfld)) && numel(stc.(fld).(subfld)) == 1 && isfield(stc.(fld).(subfld),'date') )
              %DEBUG:          disp([fld,'.',subfld]);
              grd.all_dates = [grd.all_dates,stc.(fld).(subfld).date(:)'];
            end;
          end;
        end;
      end;
    end;
    grd.all_dates = [grd.all_dates,face.NF08.dt',face.NF09.dt',sfp.latitude.date'];
    
    grd.all_lons_sorted = unique(grd.all_lons);
    grd.all_lats_sorted = unique(grd.all_lats);
    grd.all_dates_sorted = unique(grd.all_dates);
    
    disp(['Saving ',matfname]);
    save(matfname,'-struct','grd');
  end;
end; %if ( doGrid )



% load(fullfile(get_ecoforecasts_path('data'),'T_tses_1.mat'),'ts');
% ts = copy_T_ts(face,ts,'FACE_HWD_');
% ts = copy_T_ts(sefcri,ts,'SEFCRI_');
% ts = copy_T_ts(sfomc,ts,'SFOMC_');
% save(fullfile(get_ecoforecasts_path('data'),'T_tses.mat'),'ts');

clear ans ax cderf dpax fhs fld fldstub h hrs ig lh shax
