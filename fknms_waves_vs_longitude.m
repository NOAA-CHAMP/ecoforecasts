1;

if ( ~exist('doPrint','var') || isempty(doPrint) )
  doPrint = false;
end;
if ( ~exist('figspath','var') )
  %figspath = get_ecoforecasts_path('figs');
  figspath = get_coral_path('Omics');
end;

if ( ~exist('stns','var') )
  anfknms;
end;
if ( ~isfield(stns.MLRF1,'raw_nwps_sigwavehgt') )
  w = warning('OFF','Ecoforecasts:NWPS:NoFile');
  %stns = get_nwps_stations(stns,[],[],false);
  stns = get_nwps_stations(stns,[],[],true);
  warning(w);
end;

if ( ~exist('stn','var') )
  stn = get_station_from_station_name('fwyf1'); stn = load_all_ndbc_data(stn);
end;

fmg; plot_ts(stn.ndbc_wind1_speed);
xlim(stns.MLRF1.raw_nwps_sigwavehgt.date([1,end])); datetick3;
if doPrint; print('-dpng',fullfile(figspath,'fknms_waves_vs_longitude_wind_2016_Aug_Sep.png')); end;

fmg;
for stnm=sites.stnms;
  plot3(stns.(stnm{:}).raw_nwps_sigwavehgt.date,repmat(stns.(stnm{:}).lon,[numel(stns.(stnm{:}).raw_nwps_sigwavehgt.data),1]),stns.(stnm{:}).raw_nwps_sigwavehgt.data,'.');
end;
datetick3;
view(50,55);
if doPrint; print('-dpng',fullfile(figspath,'fknms_waves_vs_longitude_50_55.png')); end;
disp('Hit Enter...'); pause;

xlim(stns.MLRF1.raw_nwps_sigwavehgt.date([1,262])); datetick3;
view(90,0);
disp('Hit Enter...'); pause;

xlim(stns.MLRF1.raw_nwps_sigwavehgt.date([262,end])); datetick3;
view(90,0);
if doPrint; print('-dpng',fullfile(figspath,'fknms_waves_vs_longitude_90_0_2016_Sep.png')); end;


fmg;
%for stix=1:3:numel(sites.stnms);
for stix=1:1:numel(sites.stnms);
  stnm = sites.stnms(stix);
  stns.(stnm{:}).raw_nwps_ardhuin_surface_drift_speed.date = stns.(stnm{:}).raw_nwps_ardhuin_surface_drift_u.date;
  stns.(stnm{:}).raw_nwps_ardhuin_surface_drift_speed.data = ...
      uv_to_spd(stns.(stnm{:}).raw_nwps_ardhuin_surface_drift_u.data,...
                stns.(stnm{:}).raw_nwps_ardhuin_surface_drift_v.data);
  %lh = plot3(stns.(stnm{:}).raw_nwps_ardhuin_surface_drift_speed.date,repmat(stns.(stnm{:}).lon,[numel(stns.(stnm{:}).raw_nwps_ardhuin_surface_drift_speed.data),1]),stns.(stnm{:}).raw_nwps_ardhuin_surface_drift_speed.data,'.');
  %set(lh,'MarkerSize',10);

  stns.(stnm{:}) = verify_variable(stns.(stnm{:}),{'raw_nwps_ardhuin_surface_drift_u_5_d_sum','raw_nwps_ardhuin_surface_drift_v_5_d_sum'});
  stns.(stnm{:}).raw_nwps_ardhuin_surface_drift_speed_5_d_sum.date = stns.(stnm{:}).raw_nwps_ardhuin_surface_drift_u_5_d_sum.date;
  stns.(stnm{:}).raw_nwps_ardhuin_surface_drift_speed_5_d_sum.data = ...
      uv_to_spd(stns.(stnm{:}).raw_nwps_ardhuin_surface_drift_u_5_d_sum.data,...
                stns.(stnm{:}).raw_nwps_ardhuin_surface_drift_v_5_d_sum.data);

  ts = ts_nanify_gaps(stns.(stnm{:}).raw_nwps_ardhuin_surface_drift_speed_5_d_sum,3);
  lh = plot3(ts.date,repmat(stns.(stnm{:}).lon,[numel(ts.data),1]),ts.data.*3600/1e3,'-');
  set(lh,'LineWidth',2);
end;
%datetick3;
datetick3('x','mmmm');
view(75,66);
ylabel('Longitude'); zlabel('Travel [km]');
if doPrint; print('-dpng',fullfile(figspath,'fknms_stokes_vs_longitude_75_66.png')); end;

disp('Hit Enter...'); pause;
axis([datenum(2016,8,31),datenum(2016,10,1),-80.3,-80]); view(3);
if doPrint; print('-dpng',fullfile(figspath,'fknms_waves_vs_longitude_UK_75_66_2016_Sep.png')); end;

disp('Hit Enter...'); pause;
view(0,0);
if doPrint; print('-dpng',fullfile(figspath,'fknms_waves_vs_longitude_UK_0_0_2016_Sep.png')); end;


% %bbox = [-80.5,-80.2,25.0,25.3];
% bbox = [-80.9,-80.5,24.7,25.1];
% inix = bboxinside(sites.lons,sites.lats,bbox,false,0);
inix = 1:numel(sites.lons);

bath=[]; clear bath
%bath = read_hires_bathymetry_for_field({sites.lons,sites.lats},false,10e3);
%plot_hires_bathymetry(bath);
%axis(bbox);
bath = read_hires_bathymetry_for_field({sites.lons(inix),sites.lats(inix)},false,20e3);
bath.ngdc_hires_bathy.field(bath.ngdc_hires_bathy.field>=-0.25) = nan;
plot_hires_bathymetry(bath,-[0:10:80,100:100:1000]);
%bath = read_hires_bathymetry_for_field({sites.lons(inix),sites.lats(inix)},true,10e3);
%bath.ngdc_hires_bathy.field(bath.ngdc_hires_bathy.field>=-0.25) = nan;
%plot_hires_bathymetry(bath);
axis(axis);
disp('Drawn');

%dtix = find(stns.FKNMS_7MILE_BR.raw_nwps_ardhuin_surface_drift_u.date>=datenum(2016,8,16));
dtix = find(stns.FKNMS_7MILE_BR.raw_nwps_ardhuin_surface_drift_u.date>=datenum(2016,8,16)...
            & stns.FKNMS_7MILE_BR.raw_nwps_ardhuin_surface_drift_u.date<=datenum(2016,9,1));
% dtix = find(stns.FKNMS_7MILE_BR.raw_nwps_ardhuin_surface_drift_u.date>=datenum(2016,8,16)...
%             & stns.FKNMS_7MILE_BR.raw_nwps_ardhuin_surface_drift_u.date<=datenum(2016,8,26));
% dtix = find(stns.FKNMS_7MILE_BR.raw_nwps_ardhuin_surface_drift_u.date>=datenum(2016,8,23)...
%             & stns.FKNMS_7MILE_BR.raw_nwps_ardhuin_surface_drift_u.date<=datenum(2016,9,1));

for stix=inix(:)';
  stnm = sites.stnms(stix);
  if ( sites.depths(stix) > -2 )
    continue;
  end;
  %disp(stnm{:});

  [cu,cv] = spddir_to_uv(stns.(stnm{:}).raw_nwps_currspeed.data(dtix),...
                       stns.(stnm{:}).raw_nwps_currdir.data(dtix));
  if ( all(cu==0) || all(cv==0) )
    %disp(stnm{:});
    continue;
  end;
  lons = stns.(stnm{:}).lon + cumsum(cosd(25)*cu*3600./111e3);
  lats = stns.(stnm{:}).lat + cumsum(cv*3600./111e3);
  crlh(stix) = plot(lons,lats,'b-','LineWidth',2);
  clear lons lats

  [wu,wv] = spddir_to_uv(stns.(stnm{:}).raw_nwps_windspeed.data(dtix).*0.015,...
                       stns.(stnm{:}).raw_nwps_winddir.data(dtix));
  wu = cu + wu;
  wv = cv + wv;
  lons = stns.(stnm{:}).lon + cumsum(cosd(25)*wu*3600./111e3);
  lats = stns.(stnm{:}).lat + cumsum(wv*3600./111e3);
  wdlh(stix) = plot(lons,lats,'k-','LineWidth',2,'Color',[0.0,0.5,0.0]);
  clear lons lats

  su = cu - stns.(stnm{:}).raw_nwps_ardhuin_surface_drift_u.data(dtix);
  sv = cv - stns.(stnm{:}).raw_nwps_ardhuin_surface_drift_v.data(dtix);
  lons = stns.(stnm{:}).lon + cumsum(cosd(25)*su*3600./111e3);
  lats = stns.(stnm{:}).lat + cumsum(sv*3600./111e3);
  sdlh(stix) = plot(lons,lats,'r-','LineWidth',2);
  clear lons lats

  clear cu cv wu wv su sv
end;
legh = legend([crlh(1),wdlh(1),sdlh(1)],{'Currents','Windage','Stokes'}, 'Location','North');
titlename('Larval pathways: Wind, waves, and currents');

if doPrint; print('-dpng',fullfile(figspath,'fknms_waves_vs_longitude_pathways.png')); end;

delete(crlh); delete(wdlh); delete(sdlh); clear crlh wdlh sdlh
