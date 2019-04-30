1;

if ( ~exist('doPrint','var') )
  doPrint = false;
end;

% Read GoM HYCOM data
if ( ~exist('hflds','var') )
  hflds = read_hycom_gom_reanalysis;
end;

if ( ~exist('wflds','var') )
  w = warning('OFF','Ecoforecasts:NWPS:NoFile');
  wflds = get_nwps_fields([],[],datenum(2016,[8,9],1));
  warning(w); clear w
end;

hdts = hflds.u.date;
hlat = hflds.u.lat;
hlon = hflds.u.lon;
hu = hflds.u.field;
hv = hflds.v.field;

wdts = wflds.windspeed.date;
wlat = wflds.windspeed.lat;
wlon = wflds.windspeed.lon;

% Restrict HYCOM domain to NWPS domain
hu = hu(:,min(wlat) <= hlat & hlat <= max(wlat),min(wlon) <= hlon & hlon <= max(wlon));
hv = hv(:,min(wlat) <= hlat & hlat <= max(wlat),min(wlon) <= hlon & hlon <= max(wlon));
hlat = hlat(min(wlat) <= hlat & hlat <= max(wlat));
hlon = hlon(min(wlon) <= hlon & hlon <= max(wlon));

% Either extract drifts from netCDF, or recalculate them by hand in MATLAB

% Ardhuin et al. 2009 method
a = wflds.ardhuin_surface_drift.field;
au = wflds.ardhuin_surface_drift_u.field;
av = wflds.ardhuin_surface_drift_v.field;

% a = stokes_drift(wflds.windspeed.field(:),wflds.sigwavehgt.field(:),wflds.primwaveper.field(:),'ardhuin');
% a = reshape(a,size(wflds.windspeed.field));

% Monismith and Fong 2004 method
% m = wflds.monismith_surface_drift.field;
% mu = wflds.monismith_surface_drift_u.field;
% mv = wflds.monismith_surface_drift_v.field;

% h = 20;
h = 5;
% h = 3;
m = stokes_drift(wflds.windspeed.field(:),wflds.sigwavehgt.field(:),wflds.primwaveper.field(:),'monismith',h);
m = reshape(m,size(wflds.windspeed.field));
[mu,mv] = spddir_to_uv(m,wflds.winddir.field);


% Hourly Stokes
wu = interp3(wlat,wdts,wlon,mu,hlat,hdts,hlon);
wv = interp3(wlat,wdts,wlon,mv,hlat,hdts,hlon);

% Total surface current
u = wu + hu;
v = wv + hv;

% Speeds (HYCOM, Stokes, and total) from U and V
hspd = uv_to_spd(hu,hv);
wspd = uv_to_spd(wu,wv);
spd = uv_to_spd(u,v);

% Stokes as a percent of total
spdpct = wspd ./ spd;

%wsmooth = smooth3(wspd,'gaussian');

if 0;
  fmg;
  contour_field(wflds.monismith_surface_drift,@contourf,[],{@nanmean},[0:0.02:0.20]);
  set(gca,'CLim',[0,0.20]);
  titlename('Mean Monismith Lagrangian velocity (Python)');
  fmg;
  contour_field(wflds.monismith_surface_drift,@contourf,[],{@nanmean},[0:0.02:0.20]);
  set(gca,'CLim',[0,0.20]);
  titlename('Mean Monismith Lagrangian velocity (MATLAB)');

  fmg;
  contour_field(wflds.monismith_surface_drift,@contourf,[],{@prctile,93},[0:0.02:0.20]);
  set(gca,'CLim',[0,0.20]);
  titlename('93^{rd} percentile Monismith Lagrangian velocity (Python)');
  fmg;
  contour_field(wflds.monismith_surface_drift,@contourf,[],{@prctile,93},[0:0.02:0.20]);
  set(gca,'CLim',[0,0.20]);
  titlename('93^{rd} percentile Monismith Lagrangian velocity (MATLAB)');

  fmg;
  contour_field(wflds.ardhuin_surface_drift,@contourf,[],{@nanmean},[0:0.02:0.20]);
  set(gca,'CLim',[0,0.20]);
  titlename('Mean Ardhuin Lagrangian velocity (Python)');
  fmg;
  contour_field(wflds.ardhuin_surface_drift,@contourf,[],{@nanmean},[0:0.02:0.20]);
  set(gca,'CLim',[0,0.20]);
  titlename('Mean Ardhuin Lagrangian velocity (MATLAB)');

  fmg;
  contour_field(wflds.ardhuin_surface_drift,[],{@prctile,93},[0:0.02:0.20]);
  set(gca,'CLim',[0,0.3]);
  titlename('93^{rd} percentile Ardhuin Lagrangian velocity (MATLAB)');
  fmg;
  contour_field(wflds.ardhuin_surface_drift,[],{@prctile,93},[0:0.02:0.20]);
  set(gca,'CLim',[0,0.3]);
  titlename('93^{rd} percentile Ardhuin Lagrangian velocity (MATLAB)');
end;


if 0;
  fmg; contourf(hlon,hlat,squeeze(prctile(hspd,93)),[0:0.100:2.000]); colorbar;
  set(gca,'CLim',[0,2.0]); daspect([1,cosd(wlat(1)),1]);
  titlename('93^{rd} percentile HYCOM surface current [m/s] August 2016');
  if doPrint; print('-dpng',['figs/hycom_gom_2016_Aug_93pct.png']); end;

  %fmg; contourf(wlon,wlat,squeeze(prctile(wspd,93)),[0:0.010:0.250]); colorbar;
  fmg; contourf(wlon,wlat,squeeze(prctile(wspd,93)),[0:0.010:2.000]); colorbar;
  set(gca,'CLim',[0,2.0]); daspect([1,cosd(wlat(1)),1]);
  titlename('93^{rd} percentile Stokes drift [m/s] August 2016');
  if doPrint; print('-dpng',['figs/monismith_',num2str(h),'m_2016_Aug_93pct.png']); end;

  set(gca,'CLim',[0,0.20]); daspect([1,cosd(wlat(1)),1]);
  if doPrint; print('-dpng',['figs/monismith_',num2str(h),'m_2016_Aug_93pct_detail.png']); end;


  fmg; contourf(wlon,wlat,squeeze(nanmean(wspd)),[0:0.010:0.100]); colorbar;
  titlename('Mean Stokes drift [m/s] August 2016');

  fmg; contourf(wlon,wlat,squeeze(prctile(wspd,93)),[0:0.010:0.200]); colorbar;
  daspect([1,cosd(wlat(1)),1]);
  titlename('Stokes drift as % of HYCOM (93^{rd} p''ctile) August 2016');
  if doPrint; print('-dpng',['figs/monismith_',num2str(h),'m_hycom_pct_2016_Aug_93pct.png']); end;
end;


%% Calculate cross-shore components of Stokes and HYCOM at each gridpoint
%  during timesteps with the greatest field-wide spatial mean of Stokes drift
%  as a percentage of HYCOM currents

if 1;
  % Get bathymetry
  disp('Getting bathymetric gradients');
  %[bath,rad] = read_hires_bathymetry_for_field(wflds.windspeed,true);
  [bath,rad] = read_hires_bathymetry_for_field(wflds.windspeed,false);
  bat = bath.ngdc_hires_bathy; bath=[]; clear bath
end;
if 0;
  % Calculate bathymetric gradient (i.e., to get cross-shore direction)
  %[LON,LAT] = meshgrid(bat.lon,bat.lat);
  % [aspect_deg,slope_deg,beta_y,beta_x] = gradientm(LAT,LON,bat.field,wgs84Ellipsoid);
  % aspect_deg=[]; slope_deg=[]; clear aspect_deg slope_deg
  hx = 1e3 .* distance_wgs84(bat.lat(1),bat.lon(1),bat.lat(1),bat.lon(2));
  hy = 1e3 .* distance_wgs84(bat.lat(1),bat.lon(1),bat.lat(2),bat.lon(1));
  [beta_x,beta_y] = gradientn(bat.field,3,hx,hy);
  hx=[]; hy=[]; clear hx hy
  % We want direction of steepest DESCENT - so invert the gradient
  %beta_x = -beta_x; beta_y = -beta_y;
  beta_x = beta_x; beta_y = -beta_y;
  beta = uv_to_spd(beta_x,beta_y);
  beta_deg = uv_to_dir_curr(beta_x,beta_y);
end;

if 1;
  % Downsample bathymetry
  ulon = bat.lon(4:6:end-2);
  ulat = bat.lat(4:6:end-2);
  ubat = undersample_field(bat.field,6,@nanmax);
  % Calculate bathymetric gradient (i.e., to get cross-shore direction)
  hx = 1e3 .* distance_wgs84(ulat(1),ulon(1),ulat(1),ulon(2));
  hy = 1e3 .* distance_wgs84(ulat(1),ulon(1),ulat(2),ulon(1));
  [ubeta_x,ubeta_y] = gradientn(ubat,3,hx,hy);
  hx=[]; hy=[]; clear hx hy
  % We want direction of steepest DESCENT - so invert the gradient
  %ubeta_x = -ubeta_x; ubeta_y = -ubeta_y;
  ubeta_x = ubeta_x; ubeta_y = -ubeta_y;
  ubeta = uv_to_spd(ubeta_x,ubeta_y);
  ubeta_deg = uv_to_dir_curr(ubeta_x,ubeta_y);
end;

if 1;
  % % Field-wide spatial mean Stokes/HYCOM percentage of each time step
  % mspdpct = nanmean(spdpct(:,:),2);
  % % Time steps where Stokes was the greatest field-wide percentage of HYCOM
  % %[mx,mxix] = max(mspdpct);
  % %cutoff_mspdpct = prctile(mspdpct,93);
  % cutoff_mspdpct = prctile(mspdpct,50);
  % hiix = find(mspdpct > cutoff_mspdpct & ismember(hdts,wdts));

  %DEBUG:
  hiix = find(ismember(hdts,wdts)); % "ALL": One field every 3 h

  % Interpolate 4 km fields to 552 m bathymetry - only at WAVE time stamps
  disp('Interpolating HYCOM current fields');
  % HYCOM West-to-East
  zhu = interp3(hlat,hdts(hiix),hlon,hu(hiix,:,:),ulat,hdts(hiix),ulon);
  % HYCOM South-to-North
  zhv = interp3(hlat,hdts(hiix),hlon,hv(hiix,:,:),ulat,hdts(hiix),ulon);
  zhd = uv_to_dir_curr(zhu,zhv);

  disp('Interpolating Stokes drift fields');
%{
  % Recalculate Monismith & Fong Stokes drift with *actual* (NANMAX) bottom depths!
  [d,n,m] = size(zhu);
  zU10 = interp3(hlat,hdts(hiix),hlon,hv(hiix,:,:),ulat,hdts(hiix),ulon);
wu = interp3(wlat,wdts,wlon,mu,hlat,hdts,hlon);

  zHs  = 
  zPp  = 
  for tix=1:d
    m_tmp = stokes_drift(wflds.windspeed.field(tix,:,:), ...
                         wflds.sigwavehgt.field(tix,:,:), ...
                         wflds.primwaveper.field(tix,:,:),'monismith',ubat(:,:));
    true_m(tix,:,:) = reshape(tmp_m,size(wflds.windspeed.field)); m_tmp=[]; clear tmp_m
    [true_mu,true_mv] = spddir_to_uv(true_m,wflds.winddir.field);
  end;
%}
  % Stokes West-to-East
  zwu = interp3(hlat,hdts(hiix),hlon,wu(hiix,:,:),ulat,hdts(hiix),ulon);
  % Stokes South-to-North
  zwv = interp3(hlat,hdts(hiix),hlon,wv(hiix,:,:),ulat,hdts(hiix),ulon);
  zwd = uv_to_dir_curr(zwu,zwv);

  [d,n,m] = size(zhu);

  shelf_ix = find(-150 <= ubat & ubat < 0);
  reef_line_ix = find(-80 <= ubat & ubat < -3);

  zhx = repmat(nan,size(zhu));
  zhl = repmat(nan,size(zhu));
  zwx = repmat(nan,size(zhu));
  zwl = repmat(nan,size(zhu));
  disp(['Calculating ',num2str(d),' cross-shore current fields']);
  for tix=1:d
    if ( mod(tix,floor(d/10)) == 0 ); disp(tix); end;
    [zhx(tix,:,:),zhl(tix,:,:)] = reorient_vectors(ubeta_deg,squeeze(zhu(tix,:,:)),squeeze(zhv(tix,:,:)));
    [zwx(tix,:,:),zwl(tix,:,:)] = reorient_vectors(ubeta_deg,squeeze(zwu(tix,:,:)),squeeze(zwv(tix,:,:)));
  end;
end;

if 1;
  % Median of ALL cross-shore Stokes magnitudes as a percentage of
  % cross-shore HYCOM magnitudes: 11%
  nansummary(abs(zwx)./(abs(zhx))),
  % For ALL points, Stokes is onshore while HYCOM is offshore: 10% of the time
  numel(find(zwx(:)<-0.0&zhx(:)>0.0))./numel(zwx),

  %% RESULTS:
  % Median of SHELF cross-shore Stokes magnitude as percentage
  % of cross-shore HYCOM magnitude: 20%
  nansummary(abs(zwx(:,shelf_ix))./(abs(zhx(:,shelf_ix)))),
  % For SHELF, Stokes is onshore while HYCOM is offshore: 15%
  % of the time
  numel(find(zwx(:,shelf_ix)<-0.0&zhx(:,shelf_ix)>0.0))./numel(zwx(:,shelf_ix)),

  %% RESULTS:
  % Median of REEF-LINE cross-shore Stokes magnitude as percentage
  % of cross-shore HYCOM magnitude: 20%
  nansummary(abs(zwx(:,reef_line_ix))./(abs(zhx(:,reef_line_ix)))),
  % For REEF-LINE, Stokes is onshore while HYCOM is offshore: 18%
  % of the time
  numel(find(zwx(:,reef_line_ix)<-0.0&zhx(:,reef_line_ix)>0.0))./numel(zwx(:,reef_line_ix)),

  %% RESULTS:
  % For REEF-LINE, Stokes is SIGNIFICANTLY onshore (<-2 cm/s) while
  % HYCOM is significantly offshore (>2 cm/s): 8% of the time
  numel(find(zwx(:,reef_line_ix)<-0.02&zhx(:,reef_line_ix)>0.02))./numel(zwx(:,reef_line_ix)),

  %% RESULTS:
  % REEF-LINE-wide MEAN cross-shore Stokes 83rd %-ile is -1 cm/s
  prctile(nanmean(zwx(:,reef_line_ix),2),83),
  % REEF-LINE-wide MEAN cross-shore Stokes gets below -8 cm/s: 31 Aug 09:00
  [minzwx,tix] = nanmin(nanmean(zwx(:,reef_line_ix),2)), datestr(wdts(tix)),
end;

if 1;
  fmg; hist(nanmean(zhx(:,reef_line_ix),2),100); axis([-0.12,+0.12,0,12]); titlename('HYCOM u^.\nablah');
  if doPrint; print('-dpng',['figs/hycom_gom_2016_Aug_reef_line_hist.png']); end;

  fmg; hist(nanmean(zwx(:,reef_line_ix),2),100); axis([-0.12,+0.12,0,12]); titlename('Stokes u^.\nablah');
  if doPrint; print('-dpng',['figs/monismith_',num2str(h),'m_2016_Aug_reef_line_hist.png']); end;
end;

if 0;
  %tix = find(nanmean(zwx(:,:),2)<-0.000065,1);
  [minzwx,tix] = nanmin(nanmean(zwx(:,reef_line_ix),2));

  fmg;
  contourf(wlon,wlat,squeeze(wspd(tix,:,:))); colorbar;
  % ??? PROBLEM: ARROWS SCALE DIFFERENTLY ON SUCCESSIVE CALLS?!
  h(1) = quiver(wlon,wlat,squeeze(wu(tix,:,:)),squeeze(wv(tix,:,:))); set(h(1),'Color','r');
  h(2) = quiver(wlon,wlat,squeeze(hu(tix,:,:)),squeeze(hv(tix,:,:))); set(h(2),'Color','k');
  axis([-81.3269 -80.0855 24.5793 25.5640]); set(gca,'clim',[-0.01,0.01]);
  legend(h,'Stokes','HYCOM');
  titlename(datestr(wdts(tix)));
  if doPrint; print('-dpng',['figs/monismith_',num2str(h),'m_2016_Aug_peak_reef_line_currents.png']); end;
end;

if 1;
  [WLON,WLAT] = meshgrid(wlon,wlat);
  WLON = [WLON;WLON];
  WLAT = [WLAT;WLAT];
  l.wu = wu; l.wv = wv; l.hu = hu; l.hv = hv;
  maxspd = 0.10;
  bigix = find(abs(l.wu)>maxspd|abs(l.wv)>maxspd|abs(l.hu)>maxspd|abs(l.hv)>maxspd);
  l.wu(bigix) = nan; l.wv(bigix) = nan; l.hu(bigix) = nan; l.hv(bigix) = nan;
  U = [squeeze(l.wu(tix,:,:));squeeze(l.hu(tix,:,:))];
  V = [squeeze(l.wv(tix,:,:));squeeze(l.hv(tix,:,:))];
  fmg;
  contourf(hlon,hlat,squeeze(hspd(tix,:,:))); colorbar;
  h = quiver(WLON,WLAT,U,V,0); set(h,'Color','r');
  axis([-81.3269 -80.0855 24.5793 25.5640]);
  titlename(datestr(wdts(tix)));
end;
