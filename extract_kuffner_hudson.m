1;
%% SCRIPT to parse and process thermistor data from Kuffner and Hudson 2014
%
% CONVERTED FROM: d:/Postdoc/extract_kuffner_hudson.m (2019 Mar 01, Lew.Gramer@gmail.com)
%
% Online Resource 2 - Electronic Supplementary material (ESM_2)
% Article title: A century of warming on Florida Keys coral reefs: Historic in-situ observations
% Journal name: Estuaries and Coasts
% Author names: Ilsa B. Kuffner, Barbara H Lidz, J. Harold Hudson, Jeffery S. Anderson
% Affiliation of corresponding author: U.S. Geological Survey, St. Petersburg Coastal & Marine Geology Science Center, 600 4th Street South, St. Petersburg, FL 33701
% Email address of corresponding author: ikuffner@usgs.gov
%
% Harold Hudson and Jeffery Anderson's historic underwater temperature data from the Florida Keys, FL, U.S.A.
% All temperature data are presented in degrees Celsius
% Pre-1990 data are presented graphically in:
% Hudson, J. H., R. B. Halley, A. J. Joseph, B. H. Lidz, and D. Schroeder. 1991. Long-term thermograph records from the Upper Florida Keys. U.S. Geological Survey, Open File Report 91-344.
% "Hourly data were collected using General Oceanics model 3011 thermogrpahs which were placed at sites for a maximum period of 120 days."
% Included here are daily means of the hourly measurements (hourly data cannot be located).
%
% Snake Creek coordinates: 24° 57.1' N 80° 35.3' W
% Snake Creek station site: "Thermograph placed in shallowhole near mid-channel on ocean side of bridge. Unit partially hidden by loose rock slab with sensor exposed on channel bottom."
% Snake Creek setting: " Tidal channel - Eroded into Key Largo limestone, the channel floor is mostly bedrock with interspersed patches of sediment, algae and seagrass."
% Snake Creek water depth of sensor: 3.5 m
%
% Last Saved Time-stamp: <Fri 2019-03-01 13:48:01 Eastern Standard Time gramer>


%usgspath = get_relative_path('USGS');
usgspath = fullfile('D:','Postdoc','USGS');

if ( ~exist('doPrint','var') || isempty(doPrint) )
  %doPrint = true;
  doPrint = false;
end;
if ( ~exist('doFigs','var') || isempty(doFigs) )
  %doFigs = true;
  doFigs = false;
end;


matfname = fullfile(get_ecoforecasts_path('data'),'kuffner_hudson.mat');
xlsfname = fullfile(usgspath,'ESM_2_Hudson&Anderson_select underwater temperature data 1975 to 2008.xls');

if ( exist(matfname,'file') )
  disp(['Load ',matfname]);
  load(matfname);
  
elseif ( ~exist(xlsfname,'file') )
  error('Sorry, no MAT file and I do not know where .XLS file is!');

else
  disp(['Extracting data...']);

  snaf1.station_name = 'snaf1';
  snaf1.station_long_name = 'Snake Creek';
  snaf1.lat = 24.951667;
  snaf1.lon = -80.58833;
  snaf1.depth = 3.5;
  
  % Hen and Chickens Reef coordinates: 24° 56.0' N 80° 33.0' W
  % Hen and Chickens Reef station site: "Thermograph concealed in open-ended cavidty within star coral (Montastraea annularis) 0.5 m above seafloor."
  % Hen and Chickens Reef setting: "Inshore patch reef - A mature patch reef dominated by the star coral with colonies up to 3 m high on level rock substrate. Reef area surround by generally thin sediment cover and seagrass."
  % Hen and Chickens Reef water depth of sensor: 4 m
  
  hacf1.station_name = 'hacf1';
  hacf1.station_long_name = 'Hen and Chickens Reef';
  hacf1.lat = 24.9333;
  hacf1.lon = -80.5500;
  hacf1.depth = 4;
  
  
  % Crocker Reef coordinates: 24° 54.4' N 80° 35.5' W
  % Crocker Reef station site: "Thermograph concealed under rock outcrop so that sensor is exposed at seafloor."
  % Crocker Reef setting: "Low-profile, relict bank on outer edge of reef tract. Rocky surface populated with sea fans (Gorgonia ventilina) and other soft coral species with scattered, small colonies of hard corals (including Dichocoenia stokesii and Meandrina meandrities)."
  % Crocker Reef water depth of sensor: 5 m
  
  ccrf1.station_name = 'ccrf1';
  ccrf1.station_long_name = 'Crocker Crest Reef';
  ccrf1.lat = 24.9067;
  ccrf1.lon = -80.5917;
  ccrf1.depth = 5;
  
  
  % Cement Dome coordinates: 25° 27.8' N 80° 10.0' W	Note there must have been a typo in the USGS Open File Report, because coordinates cannot be as reported: 25° 27.8' N 82° 10.0' W
  % Cement Dome station site: "The thermograph is concealed beneath a 1 m diameter artificial coral head. Sensor on the seafloor near one of several openings at base of artificial coral head."
  % Cement Dome setting: " An artificial patch reef constructed from cement domes appearing as coral heads and seeded with transplanted stony and soft corals. Patch is situated on a rock bottom naturally populated by gorgonians and alcionarians."
  % Cement Dome water depth of sensor: 4.0 m
  
  cmdf1.station_name = 'cmdf1';
  cmdf1.station_long_name = 'Cement Dome';
  cmdf1.lat = 25.4633;
  cmdf1.lon = -80.1667;
  cmdf1.depth = 4;
  
  
  % Carysfort Reef coordinates: 25° 12.6' N 80° 12.0' W
  % Carysfort Reef station site:			"The thermograph was strapped to a cross member of the lighthouse in 1-m depth of seawater" (Gene Shinn, personal communication)
  % Carysfort Reef setting:
  % Carysfort Reef water depth of sensor: 1 m
  % Data post-1987 for Carysfort Reef are from NOAA National Oceanographic Data Center database available at:
  %  http://www.nodc.noaa.gov
  % NODC Accession 0011144
  
  cryf1.station_name = 'cryf1';
  cryf1.station_long_name = 'Carysfort Reef';
  cryf1.lat = 25.2100;
  cryf1.lon = -80.2000;
  cryf1.depth = 1;
  
  
  x = importdata(xlsfname,',',43);
  dat = x.data.Data;
  x=[]; clear x;
  
  dts = datenum(dat(:,1),dat(:,2),dat(:,3));
  
  snaf1.daily_sea_t.date = dts;
  snaf1.daily_sea_t.data = dat(:,5);
  snaf1.daily_sea_t = finite_ts(snaf1.daily_sea_t);
  
  hacf1.daily_sea_t.date = dts;
  hacf1.daily_sea_t.data = dat(:,6);
  hacf1.daily_sea_t = finite_ts(hacf1.daily_sea_t);
  
  ccrf1.daily_sea_t.date = dts;
  ccrf1.daily_sea_t.data = dat(:,7);
  ccrf1.daily_sea_t = finite_ts(ccrf1.daily_sea_t);
  
  cmdf1.daily_sea_t.date = dts;
  cmdf1.daily_sea_t.data = dat(:,8);
  cmdf1.daily_sea_t = finite_ts(cmdf1.daily_sea_t);
  
  cryf1.daily_sea_t.date = dts;
  cryf1.daily_sea_t.data = dat(:,9);
  cryf1.daily_sea_t = finite_ts(cryf1.daily_sea_t);
  
  clear ans dat dts usgspath
  
  disp(['Save ',matfname]);
  save(matfname);
end;



if ( doFigs )

if ( ~exist('mlrf1','var') )
  mlrf1 = get_station_from_station_name('mlrf1'); mlrf1=load_all_ndbc_data(mlrf1); mlrf1=verify_variable(mlrf1,'ndbc_sea_t_1_d_avg');
end;

if ( ~exist('smkf1','var') )
  smkf1 = get_station_from_station_name('smkf1'); smkf1=load_all_ndbc_data(smkf1); smkf1=verify_variable(smkf1,'ndbc_sea_t_1_d_avg');
end;

if ( ~exist('lonf1','var') )
  lonf1 = get_station_from_station_name('lonf1'); lonf1=load_all_ndbc_data(lonf1); lonf1=verify_variable(lonf1,'ndbc_sea_t_1_d_avg');
end;

mlrf1 = plot_ngdc_bathymetry(mlrf1,-[0:35],55e3,true,@contour);
lh=plot(ccrf1.lon,ccrf1.lat,'kp',cmdf1.lon,cmdf1.lat,'kp',cryf1.lon,cryf1.lat,'kp',hacf1.lon,hacf1.lat,'kp',lonf1.lon,lonf1.lat,'kp',mlrf1.lon,mlrf1.lat,'kp',smkf1.lon,smkf1.lat,'kp',snaf1.lon,snaf1.lat,'kp');
daspect([1,cosd(25),1]);
set(lh,'MarkerSize',12,'MarkerFaceColor','b');

if 0;
  fmg;
  plot_ts(ccrf1.daily_sea_t,'.',cmdf1.daily_sea_t,'.',cryf1.daily_sea_t,'.',hacf1.daily_sea_t,'.',snaf1.daily_sea_t,'.');
  legend('CCR','CMD','CRY','HAC','SNA','Location','SouthWest');
  if doPrint; print('-dpng','usgs-time-series.png'); end;
end;

if 0;
  fmg;
  plot_ts(ccrf1.daily_sea_t,'.',cmdf1.daily_sea_t,'.',cryf1.daily_sea_t,'.',hacf1.daily_sea_t,'.',snaf1.daily_sea_t,'.',lonf1.ndbc_sea_t_1_d_avg,'.',mlrf1.ndbc_sea_t_1_d_avg,'.',smkf1.ndbc_sea_t_1_d_avg,'.');
  legend('CCR','CMD','CRY','HAC','SNA','LON','MLR','SMK','Location','SouthWest');
  if doPrint; print('-dpng','usgs-time-series-all.png'); end;
end;


[cc,cm,cr,ha,lo,ml,sn] = intersect_tses(ccrf1.daily_sea_t,cmdf1.daily_sea_t,cryf1.daily_sea_t,hacf1.daily_sea_t,lonf1.ndbc_sea_t_1_d_avg,mlrf1.ndbc_sea_t_1_d_avg,snaf1.daily_sea_t);

if 1;
  fh=fmg;
  spt(3,2,1); scatter_fit_ts(ml,cc,[],[],'ML','CC',fh);
  spt(3,2,2); scatter_fit_ts(ml,cm,[],[],'ML','CM',fh);
  spt(3,2,3); scatter_fit_ts(ml,cr,[],[],'ML','CR',fh);
  spt(3,2,4); scatter_fit_ts(ml,ha,[],[],'ML','HA',fh);
  spt(3,2,5); scatter_fit_ts(ml,sn,[],[],'ML','SN',fh);
  if doPrint; print('-dpng','usgs-mlrf1-scatter-fits.png'); end;
end;

if 1;
  fh=fmg;
  spt(3,2,1); scatter_fit_ts(lo,cc,[],[],'LO','CC',fh);
  spt(3,2,2); scatter_fit_ts(lo,cm,[],[],'LO','CM',fh);
  spt(3,2,3); scatter_fit_ts(lo,cr,[],[],'LO','CR',fh);
  spt(3,2,4); scatter_fit_ts(lo,ha,[],[],'LO','HA',fh);
  spt(3,2,5); scatter_fit_ts(lo,sn,[],[],'LO','SN',fh);
  if doPrint; print('-dpng','usgs-lonf1-scatter-fits.png'); end;
end;

if 1;
  fmg;
  boxplot_tses({ml,cc,cm,cr,ha,sn},[],'LEGEND',{'ml','cc','cm','cr','ha','sn'});
  titlename('Monthly: MLRF1 vs. Hudson et al.');
end;
if 1;
  fmg;
  boxplot_tses({lo,cc,cm,cr,ha,sn},[],'LEGEND',{'lo','cc','cm','cr','ha','sn'});
  titlename('Monthly: LONF1 vs. Hudson et al.');
  if doPrint; print('-dpng','usgs-lonf1-monthly-boxplots.png'); end;
end;

% "Trends"
if 1;
  scatter_fit(cmdf1.daily_sea_t.date,cmdf1.daily_sea_t.data); datetick3;
  titlename('CMDF1 trend');
  if doPrint; print('-dpng','usgs-cmdf1-trend.png'); end;
  scatter_fit(hacf1.daily_sea_t.date,hacf1.daily_sea_t.data); datetick3;
  titlename('HACF1 trend');
  if doPrint; print('-dpng','usgs-hacf1-trend.png'); end;
  scatter_fit(smkf1.ndbc_sea_t_1_d_avg.date,smkf1.ndbc_sea_t_1_d_avg.data); datetick3;
  titlename('SMKF1 trend');
  if doPrint; print('-dpng','usgs-smkf1-trend.png'); end;
end;

% But nothing interesting here - would need hourly data?
if 0;
  station_anova_multcompare(ccrf1,'daily_sea_t',[],[],[],[1,0,0,0],{'alpha',0.37});
  if doPrint; print('-dpng','usgs-ccrf1-annual-anova.png'); end;
  station_anova_multcompare(cmdf1,'daily_sea_t',[],[],[],[1,0,0,0],{'alpha',0.37});
  if doPrint; print('-dpng','usgs-cmdf1-annual-anova.png'); end;
  station_anova_multcompare(cryf1,'daily_sea_t',[],[],[],[1,0,0,0],{'alpha',0.37});
  if doPrint; print('-dpng','usgs-cryf1-annual-anova.png'); end;
  station_anova_multcompare(hacf1,'daily_sea_t',[],[],[],[1,0,0,0],{'alpha',0.37});
  if doPrint; print('-dpng','usgs-hacf1-annual-anova.png'); end;
  station_anova_multcompare(snaf1,'daily_sea_t',[],[],[],[1,0,0,0],{'alpha',0.37});
  if doPrint; print('-dpng','usgs-snaf1-annual-anova.png'); end;
end;


%% MOST INTERESTING RESULT? (for me)

% Crocker Crest Reef vs. Hens and Chickens
% Two sites 5.13 km apart, yet some surprisingly weak seasonal correlations!
%scatter_fit_ts_seasons(hacf1.daily_sea_t,ccrf1.daily_sea_t,[],[],'HA','CC');
scatter_fit_ts_seasons(ha,cc,[],[],'HA','CC');
subplots_set(gcf,'XLim',[14,32],'YLim',[14,32]);
if doPrint; print('-dpng','usgs-hacf1-scatter-seasons-ccrf1.png'); end;

% Snake Creek turns out to be a better match for Hens and Chickens than CCR...
%scatter_fit_ts_seasons(hacf1.daily_sea_t,snaf1.daily_sea_t,[],[],'HA','SN');
scatter_fit_ts_seasons(ha,sn,[],[],'HA','SN');
subplots_set(gcf,'XLim',[14,32],'YLim',[14,32]);
if doPrint; print('-dpng','usgs-hacf1-scatter-seasons-snaf1.png'); end;

% Even Cement Dome 70 km away matches Hens and Chickens better than CCR!
%scatter_fit_ts_seasons(hacf1.daily_sea_t,cmdf1.daily_sea_t,[],[],'HA','CM');
scatter_fit_ts_seasons(ha,cm,[],[],'HA','CM');
subplots_set(gcf,'XLim',[14,32],'YLim',[14,32]);
if doPrint; print('-dpng','usgs-hacf1-scatter-seasons-cmdf1.png'); end;

% Is it outflow? If so, CCR and Snake should differ markedly...
scatter_fit_ts_seasons(cc,sn,[],[],'CC','SN');
subplots_set(gcf,'XLim',[14,32],'YLim',[14,32]);
if doPrint; print('-dpng','usgs-ccrf1-scatter-seasons-snaf1.png'); end;



% Crocker bottom sensor at 5 m, Hens and Chickens at 4 m - huh??
bathc=[]; clear bathc;
bathc.lon = mean([ccrf1.lon,hacf1.lon,snaf1.lon]); bathc.lat = mean([ccrf1.lat,hacf1.lat,snaf1.lat]);
bathc = plot_ngdc_bathymetry(bathc,-[0:0.5:8],8e3,true,@contour);
%set(gca,'clim',[-8,0]);
set(gca,'clim',[-7,-4]);
plot(hacf1.lon,hacf1.lat,'kp'); text(hacf1.lon,hacf1.lat,' HACF1');
plot(ccrf1.lon,ccrf1.lat,'kp'); text(ccrf1.lon,ccrf1.lat,' CCRF1');
plot(snaf1.lon,snaf1.lat,'kp'); text(snaf1.lon,snaf1.lat,' SNAF1');
daspect([1,cosd(25),1]);
if doPrint; print('-dpng','usgs-bathy-wide.png'); end;

% Calculate seafloor slope using 7-point centered finite differences
[bathc.ngdc_92m_bathy.dzdx,bathc.ngdc_92m_bathy.dzdy] = gradientn(bathc.ngdc_92m_bathy.field,7,84.1,92.8);
bathc.ngdc_92m_bathy.beta = sqrt( (bathc.ngdc_92m_bathy.dzdx.^2)+(bathc.ngdc_92m_bathy.dzdy.^2) );
if doPrint; print('-dpng','usgs-bathy-closeup.png'); end;

fmg; contourf(bathc.ngdc_92m_bathy.lon,bathc.ngdc_92m_bathy.lat,bathc.ngdc_92m_bathy.beta,[0.0050:0.0025:0.0300]); set(gca,'clim',[0.0050,0.030]); colorbar; titlename('\beta');
plot(hacf1.lon,hacf1.lat,'rp'); text(hacf1.lon,hacf1.lat,' HACF1');
plot(ccrf1.lon,ccrf1.lat,'rp'); text(ccrf1.lon,ccrf1.lat,' CCRF1');
plot(snaf1.lon,snaf1.lat,'rp'); text(snaf1.lon,snaf1.lat,' SNAF1');
daspect([1,cosd(25),1]);
if doPrint; print('-dpng','usgs-bathy-wide-beta.png'); end;

end;
