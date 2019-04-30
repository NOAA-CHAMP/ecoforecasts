1;
%% SCRIPT calc_spatial_dt_hc
%
% Calculate spatial pattern of cooling/warming based on Horizontal Convection
% (HC_SS, or ...) over depths H, seafloor slopes BET and HC ranges HCRNG. (If
% H,... are not found, load them from appropriately named local .MAT files.)
%
% Last Saved Time-stamp: <Thu 2017-05-18 21:58:51 Eastern Daylight Time gramer>

set_more off;

if ( ~exist('subrgn','var') )
  subrgn = '';
  %subrgn = 'SE';
  %subrgn = 'UK';
  %subrgn = 'MK';
  %subrgn = 'LK';
  %subrgn = 'DT';
end;
if ( ~exist('use_habitat_map','var') || isempty(use_habitat_map) )
  use_habitat_map = false;
  %use_habitat_map = true;
end;
if ( ~exist('allow_hard_bottom','var') || isempty(allow_hard_bottom) )
  allow_hard_bottom = true;
  %allow_hard_bottom = false;
end;

if ( ~exist('mindepth','var') || isempty(mindepth) )
  mindepth = 1;
end;
if ( ~exist('maxdepth','var') || isempty(maxdepth) )
  maxdepth = 30;
end;

dt = 3600;

R = (1-0.08);

%% Check seasonality of sea-surface radiative and turbulent stress forcing

%load('d:/thesis/data/mlrf1-heat_budget-erai-avhrr_weekly-ndbc-tpxo_tide-erai_DISSERT.mat');
%OR: load('d:/thesis/data/lonf1-heat_budget-erai-avhrr_weekly-ndbc-tpxo_tide-erai.mat');
%fmg; boxplot_ts(stn.simple_ndbc_erai_erai_30a_net_flux,[],'mean',true,'std',true);
%stn = verify_variable(stn,'simple_ndbc_erai_erai_30a_net_flux_1_d_sum');
%MLRF1: unique(get_month(stn.simple_ndbc_erai_erai_30a_net_flux_1_d_sum.date(stn.simple_ndbc_erai_erai_30a_net_flux_1_d_sum.data<-23000)))
%LONF1: unique(get_month(stn.simple_ndbc_erai_erai_30a_net_flux_1_d_sum.date(stn.simple_ndbc_erai_erai_30a_net_flux_1_d_sum.data<-18500)))
%MLRF1: unique(get_year(stn.simple_ndbc_erai_erai_30a_net_flux_1_d_sum.date(get_month(stn.simple_ndbc_erai_erai_30a_net_flux_1_d_sum.date)==6&stn.simple_ndbc_erai_erai_30a_net_flux_1_d_sum.data>6200))) % 1998
%LONF1: unique(get_year(stn.simple_ndbc_erai_erai_30a_net_flux_1_d_sum.date(get_month(stn.simple_ndbc_erai_erai_30a_net_flux_1_d_sum.date)>=6&stn.simple_ndbc_erai_erai_30a_net_flux_1_d_sum.data>7500))) % 1995

% % December/July "representative" means
% q0s = [-100,100];
% q0str={'Representative mean winter','Representative mean summer'};
% % MLRF: December/July means at MLRF1
% q0s = [-115,95];
% q0str={'Reef-crest mean winter','Reef-crest mean summer'};
% % LONF1: February 25th-percentile / May 75th-percentile
% q0s = [-140,285];
% q0str={'Back-reef flats cool winter','Back-reef flats warm summer'};
% % MLRF1: December 25th-percentile / May 75th-percentile
% q0s = [-240,320];
% q0str={'Reef-crest cool winter','Reef-crest warm summer'};

% %BASED ON 1_d_sum
% % MLRF1: December lower std. dev. and mean, June mean, June maximum
% q0s = [-240,-115,95,255];
% q0str={'Cold winter','Normal winter','Normal summer','Hot summer'};

% %BASED ON 1_d_avg
% % MLRF1: January minimum and mean, June mean, June maximum
% q0s = [-550,-115,95,255];
% q0str={'2010 cold snap','Normal winter','Normal summer','1998 bleaching'};

%{
%%%%%%%%%% HACK - DEBUGGING with BNPON and BNPMI as examples
q0s = [  -800,             -115];
Tas = [  7.5,              14.0];
q0str = {'2010 cold snap','Normal winter'};
%}


% MLRF1: January minimum and mean, June mean, June maximum
q0s = [  -800,             -115,           95,              255];
% 7th or 93rd percentile Air Temperature for the month of each event
Tas = [  7.5,              14.0,           30.0,            30.5];
q0str = {'2010 cold snap', 'Normal winter', 'Normal summer','1998 bleaching'};

nq0s = numel(q0s);


switch (subrgn),
 case '',
  error('No sub-region code (SUBRGN) was specified!');
 case 'SE',
  subrgnstr = 'SE Florida';
  subrgnbox = [-80.20,-79.95,25.45,27.00];%27.30];
 case 'UK',
  subrgnstr = 'Upper Keys';
  %subrgnbox = [-80.45,-80.05,24.95,25.55];
  subrgnbox = [-80.60,-80.05,24.90,25.55];
 case 'MK',
  subrgnstr = 'Middle Keys';
  %subrgnbox = [-81.05,-80.35,24.55,25.05];
  subrgnbox = [-81.05,-80.30,24.55,25.15];
 case 'LK',
  subrgnstr = 'Lower Keys';
  subrgnbox = [-82.05,-80.95,24.40,24.75];
 case 'DT',
  subrgnstr = 'Dry Tortugas';
  subrgnbox = [-82.95,-81.95,24.35,24.75];

 case 'FK',
  subrgnstr = 'Florida Keys';
  subrgnbox = [-82.05,-80.05,24.35,25.60];
 case 'FRT',
  subrgnstr = 'Florida Reef Tract';
  subrgnbox = [-82.95,-79.95,24.35,27.00];
 otherwise,
  error('Region %s not yet implemented!',subrgn);
end;


% Assume Horizontal Convection affects sea temperatures within about
% 270 m radius: Downsample bathymetry to the appropriate resolution
if ( ~exist('hc_radius','var') || isempty(hc_radius) )
  hc_radius = 270;
  %hc_radius = 92;
end;

% % 10-m resolution bathymetry
% bath_res = 10;
% % 30-m resolution bathymetry
% bath_res = 30;
% 92-m resolution bathymetry
bath_res = 92;

downstride = ceil(hc_radius/bath_res);
if ( downstride > 1 )
  downfunc = {@nanmean,downstride,downstride,downstride+1}; 
  downfuncstr = ['_',interpmethod_to_text(downfunc)];
else
  disp('No downsampling needed');
  downfuncstr = '';
end;

matbasename = sprintf('FRT_depth_and_beta_%dm',bath_res);
if ( use_habitat_map )
  if ( allow_hard_bottom )
    habfname = '_coral_and_hard_bottom';
  else
    habfname = '_coral';
  end;
else
  habfname = '';
end;
matbasefname = [matbasename,'_hc_range_depth',habfname,'.mat'];
matfname = [matbasename,downfuncstr,'_hc_range_depth',habfname,'.mat'];

if ( ~exist('lon','var') || ~exist('lat','var')  || ...
     ~exist('h','var') || ~exist('bet','var') || ~exist('hcrng','var') )

  if ( exist(matfname,'file') )
    disp(['Loading ',matfname]);
    load(matfname);

    % Somehow, the unaveraged bathymetry got processed differently before saving??
    if ( downstride == 1 )
      h = h';
      bet = bet';
      % No HC over land, and depths are positive
      h(h>0) = 0;
      h = abs(h);
    end;

  else
    %DEBUG:
    disp(['Not found: ',matfname]);    keyboard;
    %DEBUG:    if ( use_habitat_map );      keyboard;    end;

    lon=[]; lat=[]; h=[]; bet=[]; clear lon lat h bet
    if ( exist(matbasefname,'file') )
      disp(['Loading ',matbasefname]);
      load(matbasefname);

    elseif ( use_habitat_map )
      error(['USE_HABITAT_MAP: Unable to find ',matbasefname]);

    else
      disp('Extracting bathymetry');
      fld.lon = [-82.95,-79.95];
      fld.lat = [+24.30,+27.30];
      [x,rad] = read_hires_bathymetry_for_field(fld,false);
      bath = x.ngdc_hires_bathy; x=[]; clear x
      [ig,ig,ig,bath] = find_ngdc_slope(bath,[],[],3);
      lon = bath.lon;
      lat = bath.lat;
      h = bath.field;
      bet = bath.beta;
      disp(['Saving ',matbasefname]);
      save(matbasefname,'lon','lat','h','bet');
    end; %if ( exist(matbasefname,'file') ) else

    % No HC over land, and depths are positive
    h(h>0) = 0;
    h = abs(h);

    % Downsample to appropriate resolution for HC effects
    if ( downstride > 1 )
      disp(['Downsampling H and BETA to ~',num2str(hc_radius),' m']);
      [LAT,LON] = meshgrid(lat(1:downstride:end),lon(1:downstride:end));
      H = interp_field(lat,lon,h,LAT,LON,downfunc);
      H = reshape(H,size(LAT));
      BET = interp_field(lat,lon,bet,LAT,LON,downfunc);
      BET = reshape(BET,size(LAT));

      lon=[]; lat=[]; h=[]; bet=[]; clear lon lat h bet
      lon = unique(LON(:));
      lat = unique(LAT(:));
      h = H;
      bet = BET;
      LON=[]; LAT=[]; H=[]; BET=[]; clear LON LAT H BET

      disp(['Saving ',matfname]);
      save(matfname,'lon','lat','h','bet');
    end; %if ( downstride > 1 )

  end; %if ( exist(matfname,'file') ) else

end; %if ( ~exist('lon','var') || ...


% No HC over land, and depths are positive
hch(hch>0) = 0;
hch = abs(hch);
% HC end-point depth should never be deeper than depth ("HC height" >= "height")
hch(hch<h) = h(hch<h);


lon_rmix = find(lon<subrgnbox(1) | subrgnbox(2)<lon);
lat_rmix = find(lat<subrgnbox(3) | subrgnbox(4)<lat);
if ( ~isempty(lon_rmix) || ~isempty(lat_rmix) )
  disp(['Subsetting to sub-region ',subrgn]);
  lat(lat_rmix) = [];		lon(lon_rmix) = [];
  h(:,lat_rmix) = [];		h(lon_rmix,:) = [];
  bet(:,lat_rmix) = [];		bet(lon_rmix,:) = [];
  ang(:,lat_rmix) = [];		ang(lon_rmix,:) = [];
  hcrng(:,lat_rmix) = [];	hcrng(lon_rmix,:) = [];
  hch(:,lat_rmix) = [];		hch(lon_rmix,:) = [];
end;

if ( isempty(h) || isempty(bet) )
  error('Subsetting resulted in empty set! Try reloading...');
end;

%{
%%%%%%%%%% HACK - DEBUGGING with instrumented sites (e.g., BNPON and BNPMI) as examples
%%%%%%%%%% HACK
[LON,LAT] = meshgrid(lon,lat); 
cstnms={'BNPON','BNPMI','LONF1'}; for stix=1:numel(cstnms); stn=get_station_from_station_name(cstnms{stix}); [ig,ix]=min(abs(LON(:)-stn.lon)+abs(LAT(:)-stn.lat)); [jix(1:2,stix),iix(1:2,stix)]=ind2sub(size(LON),ix); end;
LON=[]; LAT=[]; clear LON LAT cstnms stix stn ig ix

h     = [5.45,   6.83,   2.01 ;   5.45,   6.83,   2.01];
bet   = [0.0123, 0.0070, 0.0006 ; 0.0123, 0.0070, 0.0006];
hcrng = [736,    276,    184 ;    540,    540,    540];
hch   = [6.3,    5.6,    1.8 ;    6.3,    5.6,    1.8];

% h,bet,hcrng,hch,
% 26+(squeeze((dTdt_SS(1,:,:)*14)+(dTdt_SS(2,:,:)*16))),
%%%%%%%%%% HACK
%}



if ( ~exist('rhoCp','var') )
  %DEBUG:  disp('Calculating rhoCp');
  s = repmat(35,size(h));
  t = repmat(25,size(h));
  %h = 10.96;

  g = 9.79;						%[m/s^2]
  alph = sw_alpha(s,t,h);				%[1/K]
  rho = sw_dens(s,t,h);					%[kg/m^3]
  Cp = sw_cp(s,t,h);					%[J/kg*K]
  rhoCp = rho.*Cp;					%[J/K*m^3]
  s=[]; t=[]; rho=[]; Cp=[]; clear s t Cp rho
end;

%clear dTdthc_SS dTdthc_US dTdthc_SU dTdthc_UU

disp('Calculating Horizontal Convection and dT');

[betSz1,betSz2] = size(bet);
bix1=1:betSz1;
bix2=1:betSz2;

for qix=1:nq0s;
  % Net sea-surface heat flux, Q_0
  q0 = q0s(qix);
  %DEBUG:  disp(['Q0 = ',num2str(q0)]);

  % For warming, try to account for effects of bottom reflectivity ("gamma")
  Q0 = repmat(q0,size(h));
  if ( q0 > 0 )
    Q0 = Q0 .* (tanh(2.*h./3));
  end;


  % Net sea-surface buoyancy flux B_0
  B0 = (g.*alph.*abs(Q0))./(rhoCp);
  % "Characteristic" convective velocity (Sturmann)
  uf = (B0.*h).^(1/3);
  B0=[]; clear B0

  % Volumetric flow [m^2/s]: Panel letter refers to Fig. 10, Monismith et al (2006)
  % Advective inertial (steady) momentum balance, steady thermal balance: Panel (c)
  Qv_SS = uf.*h./(bet.^(1/3));
  % Viscous (unsteady) inertial balance, balanced thermal forcing: Panel (a)
  Qv_US = sqrt((uf.^3) .* (24.*dt) .* h);
  % Advective inertial balance, unbalanced thermal forcing: Panel (f)
  Qv_SU = (bet.^(2/3)).*uf.*((uf.*(24.*dt)./h).^(3/2));
  % Viscous (unsteady) inertia, unbalanced thermal forcing: Panel (d)
  Qv_UU = bet.*(uf.^3).*((24.*dt)^2)./h;
  uf=[]; clear uf

  dTdtq0 = dt.*Q0./(rhoCp.*h);

  % Convective flow rate [m/s]
  u_SS = (5.0.*Qv_SS./h) - 0.05;	u_SS(u_SS<0) = 0;
  Qv_SS=[]; clear Qv_SS

  u_US = (3.0.*Qv_US./h) - 0.026;	u_US(u_US<0) = 0;
  Qv_US=[]; clear Qv_US

  u_SU = (2.7.*Qv_SU./h) - 0.0224;	u_SU(u_SU<0) = 0;
  Qv_SU=[]; clear Qv_SU

  u_UU = (0.1.*Qv_UU./h);		u_UU(u_UU<0) = 0;
  Qv_UU=[]; clear Qv_UU


  % Furthest distance traveled by convective flow in an hour [m]
  dx_SS = u_SS .* dt;
  dx_US = u_US .* dt;
  dx_SU = u_SU .* dt;
  dx_UU = u_UU .* dt;
  % REALITY CHECK: Convective flow is limited in spatial extent! It is
  % limited both by local topography (first MIN) and physics (second MIN).
  dx_SS = min(dx_SS,hcrng,'includenan');  dx_SS = min(dx_SS,hc_radius*2,'includenan');
  dx_US = min(dx_US,hcrng,'includenan');  dx_US = min(dx_US,hc_radius*2,'includenan');
  dx_SU = min(dx_SU,hcrng,'includenan');  dx_SU = min(dx_SU,hc_radius*2,'includenan');
  dx_UU = min(dx_UU,hcrng,'includenan');  dx_UU = min(dx_UU,hc_radius*2,'includenan');

  hx_SS = h+(bet.*dx_SS);  %hx_SS = max(hx_SS,hch,'includenan');
  hx_US = h+(bet.*dx_US);  %hx_US = max(hx_US,hch,'includenan');
  hx_SU = h+(bet.*dx_SU);  %hx_SU = max(hx_SU,hch,'includenan');
  hx_UU = h+(bet.*dx_UU);  %hx_UU = max(hx_UU,hch,'includenan');

  %%%% HACK
  %%%% DEBUGGING: 
%{
  dx_SS = hc_radius*2;
  dx_US = hc_radius*2;
  dx_SU = hc_radius*2;
  dx_UU = hc_radius*2;
  dx_SS = hc_radius;
  dx_US = hc_radius;
  dx_SU = hc_radius;
  dx_UU = hc_radius;
%}
  dx_SS = hcrng;
  dx_US = hcrng;
  dx_SU = hcrng;
  dx_UU = hcrng;

  %%%% HACK
  %%%% DEBUGGING: 
%{
%}
  hx_SS = hch;
  hx_US = hch;
  hx_SU = hch;
  hx_UU = hch;


  % Temperature change due to Q0 at furthest extent of convection
  dTdtx_SS = dt.*q0./(rhoCp.*hx_SS);
  dTdtx_US = dt.*q0./(rhoCp.*hx_US);
  dTdtx_SU = dt.*q0./(rhoCp.*hx_SU);
  dTdtx_UU = dt.*q0./(rhoCp.*hx_UU);


  % Static temperature gradient due to Q0 over depth difference between
  % observation point, and point of further extent of convection
  dTdx_SS=(dTdtq0-dTdtx_SS)./dx_SS;	dTdx_SS(~isfinite(dTdx_SS)) = 0;
  dx_SS=[]; dTdtx_SS=[]; clear dx_SS dTdtx_SS

  dTdx_US=(dTdtq0-dTdtx_US)./dx_US;	dTdx_US(~isfinite(dTdx_US)) = 0;
  dx_US=[]; dTdtx_US=[]; clear dx_US dTdtx_US

  dTdx_SU=(dTdtq0-dTdtx_SU)./dx_SU;	dTdx_SU(~isfinite(dTdx_SU)) = 0;
  dx_SU=[]; dTdtx_SU=[]; clear dx_SU dTdtx_SU

  dTdx_UU=(dTdtq0-dTdtx_UU)./dx_UU;	dTdx_UU(~isfinite(dTdx_UU)) = 0;
  dx_UU=[]; dTdtx_UU=[]; clear dx_UU dTdtx_UU


  % Rayleigh Benard instability may dampen HC during warming (Mao Lei Patterson)
  RB = repmat(1.00,size(dTdtq0));
  RB(dTdtq0 > 0) = 0.66;

  % Temperature change due to horizontal convection at observation point
  dTdthc_SS(qix,bix1,bix2) = -RB.*R.*24.*dt.*u_SS.*dTdx_SS;
  u_SS=[]; dTdx_SS=[]; clear u_SS dTdx_SS

  dTdthc_US(qix,bix1,bix2) = -RB.*R.*24.*dt.*u_US.*dTdx_US;
  u_US=[]; dTdx_US=[]; clear u_US dTdx_US

  dTdthc_SU(qix,bix1,bix2) = -RB.*R.*24.*dt.*u_SU.*dTdx_SU;
  u_SU=[]; dTdx_SU=[]; clear u_SU dTdx_SU

  dTdthc_UU(qix,bix1,bix2) = -RB.*R.*24.*dt.*u_UU.*dTdx_UU;
  u_UU=[]; dTdx_UU=[]; clear u_UU dTdx_UU
  RB=[]; clear RB

  % Net temperature change (HC+Q0) at observation point
  dTdt_SS(qix,bix1,bix2) = squeeze(dTdthc_SS(qix,:,:)) + (dTdtq0.*24);
  dTdt_US(qix,bix1,bix2) = squeeze(dTdthc_US(qix,:,:)) + (dTdtq0.*24);
  dTdt_SU(qix,bix1,bix2) = squeeze(dTdthc_SU(qix,:,:)) + (dTdtq0.*24);
  dTdt_UU(qix,bix1,bix2) = squeeze(dTdthc_UU(qix,:,:)) + (dTdtq0.*24);

  % REALITY CHECK: this convective flow model cannot reverse temperatures!
  if ( q0 > 0 )
    dTdt_SS(qix,:,:) = max(0,dTdt_SS(qix,:,:),'includenan');
    dTdt_US(qix,:,:) = max(0,dTdt_US(qix,:,:),'includenan');
    dTdt_SU(qix,:,:) = max(0,dTdt_SU(qix,:,:),'includenan');
    dTdt_UU(qix,:,:) = max(0,dTdt_UU(qix,:,:),'includenan');
  else
    dTdt_SS(qix,:,:) = min(0,dTdt_SS(qix,:,:),'includenan');
    dTdt_US(qix,:,:) = min(0,dTdt_US(qix,:,:),'includenan');
    dTdt_SU(qix,:,:) = min(0,dTdt_SU(qix,:,:),'includenan');
    dTdt_UU(qix,:,:) = min(0,dTdt_UU(qix,:,:),'includenan');
  end;


  dTdtq0_simple(qix,bix1,bix2) = dTdtq0;
  dTdtq0=[]; clear dTdtq0
  Q0=[]; clear Q0

end; %for qix=1:nq0s;


alph=[]; rhoCp=[]; clear alph rhoCp

set_more;
