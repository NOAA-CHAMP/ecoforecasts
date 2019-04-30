1;
% SCRIPT BATHYSTATS.m: Calculate spatial statistics from bathymetry. Call
% CONTOURC to select all points along isobath ISO (POSITIVE for water depth,
% DEFAULT: 20); call INTERP_FIELD to process points in a SAMPLERNG radius
% around each point (DEFAULT: 1e3 m), finding the highest (NANMAX) depth and
% steepest slope in that sampling region. Produces arrays: Ciso of [LON;LAT]
% contour points; Cmaxz of highest depths and Cmaxbeta of greatest slopes
% within SAMPLERNG; and Cmaxbeta_at of slopes *at* each point Ciso.
%
% If DOFIGS (DFLT: True) display figures; if DOPRINT (DFLT: False), print to
% PNG files in FIGSPATH (DFLT: get_relative_path('figs')); if DOSAVE (DFLT:
% False), save workspace to .MAT file before exiting.
%
% Last Saved Time-stamp: <Sun 2017-03-05 22:09:14 Eastern Standard Time gramer>

if ( ~exist('doFigs','var') )
  doFigs = true;
end;
if ( ~exist('doPrint','var') )
  doPrint = false;
end;
if ( ~exist('figspath','var') )
  figspath = get_relative_path('figs');
end;
if ( ~exist('iso','var') )
  iso = 20;
end;
if ( ~exist('samplerng','var') )
  samplerng = 1e3;
end;
if ( ~exist('doSave','var') )
  doSave = false;
end;

bath=[]; clear bath;
bath.lat = +25.00; bath.lon = -80.50;
switch (iso),
 case 20,
  bathyrng = [125e3,120e3];
 case 15,
  bathyrng = [110e3,120e3];
 case 10,
  bathyrng = [90e3,120e3];
 otherwise,
  error('Unknown isobath ISO %g',iso);
end;
bath = read_hires_bathymetry(bath,bathyrng,[],true);

% Calculate beta: seafloor slope at each point in BATH.ngdc_hires_bathy
bath.ngdc_hires_bathy = calc_ngdc_slope(bath.ngdc_hires_bathy);

dy = field_resolution_m(bath.ngdc_hires_bathy);
if ( 20 < dy && dy < 40 )
  bathysrc = 'USGS';
else
  bathysrc = 'NGDC'; %Default
end;

fbasenm = lower([mfilename,'-',bathysrc,'-',num2str(iso),'m-isobath']);


Ciso=[]; Cmaxz=[]; Cmaxzix=[]; clear Ciso Cmaxz Cmaxzix
sampleix = floor(samplerng/dy);
Ciso = contourc(bath.ngdc_hires_bathy.lon,bath.ngdc_hires_bathy.lat,bath.ngdc_hires_bathy.field,[-iso -iso]);
%goodix = find(Ciso(1,:)~=-iso);
Ciso(:,Ciso(1,:)==-iso) = [];

if ( doFigs )
  fmg; plot(Ciso(1,:),Ciso(2,:),'.');
  titlename([bathysrc,' ',num2str(iso),' m contour']);
  if doPrint; print('-dpng',fullfile(figspath,[fbasenm,'-contour.png'])); end;
end;

Cmaxz = interp_field(bath.ngdc_hires_bathy.lat,bath.ngdc_hires_bathy.lon,bath.ngdc_hires_bathy.field,Ciso(2,:),Ciso(1,:),{@nanmax,sampleix});
Cmaxbeta = interp_field(bath.ngdc_hires_bathy.lat,bath.ngdc_hires_bathy.lon,bath.ngdc_hires_bathy.beta,Ciso(2,:),Ciso(1,:),{@nanmax,sampleix});
Cmaxbeta_at = interp_field(bath.ngdc_hires_bathy.lat,bath.ngdc_hires_bathy.lon,bath.ngdc_hires_bathy.beta,Ciso(2,:),Ciso(1,:));


if ( doFigs )
  fmg; hist(Cmaxz(:),1000);
  titlename(['Highest point near ',bathysrc,' ',num2str(iso),' m isobath']);
  if doPrint; print('-dpng',fullfile(figspath,[fbasenm,'-maxz-near-HIRES.png'])); end;

  fmg; hist(Cmaxbeta(:),1000);
  xlim([0,0.20]); % There may be outliers (canyons? bullseyes in the data?)
  titlename(['Peak slope near ',bathysrc,' ',num2str(iso),' m isobath']);
  if doPrint; print('-dpng',fullfile(figspath,[fbasenm,'-maxbeta-near-HIRES.png'])); end;

  fmg; hist(Cmaxbeta_at(:),1000);
  xlim([0,0.20]);
  titlename(['Slope precisely at ',bathysrc,' ',num2str(iso),' m isobath']);
  if doPrint; print('-dpng',fullfile(figspath,[fbasenm,'-maxbeta-at-HIRES.png'])); end;
end;

if ( doSave )
  disp(['Saving ',fbasenm,'.mat']);
  save([fbasenm,'.mat'],'-v7.3');
end;
