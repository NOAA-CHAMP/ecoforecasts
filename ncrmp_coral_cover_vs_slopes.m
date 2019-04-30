1;
%% SCRIPT ncrmp_coral_cover_vs_slopes.m:
% 
% Compare NCRMP coral cover estimates for each survey site in DATASET
% (RGN_YEAR) with seafloor slope for those sites from NGDC bathymetry.
% Regions may be FGBNMS, FLK, PRICO, SEFCRI, TortugasMarq, or USVI. Years are
% available survey years in *.RDA files (see NCRMP_CORAL_COVER).
% 
% CALLS: ncrmp_coral_cover
% 
% Last Saved Time-stamp: <Sun 2019-04-14 15:49:43 Eastern Daylight Time gramer>


if ( ~exist('doPrint','var') )
  doPrint = false;
end;
if doPrint; disp('Printing figures to ./figs. Hit enter to confirm...'); pause; end;

if ( ~exist('ncrmp','var') )
  ncrmp_coral_cover;
end;

%dataset = 'FGBNMS_2013'; %N=0 (all surveys are below 20 m)
%dataset = 'FGBNMS_2015'; %N=0 (all surveys are below 20 m)
dataset = 'FLK_2014'; % SUCCESS!
%dataset = 'FLK_2016'; % Need to fix bathymetry coverage: p>0.50
%dataset = 'PRICO_2014'; %p>0.98!
%dataset = 'PRICO_2016'; % p>0.50
%dataset = 'SEFCRI_2014'; %p>0.50, B<0!
%dataset = 'SEFCRI_2016'; %p>0.35, B<0! (but Bath h ~ dive h)
%dataset = 'TortugasMarq_2014'; % Too many flat sites: p>0.50
%dataset = 'TortugasMarq_2016'; % Too many flat sites: p>0.15
%dataset = 'USVI_2013'; % p>0.90!
%dataset = 'USVI_2015'; %p>0.70
%dataset = 'USVI_2017'; %p>0.60

disp(['Finding Slopes for ',ncrmp.(dataset).basename]);

if ( exist('bath','var') )
  bbox = field_bbox(bath.ngdc_hires_bathy);
  nhits = numel(find(inside(ncrmp.(dataset).ulon,ncrmp.(dataset).ulat,bbox([1,1,2,2,1]),bbox([3,4,4,3,3]))));
  hitpct = 100*nhits/numel(ncrmp.(dataset).ulon);
  if ( hitpct < 66 )
    bath=[]; clear bath rad lh
  elseif ( hitpct < 97 )
    warning('Reusing bathymetry with only %g%% of points',hitpct);
  end;
  clear bbox nhits hitpct
end;
if ( ~exist('bath','var') )
  disp('Extracting bathymetry');
  [bath,rad] = read_hires_bathymetry_for_field({ncrmp.(dataset).ulon,ncrmp.(dataset).ulat},true,[10e3,20e3]);
  % BAD: For regions outside Florida, this gets ~1 km resolution bathymetry!
  %%[bath,rad] = read_hires_bathymetry_for_field({ncrmp.(dataset).ulon,ncrmp.(dataset).ulat},false,[10e3,20e3]);
end;

sites=[]; clear sites
if ( min(bath.ngdc_hires_bathy.yres) < 15 )
  % 10 m bathymetry: MEDIAN 9x9-point (90 m) squares with >=10 valid points
  [sites,bath] = find_ngdc_slope_sites({ncrmp.(dataset).ulon,ncrmp.(dataset).ulat},bath,9,{@nanmedian,9,9,10});
elseif ( min(bath.ngdc_hires_bathy.yres) < 50 )
  % 30 m bathymetry: MEDIAN 3x3-point (90 m) squares with >=4 valid points
  [sites,bath] = find_ngdc_slope_sites({ncrmp.(dataset).ulon,ncrmp.(dataset).ulat},bath,3,{@nanmedian,3,3,4});
else
  % 92 m bathymetry: Simple finite difference and linear interpolation
  [sites,bath] = find_ngdc_slope_sites({ncrmp.(dataset).ulon,ncrmp.(dataset).ulat},bath,2);
end;

[sites.coastdx,sites.coastaz,sites.nearestCoastCoords] = ...
    get_contour_distance({bath.ngdc_hires_bathy,0.1},[ncrmp.(dataset).ulon';ncrmp.(dataset).ulat']);

[B,Stats,fh,sh] = scatter_fit(sites.depths,-ncrmp.(dataset).ueh,...
                              ['Bath ',num2str(nanmean(bath.ngdc_hires_bathy.yres),'%3.1f'),' m'],...
                              [textize(ncrmp.(dataset).basename),' Mean Dive h']);
set(gca,'FontSize',24); set(sh,'LineWidth',3,'MarkerSize',32);
if doPrint; print('-dpng',['figs/ncrmp_survey_vs_depths_',ncrmp.(dataset).basename,'.png']); end;

% % Cut off at slopes and DEPTHS where horizontal convection and wave breaking are unlikely
% goodix = find(sites.betas>=0.0075 & 1<=ncrmp.(dataset).ueh & ncrmp.(dataset).ueh<=20.0);
% LOOSE-ish Cut off at slopes and DEPTHS where horizontal convection and wave breaking are unlikely
goodix = find(sites.betas>=0.0075 & 1<=ncrmp.(dataset).ueh & ncrmp.(dataset).ueh<=30.0);
% % LOOSE Cut off at slopes and DEPTHS where horizontal convection and wave breaking are unlikely
% goodix = find(sites.betas>=0.0025 & 1<=ncrmp.(dataset).ueh & ncrmp.(dataset).ueh<=30.0);
disp(['Good=',num2str(100*numel(goodix)/numel(sites.betas)),'%']);

% Distance has small predictive power for seafloor slope...
scatter_fit(sites.coastdx(goodix)',sites.betas(goodix)','Dx to shore','Seafloor \beta');
% But distance has no predictive power for coral cover! :) :)
scatter_fit(sites.coastdx(goodix)',ncrmp.(dataset).utc(goodix)','Dx to shore','Total Coral Cvr%');

expon = 2/3;
bets = [0.0075,0.01:0.01:0.05,0.10:0.10:0.30];


[B,Stats,fh,sh] = scatter_fit(sites.betas(goodix)'.^expon,ncrmp.(dataset).utc(goodix)',...
                              ['Seafloor Slope (\beta)^{',num2str(expon),'}'],...
                              ['Total ',textize(ncrmp.(dataset).basename),' Cvr%']);
ylim([0,100]);
set(gca,'FontSize',24); set(sh,'LineWidth',3,'MarkerSize',32);
for bet=bets(:)'
  [lh,th]=annotline(bet.^expon,[],['\beta=',num2str(bet)],'k');
  set(lh,'LineWidth',3); set(th,'FontSize',24);
end;
legend(sh,'Location','NorthEast');
if doPrint; print('-dpng',['figs/ncrmp_coral_cover_vs_slopes_',ncrmp.(dataset).basename,'.png']); end;


[B,Stats,fh,sh] = scatter_fit(sites.betas(goodix)'.^expon,ncrmp.(dataset).mtc(goodix)',...
                              ['Seafloor Slope (\beta)^{',num2str(expon),'}'],...
                              [textize(ncrmp.(dataset).basename),' Mas.Cvr%']);
ylim([0,100]);
set(gca,'FontSize',24); set(sh,'LineWidth',3,'MarkerSize',32);
for bet=bets(:)'
  [lh,th]=annotline(bet.^expon,[],['\beta=',num2str(bet)],'k');
  set(lh,'LineWidth',3); set(th,'FontSize',24);
end;
legend(sh,'Location','NorthEast');
if doPrint; print('-dpng',['figs/ncrmp_coral_cover_vs_slopes_',ncrmp.(dataset).basename,'_massive.png']); end;

clear bet lh sh th

%{
%}

disp(['Finding Waves for ',ncrmp.(dataset).basename]);

basename = 'ncrmp_coral_cover_vs_slopes';
nwpsmatfname = fullfile(get_ecoforecasts_path('data'),[basename,'_nwps.mat']);
if ( exist(nwpsmatfname,'file') )
  disp(['Loading ',nwpsmatfname]);
  load(nwpsmatfname);
else
  stns = get_nwps_mat_stations([ncrmp.(dataset).ulon,ncrmp.(dataset).ulat,ncrmp.(dataset).ustn]');
  disp(['Saving ',nwpsmatfname]);
  save(nwpsmatfname,'stns');
end;

ncrmp.(dataset).nwps_dates = stns(1).nwps_sigwavehgt.date;
flds={'nwps_sigwavehgt','nwps_primwavedir','nwps_primwaveper','nwps_swellhgt'};
for fldix=1:numel(flds)
  fld = flds{fldix};
  ncrmp.(dataset).(fld) = repmat(nan,[numel(ncrmp.(dataset).nwps_dates),numel(stns)]);
  for stix=1:numel(stns)
    ncrmp.(dataset).(fld)(:,stix) = stns(stix).(fld).data;
  end;
  ncrmp.(dataset).([fld,'_mean']) = nanmean(ncrmp.(dataset).(fld),1);
  ncrmp.(dataset).([fld,'_std']) = nanstd(ncrmp.(dataset).(fld),0,1);
end;

% Random spot-checks of results
ixen = round(rand(1,4).*numel(stns)); jxen = round(rand(1,4).*numel(stns));
fh=fmg; for ix=1:4; for jx=1:4; spt(4,4,(ix-1)*4+jx); scatter_fit_ts(stns(ixen(ix)).nwps_sigwavehgt,stns(jxen(jx)).nwps_sigwavehgt,[],[],[],[],fh,[],true,[],[],{'shortLegend','Location','North'}); axis([0,2.5,0,2.5]); end; end; subplots_set(fh,'FontSize',12);

ixen = round(rand(4,4).*numel(stns));
fh=fmg; for ix=1:16; spt(4,4,ix); rose(stns(ixen(ix)).nwps_primwavedir.data,100); end; subplots_set(fh,'FontSize',12);
fh=fmg; for ix=1:16; spt(4,4,ix); rosewgtd(stns(ixen(ix)).nwps_primwavedir,stns(ixen(ix)).nwps_sigwavehgt); end; subplots_set(fh,'FontSize',12);


%% Calculate WAVE ENERGY statistics for each site: https://cdip.ucsd.edu/?nav=documents&sub=index&xitem=waves
% . . .
