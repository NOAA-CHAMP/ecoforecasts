1;
%% SCRIPT plot_spatial_dt_hc_thermal_stress
%
% Plot spatial pattern of cooling/warming based on steady-state Horizontal
% Convection (HC_SS) over depths H and seafloor slopes BET. If dTdt_* do not
% exist, calculate them by calling script CALC_SPATIAL_DT_HC (v.)
%
% This version uses a more "realistic" pattern of sea-surface heating ("Q0")
% than CALC_SPATIAL_DT_HC (v.) And it estimates the "coral thermal stress",
% degree-heating days above an estimated maximum monthly per-pixel mean.
%
% (SEE CALC_SPATIAL_DT_HC_THERMAL_STRESS.m)
%
% SAMPLE RUN for only reef and hard-bottom habitats in Middle Keys sub-region:
% >> subrgn='MK'; use_habitat_map=true; allow_hard_bottom=true; plot_sites=true; doPrint=false; plot_spatial_dt_hc_thermal_stress
%
% TO RERUN AND PRINT ALL four SCENARIOS for ALL SUB-REGIONS:
% >> doPrint=true; for x={'DT','LK','MK','UK','SE'}; fclose('all'); ...
%  close all; clear a* b* c* d* e* f* g* h* l* m* n* o* p* q* r* s* t* u* v* w*; ...
%  subrgn=x{:}; disp(subrgn); plot_spatial_dt_hc_thermal_stress; end;
%
% Last Saved Time-stamp: <Thu 2018-12-06 21:29:11 Eastern Standard Time gramer>


set_more off;

if ( ~exist('dTdt_SS') )
  calc_spatial_dt_hc_thermal_stress;
end;
if ( ~exist('TW') )
  if ( ~exist('doReport','var') || isempty(doReport) )
    doReport = false;
  end;
  compare_spatial_dt_hc_thermal_stress;
end;

if ( ~exist('plot_coastline') || isempty(plot_coastline) )
  plot_coastline = true;
end;

if ( ~exist('plot_sites') || isempty(plot_sites) )
  %plot_sites = false;
  plot_sites = true;
end;
if ( ~exist('stnclr') || isempty(stnclr) )
  % % Greenish-gray for all scenarios
  % stnclr = [.5,.8,.5];

  % Reddish-gray for cool scenarios, greenish-gray for warm ones
  stnclr = {[.8,.5,.5],[.8,.5,.5],[.5,.8,.5],[.5,.8,.5],};
end;
if ( ~iscell(stnclr) )
  stnclr = {stnclr,stnclr,stnclr,stnclr};
end;

if ( ~exist('scenarios','var') || isempty(scenarios) )
  scenarios = 1:nq0s;
end;


if ( ~exist('doLog','var') || isempty(doLog) )
  %doLog = true;
  doLog = false;
end;
if ( ~exist('doPrint','var') )
  %doPrint = true;
  doPrint = false;
end;
if ( ~exist('figspath','var') )
  figspath = get_ecoforecasts_path('figs');
end;

figbasename = mfilename;
if ( use_habitat_map )
  figbasename = [figbasename,'_habitat'];
  if ( allow_hard_bottom )
    figbasename = [figbasename,'_rfhd'];
  end;
end;
if ( plot_sites )
  figbasename = [figbasename,'_sites'];
end;
if ( doPrint )
  disp(['Will print plots to ',fullfile(figspath,[figbasename,'_*'])]);
end;


% Load Coral and Hard-Bottom habitat map polygons
if ( ~exist('rfhd','var') )
  read_coral_and_hard_bottom;
end;


if ( any(plot_coastline) && ~exist('bath','var') )
  x = read_hires_bathymetry_for_field({lon,lat},false);
  bath = x.ngdc_hires_bathy; x=[]; clear x
end;

if (doLog) for qix=scenarios(:)';

  % Net sea-surface heat flux, Q_0
  q0 = q0s(qix);
  disp(['Log warming: Q0 = ',num2str(q0)]);
  dT = squeeze(dTdt_SS(qix,:,:));

  fh=fmg;
  dTlog = log10(abs(dT));
  dTlog(maxdepth<h | h<mindepth) = nan;
  % %contourf(lon,lat,dTlog',[-1.00:0.05:0.50]);  caxis([-1,0]);
  % contourf(lon,lat,dTlog',[-1.30:0.05:0.50]);  caxis([-1.3,0.2]);
  surf(lon,lat,zeroZ',dTlog'); caxis(log10([0.05,2.0])); shading interp
  axis(subrgnbox);
  set_pcolor_cursor(fh,@xyz_logc_id_select_cb);

  if ( q0 > 0 )
    phenomstr = 'warming';
    % Reversed HOT shows warmest spots as darkest red
    cm=hot; colormap(cm(end:-1:1,:))
  else
    phenomstr = 'cooling';
    cm=cool; colormap(cm);
  end;

  cbh=colorbar;
  % Make a logarithmic color bar
  rescale_colorbar(cbh);
  set(get(cbh,'Label'),'String','log_1_0(^oC)');

  daspect([1,cosd(lat(1)),1]); 

  % Adjust graph appearance based on shape of sub-region
  switch (subrgn),
   case 'SE',
    view(90,90);
    set(cbh,'Orientation','horizontal','Position',[0.30 0.30 0.50 0.02]);
   case 'UK',
    % view(90,90);
    % set(cbh,'Location','East');
   case 'MK',
    view(25,90);
    %set(gca,'Position',[-0.5,-0.5,2,2]);
    %set(cbh,'Orientation','horizontal','Position',[0.30 0.10 0.50 0.02]);
    set(gca,'Position',[-0.52,-0.30,2.0,1.9]);
    set(get(gca,'Title'),'Units','normalized','Position',[0.5,0.4]);
    set(cbh,'Orientation','horizontal','Position',[0.25 0.40 0.50 0.02]);
   case 'LK',
    set(cbh,'Orientation','horizontal','Position',[0.50,0.30,0.45,0.03]);
   case 'DT',
    %set(cbh,'Location','East');
  end;

  axis(axis);

  % Basic bathymetry contour lines
  if ( any(plot_coastline) )
    [cc0,ch0] = contour(bath.lon,bath.lat,-bath.field,[0 0],'-','Color',[.8,.8,.8],'LineWidth',2);
    [cc1,ch1] = contour(bath.lon,bath.lat,-bath.field,[10 30 maxdepth],'-','Color',[.8,.8,.8],'LineWidth',1.5);
    [cc2,ch2] = contour(bath.lon,bath.lat,-bath.field,[5 15],':','Color',[.8,.8,.8],'LineWidth',1.5);
  end;

  % Coral reef habitat contour lines (use REEF or RFHD)
  if ( use_habitat_map )
    % If plot is already limited to habitat regions, polygons are superfluous
  else
    if ( allow_hard_bottom )
      hhdl = plot(rfhd.polylon,rfhd.polylat,'-','Color',[.4,.4,.4],'LineWidth',0.5);
    else
      hhdl = plot(reef.polylon,reef.polylat,'-','Color',[.4,.4,.4],'LineWidth',0.5);
    end;
  end;

  if ( plot_sites && exist('stns','var') )
    stnms = fieldnames(stns);
    for stnix=1:numel(stnms)
      stnm = stnms{stnix};
      if ( isfield(stns.(stnm),'fknms_seatemp') )
        plot_station_marker(stns.(stnm),false,stnclr{qix},0);
      end;
    end;
  end;

  titlename(sprintf('%s: %s (%+d W/m^2) daily reef %s',subrgnstr,q0str{qix},q0,phenomstr));
  if ( doPrint )
    print('-dpng',fullfile(figspath,sprintf('%s_%s_%s_%dW.png',figbasename,subrgn,phenomstr,abs(q0s(qix)))));
  end;

end; end; %if (doLog); for qix=scenarios(:)';


zeroZ = repmat(0,size(dTS));

for qix=scenarios(:)';
  % Net sea-surface heat flux, Q_0
  q0 = q0s(qix);
  disp(['Extreme T: Q0 = ',num2str(q0)]);

  fh=fmg;
  cbh=colorbar;
  set(get(cbh,'Label'),'String','^oC');
  switch (qix),
   case 1,
    surf(lon,lat,zeroZ',TC'); caxis([winter_min,winter_mean]); shading interp
   case 2,
    surf(lon,lat,zeroZ',TW'); caxis([winter_min,winter_mean]); shading interp
   case 3,
    surf(lon,lat,zeroZ',TS'); caxis([summer_mean,summer_max]); shading interp
    cm=hot; colormap(cm(end:-1:1,:));
   case 4,
    surf(lon,lat,zeroZ',TB'); caxis([summer_mean,summer_max]); shading interp
    cm=hot; colormap(cm(end:-1:1,:));
  end;
  axis(subrgnbox);
  set_pcolor_cursor;
  daspect([1,cosd(lat(1)),1]); 
  % Adjust graph appearance based on shape of sub-region
  switch (subrgn),
   case 'SE',
    view(90,90);
    set(cbh,'Orientation','horizontal','Position',[0.30 0.30 0.50 0.02]);
   case 'UK',
    % view(90,90);
    % set(cbh,'Location','East');
   case 'MK',
    view(25,90);
    %set(gca,'Position',[-0.5,-0.5,2,2]);
    %set(cbh,'Orientation','horizontal','Position',[0.30 0.10 0.50 0.02]);
    set(gca,'Position',[-0.52,-0.30,2.0,1.9]);
    set(get(gca,'Title'),'Units','normalized','Position',[0.5,0.4]);
    set(cbh,'Orientation','horizontal','Position',[0.25 0.40 0.50 0.02]);
   case 'LK',
    set(cbh,'Orientation','horizontal','Position',[0.50,0.30,0.45,0.03]);
   case 'DT',
    %set(cbh,'Location','East');
  end;

  axis(axis);

  % Basic bathymetry contour lines
  if ( any(plot_coastline) )
    [cc0,ch0] = contour(bath.lon,bath.lat,-bath.field,[0 0],'-','Color',[.8,.8,.8],'LineWidth',2);
    % [cc1,ch1] = contour(bath.lon,bath.lat,-bath.field,[10 30 maxdepth],'-','Color',[.8,.8,.8],'LineWidth',1.5);
    % [cc2,ch2] = contour(bath.lon,bath.lat,-bath.field,[5 15],':','Color',[.8,.8,.8],'LineWidth',1);
    [cc1,ch1] = contour(bath.lon,bath.lat,-bath.field,[30 maxdepth],'-','Color',[.8,.8,.8],'LineWidth',1.5);
  end;

  % If we are only plotting habitat map regions, polygons are superfluous
  if ( ~use_habitat_map )
    % Coral reef habitat contour lines (REEF or ReeF-and-HarD-bottom)
    if ( allow_hard_bottom )
      hhdl = plot(rfhd.polylon,rfhd.polylat,'-','Color',[.4,.4,.4],'LineWidth',0.5);
    else
      hhdl = plot(reef.polylon,reef.polylat,'-','Color',[.4,.4,.4],'LineWidth',0.5);
    end;
  end;

  if ( plot_sites && exist('stns','var') )
    stnms = fieldnames(stns);
    for stnix=1:numel(stnms)
      stnm = stnms{stnix};
      if ( isfield(stns.(stnm),'fknms_seatemp') )
        plot_station_marker(stns.(stnm),false,stnclr{qix},0);
      end;
    end;
  end;

  titlename(sprintf('%s: %s (%+d W/m^2) extreme temperature',subrgnstr,q0str{qix},q0));
  if ( doPrint )
    print('-dpng',fullfile(figspath,sprintf('%s_%s_extreme_%dW.png',figbasename,subrgn,abs(q0s(qix)))));
  end;

  if ( 0 && strcmp(subrgn,'UK') )
    disp('Hit Enter to zoom to sub-region'); pause;
    axis([-80.21,-80.135,25.33,25.40]); daspect([1,cosd(25),1]); view(2); colorbar('Location','EastOutside');
    if ( doPrint )
      print('-dpng',fullfile(figspath,sprintf('%s_%s_extreme_%dW_BLOWUP.png',figbasename,subrgn,abs(q0s(qix)))));
    end;
  end;

end; %for qix=scenarios(:)';


set_more;
