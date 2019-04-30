1;
% SCRIPT CMPEFKEYS.m - compare contours and sections of U,V,T for output of
% model eFKEYS; if doFKEYS True, do FKEYS instead; if doGOM2, do GOM2 instead.

if ( ~exist('doPrint','var') )
  doPrint = false;
end;
if ( doPrint )
  figpath = get_ecoforecasts_path('figs');
end;

if ( ~exist('doFKEYS','var') )
  doFKEYS = false;
end;
if ( ~exist('doGOM2','var') )
  doGOM2 = false;
end;

if ( doGOM2 )
  if ( ~exist('mdl','var') || ~isfield(mdl,'model') || ~strcmpi(mdl.model,'GOM2') )
    mdl=[]; clear mdl;
    angom2;
  end;
  dts = datenum(2014,1,[1:31],12,0,0);
  mdl.model = 'GOM2';
elseif ( doFKEYS )
  if ( ~exist('mdl','var') || ~isfield(mdl,'model') || ~strcmpi(mdl.model,'FKEYS') )
    mdl=[]; clear mdl;
    anfkeys;
  end;
  dts = datenum(2008,1,1,6,0,0):(6/24):datenum(2008,1,7,18,0,0);
  mdl.model = 'FKEYS';
else
  if ( ~exist('mdl','var') || ~isfield(mdl,'model') || ~strcmpi(mdl.model,'eFKEYS') )
    mdl=[]; clear mdl;
    anefkeys;
  end;
  dts = datenum(2012,1,1,6,0,0):(6/24):datenum(2012,1,7,18,0,0);
  mdl.model = 'eFKEYS';
end;

dix = find(ismember(mdl.d,dts));
dstr = [' (',datestr(mdl.d(dix(1))),' - ',datestr(mdl.d(dix(end))),')'];


% Perform model-in situ comparisons

cstnms = {'fwyf1','mlrf1','lonf1','smkf1','looe1','sanf1'};
ncstnms = numel(cstnms);

for stix=1:ncstnms
  stnm = cstnms{stix};
  %mdl.(stnm)=[];
  %mdl = rmfield(mdl,stnm);
  if ( ~isfield(mdl,stnm) )

    mdl.(stnm) = get_station_from_station_name(stnm);
    mdl.(stnm) = station_optimal_isobath_orientation(mdl.(stnm));

    [mdl.(stnm).lonerr,mdl.(stnm).lonix] = min(abs(mdl.lon-mdl.(stnm).lon));
    [mdl.(stnm).laterr,mdl.(stnm).latix] = min(abs(mdl.lat-mdl.(stnm).lat));

    mdl.(stnm).t.z = mdl.z;
    mdl.(stnm).t.date = mdl.d;
    mdl.(stnm).t.prof = mdl.t(:,:,mdl.(stnm).latix,mdl.(stnm).lonix);
    mdl.(stnm).t.data = nanmean(mdl.(stnm).t.prof,2);
    for zix = 1:numel(mdl.z)
      fld = ['t',num2str(mdl.z(zix))];
      mdl.(stnm).(fld).z = mdl.z(1);
      mdl.(stnm).(fld).date = mdl.d;
      mdl.(stnm).(fld).data = mdl.t(:,zix,mdl.(stnm).latix,mdl.(stnm).lonix);
    end;
    clear fld zix

    if ( isfield(mdl,'s') )
      mdl.(stnm).s.z = mdl.z;
      mdl.(stnm).s.date = mdl.d;
      mdl.(stnm).s.prof = mdl.s(:,:,mdl.(stnm).latix,mdl.(stnm).lonix);
      mdl.(stnm).s.data = nanmean(mdl.(stnm).s.prof,2);
      for zix = 1:numel(mdl.z)
        fld = ['s',num2str(mdl.z(zix))];
        mdl.(stnm).(fld).z = mdl.z(1);
        mdl.(stnm).(fld).date = mdl.d;
        mdl.(stnm).(fld).data = mdl.s(:,zix,mdl.(stnm).latix,mdl.(stnm).lonix);
      end;
      clear fld zix
    end; %if ( isfield(mdl,'s') )

    mdl.(stnm).u.z = mdl.z;
    mdl.(stnm).u.date = mdl.d;
    mdl.(stnm).u.prof = mdl.u(:,:,mdl.(stnm).latix,mdl.(stnm).lonix);
    mdl.(stnm).u.data = nanmean(mdl.(stnm).u.prof,2);

    mdl.(stnm).v.z = mdl.z;
    mdl.(stnm).v.date = mdl.d;
    mdl.(stnm).v.prof = mdl.v(:,:,mdl.(stnm).latix,mdl.(stnm).lonix);
    mdl.(stnm).v.data = nanmean(mdl.(stnm).v.prof,2);

    % Reorient vector currents to local isobath: 'x'=cross-shore, 'l'=longshore
    mdl.(stnm).x.z = mdl.z;
    mdl.(stnm).x.date = mdl.d;
    mdl.(stnm).l.z = mdl.z;
    mdl.(stnm).l.date = mdl.d;
    [mdl.(stnm).x.prof,mdl.(stnm).l.prof] = ...
        reorient_vectors(mdl.(stnm).isobath_orientation,mdl.(stnm).u.prof,mdl.(stnm).v.prof);

    mdl.(stnm).x.data = nanmean(mdl.(stnm).x.prof,2);
    mdl.(stnm).l.data = nanmean(mdl.(stnm).l.prof,2);
  end; %if ( ~isfield(mdl,stnm) )

end; %for stix=1:ncstnms


min_t = +inf;
max_t = -inf;
min_s = +inf;
max_s = -inf;
min_curr = +inf;
max_curr = -inf;

for stix=1:ncstnms
  stnm = cstnms{stix};
  min_t = nanmin([min_t,nanmin(mdl.(stnm).t.prof(dix,:))]);
  max_t = nanmax([max_t,nanmax(mdl.(stnm).t.prof(dix,:))]);
  if ( isfield(mdl.(stnm),'s') )
    min_s = nanmin([min_s,nanmin(mdl.(stnm).s.prof(dix,:))]);
    max_s = nanmax([max_s,nanmax(mdl.(stnm).s.prof(dix,:))]);
  end;
  min_curr = nanmin([min_curr,nanmin(mdl.(stnm).u.prof(:)),nanmin(mdl.(stnm).v.prof(dix,:))]);
  max_curr = nanmax([max_curr,nanmax(mdl.(stnm).u.prof(:)),nanmax(mdl.(stnm).v.prof(dix,:))]);
end;

max_spd = max(abs(max_curr),abs(min_curr));


%stnixen = [1,2,4,5,6];
stnixen = [1,2,4,5,6];

if 1;
  fmg;
  ncols = numel(stnixen);
  colix = 1;
  %min_hist = floor(-max_spd); max_hist = ceil(+max_spd); d_hist = 0.025;
  min_hist = -0.2; max_hist = +0.2; d_hist = 0.01;
  for stix=stnixen(:)';
    stnm = cstnms{stix};

    %spt(2,ncols,colix);
    subplot(2,ncols,colix);
    hist(mdl.(stnm).l.data(dix),min_hist:d_hist:max_hist); 
    xlim([min_hist-d_hist,max_hist+d_hist]);
    xlabel([upper(stnm),' alongshore (m/s)']);

    %spt(2,ncols,colix+ncols);
    subplot(2,ncols,colix+ncols);
    hist(mdl.(stnm).x.data(dix),min_hist:d_hist:max_hist); 
    xlim([min_hist-d_hist,max_hist+d_hist]);
    xlabel([upper(stnm),' cross-shore (m/s)']);

    colix = colix + 1;
  end;
  if ( doPrint )
    print('-dpng',fullfile(figpath,[mdl.model,'.uv_hist.png']));
  end;
end;

%stnm = 'fwyf1';
%stnm = 'mlrf1';
%stnm = 'sanf1';
% if ( stnm )
%%for stix=1:ncstnms
%for stix=stnixen(:)'
for stix=stnixen([2,4])
  stnm = cstnms{stix};
  disp(stnm);

  % %fmg;
  % stn = []; clear stn
  % stn = get_station_from_station_name(stnm);
  % stn = station_optimal_isobath_orientation(stn);
  % stn = read_hires_bathymetry(stn,[12e3,12e3],[],false);

  % Get high-resolution bathymetry
  mdl.(stnm) = read_hires_bathymetry(mdl.(stnm),[12e3,12e3],[],false);

  % Draw bathymetry map
  mdl.(stnm) = plot_hires_bathymetry(mdl.(stnm),-[0:5:80],[12e3,12e3],true,[],false,[],false);
  fh = gcf;
  axis(axis);
  [c,h] = contour(mdl.lon,mdl.lat,squeeze(nanmean(mdl.t(dix,1,:,:))),[floor(min_t):0.25:ceil(max_t)]);
  clabel(c,h);
  c=[]; clear c h


  %fmg;
  [res,mdl.(stnm)] = plot_bathy_transect(mdl.(stnm),10,mdl.(stnm).isobath_orientation+90,'ngdc_hires_bathy');
  [dx,az] = distance_wgs84(mdl.(stnm).lat,mdl.(stnm).lon,res.lat,res.lon);
  dx(round(az) ~= mdl.(stnm).isobath_orientation+90) = -dx(round(az) ~= mdl.(stnm).isobath_orientation+90);
  res.dx = dx;

  lonix = interp1(mdl.lon,1:mdl.nlon,res.lon,'nearest');
  latix = interp1(mdl.lat,1:mdl.nlat,res.lat,'nearest');
  ix = sub2ind([size(mdl.t,3),size(mdl.t,4)],latix,lonix);
  %t = squeeze(mdl.t(1,:,ix));
  t = squeeze(nanmean(mdl.t(dix,:,ix)));
  u = squeeze(nanmean(mdl.u(dix,:,ix)));
  v = squeeze(nanmean(mdl.v(dix,:,ix)));
  % disp([nanmin(u(:)),nanmax(u(:))]);
  % disp([nanmin(v(:)),nanmax(v(:))]);
  [x,l] = reorient_vectors(mdl.(stnm).isobath_orientation,u,v);

  disp([nanmin(t(:)),nanmax(t(:))]);
  disp([nanmin(l(:)),nanmax(l(:))]);
  disp([nanmin(x(:)),nanmax(x(:))]);

  [c,h] = contourf(res.dx,-mdl.z,t,[11.0:0.5:24.5]);
  set(gca,'CLim',[11.0,24.5]);
  clabel(c,h);
  c=[]; clear c h
  plot(res.dx,res.field,'k-','LineWidth',2);
  xlim([-10,+10]); ylim([-140,0]);
  titlename([mdl.model,': ',upper(stnm),' bathymetry vs. T cross-shore profile',dstr]);
  if ( doPrint )
    print('-dpng',fullfile(figpath,[mdl.model,'.',stnm,'.section_t_mean.png']));
  end;

  %fmg;
  plot_bathy_transect(mdl.(stnm),10,mdl.(stnm).isobath_orientation+90,'ngdc_hires_bathy');
  [c,h] = contourf(res.dx,-mdl.z,l,[-0.300:0.05:+1.200]);
  set(gca,'CLim',[-0.300,+1.200]);
  clabel(c,h);
  c=[]; clear c h
  plot(res.dx,res.field,'k-','LineWidth',2);
  xlim([-10,+10]); ylim([-140,0]);
  titlename([mdl.model,': ',upper(stnm),' bathymetry vs. alongshore V',dstr]);
  if ( doPrint )
    print('-dpng',fullfile(figpath,[mdl.model,'.',stnm,'.section_ls_mean.png']));
  end;

  %fmg;
  plot_bathy_transect(mdl.(stnm),10,mdl.(stnm).isobath_orientation+90,'ngdc_hires_bathy');
  [c,h] = contourf(res.dx,-mdl.z,x,[-0.150:0.010:+0.150]);
  set(gca,'CLim',[-0.150,+0.150]);
  clabel(c,h);
  c=[]; clear c h
  plot(res.dx,res.field,'k-','LineWidth',2);
  xlim([-10,+10]); ylim([-140,0]);
  titlename([mdl.model,': ',upper(stnm),' bathymetry vs. cross-shore U',dstr]);
  if ( doPrint )
    print('-dpng',fullfile(figpath,[mdl.model,'.',stnm,'.section_xs_mean.png']));
  end;

  % Draw cross-shore section on bathymetry map
  figure(fh);
  plot(res.lon,res.lat,'k-','LineWidth',2);
  clear fh
  if ( doPrint )
    print('-dpng',fullfile(figpath,[mdl.model,'.',stnm,'.surface_t_mean.png']));
  end;
end; %for stix=stnixen([2,4])'


clear ans cstnms stix stnm
clear az dx ig ix latix lonix s t
%clear datpath hycpath matfname ncols ncstnms ans
%stn=[]; clear res stn min_t max_t min_s max_s min_curr max_curr max_spd
