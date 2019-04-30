1;

if ( ~exist('doPrint','var') )
  doPrint = false;
end;

if ( ~exist('stn','var') || ~isfield(stn,'ps') )
  stn=[]; clear stn
  stn = extract_xaymara_sites;
end;

msix=find(ismember({stn.ps.station_name},{stn.ms.station_name}));

% % Which one of these is (are) not like the other(s)??
% for ix=1:numel(stn.mlons);
%   d = distance_wgs84(stn.mlats(ix),stn.mlons(ix),stn.plats(msix),stn.plons(msix));
%   % if ( all(d > 500/111e3) ); disp({stn.ms(ix).station_name,min(d)}); end;
%   [m,mx] = min(d);
%   str=''; if ( ~strcmpi(stn.ms(ix).station_name,stn.ps(msix(mx)).station_name) ); str='  ***'; end;
%   disp({stn.ms(ix).station_name,stn.ps(msix(mx)).station_name,m,str});
% end;
% clear d ix m mx str

%for stix = 1:3
for stix = 1:stn.nps
  disp([stn.ps(stix).station_name,':',stn.ps(stix).station_desc]);
  plot_hires_bathymetry(stn.ps(stix),-[0:1:50]);
  set(gca,'CLim',[-30,0]);
  axis(axis);

  plot(stn.lons,stn.lats,'ws','MarkerFaceColor','w');

  plot(stn.mlons,stn.mlats,'k.','MarkerSize',0.5);
  plot(stn.plons,stn.plats,'r.','MarkerSize',0.5);
  for ix=1:stn.nms
    text(stn.mlons(ix),stn.mlats(ix),[' \leftarrow M:',stn.ms(ix).station_name],...
         'Color','w','HorizontalAlignment','left','Rotation',(ix*7)-90,'FontWeight','bold');
  end;
  for ix=1:stn.nps
    text(stn.plons(ix),stn.plats(ix),['P:',stn.ps(ix).station_name,' \rightarrow '],...
         'Color','w','HorizontalAlignment','right','Rotation',(ix*7)-90,'FontWeight','bold');
  end;
  if doPrint; print('-dpng',[lower(stn.ps(stix).station_name),'_bathy.png']); end;
end;


% x = get_station_from_station_name('sanf1'); x=station_optimal_isobath_orientation(x); x=station_ngdc_offshore_slope(x), interp2(stn.ngdc_hires_bathy.lon,stn.ngdc_hires_bathy.lat,stn.beta,x.lon,x.lat), [ig,xix] = min(abs(stn.ngdc_hires_bathy.lon-x.lon)), [ig,yix] = min(abs(stn.ngdc_hires_bathy.lat-x.lat)), fmg; contourf(stn.ngdc_hires_bathy.lon(xix-40:xix+40),stn.ngdc_hires_bathy.lat(yix-25:yix+25),stn.ngdc_hires_bathy.field(yix-25:yix+25,xix-40:xix+40)); colorbar;

clear ans stix ix
