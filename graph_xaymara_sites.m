1;
% Graph bathymetric properties of each of Xaymara Serrano's coral genetic
% sampling sites in the Florida Keys and Dry Tortugas

if ( ~exist('doPrint','var') )
  doPrint = false;
end;

if ( ~exist('stn','var') || ~isfield(stn,'ps') )
  stn=[]; clear stn
  stn = extract_xaymara_sites;
end;

% Are there any anomalies in site name/location between P. ast. and M. cav.?
for ix=1:numel(stn.mlons);
  msix = find(strcmp(stn.ms(ix).station_name,{stn.ps.station_name}));
  if ( isempty(msix) )
    warning('No match for M. cav. site %s',stn.ms(ix).station_name);
  else
    d = distance_wgs84(stn.mlats(ix),stn.mlons(ix),stn.plats(msix),stn.plons(msix));
    str=''; if ( d > 0.1 ); str='  ***'; end;
    disp({ix,stn.ms(ix).station_name,msix,stn.ps(msix).station_name,d,str});
  end;
end;
clear d ix m mx str


mixen = 1:numel(stn.mlons);
%mixen = 1:14;

bath = read_hires_bathymetry_for_field({stn.mlons(mixen),stn.mlats(mixen)},false);

%fmg; plot_hires_coastline(bath.ngdc_hires_bathy);
plot_hires_bathymetry(bath,-[0:5:80],[],true,@contour);
for mix=mixen(:)'
  plot(stn.ms(mix).lon,stn.ms(mix).lat,'bo');
  nm = [stn.ms(mix).station_name,' (',num2str(stn.mbetas(mix)),')'];
  text(stn.ms(mix).lon,stn.ms(mix).lat,nm,'Color','b','Rotation',((mix-1)*5));
  pix = find(strcmp(stn.ms(mix).station_name,{stn.ps.station_name}));
  plot(stn.ps(pix).lon,stn.ps(pix).lat,'ro');
  nm = [stn.ps(pix).station_name,' (',num2str(stn.pbetas(pix)),')'];
  text(stn.ps(pix).lon,stn.ps(pix).lat,nm,'Color','r','Rotation',90+((mix-1)*5));
end;
clear ans mix mixen nm pix


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
