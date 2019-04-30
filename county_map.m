function [fh,ah,lh,th,bath,isocs,isoch] = county_map(varargin)
%function [fh,ah,lh,th,bath,isocs,isoch] = county_map(['BBOX',bbox]['STATES',states],['Counties',counties],['Isobaths',isobaths],['BathymetryStruct',bath],['Figure',fh][NM1,VAL1,...])
%
% Load GADM36 shape file(s) and draw a high-resolution map showing US State and
% County boundaries. Named args: 'BBOX' (DEFAULT: [-180,180,-90,90]), 'STATES'
% (D: string({"Florida"})), 'COUNTIES' (D: string({'Palm Beach', 'Broward',
% 'Miami-Dade', 'Monroe'})), 'ISOBATHS', 'BATH' (STRUCT w/.ngdc_hires_bathy),
% 'FIG' Figure handle, 'FONTSIZE' (D: 18), 'LINEWIDTH' (D: 2.5) for isobaths,
% county 'FACECOLOR' and isobath 'LINECOLOR' (D: cmu.colors('light gray')).
% 'FOCUSCOLOR' (D: 'army green') is a separate FaceColor for COUNTIES list.
%
% CALLS: Map/SHAPEREAD, GEOSHOW; Vision/BBOXOVERLAPRATIO, TEXT. If 'ISOBATHS'
% is given, calls PLOT_HIRES_BATHYMETRY. If 'ISOBATHS' is given but 'BATH' is
% not, calls READ_HIRES_BATHYMETRY first.
%
% Last Saved Time-stamp: <Wed 2018-10-17 16:38:36 EDT lew.gramer>

  [bath,varargin] = getarg(varargin,'BATH');
  [bbox,varargin] = getarg(varargin,'BBOX','default',[-180,+180,-90,+90]);
  [fh,varargin] = getarg(varargin,'FIG');
  [isobaths,varargin] = getarg(varargin,'ISO');
  [states,varargin] = getarg(varargin,'STAT','default',string({'Florida'}));
  [counties,varargin] = getarg(varargin,'COUNTI','default',string({'Palm Beach','Broward','Miami-Dade','Monroe'}));

  [ctyclr,varargin] = getarg(varargin,'FACECO','default',cmu.colors('light gray'));
  [focclr,varargin] = getarg(varargin,'FOCUSC','default',cmu.colors('army green'));
  [fontsz,varargin] = getarg(varargin,'FONTSI','default',18);
  [linewi,varargin] = getarg(varargin,'LINEWI','default',2.5);
  [linest,varargin] = getarg(varargin,'LINEST','default','-');
  [lineco,varargin] = getarg(varargin,'LINECO','default',cmu.colors('slate gray'));

  if ( ischar(bbox) || isstring(bbox) )
    switch (lower(bbox)),
     case 'south florida',	bbox = [-88,-78,+24.5,+27.5];
    end;
  end;


  geodat = shaperead(get_ecoforecasts_path('coast/gadm36_USA_2'),'UseGeoCoords',true,'Selector',{@(name)ismember(name,states),'NAME_1'});
  
  bboxes = [geodat.BoundingBox];
  inix = find(bboxOverlapRatio(bbox2rect(bbox),bbox2rect(bboxes))>0);
  ctix = intersect(inix,find(ismember(string({geodat.NAME_2}),counties)));
  
  if ( isempty(ctix) )
    error('No Counties matching %s in States %s',join(states),join(counties));
  end;

  aLat = nanmean(geodat(ctix(1)).Lat);

  if ( isempty(fh) )
    fh = fmg;
  end;
  figure(fh);
  ah = gca;
  set(ah,'FontSize',fontsz);
  %usamap(states);
  lh = geoshow(geodat(inix),'FaceColor',ctyclr);
  if ( ~isempty(focclr) )
    lh = geoshow(geodat(ctix),'FaceColor',focclr);
  end;
  daspect([1,cosd(aLat),1]);
  for ix=1:numel(ctix)
    % %lon = nanmin(geodat(ctix(ix)).Lon)+0.05;
    % %lat = mean([nanmin(geodat(ctix(ix)).Lat),nanmax(geodat(ctix(ix)).Lat)]);
    % %th = text(lon,lat,upper(geodat(ctix(ix)).NAME_2),'HorizontalAlignment','center','FontSize',fontsz);
    % lon = mean([nanmin(geodat(ctix(ix)).Lon),nanmax(geodat(ctix(ix)).Lon)]);
    % lat = mean([nanmin(geodat(ctix(ix)).Lat),nanmax(geodat(ctix(ix)).Lat)]);
    % %th = text(lon,lat,upper(geodat(ctix(ix)).NAME_2),'HorizontalAlignment','left','FontSize',fontsz);
    % th = text(lon,lat,upper(geodat(ctix(ix)).NAME_2),'HorizontalAlignment','center','FontSize',fontsz);

    lon = nanmin(geodat(ctix(ix)).Lon)+0.05;
    lat = mean([nanmin(geodat(ctix(ix)).Lat),nanmax(geodat(ctix(ix)).Lat)]);
    th = text(lon,lat,upper(geodat(ctix(ix)).NAME_2),'HorizontalAlignment','left','FontSize',fontsz,'Color','w');
  end;

  if ( ~isempty(isobaths) )
    if ( isempty(bath) )
      bath.station_name = 'County_Map';
      bath.lon = nanmean([geodat(inix).Lon]);
      bath.lat = nanmean([geodat(inix).Lat]);
      [radx,rady] = geographic_radius_m([geodat(inix).Lon],[geodat(inix).Lat]);
      bath = read_hires_bathymetry(bath,[radx,rady],[],false);
    end;
    [bath,isocs,isoch] = plot_hires_bathymetry(bath,isobaths,[],true,@contour,false,[],[],fh);
    set(isoch,'LineWidth',linewi,'LineStyle',linest,'Color',lineco);
    colorbar off;
    title([]);
    axis(bbox);
  end;

return;
