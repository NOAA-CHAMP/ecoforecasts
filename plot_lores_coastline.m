function [ch,fh] = plot_lores_coastline(lon_range,lat_range,map,fh)
%function [ch,fh] = plot_lores_coastline(lon_range,lat_range,map,fh)
%
% Use PLOTM from the Map Toolkit (v.) to plot basin- or global-scale coast
% line within region covered by the vectors LON_RANGE and LAT_RANGE. Plots in
% current figure (FMG) unless optional arg FH is specified. Plots using map
% projection 'eqacylin', unless optional CHAR string MAP is specified.
%
% Optionally returns line handles CH returned by PLOTM, and FIGURE FH.
%
% Last Saved Time-stamp: <Sun 2018-05-27 20:12:17 Eastern Daylight Time gramer>

  persistent coastlat;
  persistent coastlon;

  if ( ~exist('map','var') || isempty(map) )
    map = 'eqacylin';
  end;
  if ( ~exist('fh','var') || isempty(fh) )
    fh = fmg;
  end;

  if ( ~exist('coastlat','var') || isempty(coastlat) )
    load coastlines
  end;

  switch (map),
   case 'ortho',
    axesm(map,'origin',[nanmean(lat_range(:)),nanmean(lon_range(:))]);
   case 'eqacylin',
    axesm(map,'maplonlimit',[min(lon_range(:)),max(lon_range(:))],...
          'maplatlimit',[min(lat_range(:)),max(lat_range(:))]);
   otherwise,
    warning('How do I convert that range to Map type "%s"?',map);
    axesm(map,'maplonlimit',[min(lon_range(:)),max(lon_range(:))],...
          'maplatlimit',[min(lat_range(:)),max(lat_range(:))]);
  end;

  axis off;
  gridm on;
  framem on;
  mlabel('equator')
  % plabel(0);
  plabel('prime');
  plabel('fontweight','bold')

  gridm('MLineLocation',10, 'MLabelLocation',30,...
        'PLineLocation',10, 'PLabelLocation',30);

  ch = plotm(coastlat,coastlon,'LineWidth',2);
  %ch = fillm(coastlat,coastlon,'FaceColor','b');

  setm(gca,'FontSize',20);
  set(gca,'FontSize',20);

return;
