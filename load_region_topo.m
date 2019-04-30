function [topo,tlat,tlon,dlat,dlon] = load_region_topo(region,bbox,doPlot)
%function [topo,tlat,tlon,dlat,dlon] = load_region_topo(region,bbox,doPlot)
%
% Load 1-minute ETOPO topography [m, <0 for seafloor] for the MISST region
% named REGION, and subset to points within BBOX ([lat1 lat2 lon1 lon2]).
% 
% Last Saved Time-stamp: <Fri 2010-01-22 15:52:53 Eastern Standard Time gramer>

  datapath = get_ecoforecasts_path('coast');

  load(fullfile(datapath, [region '.topo.mat']),'topo','tlat','tlon');

  dlat = min(diff(tlat(:)));
  dlon = min(diff(tlon(:)));

  if ( ~isempty(bbox) )

    latix = find(bbox(1) <= tlat & tlat <= bbox(2));
    tlat = tlat(latix);

    % ETOPO-01 prefers longitudes in form [-180:180]
    if (bbox(3) > 180); bbox(3) = bbox(3) - 360; end;
    if (bbox(4) > 180); bbox(4) = bbox(4) - 360; end;

    if ( bbox(3) > bbox(4) )
      % Region crosses anti-meridian
      lonix1 = find(bbox(3) <= tlon);
      lonix2 = find(tlon <= bbox(4));
      lonix = [lonix1 lonix2];
    else
      % Region is contiguous
      lonix = find(bbox(3) <= tlon & tlon <= bbox(4));
    end;
    tlon = tlon(lonix);
    tlon(tlon < 0) = 360 + tlon(tlon < 0);

    % Subset topography to just our bounding box
    topo = topo(latix,lonix);

  end;

  if ( exist('doPlot','var') && doPlot )
    figure;

    [c,h] = contour(tlon,tlat,topo,[0 -1 -3 -10 -30 -100 -300 -1000 -3000]);
    clabel(c,h);

    if ( ~isempty(bbox) )
      ttl = sprintf('%s: %s [%g %g %g %g]', mfilename, region, bbox);
    else
      ttl = sprintf('%s: %s [WHOLE REGION]', mfilename, region);
    end;
    ttl = strrep(ttl, '_', '\_');
    title(ttl);
    maximize_graph;

    if ( ~isempty(bbox) )
      plotlims = bbox;
      if (bbox(3)<0); plotlims(3)=360+bbox(3); end;
      if (bbox(4)<0); plotlims(4)=360+bbox(4); end;
      ylim([plotlims(1) plotlims(2)]); xlim([plotlims(3) plotlims(4)]);
    end;

    lontick;
  end;

return;
