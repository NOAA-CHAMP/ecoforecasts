1;

if ( ~exist('doSave','var') || isempty(doSave) )
  doSave = false;
end;
if ( ~exist('startOver','var') || isempty(startOver) )
  startOver = false;
end;
if ( ~exist('smoothPeaks','var') || isempty(smoothPeaks) )
  smoothPeaks = true;
end;

if ( ~exist('LON','var') || isempty(LON) )
  [LON,LAT] = meshgrid(lon,lat);
end;

if ( startOver )
  if ( exist('fh','var') && ishandle(fh) )
    close(fh);
  end;
  clear fh
  newdat=[]; clear newdat

  clear ans ig s nanix nanlandix yix xix previx nextix badlandix nbadland ixes Top Tops Right Rights Bottom Bottoms Left Lefts edgeix ix nanix nixen tixen;
  FLAT=[]; FLON=[]; fdat=[]; rawdat=[]; clear FLAT FLON fdat rawdat flon flat
end;

if ( ~exist('fh','var') || ~ishandle(fh) )
  if ( ~exist('newdat','var') || isempty(newdat) )
    newdat = dat;
  end;

  fh=[];
  fh = fmg; contourf(lon,lat,newdat,-80:2:6); colorbar; 
  %fh = fmg; contourf(lon,lat,newdat,-1000:50:100); colorbar;
  contour(lon,lat,newdat,[0,0],'Color','k','LineWidth',2);
  daspect([1,cosd(lat(1)),1]); set_surf_cursor;
  titlename(['Bathy/DEM at ',num2str(target_res),' m - DRAW BOX FOR SPLINE']);

  edgeix = find_field_holes(dat);
  plot(LON(edgeix),LAT(edgeix),'r.','MarkerSize',0.5);
end;

if ( startOver || ~exist('corners','var') )
  corners = [];
end;

figure(fh);
r = getrect;
b = rect2bbox(r);
ixes = bbox2ind(b,LON,LAT);
corners(end+1,1:8) = [ixes(1,1),ixes(1,2), ixes(2,1),ixes(2,2), ixes(3,1),ixes(3,2), ixes(4,1),ixes(4,2)];

rectangle('Position',r);
text(b(2),b(4),[' \leftarrow ',num2str(size(corners,1))]);

if ( exist('fh2','var') && ishandle(fh2) )
  close(fh2);
end;
clear fh2;

fill_dem_bathy;

if ( exist('fh','var') && ishandle(fh) )
  close(fh);
end;
fh=fh2;
figure(fh);

if ( startOver )
  startOver = false;
end;
