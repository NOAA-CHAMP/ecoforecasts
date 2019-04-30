1;

if ( ~exist('doSave','var') || isempty(doSave) )
  doSave = false;
  %doSave = true;
end;
if ( ~exist('startOver','var') || isempty(startOver) )
  startOver = true;
end;
if ( ~exist('smoothPeaks','var') || isempty(smoothPeaks) )
  smoothPeaks = true;
end;

%% Sub-regions that are "sufficiently" NaN-free to spline fit within

% SUB-REGION CORNER INDICES:

if ( ~exist('corners','var') || isempty(corners) )

  switch ( batnm ),
   case 'pibhmc_bathy_5m_saipan',
    corners = [ ...
        480,255  ,  779,255  ,  779,653  ,  480,653 ; ...
        872,626  ,  926,626  ,  926,673  ,  872,673 ; ...
        816,691  ,  888,691  ,  888,761  ,  816,761 ; ...
        869,611  ,  913,611  ,  913,653  ,  869,653 ; ...
        099,375  ,  116,375  ,  116,394  ,  099,394 ; ...
        357,559  ,  436,559  ,  436,604  ,  357,604 ; ...
        291,145  ,  362,145  ,  362,223  ,  291,223 ; ...
              ];

   case 'pibhmc_bathy_5m_tinian.PARTIAL',
    corners = [ ...
        264,669  ,  483,669  ,  483,916  ,  264,916 ; ...
        170,558  ,  513,558  ,  513,680  ,  170,680 ; ...
        196,495  ,  485,495  ,  485,581  ,  196,581 ; ...
        267,385  ,  472,385  ,  472,494  ,  267,494 ; ...
        280,336  ,  558,336  ,  558,384  ,  280,384 ; ...
        410,280  ,  469,280  ,  469,336  ,  410,336 ; ...
        462,439  ,  495,439  ,  495,478  ,  462,478 ; ...
              ];

   case 'pibhmc_bathy_5m_rota',
    corners = [ ...
        146,245 , 362,245 , 362,386 , 146,386 ; ...
        055,061 , 352,061 , 352,247 , 055,247 ; ...
              ];

  end;

end;

for rix=1:size(corners,1)
  Tops(rix) = corners(rix,6);
  Rights(rix) = corners(rix,3);
  Bottoms(rix) = corners(rix,2);
  Lefts(rix) = corners(rix,1);
end;

clear rix


if ( ~exist('rawdat','var') || ~all(size(dat)==size(rawdat)) )
  disp('Saving field as RAWDAT');
  rawdat = dat;
end;

if ( startOver )
  newdat=[]; clear newdat
  newdat = dat;
end;

for ix=1:numel(Tops)
  Top = Tops(ix);
  Right = Rights(ix);
  Bottom = Bottoms(ix);
  Left = Lefts(ix);

  nixen = Left:Right;
  tixen = Bottom:Top;

  flon = lon(nixen);
  flat = lat(tixen);
  fdat = newdat(tixen,nixen);

  s = warning('OFF','MATLAB:interp2:NaNstrip');
  nanix = find(isnan(fdat));
  [FLON,FLAT] = meshgrid(flon,flat);
  fdat(nanix) = interp2(flon,flat,fdat,FLON(nanix),FLAT(nanix),'spline',nan);
  warning(s); clear s

  % Try to eliminate any spurious "land" that was just introduced by splines
  if ( smoothPeaks )
    nanlandix = find(fdat(nanix)>0);
    nanlandix = nanix(nanlandix);
    if ( ~isempty(nanlandix) )
      [yix,xix] = ind2sub(size(fdat),nanlandix);
      badcix = find(1==xix | xix==size(fdat,2) | 1==yix | yix==size(fdat,1));
      if ( ~isempty(badcix) )
        warning('Peaks next to field edge!');
keyboard;
        xix(badcix) = [];
        yix(badcix) = [];
        nanlandix(badcix) = [];
      end;
      % Look for row-wise peaks
      previx = sub2ind(size(fdat),yix,xix-1);
      nextix = sub2ind(size(fdat),yix,xix+1);
      badlandix = find(fdat(previx)<0 & fdat(nextix)<0);
      fdat(nanlandix(badlandix)) = mean([fdat(previx(badlandix)),fdat(nextix(badlandix))],2);
      nbadland = numel(nanlandix(badlandix));

      % Look for column-wise peaks
      previx = sub2ind(size(fdat),yix-1,xix);
      nextix = sub2ind(size(fdat),yix+1,xix);
      badlandix = find(fdat(previx)<0 & fdat(nextix)<0);
      fdat(nanlandix(badlandix)) = mean([fdat(previx(badlandix)),fdat(nextix(badlandix))],2);
      nbadland = nbadland + numel(nanlandix(badlandix));

      if ( nbadland > 0 )
        disp(['Removed ',num2str(nbadland),' "peak" points from box #',num2str(ix)]);
      end;

      clear nanix nanlandix yix xix previx nextix badlandix nbadland
    end;
  end;

  FLON=[]; FLAT=[]; clear FLON FLAT

  newdat(tixen,nixen) = fdat;
  fdat=[]; clear fdat flat flon
end;


fh1=fmg; contourf(lon,lat,newdat,-1000:50:100); colorbar;
contour(lon,lat,newdat,[0,0],'Color','k','LineWidth',2);
daspect([1,cosd(lat(1)),1]); set_surf_cursor;
titlename(['Bathy/DEM at ',num2str(target_res),' m AFTER SPLINE']);


fh2=fmg; contourf(lon,lat,newdat,-80:2:6); colorbar; 
%fh2=fmg; contourf(lon,lat,newdat,-1000:50:100); colorbar;
contour(lon,lat,newdat,[0,0],'Color','k','LineWidth',2);
daspect([1,cosd(lat(1)),1]); set_surf_cursor;
titlename(['Bathy/DEM at ',num2str(target_res),' m AFTER SPLINE']);

edgeix = find_field_holes(rawdat);
[LON,LAT] = meshgrid(lon,lat);
plot(LON(edgeix),LAT(edgeix),'r.','MarkerSize',0.5);
LON=[]; LAT=[]; clear LON LAT

%plot(lon([Lefts,Lefts,Rights,Rights]),lat([Tops,Bottoms,Tops,Bottoms]),'k-p','MarkerFaceColor','w');
for ix=1:numel(Tops)
  Top = Tops(ix);
  Right = Rights(ix);
  Bottom = Bottoms(ix);
  Left = Lefts(ix);
  plot(lon([Left,Right,Right,Left,Left]),lat([Top,Top,Bottom,Bottom,Top]),'k-p','MarkerFaceColor','w');
  text(lon(Right),lat(Top),[' \leftarrow ',num2str(ix)]);
end;

clear Top Tops Right Rights Bottom Bottoms Left Lefts ans edgeix ix nanix nixen tixen

if ( doSave )
  disp(['Saving ',cmbpath]);
  dat = newdat;
  newdat=[]; clear newdat
  save(cmbpath,'url','demurl','lon','lat','dat','rawdat','-v7.3');
else
  disp(['NOT saving ',cmbpath]);
end;

