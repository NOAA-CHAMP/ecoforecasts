1;

if ( ~exist('lon','var') || ~exist('lat','var')  || ~exist('h','var')  || ~exist('bet','var') )
  lon=[]; lat=[]; h=[]; bet=[]; quilted_h=[]; quilted_bet=[]; clear lon lat h bet quilted_h quilted_bet
  disp('Loading FRT_depth_and_beta.mat');
  load('FRT_depth_and_beta.mat');

  % No HC over land, and depths are positive
  h(h>0) = 0;
  h = abs(h);
end;

if ( ~exist('quilted_h','var') )
  quilted_h = h;
  quilted_bet = bet;
end;

if 0;
  fmg; contourf(lon,lat,-h',-[0:10:120]); colorbar; daspect([1,cosd(lat(1)),1]); titlename('Original');
  ax=axis;
end;

if 1;
  disp('Quilting SANF1 bathymetry into FRT');
  if ( ~exist('sanf1','var') || ~isfield(sanf1,'ngdc_hires_bathy') )
    sanf1=[]; clear sanf1
    sanf1 = get_station_from_station_name('sanf1');
    sanf1 = read_hires_bathymetry(sanf1,[54e3,35e3],[],true);
    % Default F010 bathymetry resolution for south Florida is 30 m: 3 pts. ~ 90 m
    [ig,ig,ig,sanf1.ngdc_hires_bathy] = find_ngdc_slope(sanf1.ngdc_hires_bathy,[],[],3);
    ig=[]; clear ig
  end;

  patch_lon = sanf1.ngdc_hires_bathy.lon(1:29:end);
  patch_lat = sanf1.ngdc_hires_bathy.lat(1:26:end);
  [LON,LAT]=meshgrid(patch_lon,patch_lat);
  patch_h = interp_field(sanf1.ngdc_hires_bathy.lat,sanf1.ngdc_hires_bathy.lon,sanf1.ngdc_hires_bathy.field,LAT,LON,{@nanmean,26,29,90});
  patch_h = reshape(patch_h,size(LON));
  patch_bet = interp_field(sanf1.ngdc_hires_bathy.lat,sanf1.ngdc_hires_bathy.lon,sanf1.ngdc_hires_bathy.beta,LAT,LON,{@nanmean,26,29,90});
  patch_bet = reshape(patch_bet,size(LON));
  LON=[]; LAT=[]; clear LON LAT
  patch_bet(patch_h>0) = nan;
  %patch_h(patch_h>0) = 0;
  patch_h(patch_h>0) = nan;
  patch_h = abs(patch_h);

  badix = find(isnan(quilted_h));
  [badlatix,badlonix] = ind2sub(size(quilted_h),badix);
  %quilted_h(badix) = interp_field(patch_lat,patch_lon,patch_h,lat(badlatix),lon(badlonix),'cubic');
  quilted_h(badix) = interp_field(patch_lat,patch_lon,patch_h,lat(badlatix),lon(badlonix));
  badix=[]; badlatix=[]; badlonix=[]; clear badix badlatix badlonix

  badix = find(isnan(quilted_bet));
  [badlatix,badlonix] = ind2sub(size(quilted_h),badix);
  quilted_bet(badix) = interp_field(patch_lat,patch_lon,patch_bet,lat(badlatix),lon(badlonix));
  badix=[]; badlatix=[]; badlonix=[]; clear badix badlatix badlonix

  if 0;
    fmg; contourf(patch_lon,patch_lat,-patch_h',-[0:10:120]); colorbar; daspect([1,cosd(lat(1)),1]); titlename('SANF1 Patch');
    axis(ax);
    fmg; contourf(lon,lat,-quilted_h',-[0:10:120]); colorbar; daspect([1,cosd(lat(1)),1]); titlename('SANF1 Quilt');
  end;
  patch_h=[]; patch_lon=[]; patch_lat=[]; clear patch_h patch_lon patch_lat
  sanf1=[]; clear sanf1
end;

if 1;
  disp('Quilting PVGF1 bathymetry into FRT');
  if ( ~exist('pvgf1','var') || ~isfield(pvgf1,'ngdc_hires_bathy') )
    pvgf1=[]; clear pvgf1
    pvgf1 = get_station_from_station_name('pvgf1');
    pvgf1 = read_hires_bathymetry(pvgf1,[25e3,40e3],[],true);
    % Default F010 bathymetry resolution for south Florida is 30 m: 3 pts. ~ 90 m
    [ig,ig,ig,pvgf1.ngdc_hires_bathy] = find_ngdc_slope(pvgf1.ngdc_hires_bathy,[],[],3);
    ig=[]; clear ig
  end;

  patch_lon = pvgf1.ngdc_hires_bathy.lon(1:29:end);
  patch_lat = pvgf1.ngdc_hires_bathy.lat(1:26:end);
  [LON,LAT]=meshgrid(patch_lon,patch_lat);
  patch_h = interp_field(pvgf1.ngdc_hires_bathy.lat,pvgf1.ngdc_hires_bathy.lon,pvgf1.ngdc_hires_bathy.field,LAT,LON,{@nanmean,26,29});
  patch_h = reshape(patch_h,size(LON));
  patch_bet = interp_field(pvgf1.ngdc_hires_bathy.lat,pvgf1.ngdc_hires_bathy.lon,pvgf1.ngdc_hires_bathy.beta,LAT,LON,{@nanmean,26,29,90});
  patch_bet = reshape(patch_bet,size(LON));
  LON=[]; LAT=[]; clear LON LAT
  patch_bet(patch_h>0) = nan;
  %patch_h(patch_h>0) = 0;
  patch_h(patch_h>0) = nan;
  patch_h = abs(patch_h);

  badix = find(isnan(quilted_h));
  [badlatix,badlonix] = ind2sub(size(quilted_h),badix);
  %quilted_h(badix) = interp_field(patch_lat,patch_lon,patch_h,lat(badlatix),lon(badlonix),'cubic');
  quilted_h(badix) = interp_field(patch_lat,patch_lon,patch_h,lat(badlatix),lon(badlonix));
  badix=[]; badlatix=[]; badlonix=[]; clear badix badlatix badlonix

  badix = find(isnan(quilted_bet));
  [badlatix,badlonix] = ind2sub(size(quilted_h),badix);
  quilted_bet(badix) = interp_field(patch_lat,patch_lon,patch_bet,lat(badlatix),lon(badlonix));
  badix=[]; badlatix=[]; badlonix=[]; clear badix badlatix badlonix

  if 0;
    fmg; contourf(patch_lon,patch_lat,-patch_h',-[0:10:120]); colorbar; daspect([1,cosd(lat(1)),1]); titlename('PVGF1 Patch');
    axis(ax);
    fmg; contourf(lon,lat,-quilted_h',-[0:10:120]); colorbar; daspect([1,cosd(lat(1)),1]); titlename('PVGF1 Quilt');
  end;
  patch_h=[]; patch_lon=[]; patch_lat=[]; clear patch_h patch_lon patch_lat
  pvgf1=[]; clear pvgf1
end;

if 0;
  disp('Plotting Original H');
  fmg; contourf(lon,lat,-h,-[0:10:200]); colorbar; daspect([1,cosd(lat(1)),1]); titlename('Original');
  ax=axis;
  disp('Plotting Quilted H result');
  fmg; contourf(lon,lat,-quilted_h,-[0:10:200]); colorbar; daspect([1,cosd(lat(1)),1]); titlename('Quilted');
  %axis(ax);
end;
if 1;
  disp('Plotting Original BET');
  fmg; contourf(lon,lat,bet,[0:0.005:0.10]); colorbar; daspect([1,cosd(lat(1)),1]); titlename('Original BET');
  ax=axis;
  disp('Plotting Quilted BET result');
  fmg; contourf(lon,lat,quilted_bet,[0:0.005:0.10]); colorbar; daspect([1,cosd(lat(1)),1]); titlename('Quilted BET');
  %axis(ax);
end;


if 1;
  disp('Saving FRT_complete_depth_and_beta.mat');
  h = quilted_h;
  bet = quilted_bet;
  quilted_h=[]; quilted_bet=[]; clear quilted_h quilted_bet
  save('FRT_complete_depth_and_beta.mat','lon','lat','h','bet');
end;

% >> field_bbox({lon,lat})
% ans =
%   -82.1389  -79.9876   24.3906   26.2517
