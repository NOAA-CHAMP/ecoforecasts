1;

set_more off;

if ( ~exist('rfhd','var') )
  read_coral_and_hard_bottom;
end;

checkin_count = 6000;
%checkin_count = +inf;

if ( ~exist('FRT_depth_and_beta_92m_NANMEAN_3_3_4_coral_and_hard_bottom.mat','file') )
disp('NANMEAN - RfHd - Hit Enter to continue...'); pause;
tic,
  h=[]; bet=[]; LON=[]; LAT=[]; clear lon lat h bet LON LAT
  load('FRT_depth_and_beta_92m_NANMEAN_3_3_4.mat');
  [LON,LAT] = meshgrid(lon,lat);
  inix=[];
  disp(numel(shps));
  dlat = 1.1*abs(lat(1)-lat(2));
  for shpix=1:numel(shps)
    if (mod(shpix,checkin_count) == 0); disp(shpix); toc; tic; end;
    [lonix,latix] = fieldfind(lon,lat,shps(shpix).X(1:end-1),shps(shpix).Y(1:end-1));
    goodix = find(~isnan(lonix) & ~isnan(latix));
    lonix = lonix(goodix);
    latix = latix(goodix);
    if ( ~isempty(lonix) && ~isempty(latix) )
      tr = delaunayTriangulation(lonix',latix');
      ix = sub2ind(size(LON),tr.Points(:,2),tr.Points(:,1));
      ix = find(inpolygon(LON,LAT,LON(ix),LAT(ix))==1);
      inix = union(inix,ix);
    end;
  end;
  newht = repmat(nan,size(LON));
  ht=h'; newht(inix)=ht(inix); ht=[]; clear ht
  h=[]; h=newht'; newht=[]; clear newht
  newbett = repmat(nan,size(LON));
  bett=bet'; newbett(inix)=bett(inix); bett=[]; clear bett
  bet=[]; bet=newbett'; newbett=[]; clear newbett
  disp(numel(find(~isnan(h))));
  save('FRT_depth_and_beta_92m_NANMEAN_3_3_4_coral_and_hard_bottom.mat','lon','lat','h','bet');
  clear inix ix shpix;
toc,
end;

if ( ~exist('FRT_depth_and_beta_92m_NANMEAN_3_3_4_coral.mat','file') )
disp('NANMEAN - Reef - Hit Enter to continue...'); pause;
tic,
  h=[]; bet=[]; LON=[]; LAT=[]; clear lon lat h bet LON LAT
  load('FRT_depth_and_beta_92m_NANMEAN_3_3_4.mat');
  [LON,LAT] = meshgrid(lon,lat);
  inix=[];
  disp(numel(shps));
  dlat = 1.1*abs(lat(1)-lat(2));
  for shpix=1:numel(shps)
    if (mod(shpix,checkin_count) == 0); disp(shpix); toc; tic; end;
    if ( strcmp(shps(shpix).DESCRIPT,'Coral Reef') )
      [lonix,latix] = fieldfind(lon,lat,shps(shpix).X(1:end-1),shps(shpix).Y(1:end-1));
      goodix = find(~isnan(lonix) & ~isnan(latix));
      lonix = lonix(goodix);
      latix = latix(goodix);
      if ( ~isempty(lonix) && ~isempty(latix) )
        tr = delaunayTriangulation(lonix',latix');
        ix = sub2ind(size(LON),tr.Points(:,2),tr.Points(:,1));
        ix = find(inpolygon(LON,LAT,LON(ix),LAT(ix))==1);
        inix = union(inix,ix);
      end;
    end;
  end;
  newht = repmat(nan,size(LON));
  ht=h'; newht(inix)=ht(inix); ht=[]; clear ht
  h=[]; h=newht'; newht=[]; clear newht
  newbett = repmat(nan,size(LON));
  bett=bet'; newbett(inix)=bett(inix); bett=[]; clear bett
  bet=[]; bet=newbett'; newbett=[]; clear newbett
  disp(numel(find(~isnan(h))));
  save('FRT_depth_and_beta_92m_NANMEAN_3_3_4_coral.mat','lon','lat','h','bet');
  clear inix ix shpix;
toc,
end;



if ( ~exist('FRT_depth_and_beta_92m_coral_and_hard_bottom.mat','file') )
disp('Full - RfHd - Hit Enter to continue...'); pause;
tic,
  h=[]; bet=[]; LON=[]; LAT=[]; clear lon lat h bet LON LAT
  load('FRT_depth_and_beta_92m.mat');
  [LON,LAT] = meshgrid(lon,lat);
  inix=[];
  disp(numel(shps));
  dlat = 1.1*abs(lat(1)-lat(2));
  for shpix=1:numel(shps)
    if (mod(shpix,checkin_count) == 0); disp(shpix); toc; tic; end;
    [lonix,latix] = fieldfind(lon,lat,shps(shpix).X(1:end-1),shps(shpix).Y(1:end-1));
    goodix = find(~isnan(lonix) & ~isnan(latix));
    lonix = lonix(goodix);
    latix = latix(goodix);
    if ( ~isempty(lonix) && ~isempty(latix) )
      tr = delaunayTriangulation(lonix',latix');
      ix = sub2ind(size(LON),tr.Points(:,2),tr.Points(:,1));
      ix = find(inpolygon(LON,LAT,LON(ix),LAT(ix))==1);
      inix = union(inix,ix);
    end;
  end;
  newht = repmat(nan,size(LON));
  ht=h; newht(inix)=ht(inix); ht=[]; clear ht
  h=[]; h=newht; newht=[]; clear newht
  newbett = repmat(nan,size(LON));
  bett=bet; newbett(inix)=bett(inix); bett=[]; clear bett
  bet=[]; bet=newbett; newbett=[]; clear newbett
  disp(numel(find(~isnan(h))));
  save('FRT_depth_and_beta_92m_coral_and_hard_bottom.mat','lon','lat','h','bet');
  clear inix ix shpix;
toc,
end;

if ( ~exist('FRT_depth_and_beta_92m_coral.mat','file') )
disp('Full - Reef - Hit Enter to continue...'); pause;
tic,
  h=[]; bet=[]; LON=[]; LAT=[]; clear lon lat h bet LON LAT
  load('FRT_depth_and_beta_92m.mat');
  [LON,LAT] = meshgrid(lon,lat);
  inix=[];
  disp(numel(shps));
  dlat = 1.1*abs(lat(1)-lat(2));
  for shpix=1:numel(shps)
    if (mod(shpix,checkin_count) == 0); disp(shpix); toc; tic; end;
    if ( strcmp(shps(shpix).DESCRIPT,'Coral Reef') )
      [lonix,latix] = fieldfind(lon,lat,shps(shpix).X(1:end-1),shps(shpix).Y(1:end-1));
      goodix = find(~isnan(lonix) & ~isnan(latix));
      lonix = lonix(goodix);
      latix = latix(goodix);
      if ( ~isempty(lonix) && ~isempty(latix) )
        tr = delaunayTriangulation(lonix',latix');
        ix = sub2ind(size(LON),tr.Points(:,2),tr.Points(:,1));
        ix = find(inpolygon(LON,LAT,LON(ix),LAT(ix))==1);
        inix = union(inix,ix);
      end;
    end;
  end;
  newht = repmat(nan,size(LON));
  ht=h; newht(inix)=ht(inix); ht=[]; clear ht
  h=[]; h=newht; newht=[]; clear newht
  newbett = repmat(nan,size(LON));
  bett=bet; newbett(inix)=bett(inix); bett=[]; clear bett
  bet=[]; bet=newbett; newbett=[]; clear newbett
  disp(numel(find(~isnan(h))));
  save('FRT_depth_and_beta_92m_coral.mat','lon','lat','h','bet');
  clear inix ix shpix;
toc,
end;

set_more;
