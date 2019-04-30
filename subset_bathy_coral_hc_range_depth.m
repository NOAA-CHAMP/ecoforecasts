1;
%% SCRIPT subset_bathy_coral_hc_range_depth.m
%
% Subset contents of several .MAT files to just the areas which correspond to
% coral reef polygons (REEF) and coral-and-hard-bottom combined polygons
% (RFHD) loaded with READ_CORAL_AND_HARD_BOTTOM (v.)
%
% Last Saved Time-stamp: <Mon 2017-05-08 09:29:06 Eastern Daylight Time gramer>

set_more off;

if ( ~exist('rfhd','var') )
  read_coral_and_hard_bottom;
end;

checkin_count = 6000;
%checkin_count = +inf;

if ( ~exist('FRT_depth_and_beta_92m_NANMEAN_3_3_4_hc_range_depth_coral_and_hard_bottom.mat','file') )
disp('NANMEAN - RfHd - Hit Enter to continue...'); pause;
tic,
  h=[]; bet=[]; LON=[]; LAT=[]; clear lon lat h bet LON LAT
  load('FRT_depth_and_beta_92m_NANMEAN_3_3_4_hc_range_depth.mat');
  [LON,LAT] = meshgrid(lon,lat);
  inix=[];
  disp(numel(shps));
  dlat = 1.1*abs(lat(1)-lat(2));

  w = warning('OFF','MATLAB:triangulation:EmptyTri2DWarnId'); % We don't care about these
  for shpix=1:numel(shps)
    if (mod(shpix,checkin_count) == 0); disp(shpix); toc; tic; end;
    [lonix,latix] = fieldfind(lon,lat,shps(shpix).X(1:end-1),shps(shpix).Y(1:end-1));
    goodix = find(~isnan(lonix) & ~isnan(latix));
    lonix = lonix(goodix);    latix = latix(goodix);
    lix = unique([lonix;latix]','rows')';
    if ( ~isempty(lix) )
      lonix = lix(1,:);
      latix = lix(2,:);
      if ( numel(lonix) < 3 || shps(shpix).Shapearea < 2e5 )
        ix = sub2ind(size(LON),latix,lonix);
        inix = union(inix,ix);
      else
        tr = delaunayTriangulation(lonix',latix');
        %DEBUG:        [w,d]=lastwarn; if strcmp(d,'MATLAB:delaunayTriangulation:DupPtsWarnId'); keyboard; lastwarn(''); end;
        ix = sub2ind(size(LON),tr.Points(:,2),tr.Points(:,1));
        ix = find(inpolygon(LON,LAT,LON(ix),LAT(ix))==1);
        %DEBUG:        [w,d]=lastwarn; if length(w)>0; keyboard; lastwarn(''); end;
        inix = union(inix,ix);
      end; %if ( numel(lonix) < 3 || shps(shpix).Shapearea < 2e5 )
    end; %if ( ~isempty(lix) )
  end; %for shpix=1:numel(shps)
  warning(w);

  newht = repmat(nan,size(LON));
  ht=h'; newht(inix)=ht(inix); ht=[]; clear ht
  h=[]; h=newht'; newht=[]; clear newht

  newbett = repmat(nan,size(LON));
  bett=bet'; newbett(inix)=bett(inix); bett=[]; clear bett
  bet=[]; bet=newbett'; newbett=[]; clear newbett

  newangt = repmat(nan,size(LON));
  angt=ang'; newangt(inix)=angt(inix); angt=[]; clear angt
  ang=[]; ang=newangt'; newangt=[]; clear newangt

  newhcrngt = repmat(nan,size(LON));
  hcrngt=hcrng'; newhcrngt(inix)=hcrngt(inix); hcrngt=[]; clear hcrngt
  hcrng=[]; hcrng=newhcrngt'; newhcrngt=[]; clear newhcrngt

  newhcht = repmat(nan,size(LON));
  hcht=hch'; newhcht(inix)=hcht(inix); hcht=[]; clear hcht
  hch=[]; hch=newhcht'; newhcht=[]; clear newhcht

  disp(numel(find(~isnan(h))));

  save('FRT_depth_and_beta_92m_NANMEAN_3_3_4_hc_range_depth_coral_and_hard_bottom.mat',...
       'lon','lat','h','bet','ang','hcrng','hch');
  clear inix ix shpix;
toc,
end;

if ( ~exist('FRT_depth_and_beta_92m_NANMEAN_3_3_4_hc_range_depth_coral.mat','file') )
disp('NANMEAN - Reef - Hit Enter to continue...'); pause;
tic,
  h=[]; bet=[]; LON=[]; LAT=[]; clear lon lat h bet LON LAT
  load('FRT_depth_and_beta_92m_NANMEAN_3_3_4_hc_range_depth.mat');
  [LON,LAT] = meshgrid(lon,lat);
  inix=[];
  disp(numel(shps));
  dlat = 1.1*abs(lat(1)-lat(2));

  w = warning('OFF','MATLAB:triangulation:EmptyTri2DWarnId'); % We don't care about these
  for shpix=1:numel(shps)
    if (mod(shpix,checkin_count) == 0); disp(shpix); toc; tic; end;
    if ( strcmp(shps(shpix).DESCRIPT,'Coral Reef') )
      [lonix,latix] = fieldfind(lon,lat,shps(shpix).X(1:end-1),shps(shpix).Y(1:end-1));
      goodix = find(~isnan(lonix) & ~isnan(latix));
      lonix = lonix(goodix);
      latix = latix(goodix);
      lix = unique([lonix;latix]','rows')';
      if ( ~isempty(lix) )
        lonix = lix(1,:);
        latix = lix(2,:);
        if ( numel(lonix) < 3 || shps(shpix).Shapearea < 2e5 )
          ix = sub2ind(size(LON),latix,lonix);
          inix = union(inix,ix);
        else
          tr = delaunayTriangulation(lonix',latix');
          %DEBUG:        [w,d]=lastwarn; if strcmp(d,'MATLAB:delaunayTriangulation:DupPtsWarnId'); keyboard; lastwarn(''); end;
          ix = sub2ind(size(LON),tr.Points(:,2),tr.Points(:,1));
          ix = find(inpolygon(LON,LAT,LON(ix),LAT(ix))==1);
          %DEBUG:        [w,d]=lastwarn; if length(w)>0; keyboard; lastwarn(''); end;
          inix = union(inix,ix);
        end; %if ( numel(lonix) < 3 || shps(shpix).Shapearea < 2e5 )
      end; %if ( ~isempty(lix) )
    end; %if ( strcmp(shps(shpix).DESCRIPT,'Coral Reef') )
  end; %for shpix=1:numel(shps)
  warning(w);

  newht = repmat(nan,size(LON));
  ht=h'; newht(inix)=ht(inix); ht=[]; clear ht
  h=[]; h=newht'; newht=[]; clear newht

  newbett = repmat(nan,size(LON));
  bett=bet'; newbett(inix)=bett(inix); bett=[]; clear bett
  bet=[]; bet=newbett'; newbett=[]; clear newbett

  newangt = repmat(nan,size(LON));
  angt=ang'; newangt(inix)=angt(inix); angt=[]; clear angt
  ang=[]; ang=newangt'; newangt=[]; clear newangt

  newhcrngt = repmat(nan,size(LON));
  hcrngt=hcrng'; newhcrngt(inix)=hcrngt(inix); hcrngt=[]; clear hcrngt
  hcrng=[]; hcrng=newhcrngt'; newhcrngt=[]; clear newhcrngt

  newhcht = repmat(nan,size(LON));
  hcht=hch'; newhcht(inix)=hcht(inix); hcht=[]; clear hcht
  hch=[]; hch=newhcht'; newhcht=[]; clear newhcht

  disp(numel(find(~isnan(h))));
  save('FRT_depth_and_beta_92m_NANMEAN_3_3_4_hc_range_depth_coral.mat',...
       'lon','lat','h','bet','ang','hcrng','hch');
  clear inix ix shpix;
toc,
end;



if ( ~exist('FRT_depth_and_beta_92m_hc_range_depth_coral_and_hard_bottom.mat','file') )
disp('Full - RfHd - Hit Enter to continue...'); pause;
tic,
  h=[]; bet=[]; LON=[]; LAT=[]; clear lon lat h bet LON LAT
  load('FRT_depth_and_beta_92m_hc_range_depth.mat');
  [LON,LAT] = meshgrid(lon,lat);
  inix=[];
  disp(numel(shps));
  dlat = 1.1*abs(lat(1)-lat(2));

  w = warning('OFF','MATLAB:triangulation:EmptyTri2DWarnId'); % We don't care about these
  for shpix=1:numel(shps)
    if (mod(shpix,checkin_count) == 0); disp(shpix); toc; tic; end;
    [lonix,latix] = fieldfind(lon,lat,shps(shpix).X(1:end-1),shps(shpix).Y(1:end-1));
    goodix = find(~isnan(lonix) & ~isnan(latix));
    lonix = lonix(goodix);
    latix = latix(goodix);
    lix = unique([lonix;latix]','rows')';
    if ( ~isempty(lix) )
      lonix = lix(1,:);
      latix = lix(2,:);
      if ( numel(lonix) < 3 || shps(shpix).Shapearea < 2e5 )
        ix = sub2ind(size(LON),latix,lonix);
        inix = union(inix,ix);
      else
        tr = delaunayTriangulation(lonix',latix');
        %DEBUG:        [w,d]=lastwarn; if strcmp(d,'MATLAB:delaunayTriangulation:DupPtsWarnId'); keyboard; lastwarn(''); end;
        ix = sub2ind(size(LON),tr.Points(:,2),tr.Points(:,1));
        ix = find(inpolygon(LON,LAT,LON(ix),LAT(ix))==1);
        %DEBUG:        [w,d]=lastwarn; if length(w)>0; keyboard; lastwarn(''); end;
        inix = union(inix,ix);
      end; %if ( numel(lonix) < 3 || shps(shpix).Shapearea < 2e5 )
    end; %if ( ~isempty(lix) )
  end; %for shpix=1:numel(shps)
  warning(w);

  newht = repmat(nan,size(LON));
  ht=h; newht(inix)=ht(inix); ht=[]; clear ht
  h=[]; h=newht; newht=[]; clear newht

  newbett = repmat(nan,size(LON));
  bett=bet; newbett(inix)=bett(inix); bett=[]; clear bett
  bet=[]; bet=newbett; newbett=[]; clear newbett

  newangt = repmat(nan,size(LON));
  angt=ang; newangt(inix)=angt(inix); angt=[]; clear angt
  ang=[]; ang=newangt; newangt=[]; clear newangt

  newhcrngt = repmat(nan,size(LON));
  hcrngt=hcrng; newhcrngt(inix)=hcrngt(inix); hcrngt=[]; clear hcrngt
  hcrng=[]; hcrng=newhcrngt; newhcrngt=[]; clear newhcrngt

  newhcht = repmat(nan,size(LON));
  hcht=hch; newhcht(inix)=hcht(inix); hcht=[]; clear hcht
  hch=[]; hch=newhcht; newhcht=[]; clear newhcht

  disp(numel(find(~isnan(h))));
  save('FRT_depth_and_beta_92m_hc_range_depth_coral_and_hard_bottom.mat',...
       'lon','lat','h','bet','ang','hcrng','hch');
  clear inix ix shpix;
toc,
end;

if ( ~exist('FRT_depth_and_beta_92m_hc_range_depth_coral.mat','file') )
disp('Full - Reef - Hit Enter to continue...'); pause;
tic,
  h=[]; bet=[]; LON=[]; LAT=[]; clear lon lat h bet LON LAT
  load('FRT_depth_and_beta_92m_hc_range_depth.mat');
  [LON,LAT] = meshgrid(lon,lat);
  inix=[];
  disp(numel(shps));
  dlat = 1.1*abs(lat(1)-lat(2));

  w = warning('OFF','MATLAB:triangulation:EmptyTri2DWarnId'); % We don't care about these
  for shpix=1:numel(shps)
    if (mod(shpix,checkin_count) == 0); disp(shpix); toc; tic; end;
    if ( strcmp(shps(shpix).DESCRIPT,'Coral Reef') )
      [lonix,latix] = fieldfind(lon,lat,shps(shpix).X(1:end-1),shps(shpix).Y(1:end-1));
      goodix = find(~isnan(lonix) & ~isnan(latix));
      lonix = lonix(goodix);
      latix = latix(goodix);
      lix = unique([lonix;latix]','rows')';
      if ( ~isempty(lix) )
        lonix = lix(1,:);
        latix = lix(2,:);
        if ( numel(lonix) < 3 || shps(shpix).Shapearea < 2e5 )
          ix = sub2ind(size(LON),latix,lonix);
          inix = union(inix,ix);
        else
          tr = delaunayTriangulation(lonix',latix');
          %DEBUG:        [w,d]=lastwarn; if strcmp(d,'MATLAB:delaunayTriangulation:DupPtsWarnId'); keyboard; lastwarn(''); end;
          ix = sub2ind(size(LON),tr.Points(:,2),tr.Points(:,1));
          ix = find(inpolygon(LON,LAT,LON(ix),LAT(ix))==1);
          %DEBUG:        [w,d]=lastwarn; if length(w)>0; keyboard; lastwarn(''); end;
          inix = union(inix,ix);
        end; %if ( numel(lonix) < 3 || shps(shpix).Shapearea < 2e5 )
      end; %if ( ~isempty(lix) )
    end; %if ( strcmp(shps(shpix).DESCRIPT,'Coral Reef') )
  end; %for shpix=1:numel(shps)
  warning(w);

  newht = repmat(nan,size(LON));
  ht=h; newht(inix)=ht(inix); ht=[]; clear ht
  h=[]; h=newht; newht=[]; clear newht

  newbett = repmat(nan,size(LON));
  bett=bet; newbett(inix)=bett(inix); bett=[]; clear bett
  bet=[]; bet=newbett; newbett=[]; clear newbett

  newangt = repmat(nan,size(LON));
  angt=ang; newangt(inix)=angt(inix); angt=[]; clear angt
  ang=[]; ang=newangt; newangt=[]; clear newangt

  newhcrngt = repmat(nan,size(LON));
  hcrngt=hcrng; newhcrngt(inix)=hcrngt(inix); hcrngt=[]; clear hcrngt
  hcrng=[]; hcrng=newhcrngt; newhcrngt=[]; clear newhcrngt

  newhcht = repmat(nan,size(LON));
  hcht=hch; newhcht(inix)=hcht(inix); hcht=[]; clear hcht
  hch=[]; hch=newhcht; newhcht=[]; clear newhcht

  disp(numel(find(~isnan(h))));
  save('FRT_depth_and_beta_92m_hc_range_depth_coral.mat',...
       'lon','lat','h','bet','ang','hcrng','hch');
  clear inix ix shpix;
toc,
end;

set_more;
