1;
%% SCRIPT process_F010.m
% Process "F010" USGS DEM high-res. (20 m, UTM grid) coastal bathymetry file
% for south  Florida: loads the ASCII file using the confusiingly-named
% USGS24KDEM (MATLAB Map Toolbox), interpolates onto a latitude-longitude
% grid, and saves the resulting variables (lon,lat,dat) in a .MAT file.
%
% Data is documented in file "F010_25081C1_BIG.dem.IS.FLORIDA.30m.bathy.txt"
%
% Data source URL:
%  http://estuarinebathymetry.noaa.gov/bathy_htmls/F010.html
% USGS DEM file format reference URL:
%  https://en.wikipedia.org/wiki/USGS_DEM
%
% Last Saved Time-stamp: <Mon 2016-07-11 11:49:37 Eastern Daylight Time gramer>

set_more off

coastpath = get_ecoforecasts_path('coast');
fname = 'F010_25081C1_BIG.dem';

if ( ~exist('LAT','var') ||~exist('LON','var') ||~exist('DAT','var') )
  DAT=[]; LAT=[]; LON=[]; clear DAT LAT LON
  intermed_matfname = fullfile(coastpath,[fname,'-intermed.mat']);
  if ( exist(intermed_matfname,'file') )
    disp(['Loading ',intermed_matfname]);
    load(intermed_matfname);
  else
    fpath = fullfile(coastpath,fname);
    disp(['USGS24KDEM(',fpath,',1)']);
    [LAT,LON,DAT] = usgs24kdem(fpath,1);
    disp(['Saving ',intermed_matfname]);
    save(intermed_matfname,'LAT','LON','DAT','fpath','-v7.3');
  end;

  % Data is loaded by USGS24KDEM with highest latitude in first row: reverse it
  LAT = LAT(end:-1:1,:);
  LON = LON(end:-1:1,:);
  DAT = DAT(end:-1:1,:);
end;

% Linearly interpolate back onto a latitude-longitude grid
lon = linspace(min(LON(:)),max(LON(:)),size(DAT,2));
lat = linspace(min(LAT(:)),max(LAT(:)),size(DAT,1));

[MLON,MLAT] = meshgrid(lon,lat);

dat=[]; clear dat
dat = repmat(nan,size(DAT));


% Grid is too big to call GRIDDATA all at once: parcel up into more or less
% equal-sized tiles, with a significant-enough overlap between them

roverlap = floor(size(DAT,1)/10);
coverlap = floor(size(DAT,2)/10);

rstrd = ceil(size(DAT,1)/4);
cstrd = ceil(size(DAT,2)/4);
begrixen = 1:rstrd:size(DAT,1);
begcixen = 1:cstrd:size(DAT,2);

for begrixix=1:numel(begrixen);
  fprintf(1,'%g of %g rows\n',begrixix,numel(begrixen)); 
  tic,

  begrix = begrixen(begrixix);
  rixen = begrix-roverlap:begrix+rstrd-1+roverlap;
  rixen = rixen(1 <= rixen & rixen <= size(DAT,1));
  nestrix = begrix:begrix+rstrd-1;
  nestrix = nestrix(1 <= nestrix & nestrix <= size(DAT,1));
  gridrix = find(ismember(rixen,nestrix));

  for begcixix=1:numel(begcixen); 
    begcix=begcixen(begcixix);
    cixen = begcix-coverlap:begcix+cstrd-1+coverlap;
    cixen = cixen(1 <= cixen & cixen <= size(DAT,2));
    nestcix = begcix:begcix+cstrd-1;
    nestcix = nestcix(1 <= nestcix & nestcix <= size(DAT,2));
    gridcix = find(ismember(cixen,nestcix));

    disp({min(nestrix),max(nestrix),min(nestcix),max(nestcix)});
    d = griddata(LON(rixen,cixen),LAT(rixen,cixen),DAT(rixen,cixen),MLON(rixen,cixen),MLAT(rixen,cixen),'nearest');
    dat(nestrix,nestcix) = d(gridrix,gridcix);
    d=[]; clear d
  end; 
  clear begcixix begcix cixen

  toc,
end;
clear begrixix begrix rixen

MLAT=[]; MLON=[]; clear MLAT MLON

if 1;
  matfname = fullfile(coastpath,[fname,'.mat']);
  disp(['Saving ',matfname]);
  save(matfname,'lat','lon','dat','fname','-v7.3');
end;

if 1;
  all_lat = lat;
  all_lon = lon;
  all_dat = dat;
  lat = []; lon = []; dat = []; clear lat lon dat

  latquadir = 'SN';
  lonquadir = 'WE';
  latix = [1,ceil(numel(all_lat)/2),floor(numel(all_lat)/2),numel(all_lat)];
  lonix = [1,ceil(numel(all_lon)/2),floor(numel(all_lon)/2),numel(all_lon)];
  for dlat=1:2
    lat = all_lat(latix((dlat*2)-1):latix((dlat*2)));
    for dlon=1:2
      lon = all_lon(lonix((dlon*2)-1):lonix((dlon*2)));
      dat = all_dat(latix((dlat*2)-1):latix((dlat*2)), lonix((dlon*2)-1):lonix((dlon*2)));

      quadir = [latquadir(dlat),lonquadir(dlon)];
      matfbasename = sprintf('%s_%s.mat',fname,quadir);
      matfname = fullfile(coastpath,matfbasename);
      disp(['Saving ',matfname]);
      save(matfname,'lat','lon','dat','fname','-v7.3');
      disp({matfbasename,min(lon),max(lon),min(lat),max(lat)});
      lon = []; dat = []; clear lon dat
    end;
    lat = []; clear lat
  end;

  lat = all_lat;
  lon = all_lon;
  dat = all_dat;
  all_lat = []; all_lon = []; all_dat = []; clear all_lat all_lon all_dat

  clear ans dlat dlon latix lonix latquadir lonquadir matfbasename matfname quadir
end;

set_more

%dat=[]; clear dat; clear ans begc* begr* coastpath coverlap roverlap cstrd rstrd fname fpath gridcix gridrix intermed_matfname lat lon latquadir lonquadir matfname nestcix nestrix quadir
