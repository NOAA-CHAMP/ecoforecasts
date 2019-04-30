function [lon,lat,dat,matfname,utmzone] = asc2mat_tithe(dataset,coastpath,matpath,reconvert,utmzone)
%function [lon,lat,dat,matfname,utmzone] = asc2mat_tithe(dataset,coastpath,matpath,reconvert,utmzone)
%
% **HACK**: This version of ASC2MAT breaks a massive file up into 10 parts.
%
% Convert data in an ASC file (common ArcGIS ASCII format with 6-line header)
% into a MAT file. Return N-vector of LONgitudes, M-vector of LATitudes, and
% MxN matrix of DATa. Input file is FULLFILE(COASTPATH,DATASET,'.asc') and
% output file FULLFILE(MATPATH,DATASET,'.mat'). By default, COASTPATH is the
% 'coast' subdirectory of Ecoforecasts toolbox, and MATPATH==COASTPATH.
%
% NOTE: If MAT file exists, simply load it, unless RECONVERT is True. 
%
% If no MAT file exists yet (or RECONVERT is True) and optional UTMZONE is
% given, ASC file coordinates are assumed to be in UTM (Easting and Northing)
% rather than lat/lon; pass UTMZONE as a third arg to UTM2LATLON (v.), to
% convert all coordinates before saving to MAT file. Zone should be "17R" for
% areas of Florida. Other zones of interest for coral reefs: 55P Saipan/CNMI,
% 2L Anerican Samoa, 19Q for western Puerto Rico and Dominican Republic, and
% 20Q for USVI. 
%
% Last Saved Time-stamp: <Sat 2016-06-04 16:44:47 Eastern Daylight Time gramer>

  %DEBUG:  tic,

  if ( ~exist('coastpath','var') || isempty(coastpath) )
    coastpath = get_ecoforecasts_path('coast');
  end;
  if ( ~exist('matpath','var') || isempty(matpath) )
    matpath = coastpath;
  end;
  if ( ~exist('reconvert','var') || isempty(reconvert) )
    reconvert = false;
  end;
  if ( ~exist('utmzone','var') || isempty(utmzone) )
    utmzone = [];
  end;

  ascfname = fullfile(coastpath,[dataset,'.asc']);
  matfname = fullfile(coastpath,[dataset,'.mat']);

  if ( ~reconvert && exist(matfname,'file') )
    disp(['Loading ',matfname]);
    load(matfname);

  else
    fid = fopen(ascfname,'r');
    if ( fid < 0 )
      error('Unable to open %s',ascfname);
    end;
    A = fscanf(fid,'%*s %g\n',6);
    ncols = A(1);
    nrows = A(2);
    %dat = fscanf(fid,'%g',[nrows,ncols]);
    dat = fscanf(fid,'%f',[ncols,nrows]);
    % Columns first and latitude in reverse - ouch my dyslexia!
    dat = dat(:,end:-1:1)';
    fclose(fid);

    xllc = A(3);
    yllc = A(4);
    cellsz = A(5);
    badval = A(6);

    lon = xllc + ([0:ncols-1].*cellsz);
    lat = yllc + ([0:nrows-1].*cellsz);

    dat(dat==badval) = nan;

    if ( ~isempty(utmzone) )
      x = lon;
      y = lat;
      ws = warning('OFF','Ecoforecasts:UTMgrid');
      [LAT,LON] = utm2latlon(x,y,utmzone);
      warning(ws);
      % Linearly interpolate back onto a grid
      lon = linspace(min(LON(:)),max(LON(:)),numel(unique(x)));
      lat = linspace(min(LAT(:)),max(LAT(:)),numel(unique(y)));
      dat = griddata(LON,LAT,dat,lon,lat');
      LAT=[]; LON=[]; clear LAT LON
    end;

    clear A fid
    disp(['Saving ',matfname]);
    save(matfname);

  end;

  %DEBUG:  toc,

return;
