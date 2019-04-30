function stn = get_avhrr_field(stn_or_stnm,calcFieldTerms,region,interpMethod,minPct,xradius,yradius,doUpdate)
%function stn = get_avhrr_field(stn_or_stnm,calcFieldTerms,region,interpMethod,minPct,xradius,yradius,doUpdate)
%
% Query synoptic AVHRR SST imagery from USF online archive: Download images
% in USF area REGION (DEFAULT: 'florida'); subset each image to a 19x17-pixel
% box surrounding the station in or named by STN_OR_STNM (a STRUCT or CHAR).
%
% Returns STN struct with field .avhrr_sst_field, containing the subfields
% .lon, .lat, .date, .field. If optional CALCFIELDTERMS==True, calculate
% gradients and field Laplacian of SST field, .gradient_x, .gradient_y,
% .laplacian. Interpolates fields using INTERPMETHOD (DEFAULT: 'linear') to
% create time series .avhrr_sst, and .avhrr_sst_x,.avhrr_sst_y,.avhrr_sst_l.
%
% NOTE: For hourly data (see CALC_FIELD_TERMS), time series fields for _x,
% _y, _l (and _[dxy][lxy]) may each contain different numbers of points due
% to clouds/land. You may correct for this by calling INTERSECT_TSES.
%
% NOTE: Only returns SST fields where at least MINPCT (DEFAULT: 0.50) of the
% pixels in the field are non-NaN (not cloud, land mask, or error).
%
% Optional XRADIUS, YRADIUS determine the pixel size of the extracted SST
% fields (DEFAULT is as indicated above: XRADIUS==8, YRADIUS==9).
%
% A MAT file is created after the first call for a given station and x- and
% y-radii. If optional DOUPDATE (DEFAULT: false), the existing MAT file is
% updated starting from the last month for which valid data is found.
%
% Last Saved Time-stamp: <Fri 2012-05-11 17:27:05  Lew.Gramer>

  set_more off;

  %DEBUG:
  tic,

  datapath = get_ecoforecasts_path('data');
  avhrrpath = fullfile(datapath,'avhrr');

  stn = get_station_from_station_name(stn_or_stnm);

  if ( ~exist('calcFieldTerms','var') || isempty(calcFieldTerms) )
    calcFieldTerms = false;
  end;
  if ( ~exist('region','var') || isempty(region) )
    region = 'florida';
  end;
  if ( ~exist('interpMethod','var') || isempty(interpMethod) )
    interpMethod = 'linear';
  end;
  if ( ~exist('minPct','var') || isempty(minPct) )
    minPct = 0.50;
  end;
  if ( ~exist('xradius','var') || isempty(xradius) )
    xradius = 8;
  end;
  if ( ~exist('yradius','var') || isempty(yradius) )
    yradius = 9;
  end;
  if ( ~exist('doUpdate','var') || isempty(doUpdate) )
    doUpdate = false;
  end;

  % If user wants CALC_FIELD_TERMS, how many points to use in TEMPLATE
  npts = 5;

  switch ( region )
   case 'florida',
    minlon = -91; maxlon = -79; dlon = (1/110);
    minlat =  22; maxlat =  31; dlat = (1/110);
   otherwise,
    error('Do not yet know how to extract region "%s"!', region);
  end;

  rgnlon = [minlon:dlon:maxlon]';
  rgnlat = [maxlat:-dlat:minlat]';

  matfname = fullfile(datapath,sprintf('avhrr_sst_%s_%02dx%02d.mat',...
                                       stn.station_name,...
                                       ((yradius*2)+1),((xradius*2)+1)));


  if ( exist(matfname,'file') )
    disp(['Reloading ' matfname]);
    load(matfname,'station');
  end;

  if ( ~exist(matfname,'file') || doUpdate )

    disp('Extracting from raw PNG files...');

    % Start of USF AVHRR coverage with 'n11.19930825.2209.florida.png'
    zerodt = datenum(1993,8,25);

    if ( doUpdate && exist('station','var') && isfield(station,'avhrr_sst_field') )
      lastdt = station.avhrr_sst_field.date(end);
      begdt = datenum(get_year(lastdt),get_month(lastdt),1);
      disp(['(Updating from ',datestr(begdt),' onward)']);
      alldts = begdt:now;

      overlapix = find(station.avhrr_sst_field.date>=begdt);
      station.avhrr_sst_field.date(overlapix) = [];
      station.avhrr_sst_field.field(overlapix,:,:) = [];

      starting_index = length(station.avhrr_sst_field.date);

    else
      alldts = zerodt:now;

      starting_index = 0;

      station = get_station_from_station_name(stn.station_name);
      station.avhrr_sst_midx = interp1(rgnlon,1:length(rgnlon),stn.lon,'nearest');
      station.avhrr_sst_midy = interp1(rgnlat,1:length(rgnlat),stn.lat,'nearest');

      station.avhrr_sst_xix = [station.avhrr_sst_midx-xradius:station.avhrr_sst_midx+xradius]';
      station.avhrr_sst_xix(~ismember(station.avhrr_sst_xix,1:length(rgnlon)))=[];
      station.avhrr_sst_yix = [station.avhrr_sst_midy-yradius:station.avhrr_sst_midy+yradius]';
      station.avhrr_sst_yix(~ismember(station.avhrr_sst_yix,1:length(rgnlat)))=[];

      station.avhrr_sst_field.lon = rgnlon(station.avhrr_sst_xix);
      station.avhrr_sst_field.lat = rgnlat(station.avhrr_sst_yix);
      station.avhrr_sst_field.date = [];
      station.avhrr_sst_field.field = repmat(nan,[0,numel(station.avhrr_sst_yix),numel(station.avhrr_sst_xix)]);
    end;

    %DEBUG:    alldts = datenum(1993,8,25):datenum(1994,1,31);
    %DEBUG:    alldts = datenum(2011,7,1):now;

    [yrs,mos,dys] = datevec(alldts);
    yrmos = unique((yrs*100) + mos);


    N = 0;
    skipped = 0;

    for yrmoix = 1:length(yrmos)

      yrmo = yrmos(yrmoix);

      yr = floor(yrmo / 100);
      mo = rem(yrmo, 100);
      % fprintf('%04d-%02d', yr, mo);
      % fprintf('\n');

      % % *OLD* Sample URL (synoptic images) as of 2011 July:
      % % http://www.imars.usf.edu/husf_avhrr/products/images/florida/2009.03/n15.20090303.2211.florida.true.png
      % baseurl = sprintf('http://www.imars.usf.edu/husf_avhrr/products/images/%s/%04d.%02d', ...
      %                   region,yr,mo);

      % *NEW* Sample URL (synoptic images) for all years/months:
      % http://www.imars.usf.edu/husf_avhrr/daily/1km/sst/florida/2009.03/n15.20090303.2211.florida.true.png
      baseurl = sprintf('http://www.imars.usf.edu/husf_avhrr/daily/1km/sst/%s/%04d.%02d', ...
                        region,yr,mo);


      % Reading each monthly dir each time we load a new station is too slow
      %% fnames = urlread(baseurl);

      fnamespath = fullfile(avhrrpath,sprintf('images_%s_%04d_%02d.html',region,yr,mo));
      if ( ~exist(fnamespath,'file') )
        try,
          %DEBUG:
          disp(['URLWRITE ' baseurl]);
          urlwrite(baseurl,fnamespath);
          fprintf('%04d-%02d', yr, mo);
          fprintf('\n');
        catch
          warning('avhrrsst:BadDirURL','Unable to download "%s"!',baseurl);
          continue;
        end;
      end; %if ( ~exist(fnamespath,'file') )

      fnames = urlread(['file:' fnamespath]);


      % Go back and find all dates within this year-month
      dtixes = find((yrs == yr) & (mos == mo));

      for dtix = dtixes(:)'

        dy = dys(dtix);
        %DEBUG:        fprintf(' %02d', dy);

        % Extract matching filenames from USF's directory listing
        patt = sprintf('>n[0-9][0-9]*[.]%04d%02d%02d[.][0-9][0-9][0-9][0-9][.]%s[.]true.png<', ...
                       yr, mo, dy, region);
        begix = regexp(fnames, patt);

        if ( isempty(begix) )
          warning('avhrrsst:NoMatches','No match: "%s"',patt);
          continue;
        end;

        % Figure out which filename is which
        for ix = begix(:)'
          endix = strfind(fnames(ix:end), '<');
          fname = fnames(ix+1:ix+endix(1)-2);

          fpath = fullfile(avhrrpath,fname);
          url = sprintf('%s/%s', baseurl, fname);

          if ( ~exist(fpath,'file') )
            try,
              %DEBUG:
              disp(['URLWRITE ' url]);
              urlwrite(url,fpath);
            catch
              warning('avhrrsst:BadURL','Unable to download "%s"!',url);
              continue;
            end;
          end; %if ( ~exist(fpath,'file') )

          sstbytes = imread(fpath);

          [ig,ig,YYYY,MM,DD,hh,mm] = parseusfurl(url);

          sst = read_avhrr_subset(sstbytes);
          sst = sst(station.avhrr_sst_yix,station.avhrr_sst_xix);

          % Only take images with some cloud-free data surrounding our site!
          if ( length(find(isfinite(sst(:)))) >= (minPct*numel(sst(:))) )
            N = N + 1;
            station.avhrr_sst_field.date(starting_index+N,1) = datenum(YYYY,MM,DD,hh,mm,0);
            station.avhrr_sst_field.field(starting_index+N,1:size(sst,1),1:size(sst,2)) = sst;
          else
            skipped = skipped + 1;
          end;

          clear sst;
          clear sstbytes;

        end; %for ix = begix(:)'

      end; %for dtix = dtixes(:)'

    end; %for yrmoix = 1:length(yrmos)

    fprintf('\n');
    % DEBUG:
    if ( skipped > 0 )
      fprintf('\nSkipped %d cloudy or bad images!\n', skipped);
      fprintf('%d (%g%%) useable images\n', N, (100*N)/(N+skipped));
    end;

    % Field .date should be a column vector (Nx1)
    station.avhrr_sst_field.date = station.avhrr_sst_field.date(:);

    % Make sure our dates increase monotonically
    [station.avhrr_sst_field.date,sortix] = sort(station.avhrr_sst_field.date);
    station.avhrr_sst_field.field = station.avhrr_sst_field.field(sortix,:,:);


    %DEBUG:    disp(['NOT Saving ' matfname]);
    disp(['Saving ' matfname]);
    save(matfname,'station');

  end; %if ( exist(matfname,'file') )


  flds = fieldnames(station);
  for fldix = 1:length(flds)
    fld = flds{fldix};
    stn.(fld) = station.(fld);
  end;
  station = []; clear station;


  if ( isfield(stn,'lat') && isfield(stn,'lon') )
    stn.avhrr_sst.date = stn.avhrr_sst_field.date;
    stn.avhrr_sst.data = ...
        interp_field(stn.avhrr_sst_field.lat,stn.avhrr_sst_field.lon,...
                     stn.avhrr_sst_field.field,stn.lat,stn.lon,interpMethod);

    % Interpolate hourly time series from synoptic (intermittent) time
    % series: normally we hope for at least one SST value every 3.5 days.
    % Also note default 'spline' interpolation can give bizarre results!
    stn.hourly_avhrr_sst = interp_ts(stn.avhrr_sst,3.5,[],'linear');
  end;

  if ( calcFieldTerms )
    % Get hourly time series from synoptic fields with appropriate INTERP_TS args
    if ( isfield(stn,'lat') && isfield(stn,'lon') )
      % Interpolate hourly time series from synoptic (intermittent) time
      % series: normally we hope for at least one gradient value every week.
      stn = calc_field_terms(stn,'avhrr_sst_field','avhrr_sst',interpMethod,stn.lat,stn.lon,npts,[],true,7,[],'linear');
    else
      stn = calc_field_terms(stn,'avhrr_sst_field',[],[],[],[],npts,[],true,7,[],'linear');
    end;
  end;

  %DEBUG:
  toc,

  set_more;

return;
