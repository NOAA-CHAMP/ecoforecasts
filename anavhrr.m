function [dts,ssts,urls] = anavhrr(indts, bbox_or_stanm, region, acceptable_bad_pct)
%function [dts,ssts,urls] = anavhrr(indts, bbox_or_stanm, region, acceptable_bad_pct)
%
% Analyze synoptic AVHRR SST imagery from USF online archive: Get all images
% for DATES (not specific timestamps) matching any element of 'indts'; subset
% each image to bounding box specified in, or a 81x81-pixel box surrounding
% SEAKEYS station named by, 'bbox_or_stanm' (DEFAULT use bounding box showing
% whole FKNMS). Use dataset for USF area 'region' (DEFAULT region='florida'.)
%
% Returns a vector of actual image-timestamps, and an NxMxP matrix of SST
% values. (N is the number of useable images found; MxP is the number of
% pixels in the area specified by 'bbox_or_stanm'.) Also optionally returns
% a cell array of the URL strings (full web addresses) for each USF AVHRR SST
% image successfully processed into the 'dts' and 'ssts' return values.
% 
% NOTE: Only returns SST fields where at least 50% of the 16x16-pixel area
% surrounding the center of the bounding box (or station site) are not NaN.
% Optionally, fourth arg 'acceptable_bad_pct' can specify a number other 
% than the default "0.50" for this acceptable percentage of bad pixels.
%
% Last Saved Time-stamp: <Fri 2011-08-19 11:29:09  Lew.Gramer>


  dts = [];
  ssts = [];
  urls = {};

  sstix = 1;
  skipped = 0;

  % Store/retrieve data in this M-file's local directory
  if ( ~exist('datapath', 'var') || isempty(datapath) )
    datapath = get_ecoforecasts_path('data');
  end;
  avhrrpath = fullfile(datapath,'avhrr');

  if ( ~exist('bbox_or_stanm', 'var') || isempty(bbox_or_stanm) )
    % DEFAULT: Show whole Florida Reef Tract and Straits of Florida
    bbox_or_stanm = [ -84 -79 23 26 ];
  end;
  if ( ~exist('region','var') || isempty(region) )
    region = 'florida';
  end;
  if ( ~exist('acceptable_bad_pct','var') || isempty(acceptable_bad_pct) )
    acceptable_bad_pct = 0.50;
  end;


  minval = +0.0; maxval = +40.0;


  switch ( region )
   case 'florida',
    minlon = -91; maxlon = -79; dlon = (1/110);
    minlat =  22; maxlat =  31; dlat = (1/110);
   otherwise,
    error('Unrecognized region "%s"!', region);
  end;

  boxradius = 40;

  if ( ~ischar(bbox_or_stanm) )
    stanm = 'fknms';
    bbox = bbox_or_stanm;
    miny = (bbox(1) - minlon) / dlon;
    maxy = (bbox(2) - minlon) / dlon;
    maxx = (maxlat - bbox(3)) / dlat;
    minx = (maxlat - bbox(4)) / dlat;
    stn_x = minx + round((maxx+minx)/2);
    stn_y = miny + round((maxy+miny)/2);

  else
    stanm = bbox_or_stanm;

    [stn_lat, stn_lon] = station_latlon(stanm);

    stn_x = round((maxlat - stn_lat) / dlat);
    stn_y = round((stn_lon - minlon) / dlon);

    minx = stn_x - boxradius; maxx = stn_x + boxradius;
    miny = stn_y - boxradius; maxy = stn_y + boxradius;
  end;

  sstdims = [maxx-minx+1, maxy-miny+1];


  [yrs,mos,dys] = datevec(indts);


  yrmos = unique((yrs*100) + mos);

  for yrmoix = 1:length(yrmos)

    yrmo = yrmos(yrmoix);

    yr = floor(yrmo / 100);
    mo = rem(yrmo, 100);
    fprintf('\n%04d %02d: ', yr, mo);

    % Sample URL (synoptic images):
    % http://www.imars.usf.edu/husf_avhrr/products/images/florida/2009.03/n15.20090303.2211.florida.true.png

    baseurl = sprintf('http://www.imars.usf.edu/husf_avhrr/products/images/%s/%04d.%02d', ...
                      region, yr, mo);

    fnames = urlread(baseurl);

    % Go back and find all dates within this year-month
    dtixes = find((yrs == yr) & (mos == mo));

    %%%% ??? NEED TO UNIQUE-IFY DATES BEFORE THIS LOOP!

    for dtix = dtixes

      dy = dys(dtix);
      fprintf('%02d... ', dy);

      % dailysst = repmat(inf, [maxx-minx+1, maxy-miny+1]);
      dailysst = [];
      N = 0;

      % Extract matching filenames from USF's directory listing
      patt = sprintf('>n[0-9][0-9]*[.]%04d%02d%02d[.][0-9][0-9][0-9][0-9][.]%s[.]true.png<', ...
                     yr, mo, dy, region);
      begix = regexp(fnames, patt);

      if ( isempty(begix) )
        warning('querysst:NoMatches', ...
                'Found no matches for "%s/%s"!', ...
                baseurl, patt);
        continue;
      end;

      % Figure out which filename is which
      for ix = begix
        endix = strfind(fnames(ix:end), '<');
        fname = fnames(ix+1:ix+endix(1)-2);

        url = sprintf('%s/%s', baseurl, fname);

        pname = fullfile(avhrrpath, [stanm '-' fname]);
        if ( exist(pname, 'file') )
          try, sstbytes = imread(pname);
          catch, delete(pname);
          end;
        end;

        if ( ~exist(pname, 'file') )
          try, sstbytes = imread(url);
          catch
            warning('querysst:BadURL', ...
                    'Unable to download "%s/%s"!', ...
                    url);
            continue;
          end;
          % Subset to user-specified bounding box - default is Straits of Florida
          sstbytes = sstbytes(minx:maxx, miny:maxy);
        end;

        if ( length(sstbytes) < 2 )
          % A one-byte image file means there was NO data from that pass
          skipped = skipped + 1;
          continue;
        end;

        sst = read_avhrr_subset(sstbytes);

        % Only take images with some cloud-free data surrounding our site!
        ctr = round(size(sst)./2);
        ctrsst = sst((ctr(1)-8+1):(ctr(1)+8),(ctr(2)-8+1):(ctr(2)+8));
        % if ( length(find(isnan(ctrsst))) < (0.50*numel(ctrsst)) && ...
        %      ~isnan(sst(ctr(1), ctr(2))) )
        if ( length(find(isnan(ctrsst))) < (acceptable_bad_pct*numel(ctrsst)) )

          % Correct for diurnal cycle by removing synoptic median
          sst = sst - nanmean(sst(:));

          [ig,ig,YYYY,MM,DD,hh,mm] = parseusfurl(url);

          dts(sstix) = datenum(YYYY,MM,DD,hh,mm,0);
          if ( nargout > 1 )
            ssts(sstix, 1:size(sst,1), 1:size(sst,2)) = sst;
          end;
          if ( nargout > 2 )
            urls{sstix} = url;
          end;
          sstix = sstix + 1;

        else
          skipped = skipped + 1;
          sstbytes = cast(0, 'uint8');
        end;

        % Try to save STATION-specific image subset for future reference.
        % NOTE: If data in image was unuseable, 'sstbytes' will just be [0].
        if ( numel(sstbytes) <= 1 || ~exist(pname, 'file') )
          imwrite(sstbytes, pname);
        end;

        clear sst;
        clear ctrsst;
        clear sstbytes;

      end;

    end;

  end;

  % DEBUG:
  sstix = sstix - 1;
  fprintf('\nSkipped %d cloudy or bad images!\n', skipped);
  fprintf('%d (%g%%) useable images\n', sstix, ...
          (sstix/(sstix+skipped)));

return;
