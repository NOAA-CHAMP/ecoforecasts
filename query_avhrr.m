function dts = query_avhrr(indts, bbox_or_stanm_list, region)
%function dts = query_avhrr(indts, bbox_or_stanm_list, region)
%
% Query synoptic AVHRR SST imagery from USF online archive: Check all images
% for DATES (not specific timestamps) matching any element of INDTS; subset
% each image to bounding box specified in, or a 81x81-pixel box surrounding
% SEAKEYS station named by, each element of cell array BBOX_OR_STANM_LIST.
% Use dataset for USF area REGION (DEFAULT: 'florida').
%
% Returns a cell array 'dts' the same size as 'bbox_or_stanm_list'; each
% element of array is a vector of useable image-timestamps for that box.
%
% NOTE: Only returns SST fields where at least 50% of the 16x16-pixel area
% surrounding the center of the bounding box (or station site) are not NaN.
%
% Last Saved Time-stamp: <Mon 2016-02-01 16:12:35 Eastern Standard Time gramer>

  dts = {};
  sstix = {};
  skipped = {};

  % Store/retrieve data in this M-file's local directory
  %[pathroot, ig, ig, ig] = fileparts(mfilename('fullpath'));
  [pathroot, ig, ig] = fileparts(mfilename('fullpath'));
  if ( ~exist('datapath', 'var') || isempty(datapath) )
    datapath = get_ecoforecasts_path('data');
  end;
  avhrrpath = fullfile(datapath,'avhrr');

  if ( ~exist('bbox_or_stanm_list', 'var') || isempty(bbox_or_stanm_list) )
    bbox_or_stanm_list = {};
  end;
  if ( ~exist('region','var') || isempty(region) )
    region = 'florida';
  end;


  switch ( region )
   case 'florida',
    minlon = -91; maxlon = -79; dlon = (1/110);
    minlat =  22; maxlat =  31; dlat = (1/110);
   otherwise,
    error('Unrecognized region "%s"!', region);
  end;


  boxradius = 40;

  for boxix = 1:length(bbox_or_stanm_list)
    dts{boxix} = [];
    sstix{boxix} = 1;
    skipped{boxix} = 0;
  end;


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


    % yrmo_ext = fullfile(datapath, ...
    %                     sprintf( '%06d-%03dx%03d.mat', yrmo, ...
    %                              ((boxradius*2)+1), ((boxradius*2)+1) ));


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

        try, sstbytes_all = imread(url);
        catch
          warning('querysst:BadURL', ...
                  'Unable to download "%s/%s"!', ...
                  url);
          continue;
        end;

        [ig,ig,YYYY,MM,DD,hh,mm] = parseusfurl(url);

        for boxix = 1:length(bbox_or_stanm_list)
          bbox_or_stanm = bbox_or_stanm_list{boxix};
          if ( isempty(bbox_or_stanm) )
            % DEFAULT: Show whole Florida Reef Tract and Straits of Florida
            stanm{boxix} = 'fknms';
            bbox_or_stanm = [ -83.5 -79.5 24 26 ];
          end;
          if ( ~ischar(bbox_or_stanm) )
            stanm{boxix} = 'fknms';
            bbox = bbox_or_stanm;
            miny = (bbox(1) - minlon) / dlon;
            maxy = (bbox(2) - minlon) / dlon;
            maxx = (maxlat - bbox(3)) / dlat;
            minx = (maxlat - bbox(4)) / dlat;
            stn_x = minx + round((maxx+minx)/2);
            stn_y = miny + round((maxy+miny)/2);
          else
            stanm{boxix} = bbox_or_stanm;
            [stn_lat, stn_lon] = station_latlon(stanm{boxix});
            stn_x = round((maxlat - stn_lat) / dlat);
            stn_y = round((stn_lon - minlon) / dlon);
            minx = stn_x - boxradius; maxx = stn_x + boxradius;
            miny = stn_y - boxradius; maxy = stn_y + boxradius;
          end;
          sstdims = [maxx-minx+1, maxy-miny+1];

          pname = fullfile(avhrrpath, [stanm{boxix} '-' fname]);
%           if ( exist(pname, 'file') )
%             dts{boxix}{sstix{boxix}} = datenum(YYYY,MM,DD,hh,mm,0);
%             sstix{boxix} = sstix{boxix} + 1;
%             continue;
%           end;

          sstbytes = sstbytes_all(minx:maxx, miny:maxy);

          sst = read_avhrr_subset(sstbytes);

          % Only take images with some cloud-free data surrounding our site!
          ctr = round(size(sst)./2);
          ctrsst = sst((ctr(1)-8+1):(ctr(1)+8),(ctr(2)-8+1):(ctr(2)+8));
          % if ( length(find(isnan(ctrsst))) < (0.50*numel(ctrsst)) && ...
          %      ~isnan(sst(ctr(1), ctr(2))) )
          if ( length(find(isnan(ctrsst))) < (0.50*numel(ctrsst)) )
            dts{boxix}{sstix{boxix}} = datenum(YYYY,MM,DD,hh,mm,0);
            sstix{boxix} = sstix{boxix} + 1;
          else
            skipped{boxix} = skipped{boxix} + 1;
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

  end;

  % DEBUG:
  for boxix = 1:length(bbox_or_stanm_list)
    sstix{boxix} = sstix{boxix} - 1;
    fprintf('\nBOX "%s"\nSkipped %d cloudy or bad images!\n', stanm{boxix}, skipped{boxix});
    fprintf('%d (%g%%) useable images\n', sstix{boxix}, ...
            (100*sstix{boxix}/(sstix{boxix}+skipped{boxix})));
  end;

return;
