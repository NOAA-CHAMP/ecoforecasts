function stn = extract_avhrr_field(stn_or_stnm,indts,interpMethod)
%function stn = extract_avhrr_field(stn_or_stnm,indts,interpMethod)
%
% Query synoptic AVHRR SST imagery from USF online archive: Download images
% in USF area REGION (DEFAULT: 'florida'), matching *DATES* (*NOT* specific
% timestamps) INDTS; subset each image to a 17x19-pixel box surrounding the
% station in or named by STN_OR_STNM (must be a STRUCT or CHAR).
%
% Returns STN struct with field .avhrr_sst_field, containing the subfields
% '.lon,'.lat','.date','.field'; also interpolates field using INTERPMETHOD
% (DEFAULT: 'linear') to produce time series field STATION.avhrr_sst.
%
% NOTE: Only returns SST fields where >=50% of the pixels surrounding the
% center of the field are valid (not NaN).
%
% Last Saved Time-stamp: <Fri 2011-08-19 11:29:04  Lew.Gramer>

  set_more off;

  % Store/retrieve data in this M-file's local directory
  datapath = get_ecoforecasts_path('data');
  avhrrpath = fullfile(datapath,'avhrr');

  stn = get_station_from_station_name(stn_or_stnm);

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

  rgnlon = [minlon:dlon:maxlon]';
  rgnlat = [minlat:dlat:maxlat]';

  xradius = 8;
  yradius = 9;

  matfname = fullfile(datapath,sprintf('%s-%03dx%03d.mat',...
                                       station.station_name,...
                                       ((xradius*2)+1),((yradius*2)+1)));


  if ( exist(matfname,'file') )
    disp(['Reloading ' matfname]);
    load(matfname,'station');

  else
    disp('Extracting from raw PNG files...');

    station = get_station_from_station_name(stn.station_name);

    station = safe_rmfield(station,{'avhrr_sst_field','avhrr_sst'});
    station.avhrr_sst_midx = interp1([1:length(rgnlon)]',rgnlon,stn.lon);
    station.avhrr_sst_midy = interp1([1:length(rgnlat)]',rgnlat,stn.lat);

    station.avhrr_sst_xix = [station.avhrr_sst_midx-xradius:station.avhrr_sst_midx+xradius]';
    station.avhrr_sst_yix = [station.avhrr_sst_midy-yradius:station.avhrr_sst_midy+yradius]';

    % Start of USF AVHRR coverage with 'n11.19930825.2209.florida.png'
    alldts = datenum(1993,8,25):now;
    [yrs,mos,dys] = datevec(alldts);
    yrmos = unique((yrs*100) + mos);
    skipped = {};

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
            minx = stn_x - xradius; maxx = stn_x + xradius;
            miny = stn_y - yradius; maxy = stn_y + yradius;
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

    disp(['Saving ' matfname]);
    save(matfname,'station');

  end;

  % DEBUG:
  for boxix = 1:length(bbox_or_stanm_list)
    sstix{boxix} = sstix{boxix} - 1;
    fprintf('\nBOX "%s"\nSkipped %d cloudy or bad images!\n', stanm{boxix}, skipped{boxix});
    fprintf('%d (%g%%) useable images\n', sstix{boxix}, ...
            (100*sstix{boxix}/(sstix{boxix}+skipped{boxix})));
  end;

  set_more;

return;
