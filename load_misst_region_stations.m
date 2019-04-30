function stns = load_misst_region_stations(region,stns_or_stnms)
%function stns = load_misst_region_stations(region,stns_or_stnms)
%
% Calculate MISST indices and load MISST time series for multiple sites
%
% Last Saved Time-stamp: <Tue 2012-07-17 22:19:03  lew.gramer>

  set_more off;
  %DEBUG:
  tic,

  datapath = get_ecoforecasts_path('data');

  if ( ~exist('region','var') || ~ischar(region) )
    error('First arg must be a MISST region name string!');
  end;

  if ( ~exist('stns_or_stnms','var') || isempty(stns_or_stnms) )
    stns = load_region_metadata(region);
  elseif ( ischar(stns_or_stnms) )
    stns.station_name = stns_or_stnms;
  elseif ( iscellstr(stns_or_stnms) )
    for ix = 1:length(stns_or_stnms)
      stns(ix).station_name = stns_or_stnms{ix};
    end;
  elseif ( isstruct(stns_or_stnms) )
    stns = stns_or_stnms;
  else
    error('Second arg if given must be station name(s) or STN struct(s)!');
  end;

  if ( ~isfield(stns,'station_name') )
    error('Found no .station_name for stations!');
  elseif ( ~isfield(stns,'lon') || ~isfield(stns,'lat') )
    % Make every effort possible to look up lat/lon for each station
    stns(1).lon = []; stns(1).lat = [];
    for ix = 1:length(stns)
      try
        [stns(ix).lon,stns(ix).lat,stns(ix).depth] = get_station_coords(stns(ix).station_name);
      catch
      end;
    end;
    if ( any(isempty({stns.lon})) || any(isempty({stns.lat})) )
      rgnstns = load_region_metadata(region);
      rgnflds = fieldnames(rgnstns);
      badix = [];
      for ix = 1:length(stns)
        rgnix = find(strcmpi(stns(ix).station_name,{rgnstns.station_name}),1);
        if ( isempty(rgnix) )
          warning('MISST:NoCoords',...
                  'No coordinates could be found for "%s"!',...
                  stns(ix).station_name);
          badix(end+1) = ix;
        else
          for cfld = rgnflds{:}
            fld = cfld{:};
            stns(ix).(fld) = rgnstns(rgnix).(fld);
          end;
        end;
      end;
      clear rgnstns;
      stns(badix) = [];
      if ( isempty(stns) )
        error('No valid stations found!');
      end;
    end;
  end;


  region = lower(region);

  worldOnly = ( strcmpi(region,'world') );

  disp(worldOnly);

  needix = 1:numel(stns);
  for ix = 1:numel(stns)
    matfname = fullfile(datapath,[stns(ix).station_name '_misst_' region '.mat']);
    % If we already loaded data one MISST file at a time, don't do it again
    if ( exist(matfname, 'file') )
      disp(['Reloading ' matfname]);
      load(matfname,'station');
      stations(ix) = station;
      station = []; clear station
      needix(needix == ix) = [];
    else
      flds = {'station_name','name','fname','lon','lat','depth','misst_region','misst_lonix','misst_latix'};
      for fldix = 1:length(flds)
        fld = flds{fldix};
        if ( isfield(stns,fld) )
          stations(ix).(fld) = stns(ix).(fld);
        end;
      end;
      if ( isfield(stations(ix),'station_name') && ~isfield(stations(ix),'name') )
        stations(ix).name = stations(ix).station_name;
      end;
    end;
  end;


  % Extract MISST SSTs as needed
  if ( ~isempty(needix) )

    disp('Extracting SST from original MISST binary files...');

    % Find lat/lon indices for each station in MISST dataset for REGION
    if ( worldOnly )
      lonfld = 'misst_lonix';
      latfld = 'misst_latix';
    else
      lonfld = ['misst_' region '_lonix'];
      latfld = ['misst_' region '_latix'];
    end;
    if ( ~isfield(stations,lonfld) )
      stations(1).(lonfld) = [];
    end;
    if ( ~isfield(stations,latfld) )
      stations(1).(latfld) = [];
    end;
    [lon,lat,ig,dx] = read_misst_region(region);
    stations(needix) = load_misst_cfg(stations(needix),region);
    for ix=needix(:)'
      if ( isempty(stations(ix).(lonfld)) || isempty(stations(ix).(latfld)) )
        stations(ix).(lonfld) = interp1(lon,1:length(lon),stations(ix).lon,'nearest');
        stations(ix).(latfld) = interp1(lat,1:length(lat),stations(ix).lat,'nearest');
      end;
    end;

    % If region is not 'world', find corresponding indices for global dataset also
    if ( ~worldOnly )
      [wlon,wlat,ig,wdx] = read_misst_region('world');
      rgnlonix = interp1(wlon,1:length(wlon),lon,'nearest');
      rgnlatix = interp1(wlat,1:length(wlat),lat,'nearest');

      % stations(needix) = load_misst_cfg(stations(needix),'world');
      for ix=needix(:)'
        stations(ix).misst_lonix = interp1(wlon,1:length(wlon),lon(stations(ix).(lonfld)),'nearest');
        stations(ix).misst_latix = interp1(wlat,1:length(wlat),lat(stations(ix).(latfld)),'nearest');
      end;
    end;


    thisyr = get_year(now - 1);
    thisjd = get_jday(now - 1);
    if ( worldOnly )
      % Once per day, Julian Day 2005-235 to "the present"
      begyear = 2005;  begjday = 235;
    else
      % Once per day, Julian Day 2002-185 to "the present"
      begyear = 2002;  begjday = 185;
    end;
    %endyear = 2011;  endjday = 61;
    endyear = thisyr;  endjday = thisjd;
    alldts = datenum(begyear,1,begjday):datenum(endyear,1,endjday);

    % Make subsets rectangular, to avoid dyslexic errors
    xrad = 4;
    yrad = 3;

    for ix = needix(:)'
      stations(ix).misst_sst.date = repmat(nan,[length(alldts) 1]);
      stations(ix).misst_sst.data = repmat(nan,[length(alldts) 1]);

      stations(ix).misst_sst_field.date = repmat(nan,[length(alldts) 1]);

      xix = min(stations(ix).(lonfld))-xrad:max(stations(ix).(lonfld))+xrad;
      yix = min(stations(ix).(latfld))-yrad:max(stations(ix).(latfld))+yrad;
      stations(ix).misst_sst_field.lon = lon(xix);
      stations(ix).misst_sst_field.lat = lat(yix);

      stations(ix).misst_sst_field.field = repmat(nan,[length(alldts),length(lat(yix)),length(lon(xix))]);
    end;


    warning('off','MISST:NoRawFile');

    yrs = get_year(alldts(1)):get_year(alldts(end));

    for yr = yrs(:)'

      %DEBUG:
      disp(yr);

      switch (yr),
       case begyear,
        jds = begjday:365;
       case endyear,
        jds = 1:endjday;
       otherwise,
        jds = 1:365;
        if ( mod(yr,4) == 0 )
          jds = 1:366;
        end;
      end;

      for jd = jds(:)'
        dt = datenum(yr,1,1) + jd - 1;
        %DEBUG:        disp(datestr(dt));

        dtix = find(alldts == dt);
        if ( isempty(dtix) )
          error('Found no matching date for %g ("%s")?!',dt,datestr(dt));
        end;

        [ig,ig,sst] = read_misst_region(region,yr,jd);
        if ( isempty(sst) && ~worldOnly )
          [ig,ig,wsst] = read_misst_region('world',yr,jd);
          % Subset so each site's LATFLD and LONFLD will still index properly!
          if ( ~isempty(wsst) )
            if ( yr < 2010 )
              warning('MISST:DefaultToWorld','Using WORLD file "%s"',datestr(dt));
            end;
            sst = wsst(rgnlatix,rgnlonix);
          end;
        end;

        if ( isempty(sst) )
          warning('MISST:MissingDay','No MISST data for year %g day %g',yr,jd);

        else
          for ix = needix(:)'
            dat = sst(stations(ix).(latfld),stations(ix).(lonfld));

            stations(ix).misst_sst.date(dtix,1) = dt;
            stations(ix).misst_sst.data(dtix,1) = nanmedian(dat(:));

            xix = min(stations(ix).(lonfld))-xrad:max(stations(ix).(lonfld))+xrad;
            yix = min(stations(ix).(latfld))-yrad:max(stations(ix).(latfld))+yrad;
            stations(ix).misst_sst_field.date(dtix,1) = dt;
            stations(ix).misst_sst_field.field(dtix,:,:) = sst(yix,xix);
          end; %for ix = needix(:)'

        end; %if ( isempty(sst) ) else

      end; %for jd = jds(:)'

    end; %for yr = yrs(:)'

    warning('on','MISST:NoRawFile');


    for ix = needix(:)'
      % Basic QA - remove days with missing or bad files
      badix = find( ~isfinite(stations(ix).misst_sst.date) | ...
                    ~isfinite(stations(ix).misst_sst_field.date) );
      %             ~isfinite(stations(ix).misst_sst.data) | ...
      %             all(~isfinite(stations(ix).misst_sst_field.field(:))) | ...
      stations(ix).misst_sst.date(badix) = [];
      stations(ix).misst_sst.data(badix) = [];
      stations(ix).misst_sst_field.date(badix) = [];
      stations(ix).misst_sst_field.field(badix,:,:) = [];

      matfname = fullfile(datapath,[stations(ix).station_name '_misst_' region '.mat']);
      disp(['Saving ' matfname]);
      station = stations(ix);
      save(matfname,'station');
      station = []; clear station
    end; %for ix = needix(:)'

  end; %if ( ~isempty(needix) )


  flds = fieldnames(stations);
  for fldix = 1:length(flds)
    fld = flds{fldix};
    for ix = 1:length(stations)
      stns(ix).(fld) = stations(ix).(fld);
      stations(ix).(fld) = [];
    end;
  end;
  stations = []; clear stations

  %DEBUG:
  toc,
  set_more;

return;
