function stn = read_fdep_stevens_data(stn_or_stcd,dts,forceReload,doQC)
%function stn = read_fdep_stevens_data(stn_or_stcd,dts,forceReload,doQC)
% 
% Read CSV files(s) previously extracted (v. EXTRACT_FDEP_STEVENS_DATA) for
% individual month(s) of environmental data for a Florida Dept. of Envir.
% Protection (FDEP) Stevens Water Monitoring Systems coastal monitoring
% station STN (e.g., with .station_name='fdepk') or STCD (DEFAULT: 'k' for
% their "Station K" near St. Lucie Inlet). See: http://www.fldep-stevens.com
%
% If FORCERELOAD (Default: false), ignore any existing 'FDEP?.MAT' file.
% If DOQC (Default: true), remove timestamps with *any* "obviously bad" data.
%
% Returns STN with time series fields added or overwritten. If DTS a vector
% of DATENUM, time series are restricted to the range [MIN(DTS),MAX(DTS)].
% 
% NOTE: If STN is passed in with FDEP data outside the range of DTS, those
% data will remain in the STN struct when it is returned! If you do not want
% that, use RMFIELD to remove those fields from STN first.
% 
% Last Saved Time-stamp: <Wed 2016-08-24 01:03:38 Eastern Daylight Time gramer>

  datpath = get_ecoforecasts_path('data');

  if ( ~exist('stn_or_stcd','var') || isempty(stn_or_stcd) )
    stn_or_stcd = 'k';
  end;
  if ( ~exist('dts','var') || isempty(dts) )
    [stn,dts] = get_fdep_stevens_metadata(stn_or_stcd);
  else
    [stn,ig] = get_fdep_stevens_metadata(stn_or_stcd);
    ig=[]; clear ig
  end;
  clear stn_or_stcd

  if ( ~exist('forceReload','var') || isempty(forceReload) )
    forceReload = false;
  end;
  if ( ~exist('doQC','var') || isempty(doQC) )
    doQC = true;
  end;

  matfname = fullfile(datpath,sprintf('%s.mat',stn.station_name));
  if ( ~forceReload && exist(matfname,'file') )
    disp(['Loading ',matfname]);
    load(matfname);

  else
    mindt = min(dts);
    maxdt = max(dts);
    yrmos = unique(get_yearmonth(dts));

    curyr = 0;

    for begdt = yrmos(:)';
      yr = get_year(begdt);
      mo = get_month(begdt);
      if ( curyr ~= yr )
        curyr = yr;
        disp(curyr);
      end;

      monmatfname = fullfile(datpath,sprintf('%s-%04d%02d.mat',stn.station_name,yr,mo));
      if ( exist(monmatfname,'file') )
        %disp(['Loading ',monmatfname]);
        load(monmatfname);

      else

        enddt = datenum(yr,mo+1,0,23,59,59);

        fname = fullfile(datpath,sprintf('%s-%04d%02d.csv',stn.station_name,yr,mo));

        s = fileread(fname);
        if ( isempty(regexp(s,'\w')) )
          warning('Skipping empty %s',fname);
          continue;
        end;
        [chdrln,begdat] = textscan(s,'%[^\n]',1);
        hdrln = chdrln{1}{1};
        chdrs = textscan(hdrln,'%s','Delimiter',',');
        hdrs = chdrs{1};

        % SAMPLE HEADER:
        % Date / Time (UTC),Battery (v),Water Temp (C) Lower,Water Conductivity (mS/cm) Lower,Water Salinity (ppt) Lower,Water Temp (C) Upper,Water Conductivity (mS/cm) Upper,Water Salinity (ppt) Upper,Wind Direction (Deg),Wind Speed (knots),Air Temp (C),Relative Humidity (%),Barometric  Pressure (hBa),Rainfall (mm),NAVD88 Water Level (m)

        fmt = '"%[^"]"';
        flds = {'datetime'};
        for hdrix = 2:numel(hdrs)
          % If we want to skip a field, make its FLDS elt. an empty string
          switch (hdrs{hdrix})
           case 'Battery (v)',			fmt = [fmt,',"%f"'];	flds{end+1} = 'fdep_battv';
           case 'Water Temp (C) Lower',		fmt = [fmt,',"%f"'];	flds{end+1} = 'fdep_seatemp_deep';
           case 'Water Conductivity (mS/cm) Lower', fmt = [fmt,',"%f"']; flds{end+1} = 'fdep_seacond_deep';
           case 'Water Salinity (ppt) Lower',	fmt = [fmt,',"%f"'];	flds{end+1} = 'fdep_salin_deep';
           case 'Water Temp (C) Upper',		fmt = [fmt,',"%f"'];	flds{end+1} = 'fdep_seatemp_shallow';
           case 'Water Conductivity (mS/cm) Upper', fmt = [fmt,',"%f"']; flds{end+1} = 'fdep_seacond_shallow';
           case 'Water Salinity (ppt) Upper',	fmt = [fmt,',"%f"'];	flds{end+1} = 'fdep_salin_shallow';
           case 'Wind Direction (Deg)',		fmt = [fmt,',"%f"'];	flds{end+1} = 'fdep_wind_dir';
           case 'Wind Speed (knots)',		fmt = [fmt,',"%f"'];	flds{end+1} = 'fdep_wind_speed';
           case 'Air Temp (C)',			fmt = [fmt,',"%f"'];	flds{end+1} = 'fdep_airtemp';
           case 'Relative Humidity (%)',	fmt = [fmt,',"%f"'];	flds{end+1} = 'fdep_relhumid';
           case 'Barometric  Pressure (hBa)',	fmt = [fmt,',"%f"'];	flds{end+1} = 'fdep_airpres';
           case 'Rainfall (mm)',		fmt = [fmt,',"%f"'];	flds{end+1} = 'fdep_rainfall';
           case 'NAVD88 Water Level (m)',	fmt = [fmt,',"%f"'];	flds{end+1} = 'fdep_reltide';
           otherwise,
            error('Unrecognized header column #%g (%s) in %s',hdrix,hdrs{hdrix},fname);
          end; %switch (hdrs{hdrix})
        end; %for hdrix = 2:numel(hdrs)

        % Stupid comma in pressure (e.g., "1,024.3") needs a more elegant fix!
        s = regexprep(s,',"([0-9][0-9]*),',',"$1');

        C = textscan(s,fmt,'Delimiter',',','EndOfLine','\r\n','Headerlines',1,'TreatAsEmpty','-');
        if ( isempty(C) || isempty(C{1}) )
          warning('Skipping dateless %s',fname);
          continue;
        end;

        dts = datenum(C{1});

        res = [];

        % These query results are generally time-inverted
        [dts,sortix] = sort(dts);
        for fldix=2:numel(flds)
          fld = flds{fldix};
          if ( ~isempty(fld) )
            res.(fld).date = dts;
            res.(fld).data = C{fldix}(sortix);
          end;
        end; %for fldix=2:numel(flds)

        if ( isempty(res) )
          warning('Skipping data-free %s',fname);
          continue;
        end;

        % Save raw data for this year-month to its own MAT file
        disp(['Saving ',monmatfname]);
        save(monmatfname,'res');

      end; %if ( exist(monmatfname,'file') )

      % Limit output to the requested dates, and to good data...
      % (Easier to do these checks for each month individually)

      flds = fieldnames(res);
      for fldix=1:numel(flds)
        fld = flds{fldix};
        badix = find(mindt > res.(fld).date | res.(fld).date > maxdt);
        res.(fld).date(badix) = [];
        res.(fld).data(badix) = [];
      end; %for fldix=1:numel(flds)

      if ( doQC )
        baddts = [];
        if ( isfield(stn,'fdep_battv') )
          baddts = union(baddts,res.fdep_battv.date(11>res.fdep_battv.data | res.fdep_battv.data>14));
        end;
        if ( isfield(stn,'fdep_airtemp') )
          baddts = union(baddts,res.fdep_airtemp.date(res.fdep_airtemp.data<=-5));
        end;
        if ( isfield(stn,'fdep_airpres') )
          baddts = union(baddts,res.fdep_airpres.date(res.fdep_airpres.data<=980));
        end;
        if ( isfield(stn,'fdep_seatemp_shallow') )
          baddts = union(baddts,res.fdep_seatemp_shallow.date(res.fdep_seatemp_shallow.data<=0));
        end;
        if ( isfield(stn,'fdep_seatemp_deep') )
          baddts = union(baddts,res.fdep_seatemp_deep.date(res.fdep_seatemp_deep.data<=0));
        end;
        if ( isfield(stn,'fdep_reltide') )
          baddts = union(baddts,res.fdep_reltide.date(-3>res.fdep_reltide.data | res.fdep_reltide.data>+3));
        end;
        if ( numel(baddts) > 0 )
          baddts = unique(floor(baddts));

          %DEBUG:
          disp(['Removing ',num2str(numel(baddts)),' bad dates: ',datestr(begdt)]);
          flds = fieldnames(res);
          for fldix=1:numel(flds)
            fld = flds{fldix};
            rawfld = ['raw_',fld];
            badix = find(ismember(floor(res.(fld).date),baddts));
            res.(rawfld) = res.(fld);
            res.(fld).date(badix) = [];
            res.(fld).data(badix) = [];
          end; %for fldix=1:numel(flds)
        end;
      end;

      w = warning('OFF','Ecoforecasts:mergedNonTS');
      stn = merge_station_data(stn,res);
      warning(w);

      res = []; clear res

    end; %for begdt = yrmos(:)';

    disp(['Saving ',matfname]);
    save(matfname,'stn','-v7.3');

  end; %if ( ~forceReload && exist(matfname,'file') ) else

return;
