function stns = get_hrly_gom_hycom(stns_or_stnms,mindt,maxdt,vars,flds,interpMethod,baseurl,xrad,yrad,datapath)
%function stns = get_hrly_gom_hycom(stns_or_stnms,mindt,maxdt,vars,flds,interpMethod,baseurl,xrad,yrad,datapath)
%
% Add fields for ocean surface U and V current components, sea temperature,
% salinity, mixed-layer depth, and 9x9 U, V, and sea temperature fields from
% assimilative NRL HYCOM+NCODA Gulf of Mexico 1/25 Degree Analysis (Prasad
% and Hogan 2007), to each member of struct array (or scalar) STNS. If a
% station name string or cell array of strings STNMS (five chars., e.g.,
% 'mlrf1', or {'smkf1','tnrf1'}) is given instead of structs, or if structs
% have no valid .lon and .lat fields (e.g., -80.38, 25.01), but all do have
% .station_name fields, trie to retrieve locations with GET_STATION_COORDS.
%
% Optional VARS is a string or cellstr specifying which variables to extract.
% No online catalog seems to be available, but see the URL templates below
% for a hard-won list of all the variables available from this dataset. If
% VARS is specified, optional FLDS may also be specified as the corresp.
% list of field names to be added to structs in STNS, e.g., 'seatemp'.
%
% INTERP_FIELD is passed INTERPMETHOD (DEFAULT 'linear') to subset points.
%
% DEFAULT VARS and corresponding FLDS cell arrays are as follows:
%      vars = { 'u', 'v', 'mld', 'qtot',           'temperature',   'salinity', ...
%               'u',       'v',       'mld',       'temperature',   'salinity'       };
%      flds = { 'u', 'v', 'mld', 'net_heat_flux',  'seatemp',       'salinity', ...
%               'u_field', 'v_field', 'mld_field', 'seatemp_field', 'salinity_field' };
%
% NOTE: The string 'hrly_gom_hycom_' is prepended to each of the names in FLDS.
%
% DEFAULT for BASEURL:
%  'http://hycom-tds1.coaps.fsu.edu/thredds/dodsC/GOMl0.04'
% Information here: http://hycom-tds1.coaps.fsu.edu/thredds/dodsC/GOMl0.04/expt_31.0.html
% If you have another ocean model with a THREDDS interface, specify it here.
%
% (OLD GoM HYCOM BASEURL:
%  'http://tds.hycom.org/thredds/dodsC/GOMl0.04')
%
% SAMPLE URLs (files are stored by full year, experiment 20.1 ended July 2010):
%  http://tds.hycom.org/thredds/dodsC/GOMl0.04/expt_20.1/2003
%  http://tds.hycom.org/thredds/dodsC/GOMl0.04/expt_20.1/2010
%  http://tds.hycom.org/thredds/dodsC/GOMl0.04/expt_30.1/2010
%  http://tds.hycom.org/thredds/dodsC/GOMl0.04/expt_30.1/2011
%
% CALLS: MDATASET (netCDF-Java); INTERP_FIELD, GET_STATION_COORDS (Ecoforecasts).
%
% Last Saved Time-stamp: <Fri 2017-12-01 11:25:37 Eastern Standard Time gramer>

  set_more off;

  % Grid-point radii for time-series fields (e.g., seatemp_field)
  if ( ~exist('xrad','var') || isempty(xrad) )
    %xrad = 4;
    xrad = 1;
  end;
  if ( ~exist('yrad','var') || isempty(yrad) )
    %% Choose one axis slightly larger to guard against transposition errors
    %yrad = 6;
    yrad = 1;
  end;


  if ( ~exist('datapath','var') || isempty(datapath) )
    datapath = get_ecoforecasts_path('data');
  end;

  % % Not currently used - HDF files all accessed remotely
  % gompath = fullfile(datapath,'hycom','GOM');

  if ( isempty(stns_or_stnms) )
    error('First arg was empty!');
  elseif ( ischar(stns_or_stnms) )
    for ix=1:size(stns_or_stnms,1)
      stns(ix).station_name = stns_or_stnms(ix,:);
    end;
  elseif ( iscellstr(stns_or_stnms) )
    for ix=1:numel(stns_or_stnms)
      stns(ix).station_name = stns_or_stnms{ix};
    end;
  elseif ( isstruct(stns_or_stnms) )
    stns = stns_or_stnms;
    if ( ~isfield(stns,'station_name') )
      error('Station STRUCT(s) without station_name field(s)!');
    end;
  else
    error('First arg must either be station STRUCT(S) or station name string(s)!');
  end;
  clear stns_or_stnms;

  if ( ~exist('mindt','var') || isempty(mindt) )
    mindt = -Inf;
  end;
  if ( ~exist('maxdt','var') || isempty(maxdt) )
    maxdt = +Inf;
  end;

  if ( ~exist('vars','var') || isempty(vars) )
    % vars = { 'u', 'v', 'mld', 'qtot',           'temperature',   'salinity', ...
    %          'u',       'v',       'mld',       'temperature',   'salinity'       };
    % flds = { 'u', 'v', 'mld', 'net_heat_flux',  'seatemp',       'salinity', ...
    %          'u_field', 'v_field', 'mld_field', 'seatemp_field', 'salinity_field' };
    vars = { 'u', 'v', 'mld', 'qtot',           'temperature',   'salinity', };
    flds = { 'u', 'v', 'mld', 'net_heat_flux',  'seatemp',       'salinity', };
  end;
  if ( ~exist('flds','var') || isempty(flds) )
    flds = lower(vars);
  end;

  if ( ~exist('interpMethod','var') || isempty(interpMethod) )
    interpMethod = 'linear';
  end;
  if ( interpMethod(1) == '*' )
    interpMethod = interpMethod(2:end);
  end;

  if ( ~exist('baseurl','var') || isempty(baseurl) )
    %% Ancient "Legacy" THREDDS interface
    %% baseurl = 'http://tds.hycom.org/opendap/nph-dods/datasets/hycom/GOMl0.04/expt_20.1';

    %% Recent "Legacy" THREDDS interface
    %% baseurl = 'http://tds.hycom.org/thredds/dodsC/GOMl0.04/expt_20.1';
    %baseurl = 'http://tds.hycom.org/thredds/dodsC/GOMl0.04';

    baseurl = 'http://tds.hycom.org/thredds/dodsC/GOMl0.04/expt_32.5/hrly';
  end;

  % Maximum number of times to reopen and process a year-file with MDATASET
  if ( ~exist('maxNCRetries','var') || isempty(maxNCRetries) )
    maxNCRetries = 3;
  end;


  if ( ~isfield(stns,'lon') || any(cellfun(@isempty,{stns.lon})) || ~isfield(stns,'lat') || any(cellfun(@isempty,{stns.lat})) )
    badix = [];
    for ix = 1:numel(stns)
      if ( ~isfield(stns(ix),'lon') || isempty(stns(ix).lon) || ~isfield(stns(ix),'lat') || isempty(stns(ix).lat) )
        if ( isfield(stns(ix),'station_name') && ~isempty(stns(ix).station_name) )
          try
            [stns(ix).lon,stns(ix).lat,stns(ix).depth] = get_station_coords(stns(ix).station_name);
          catch
          end;
        end;
      end;
      if ( ~isfield(stns(ix),'lon') || isempty(stns(ix).lon) || ~isfield(stns(ix),'lat') || isempty(stns(ix).lat) )
        badix = [badix ix];
      end;
    end;
    if ( ~isempty(badix) )
      stns(badix) = [];
      if ( isempty(stns) )
        error('No station found with coordinates or name!');
      else
        warning('Ecoforecasts:HrlyGomHycom:MissingCoordinates',...
                '%g station(s) specified without coordinates or name!',...
                length(badix));
      end;
    end;
  end;


  needix = 1:numel(stns);

  for ix = 1:numel(stns)
    matfname = fullfile(datapath,[lower(stns(ix).station_name) '_hrly_gom_hycom.mat']);
    % If we already did this before - do not need to load again.
    if ( exist(matfname,'file') )
      disp(['Reloading MAT file ' matfname]);
      load(matfname,'station');
      stations(ix) = station;
      station = []; clear station;
      needix(needix == ix) = [];
    else
      if ( isfield(stns(ix),'station_name') )
        stations(ix).station_name = stns(ix).station_name;
      end;
      stations(ix).lon = stns(ix).lon; stations(ix).lat = stns(ix).lat;
      stations(ix).hrly_gom_hycom_interp_method = lower(interpMethod);
    end;
  end;

  if ( ~isempty(needix) )

    disp('Loading original data from OPeNDAP...');

    % NJTBX-2.0 Toolbox does not handle unclosed Datasets very well
    try
      nc = mDataset(baseurl);
      if ( ~isempty(nc) )
        dts = cast(getTimes(nc{vars{1}}),'double');
        % WORKAROUND for weird sub-second drift in netCDF Java Toolbox
        dts = round(dts);
        zerodt = dts(1);
        lastdt = dts(end);
      end;
      close(nc); clear nc;
    catch
      if ( exist('nc','var') && ~isempty(nc) )
        close(nc);
      end;
      clear nc;
    end;
    if ( ~exist('zerodt','var') || isempty(zerodt) || zerodt > maxdt || ...
         ~exist('lastdt','var') || isempty(lastdt) || lastdt < mindt )
      error('Ecoforecasts:HrlyGomHycom:BadURL',...
            'Model timestamps from "%s" outside range',baseurl);
    end;

    mindt = max(mindt,zerodt);
    maxdt = min(maxdt,lastdt);

    dtix = find(mindt <= dts & dts <= maxdt);
    alldts = dts(dtix);

    % Preallocate for fast execution
    for vix = 1:length(vars)
      fld = ['hrly_gom_hycom_' flds{vix}];
      if ( ~isfield(stations,fld) )
        for ix = needix(:)'
          stations(ix).(fld) = [];
          if ( ~isempty(strfind(fld,'_field')) )
            stations(ix).(fld).date = repmat(nan,[length(alldts),1]);
            stations(ix).(fld).field = repmat(nan,[length(alldts),(2*yrad)+1,(2*xrad)+1]);
          else
            stations(ix).(fld).date = repmat(nan,[length(alldts),1]);
            stations(ix).(fld).data = repmat(nan,[length(alldts),1]);
          end;
        end; %for ix = needix(:)'
      end;
    end;

    lats = [];
    lons = [];

    nc = mDataset(baseurl);
    if ( isempty(nc) )
      error('Ecoforecasts:HrlyGomHycom:InvalidURL',...
            'Invalid URL "%s"??!',baseurl);
    end;

    % NJTBX-2.0 Toolbox does not handle unclosed Datasets very well
    try

      if ( isempty(lats) || isempty(lons) )
        lats = cast(nc{'Latitude'}(:),'double');
        lons = cast(nc{'Longitude'}(:),'double');

        for ix = needix(:)'
          % [yix(ix),xix(ix)] = query_hrly_gom_hycom_indices(stations(ix).lon,stations(ix).lat);
          [yerr,yix(ix)] = min(abs(lats - stations(ix).lat));
          [xerr,xix(ix)] = min(abs(lons - stations(ix).lon));
          if ( yerr > max(diff(unique(lats))) || xerr > max(diff(unique(lons))) )
            warning('Ecoforecasts:HrlyGomHycom:StationOutsideDomain',...
                    'Station "%s" at %g,%g is outside GoM HYCOM domain!',...
                    stations(ix).station_name,stations(ix).lat,stations(ix).lon);
            needix(needix == ix) = [];
          else
            lat{ix} = lats(yix(ix)-2:yix(ix)+2);
            lon{ix} = lons(xix(ix)-3:xix(ix)+3);
            stations(ix).hrly_gom_hycom_latix = yix(ix);
            stations(ix).hrly_gom_hycom_lonix = xix(ix);
          end;
        end;
      end; %if ( isempty(lats) || isempty(lons) )

      for vix = 1:length(vars)
        var = vars{vix};
        fld = ['hrly_gom_hycom_' flds{vix}];
        %DEBUG:        disp(fld);

        datsz = getShape(nc{var});

        for ix = needix(:)'
          %DEBUG:
          disp([stations(ix).station_name,'.',fld]);
          %DEBUG:tic,
          stations(ix).(fld).date(:,1) = alldts';

          if ( ~isempty(strfind(fld,'_field')) )
            if ( ~isfield(stations(ix).(fld),'lat') || isempty(stations(ix).(fld).lat) || ...
                 ~isfield(stations(ix).(fld),'lon') || isempty(stations(ix).(fld).lon) )
              stations(ix).(fld).lat = lats(yix(ix)-yrad:yix(ix)+yrad);
              stations(ix).(fld).lon = lons(xix(ix)-xrad:xix(ix)+xrad);
            end;

            %DEBUG:
keyboard;

            % 2-D data elements
            if ( length(datsz) == 3 )
              stations(ix).(fld).field(:,:,:) = ...
                  squeeze(cast(nc{var}(dtix,yix(ix)-yrad:yix(ix)+yrad,xix(ix)-xrad:xix(ix)+xrad),'double'));
            % 3-D data elements - get the first (highest) Z index
            else
              stations(ix).(fld).field(:,:,:) = ...
                  squeeze(cast(nc{var}(dtix,1,yix(ix)-yrad:yix(ix)+yrad,xix(ix)-xrad:xix(ix)+xrad),'double'));
            end;

          else

            %DEBUG:
keyboard;

            % 2-D data elements
            if ( length(datsz) == 3 )
              dat = squeeze(cast(nc{var}(dtix,yix(ix)-2:yix(ix)+2,xix(ix)-3:xix(ix)+3),'double'));
            % 3-D data elements - get the first (highest) Z index
            else
              dat = squeeze(cast(nc{var}(dtix,1,yix(ix)-2:yix(ix)+2,xix(ix)-3:xix(ix)+3),'double'));
            end;
            if ( isempty(dat) )
              warning('Ecoforecasts:HrlyGomHycom:EmptyArray',...
                      'Received empty array! Station %s, field %s, var %s, URL "%s"',...
                      stations(ix).station_name,fld,var,baseurl);
              %DEBUG:
keyboard;
            else
              % Avoid more X-Y vs. row-column confusion: stupid MATLAB...
              stations(ix).(fld).data(:,1) = ...
                  interp_field(lat{ix},lon{ix},dat,stations(ix).lat,stations(ix).lon,interpMethod,'warn');
            end;
            dat = []; clear dat;

          end;

          %DEBUG:toc,
          % Give the THREDDS server some time to chillax
          % pause(2);
          pause(0.1);
        end; %for ix = needix(:)'

      end; %for vix = 1:length(vars)

    catch
      if ( exist('nc','var') && ~isempty(nc) )
        close(nc);
      end;
      clear nc;
      rethrow(lasterror);
    end;

    if ( exist('nc','var') && ~isempty(nc) )
      close(nc);
    end;
    clear nc;

    %DEBUG:
    toc,


    for ix = needix(:)'
      % Calculate speed and direction from model U and V currents
      if ( isfield(stations(ix),'hrly_gom_hycom_u') && isfield(stations(ix),'hrly_gom_hycom_v') )
        stations(ix).hrly_gom_hycom_speed.date = stations(ix).hrly_gom_hycom_u.date;
        stations(ix).hrly_gom_hycom_speed.data = uv_to_spd(stations(ix).hrly_gom_hycom_u.data,stations(ix).hrly_gom_hycom_v.data);
        stations(ix).hrly_gom_hycom_dir.date = stations(ix).hrly_gom_hycom_u.date;
        stations(ix).hrly_gom_hycom_dir.data = uv_to_dir_curr(stations(ix).hrly_gom_hycom_u.data,stations(ix).hrly_gom_hycom_v.data);
      end;

      matfname = fullfile(datapath,[lower(stations(ix).station_name) '_hrly_gom_hycom.mat']);
      disp(['Saving to MAT file ' matfname]);
      station = stations(ix);
      save(matfname,'station');
      station = []; clear station;
    end; %for ix = needix(:)'

  end; %if ( ~isempty(needix) )


  % Gross quality control - did we get all the fields we expected?
  for vix = 1:length(vars)
    fld = ['hrly_gom_hycom_' flds{vix}];
    for ix = 1:numel(stations)
      if ( ~isfield(stations(ix),fld) || isempty(stations(ix).(fld)) )
        warning('Ecoforecasts:HrlyGomHycom:MissingField',...
                'No field "%s" found after load!',fld);
      end;
    end; %for ix = 1:numel(stations)
  end; %for vix = 1:length(vars)

  % Whether we loaded data, extracted data from netCDF files, or both, limit
  % our return to only dates the user requested - that have reasonable data
  allflds = fieldnames(stations);
  for fldix = 1:length(allflds)
    fld = allflds{fldix};

    for ix = 1:numel(stations)
      % Slightly finer quality control
      if ( isfield(stations(ix).(fld),'field') )
        % Make sure lon and lat are row vectors as for other models
        stations(ix).(fld).lon = stations(ix).(fld).lon(:);
        stations(ix).(fld).lat = stations(ix).(fld).lat(:);

        % With qtot in the mix, we are not enforcing tight bounds on any var
        ixes = find(-3000 >= stations(ix).(fld).field | stations(ix).(fld).field >= 3000);
        stations(ix).(fld).field(ixes) = nan;

        % ixes = find(mindt <= stations(ix).(fld).date & stations(ix).(fld).date <= maxdt);
        % stations(ix).(fld).date = stations(ix).(fld).date(ixes);
        % stations(ix).(fld).field = stations(ix).(fld).field(ixes,:,:);

      elseif ( isfield(stations(ix).(fld),'data') )
        % With qtot in the mix, we are not enforcing tight bounds on any var
        ixes = find(-3000 <= stations(ix).(fld).data & stations(ix).(fld).data <= 3000);
        stations(ix).(fld).date = stations(ix).(fld).date(ixes);
        stations(ix).(fld).data = stations(ix).(fld).data(ixes);

        % ixes = find(mindt <= stations(ix).(fld).date & stations(ix).(fld).date <= maxdt);
        % stations(ix).(fld).date = stations(ix).(fld).date(ixes);
        % stations(ix).(fld).data = stations(ix).(fld).data(ixes);
      end;
    end; %for ix = 1:numel(stations)
  end; %for fldix = 1:length(allflds)

  % If caller wants different interpolation method than we saved in MAT file,
  % re-call INTERP_FIELD on all the appropriate "_field" fields (if present)
  for ix = 1:numel(stations)
    if ( ~isfield(stations(ix),'hrly_gom_hycom_interp_method') || ...
         ~strcmpi(stations(ix).hrly_gom_hycom_interp_method,interpMethod) )
      disp(['Reinterpolating time series for ' stations(ix).station_name ...
            ' using ' upper(interpMethod)]);

      for vix = 1:length(vars)
        fld = ['hrly_gom_hycom_' flds{vix}];
        if ( isfield(stations(ix).(fld),'data') )
          fnm = [fld '_field'];
          if ( ~isfield(stations(ix),fnm) )
            warning('Ecoforecasts:HrlyGomHycom:NoFieldToReinterp',...
                    '%s: No field available to reinterpolate "%s"',...
                    stations(ix).station_name,fld);
          else
            f = stations(ix).(fnm);
            stations(ix).(fld).data = ...
                interp_field(f.lat,f.lon,f.field,stations(ix).lat,stations(ix).lon,interpMethod,'warn');
          end;
        end;
      end; %for vix = 1:length(vars)

      % Recalculate speed/dir time series after (presumably) reinterpolating U/V
      if ( isfield(stations(ix),'hrly_gom_hycom_u') && isfield(stations(ix),'hrly_gom_hycom_v') )
        stations(ix).hrly_gom_hycom_speed.data = uv_to_spd(stations(ix).hrly_gom_hycom_u.data,stations(ix).hrly_gom_hycom_v.data);
        stations(ix).hrly_gom_hycom_dir.data = uv_to_dir_curr(stations(ix).hrly_gom_hycom_u.data,stations(ix).hrly_gom_hycom_v.data);
      end;

    end; %if ~strcmpi(stations(ix).hrly_gom_hycom_interp_method,interpMethod)
  end; %for ix = 1:numel(stations)

  for ix = 1:numel(stations)
    allflds = grepstruct(stations(ix),'hrly_gom_hycom');
    for fldix = 1:length(allflds)
      fld = allflds{fldix};
      stns(ix).(fld) = stations(ix).(fld);
      % Keep memory growth in check for a big multi-station request
      stations(ix).(fld) = [];
    end;
  end;
  stations = []; clear stations;

  set_more;

return;
