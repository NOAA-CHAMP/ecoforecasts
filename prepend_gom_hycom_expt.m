function stns = prepend_gom_hycom_expt(stns_or_stnms,mindt,maxdt,expt,vars,flds,interpMethod,baseurl)
error('NOT YET READY FOR PRIME-TIME!');
%function stns = prepend_gom_hycom_expt(stns_or_stnms,mindt,maxdt,expt,vars,flds,interpMethod,baseurl)
%
% NOTE: This version of 'GET_GOM_HYCOM' (v.) is essentially a code fork that
%  is simplified by only having to look at the MOST RECENT online experiment.
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
% NOTE: The string 'gom_hycom_' is prepended to each of the names in FLDS.
%
% DEFAULT for BASEURL:
%  baseurl = 'http://tds.hycom.org/thredds/dodsC/GOMl0.04';
% If you have another ocean model with a THREDDS interface, specify it here.
%
% DEFAULT for EXPT:
%  expt = 'expt_32.5';
%
% Information here: http://tds.hycom.org/thredds/dodsC/GOMl0.04/expt_32.5.html
%
% HISTORICAL NOTE: Files are stored by full year, experiment 20.1 ended July 2010:
%  http://tds.hycom.org/thredds/dodsC/GOMl0.04/expt_20.1/2003
%   . . .
%  http://tds.hycom.org/thredds/dodsC/GOMl0.04/expt_20.1/2010
%  http://tds.hycom.org/thredds/dodsC/GOMl0.04/expt_30.1/2010
%   . . .
%  http://tds.hycom.org/thredds/dodsC/GOMl0.04/expt_30.1/2014
%  http://tds.hycom.org/thredds/dodsC/GOMl0.04/expt_32.5/2014
%   . . .
%
%
% CALLS: MDATASET (netCDF-Java); INTERP_FIELD, GET_STATION_COORDS (Ecoforecasts).
%
% Last Saved Time-stamp: <Mon 2017-12-04 15:04:32 Eastern Standard Time gramer>

  set_more off;

  datapath = get_thesis_path('../data');

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
    vars = { 'u', 'v', 'mld', 'qtot',           'temperature',   'salinity', ...
             'u',       'v',       'mld',       'temperature',   'salinity'       };
    flds = { 'u', 'v', 'mld', 'net_heat_flux',  'seatemp',       'salinity', ...
             'u_field', 'v_field', 'mld_field', 'seatemp_field', 'salinity_field' };
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
    baseurl = 'http://tds.hycom.org/thredds/dodsC/GOMl0.04';
  end;

  if ( ~exist('expt','var') || isempty(expt) )
    expt = 'expt_32.5';
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
        warning('Ecoforecasts:GomHycom:MissingCoordinates',...
                '%g station(s) specified without coordinates or name!',...
                length(badix));
      end;
    end;
  end;


  needix = 1:numel(stns);

  for ix = 1:numel(stns)
    stns(ix).matfname = fullfile(datapath,[lower(stns(ix).station_name) '_gom_hycom_',expt,'.mat']);

    % If we already did this before - do not need to load again.
    if ( exist(stns(ix).matfname,'file') )
      disp(['Reloading MAT file ',stns(ix).matfname]);
      load(stns(ix).matfname,'station');
      stations(ix) = station;
      station = []; clear station;
%      needix(needix == ix) = [];
    else
      error('Nothing to prepend to??? MISSING "%s"!',stns(ix).matfname);
      % if ( isfield(stns(ix),'station_name') )
      %   stations(ix).station_name = stns(ix).station_name;
      % end;
      % stations(ix).lon = stns(ix).lon; stations(ix).lat = stns(ix).lat;
      % stations(ix).gom_hycom_interp_method = lower(interpMethod);
    end;
  end;

  if ( ~isempty(needix) )

    disp(['Loading original data from ',baseurl]);

    % Grid-point radii for time-series fields (e.g., seatemp_field)
    % Choose one axis slightly larger to guard against transposition errors
    yrad = 4;
    xrad = 5;

    % What are the first and last days with available model data?
    %zerodt = datenum(2003,1,2);
    zerodt = datenum(2014,4,2);

    % lastdt = floor(now);		%%% Most recent day of experiment 
    %%% Nope - we actually have to query the most recent file
    lastdt = [];
    yearTried = 0;
    while ( isempty(lastdt) )
      % If we're running this in January - use last year instead
      url = sprintf('%s/%s/%04d',baseurl,expt,get_year(now-32));
      % NJTBX-2.0 Toolbox does not handle unclosed Datasets very well
      try
        nc = mDataset(url);
        if ( ~isempty(nc) )
          dts = cast(getTimes(nc{vars{1}}),'double');
          % WORKAROUND for weird sub-second drift in netCDF Java Toolbox
          dts = round(dts);
          lastdt = dts(end);
        end;
        close(nc); clear nc;
      catch
        if ( exist('nc','var') && ~isempty(nc) )
          close(nc);
        end;
        clear nc;
      end;

      if ( isempty(lastdt) )
        yearTried = yearTried + 1;
        if ( yearTried < maxNCRetries )
          warning('Ecoforecasts:GomHycom:RetryingLatestYear',...
                  'Unable to extract latest model date: retrying URL "%s"...',url);
          pause(10);
        else
          warning('Ecoforecasts:GomHycom:InvalidURL',...
                  'Retried URL %g times - ending on previous year instead!',...
                  maxNCRetries);
          lastdt = datenum(get_year(now)-1,12,31);
        end;
      end;
    end;

    % Model outputs are once-daily time series - assume midnight GMT
    alldts = zerodt:lastdt;

    % First time ever calling this function, we want to get ALL the data
    [begyr,ig,ig] = datevec(zerodt);
    [endyr,ig,ig] = datevec(lastdt);
    allyrs = [begyr:endyr];

    % % Preallocate for fast execution
    % for vix = 1:length(vars)
    %   fld = ['gom_hycom_' flds{vix}];
    %   if ( ~isfield(stations,fld) )
    %     for ix = needix(:)'
    %       stations(ix).(fld) = [];
    %       if ( ~isempty(strfind(fld,'_field')) )
    %         stations(ix).(fld).date = repmat(nan,[length(alldts),1]);
    %         stations(ix).(fld).field = repmat(nan,[length(alldts),(2*yrad)+1,(2*xrad)+1]);
    %       else
    %         stations(ix).(fld).date = repmat(nan,[length(alldts),1]);
    %         stations(ix).(fld).data = repmat(nan,[length(alldts),1]);
    %       end;
    %     end; %for ix = needix(:)'
    %   end;
    % end;

    lats = [];
    lons = [];

    for yr=allyrs(:)'
     %DEBUG:
     disp(yr);
     %DEBUG:
     tic,

     dts = [];

     % expt = 'expt_32.5';
     if ( yr < 2014 )
       error('NO DATA IN "%s/%s" before year 2014!',baseurl,expt);
     end;

     yearDone = false;
     yearTried = 0;

     while ( ~yearDone )

      url = sprintf('%s/%s/%04d',baseurl,expt,yr);
      %DEBUG:
      disp(url);
      nc = mDataset(url);
      if ( isempty(nc) )
        warning('Ecoforecasts:GomHycom:InvalidURL',...
                'Invalid URL "%s"!',url);
        continue;
      end;

      % NJTBX-2.0 Toolbox does not handle unclosed Datasets very well
      try

        if ( isempty(lats) || isempty(lons) )
          lats = cast(nc{'Latitude'}(:),'double');
          lons = cast(nc{'Longitude'}(:),'double');

          for ix = needix(:)'
            % [yix(ix),xix(ix)] = query_gom_hycom_indices(stations(ix).lon,stations(ix).lat);
            [yerr,yix(ix)] = min(abs(lats - stations(ix).lat));
            [xerr,xix(ix)] = min(abs(lons - stations(ix).lon));
            if ( yerr > max(diff(unique(lats))) || xerr > max(diff(unique(lons))) )
              warning('Ecoforecasts:GomHycom:StationOutsideDomain',...
                      'Station "%s" at %g,%g is outside GoM HYCOM domain!',...
                      stations(ix).station_name,stations(ix).lat,stations(ix).lon);
              needix(needix == ix) = [];
            else
              lat{ix} = lats(yix(ix)-2:yix(ix)+2);
              lon{ix} = lons(xix(ix)-3:xix(ix)+3);
              stations(ix).gom_hycom_latix = yix(ix);
              stations(ix).gom_hycom_lonix = xix(ix);
            end;
          end;
        end; %if isempty(lats)


        for vix = 1:length(vars)
          var = vars{vix};
          fld = ['gom_hycom_' flds{vix}];
          %DEBUG:
          disp(fld);

          if ( isempty(dts) )
            dts = cast(getTimes(nc{var}),'double');
            % WORKAROUND for weird sub-second drift in netCDF Java Toolbox
            dts = round(dts);
            ndts = numel(dts);

            dtix = find(ismember(alldts,dts));
            if ( length(dtix) ~= ndts )
%keyboard;
              close(nc); clear nc;
              error('Date count mismatch %d vs %d! URL "%s"',length(dtix),ndts,url);
            end;
          end;

          datsz = getShape(nc{var});

          for ix = needix(:)'
            % PREPEND!
            stations(ix).(fld).date = [dts' ; stations(ix).(fld).date];

            %DEBUG:
            disp([stations(ix).station_name,' <-- ',var]);

            if ( ~isempty(strfind(fld,'_field')) )
              if ( ~isfield(stations(ix).(fld),'lat') || isempty(stations(ix).(fld).lat) || ...
                   ~isfield(stations(ix).(fld),'lon') || isempty(stations(ix).(fld).lon) )
                stations(ix).(fld).lat = lats(yix(ix)-yrad:yix(ix)+yrad);
                stations(ix).(fld).lon = lons(xix(ix)-xrad:xix(ix)+xrad);
              end;

              % 2-D data elements
              if ( length(datsz) == 3 )
                newfld(dtix,:,:) = ...
                    squeeze(cast(nc{var}(:,yix(ix)-yrad:yix(ix)+yrad,xix(ix)-xrad:xix(ix)+xrad),'double'));
              % 3-D data elements - get the first (highest) Z index
              else
                newfld(dtix,:,:) = ...
                    squeeze(cast(nc{var}(:,1,yix(ix)-yrad:yix(ix)+yrad,xix(ix)-xrad:xix(ix)+xrad),'double'));
              end;
              % PREPEND!
              olddat = stations(ix).(fld).field;
              stations(ix).(fld).field(numel(dtix)+1:numel(dtix)+size( = 

            else

              % 2-D data elements
              if ( length(datsz) == 3 )
                dat = squeeze(cast(nc{var}(:,yix(ix)-2:yix(ix)+2,xix(ix)-3:xix(ix)+3),'double'));
                % 3-D data elements - get the first (highest) Z index
              else
                dat = squeeze(cast(nc{var}(:,1,yix(ix)-2:yix(ix)+2,xix(ix)-3:xix(ix)+3),'double'));
              end;
              if ( isempty(dat) )
%keyboard;
                close(nc); clear nc;
                error('Received empty array! Field %s, var %s, file "%s"',fld,var,url);
              end;
              % Avoid more X-Y vs. row-column confusion: stupid MATLAB...
              stations(ix).(fld).data(dtix,1) = ...
              newdat = interp_field(lat{ix},lon{ix},dat,stations(ix).lat,stations(ix).lon,interpMethod,'warn');
              % PREPEND!
              stations(ix).(fld).data = [newdat ; stations(ix).(fld).data];
              dat = []; clear dat;

            end;

            % Give the THREDDS server some time to chillax
            pause(2);
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

      % Give the THREDDS server a second to chillax
      pause(1);

      % Make sure we got >= 90% of the data we expected, or else try again!
      N = +Inf;
      targetN = length(dtix);
      for vix = 1:length(vars)
        var = vars{vix};
        fld = ['gom_hycom_' flds{vix}];
        if ( ~isfield(stations(1),fld) )
          N = 0;
        else
        for ix = needix(:)'
          if ( isfield(stations(ix).(fld),'data') )
            N = min(N,length(isfinite(stations(ix).(fld).data(dtix))));
          elseif ( isfield(stations(ix).(fld),'field') )
            dat = nanmean(nanmean(stations(ix).(fld).field(dtix,:,:),3),2);
            N = min(N,length(isfinite(dat)));
          end;
        end;
      end;
      end;
      if ( N >= 0.90*targetN )
        yearDone = true;
      else
        yearTried = yearTried + 1;
        if ( yearTried < maxNCRetries )
          warning('Ecoforecasts:GomHycom:RetryingYear',...
                  'Expected %g values for %s, got %g: retrying year %g...',...
                  targetN,var,N,yr);
          pause(60);
        else
          error('Retried %g times - giving up! "%s"',maxNCRetries,url);
        end;
      end; %if ( N >= 0.90*targetN ) else

     end; %while ( ~yearDone )

     % for ix = needix(:)'
     %   % Calculate speed and direction from model U and V currents
     %   if ( isfield(stations(ix),'gom_hycom_u') && isfield(stations(ix),'gom_hycom_v') )
     %     stations(ix).gom_hycom_speed.date = stations(ix).gom_hycom_u.date;
     %     stations(ix).gom_hycom_speed.data = uv_to_spd(stations(ix).gom_hycom_u.data,stations(ix).gom_hycom_v.data);
     %     stations(ix).gom_hycom_dir.date = stations(ix).gom_hycom_u.date;
     %     stations(ix).gom_hycom_dir.data = uv_to_dir_curr(stations(ix).gom_hycom_u.data,stations(ix).gom_hycom_v.data);
     %   end;
     %
     %   disp(['Precautionary PRE-SAVE to MAT file ',stns(ix).matfname,' after year ',num2str(yr)]);
     %   station = stations(ix);
     %   save(stns(ix).matfname,'station');
     %   station = []; clear station;
     % end; %for ix = needix(:)'

    end; %for yr=allyrs(:)'

    % ALSO COPYING THIS WITHIN EACH YEAR LOOP ABOVE AS A PRECAUTION
    for ix = needix(:)'
      % Calculate speed and direction from model U and V currents
      if ( isfield(stations(ix),'gom_hycom_u') && isfield(stations(ix),'gom_hycom_v') )
        stations(ix).gom_hycom_speed.date = stations(ix).gom_hycom_u.date;
        stations(ix).gom_hycom_speed.data = uv_to_spd(stations(ix).gom_hycom_u.data,stations(ix).gom_hycom_v.data);
        stations(ix).gom_hycom_dir.date = stations(ix).gom_hycom_u.date;
        stations(ix).gom_hycom_dir.data = uv_to_dir_curr(stations(ix).gom_hycom_u.data,stations(ix).gom_hycom_v.data);
      end;

      disp(['Saving to MAT file ',stns(ix).matfname]);
      station = stations(ix);
      save(stns(ix).matfname,'station');
      station = []; clear station;
    end; %for ix = needix(:)'

  end; %if ( ~isempty(needix) )


  % Gross quality control - did we get all the fields we expected?
  for vix = 1:length(vars)
    fld = ['gom_hycom_' flds{vix}];
    for ix = 1:numel(stations)
      if ( ~isfield(stations(ix),fld) || isempty(stations(ix).(fld)) )
        warning('Ecoforecasts:GomHycom:MissingField',...
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

        ixes = find(mindt <= stations(ix).(fld).date & stations(ix).(fld).date <= maxdt);
        stations(ix).(fld).date = stations(ix).(fld).date(ixes);
        stations(ix).(fld).field = stations(ix).(fld).field(ixes,:,:);

      elseif ( isfield(stations(ix).(fld),'data') )
        % With qtot in the mix, we are not enforcing tight bounds on any var
        ixes = find(-3000 <= stations(ix).(fld).data & stations(ix).(fld).data <= 3000);
        stations(ix).(fld).date = stations(ix).(fld).date(ixes);
        stations(ix).(fld).data = stations(ix).(fld).data(ixes);

        ixes = find(mindt <= stations(ix).(fld).date & stations(ix).(fld).date <= maxdt);
        stations(ix).(fld).date = stations(ix).(fld).date(ixes);
        stations(ix).(fld).data = stations(ix).(fld).data(ixes);
      end;
    end; %for ix = 1:numel(stations)
  end; %for fldix = 1:length(allflds)

  % If caller wants different interpolation method than we saved in MAT file,
  % re-call INTERP_FIELD on all the appropriate "_field" fields (if present)
  for ix = 1:numel(stations)
    if ( ~isfield(stations(ix),'gom_hycom_interp_method') || ...
         ~strcmpi(stations(ix).gom_hycom_interp_method,interpMethod) )
      disp(['Reinterpolating time series for ' stations(ix).station_name ...
            ' using ' upper(interpMethod)]);

      for vix = 1:length(vars)
        fld = ['gom_hycom_' flds{vix}];
        if ( isfield(stations(ix).(fld),'data') )
          fnm = [fld '_field'];
          if ( ~isfield(stations(ix),fnm) )
            warning('Ecoforecasts:GomHycom:NoFieldToReinterp',...
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
      if ( isfield(stations(ix),'gom_hycom_u') && isfield(stations(ix),'gom_hycom_v') )
        stations(ix).gom_hycom_speed.data = uv_to_spd(stations(ix).gom_hycom_u.data,stations(ix).gom_hycom_v.data);
        stations(ix).gom_hycom_dir.data = uv_to_dir_curr(stations(ix).gom_hycom_u.data,stations(ix).gom_hycom_v.data);
      end;

    end; %if ~strcmpi(stations(ix).gom_hycom_interp_method,interpMethod)
  end; %for ix = 1:numel(stations)

  for ix = 1:numel(stations)
    allflds = grepstruct(stations(ix),'gom_hycom');
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
