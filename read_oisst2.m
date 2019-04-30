function stns = read_oisst2(stns_or_stnms,mindt,maxdt,vars,flds,interpMethod,baseurl)
%function stns = read_oisst2(stns_or_stnms,mindt,maxdt,vars,flds,interpMethod,baseurl)
%
% Add fields from Optimal Interpolated Sea Surface Temperature V2 ("OISST2")
% to each member of struct array (or scalar) STNS. If a station name string
% or cell array of strings STNMS (five characters, e.g., 'mlrf1', or
% {'smkf1','tnrf1'}) is given instead of structs, or if any struct lacks a
% valid .lon or .lat field (e.g., -80.38, 25.01) but all do have valid
% .station_name fields, tries to retrieve locations with GET_STATION_COORDS.
%
% Documentation for OISST2 was accessed online as of 2011 Apr 04 from:
%  http://www.ncdc.noaa.gov/oa/climate/research/sst/oi-daily.php
%
% Optional VARS is a string or cellstr specifying which variables to extract.
% If VARS is specified, optional FLDS may also be specified as the corresp.
% list of field names to be added to structs in STNS, e.g., 'sst'.
%
% INTERP_FIELD is passed INTERPMETHOD (DEFAULT 'linear') to subset points.
%
% DEFAULT VARS and corresponding FLDS cell arrays are as follows:
%    vars = { 'sst',      'anom',      'err', ...
%             'sst', };
%    flds = { 'sst',       'sst_anom', 'sst_err', ...
%             'sst_field', };
%
% NOTE: String 'oisst2_' is prepended to each of the names in FLDS.
%
% Sample URL:
% http://nomads.ncdc.noaa.gov/thredds/dodsC/oisst2/2011/AVHRR-AMSR/amsr-avhrr-v2.20110101.nc
%
% DEFAULT for BASEURL:
%  'http://nomads.ncdc.noaa.gov/thredds/dodsC/oisst2'
%
% CALLS: MDATASET (netCDF-Java); INTERP_FIELD, GET_STATION_COORDS (Ecoforecasts).
%
% Last Saved Time-stamp: <Thu 2011-04-14 10:23:39  Lew.Gramer>

  set_more off;

  datapath = get_ecoforecasts_path('data');

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
  else
    error('First arg must be station STRUCT(s) or station name string(s)');
  end;
  clear stns_or_stnms;

  if ( ~isfield(stns,'station_name') || any(cellfun(@isempty,{stns.station_name})) )
    if ( ~(~isfield(stns,'name') || any(cellfun(@isempty,{stns.name}))) )
      for ix=1:numel(stns)
        stns(ix).station_name = stns(ix).name;
      end;
    else
      error('Station struct(s) passed in without .station_name field!');
    end;
  end;
  if ( ~isfield(stns,'lon') || ~isfield(stns,'lat') )
    for ix=1:numel(stns)
      stns(ix).lon = [];
      stns(ix).lat = [];
    end;
  end;


  if ( ~exist('mindt','var') || isempty(mindt) )
    mindt = -Inf;
  end;
  if ( ~exist('maxdt','var') || isempty(maxdt) )
    maxdt = +Inf;
  end;


  % Variable name prefix: Daily Optimal Interpolated SST v2 (Reynolds et al)
  PFX = 'oisst2_';

  if ( ~exist('vars','var') || isempty(vars) )
    vars = { 'sst',      'anom',      'err', ...
             'sst', };
  end;
  if ( ~exist('flds','var') || isempty(flds) )
    flds = { 'sst',       'sst_anom', 'sst_err', ...
             'sst_field', };
  end;

  if ( ~exist('interpMethod','var') || isempty(interpMethod) )
    interpMethod = 'linear';
  end;
  if ( ischar(interpMethod) && interpMethod(1) == '*' )
    interpMethod = interpMethod(2:end);
  end;


  if ( ~exist('baseurl','var') || isempty(baseurl) )
    baseurl = 'http://nomads.ncdc.noaa.gov/thredds/dodsC/oisst2';
  end;

  % % Full dataset 1981-present
  % yrs = 1981:get_year(now);

  % yrs = 1988:get_year(now);
  %%%% MODIFY THE FOLLOWING LINE TO CONTINUE A PARTIAL EXTRACTION
  yrs = 1988:get_year(now);

  %DEBUG:  yrs = 2009:2010;

  alldts = datenum(yrs(1),1,1):floor(now);

  % Grid-point radii for time-series fields (e.g., seatemp_field)
  % Choose one axis slightly larger to guard against transposition errors
  % Dyslexics of the world untie
  yrad = 2;
  xrad = 3;


  needix = 1:numel(stns);

  for ix=1:numel(stns)
    matfname = fullfile(datapath,[lower(stns(ix).station_name) '_oisst2.mat']);
    % If we already did this before - do not need to load again.
    if ( exist(matfname,'file') )
      disp(['Reloading MAT file ' matfname]);
      load(matfname,'station');
      stations(ix) = station;
      station = []; clear station;
      if ( isempty(stns(ix).lon) || isempty(stns(ix).lat) )
        stns(ix).lon = stations(ix).lon; stns(ix).lat = stations(ix).lat;
      end;
      % Take this station out of the firing line.
      %%%% COMMENT OUT THE FOLLOWING LINE TO CONTINUE A PARTIAL EXTRACTION
      needix(needix == ix) = [];

    else
      if ( isempty(stns(ix).lon) || isempty(stns(ix).lat) )
        [stns(ix).lon,stns(ix).lat,stns(ix).depth] = get_station_coords(stns(ix).station_name);
      end;

      stations(ix).station_name = stns(ix).station_name;
      stations(ix).lon = stns(ix).lon; stations(ix).lat = stns(ix).lat;
      stations(ix).([PFX 'interp_method']) = lower(char(interpMethod));
      for vix=1:numel(vars)
        fld = [PFX flds{vix}];
        % Preallocate for fast execution
        stations(ix).(fld) = [];
        if ( ~isempty(strfind(fld,'_field')) )
          stations(ix).(fld).date = repmat(nan,[length(alldts),1]);
          stations(ix).(fld).field = repmat(nan,[length(alldts),(2*yrad)+1,(2*xrad)+1]);
        else
          stations(ix).(fld).date = repmat(nan,[length(alldts),1]);
          stations(ix).(fld).data = repmat(nan,[length(alldts),1]);
        end;
      end;
    end;
  end;


  if ( ~isempty(needix) )

    disp('Loading original data from online HDF files...');

    lats = [];
    lons = [];

    maxRetries = 3;

    for yr=yrs(:)'
      %DEBUG:
      disp(yr);      tic,

      jds = 1:365;
      if ( mod(yr,4) == 0 )
        jds = 1:366;
      end;

      for jd=jds(:)'

        tryCount = 1;
        dayDone = false;

        while ( ~dayDone )
          
          dt = datenum(yr,1,0) + jd;
          [ig,mn,dy] = datevec(dt);
          [ig,dtix] = min(abs(alldts - dt));

          %Sample URL:
          %http://nomads.ncdc.noaa.gov/thredds/dodsC/oisst2/2002/AVHRR/avhrr-only-v2.20021231.nc
          %http://nomads.ncdc.noaa.gov/thredds/dodsC/oisst2/2011/AVHRR-AMSR/amsr-avhrr-v2.20110101.nc
          dataset = 'AVHRR-AMSR'; databasebasename = 'amsr-avhrr-v2'; 
          if ( dt < datenum(2002,6,1) )
            dataset = 'AVHRR'; databasebasename = 'avhrr-only-v2'; 
          end;

          url = sprintf('%s/%04d/%s/%s.%04d%02d%02d.nc',...
                        baseurl,yr,dataset,databasebasename,yr,mn,dy);
          %DEBUG:          disp(url);
          dat = [];
          nc = mDataset(url);
          if ( isempty(nc) )
            warning('Ecoforecasts:OISST:BadURL',...
                    'Bad URL "%s"',url);
          else
            try
              if ( isempty(lats) || isempty(lons) )
                lons = cast(nc{'lon'}(:),'double');
                lats = cast(nc{'lat'}(:),'double');
                lons(lons>180) = lons(lons>180) - 360;

                for ix = needix(:)'
                  [yerr,yix(ix)] = min(abs(lats - stations(ix).lat));
                  [xerr,xix(ix)] = min(abs(lons - stations(ix).lon));
                  if ( yerr > max(diff(unique(lats))) || xerr > max(diff(unique(lons))) )
                    % Global product - this should never happen!
                    warning('Ecoforecasts:OISST:StationOutsideDomain',...
                            'Station "%s" at %g,%g is outside OISST domain?!',...
                            stations(ix).station_name,stations(ix).lat,stations(ix).lon);
                    needix(needix == ix) = [];
                  else
                    lat{ix} = lats(yix(ix)-1:yix(ix)+1);
                    lon{ix} = lons(xix(ix)-2:xix(ix)+2);

                    stations(ix).([PFX 'lonix']) = xix(ix);
                    stations(ix).([PFX 'latix']) = yix(ix);
                    for vix=1:numel(vars)
                      var = vars{vix};
                      fld = [PFX flds{vix}];
                      if ( ~isempty(strfind(fld,'_field')) )
                        stations(ix).(fld).lon = lons(xix(ix)-xrad:xix(ix)+xrad);
                        stations(ix).(fld).lat = lats(yix(ix)-yrad:yix(ix)+yrad);
                      end;
                    end;
                  end;
                end;
              end; %if ( isempty(lats) || isempty(lons) )

              for vix=1:numel(vars)
                var = vars{vix};
                fld = [PFX flds{vix}];
                %DEBUG:            disp(fld);
                for ix=needix(:)'
                  dat = [];
                  stations(ix).(fld).date(dtix,1) = dt;
                  if ( ~isempty(strfind(fld,'_field')) )
                    dat = cast(nc{var}(:,:,yix(ix)-yrad:yix(ix)+yrad,xix(ix)-xrad:xix(ix)+xrad),'double');
                    stations(ix).(fld).field(dtix,:,:) = dat;
                  else
                    dat = cast(nc{var}(:,:,yix(ix)-1:yix(ix)+1,xix(ix)-2:xix(ix)+2),'double');
                    % stations(ix).(fld).data(dtix,1) = ...
                    %     interp2(lon{ix},lat{ix},dat,stations(ix).lon,stations(ix).lat,interpMethod);
                    stations(ix).(fld).data(dtix,1) = ...
                        interp_field(lat{ix},lon{ix},dat,stations(ix).lat,stations(ix).lon,interpMethod);
                  end; %if ( ~isempty(strfind(fld,'_field')) )
                end; %for ix=needix(:)'
              end; %for vix=1:numel(vars)
            catch
              if ( exist('nc','var') && ~isempty(nc) ); close(nc); end;
              clear nc;
              % rethrow(lasterror);
              lerr = lasterror;
              msg = 'CAUGHT ERROR!';
              if ( isfield(lerr,'identifier') )
                msg = [ msg ' ' lerr.identifier ];
              end;
              if ( isfield(lerr,'message') )
                msg = [ msg ' : ' lerr.message ];
              end;
              warning(msg);
              dat = [];
            end; %try catch
          end; %if ( isempty(nc) ) else

          if ( exist('nc','var') && ~isempty(nc) ); close(nc); end;
          clear nc;

          % Give remote server a breather
          pause(0.1);

          dayDone = true;
          % If something went wrong with extracting this file
          if ( isempty(dat) )
            dayDone = false;
            tryCount = tryCount + 1;
            if ( tryCount <= maxRetries )
              warning('Ecoforecasts:OISST:NCError',...
                      'Failure to extract data! Will retry "%s"',url);
              pause(5);
            else
              % error('Ecoforecasts:OISST:NCError',...
              %       'Max retries exceeded! Gave up at "%s"',url);
              warning('Ecoforecasts:OISST:NCError',...
                      'Max retries exceeded! Gave up at "%s"',url);
              dayDone = true;
              pause(20);
            end; %if ( tryCount <= maxRetries ) else
          end; %if ( isempty(dat) )

        end; %while ( ~dayDone )

      end; %for jd

      %DEBUG:
      toc,

      % FAILSAFE: Resave at the end of each year
      for ix=needix(:)'
        matfname = fullfile(datapath,[lower(stations(ix).station_name) '_oisst2.mat']);
        disp(['Saving MAT file ' matfname]);
        station = stations(ix);
        save(matfname,'station');
        station=[]; clear station;
      end; %for ix=needix(:)'

    end; %for yr

    for ix=needix(:)'
      matfname = fullfile(datapath,[lower(stations(ix).station_name) '_oisst2.mat']);
      disp(['Saving MAT file ' matfname]);
      station = stations(ix);
      save(matfname,'station');
      station=[]; clear station;
    end; %for ix=needix(:)'

  end; %if ( ~isempty(needix) )

  % Whether we loaded data, extracted data from netCDF files, or both, limit
  % our return to only dates the user requested - that have reasonable data
  % NOTE this also removes any time series points that still have NaN dates.
  for vix = 1:length(vars)
    fld = [PFX flds{vix}];

    for ix = 1:numel(stations)
      % Gross quality control
      if ( ~isfield(stations(ix),fld) || isempty(stations(ix).(fld)) )
        warning('Ecoforecasts:OISST:MissingField',...
                'No field "%s" found after load!',fld);
      elseif ( ~isempty(strfind(fld,'_field')) )
        % Make sure lon and lat are row vectors as for other products
        stations(ix).(fld).lon = stations(ix).(fld).lon(:);
        stations(ix).(fld).lat = stations(ix).(fld).lat(:);

        ixes = find(0 >= stations(ix).(fld).field | stations(ix).(fld).field >= 40);
        stations(ix).(fld).field(ixes) = nan;

        ixes = find(mindt <= stations(ix).(fld).date & stations(ix).(fld).date <= maxdt);
        stations(ix).(fld).date = stations(ix).(fld).date(ixes);
        stations(ix).(fld).field = stations(ix).(fld).field(ixes,:,:);
      else
        ixes = find(0 <= stations(ix).(fld).data & stations(ix).(fld).data <= 40);
        stations(ix).(fld).date = stations(ix).(fld).date(ixes);
        stations(ix).(fld).data = stations(ix).(fld).data(ixes);

        ixes = find(mindt <= stations(ix).(fld).date & stations(ix).(fld).date <= maxdt);
        stations(ix).(fld).date = stations(ix).(fld).date(ixes);
        stations(ix).(fld).data = stations(ix).(fld).data(ixes);
      end;
    end; %for ix = 1:numel(stations)
  end; %for vix = 1:length(vars)


  % If caller wants different interpolation method than we saved in MAT file,
  % re-call INTERP_FIELD on all the appropriate "_field" fields (if present)
  for ix = 1:numel(stations)

    if ( ~isfield(stations(ix),[PFX 'interp_method']) || ...
         ~strcmpi(stations(ix).([PFX 'interp_method']),char(interpMethod)) )
      disp(['Reinterpolating time series for ' stations(ix).station_name ...
            ' using ' upper(char(interpMethod))]);

      for vix = 1:length(vars)
        fld = [PFX flds{vix}];
        if ( isfield(stations(ix).(fld),'data') )
          fnm = [fld '_field'];
          if ( ~isfield(stations(ix),fnm) )
            warning('Ecoforecasts:Pathfinder:NoFieldToReinterp',...
                    '%s: No field available to reinterpolate "%s"',...
                    stations(ix).station_name,fld);
          else
            f = stations(ix).(fnm);
            stations(ix).(fld).data = ...
                interp_field(f.lat,f.lon,f.field,stations(ix).lat,stations(ix).lon,interpMethod);
          end;
        end;
      end; %for vix = 1:length(vars)

    end; %if ~strcmpi(stations(ix).([PFX 'interp_method']),interpMethod)

  end; %for ix = 1:numel(stations)


  for ix = 1:numel(stations)
    allflds = grepstruct(stations(ix),PFX);
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
