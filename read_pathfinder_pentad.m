function stns = read_pathfinder_pentad(stns_or_stnms,mindt,maxdt,vars,flds,interpMethod,baseurl)
%function stns = read_pathfinder_pentad(stns_or_stnms,mindt,maxdt,vars,flds,interpMethod,baseurl)
%
% Add fields for RSMAS/NODC Pathfinder v5.0 Cloud-filtered 5-day mean ("Pentad")
% Sea Surface Temperature to each member of struct array (or scalar) STNS. If a
% station name string or cell array of strings STNMS (five characters, e.g.,
% 'mlrf1', or {'smkf1','tnrf1'}) is given instead of structs, or if any struct
% lacks a valid .lon or .lat field (e.g., -80.38, 25.01) but all do have valid
% .station_name fields, tries to retrieve locations with GET_STATION_COORDS.
%
% Documentation for Pathfinder was accessed online as of 2011 Mar 24 from:
%  http://www.nodc.noaa.gov/SatelliteData/pathfinder4km/
%
% Optional VARS is a string or cellstr specifying which variables to extract.
% If VARS is specified, optional FLDS may also be specified as the corresp.
% list of field names to be added to structs in STNS, e.g., 'sst'.
%
% INTERPMETHOD is passed to INTERP_FIELD (DEFAULT 'linear') to do subsets.
% NOTE: Regardless INTERPMETHOD, INTER_FIELD's TRIANG arg is set to True.
%
% DEFAULT VARS and corresponding FLDS cell arrays are as follows:
%     vars = { 'SST4_daynitavg', 'STDV4_daynitavg', 'COUNTS4_daynitavg', ...
%              'SST4_daynitavg', };
%     flds = { 'sst',            'sst_std',         'sst_n', ...
%              'sst_field',      };
%
% NOTE: String 'pfv5c_pentad_' is prepended to each of the names in FLDS.
%
% DEFAULT for BASEURL:
%  'http://data.nodc.noaa.gov/thredds/dodsC/pathfinder/Version5.0_CloudScreened/5day/FullRes'
%
% CALLS: MDATASET (netCDF-Java); INTERP_FIELD, GET_STATION_COORDS (Ecoforecasts).
%
% Last Saved Time-stamp: <Thu 2011-04-14 10:13:34  Lew.Gramer>

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


  % Variable name prefix: Path Finder v5.0 Cloud-filtered 5-day mean
  PFX = 'pfv5c_pentad_';

  if ( ~exist('vars','var') || isempty(vars) )
    vars = { 'SST4_daynitavg', 'STDV4_daynitavg', 'COUNTS4_daynitavg', ...
             'SST4_daynitavg', };
  end;
  if ( ~exist('flds','var') || isempty(flds) )
    flds = { 'sst',            'sst_std',         'sst_n', ...
             'sst_field',      };
  end;

  if ( ~exist('interpMethod','var') || isempty(interpMethod) )
    interpMethod = 'linear';
  end;
  if ( ischar(interpMethod) && interpMethod(1) == '*' )
    interpMethod = interpMethod(2:end);
  end;


  if ( ~exist('baseurl','var') || isempty(baseurl) )
    %Sample URL:
    %http://data.nodc.noaa.gov/thredds/dodsC/pathfinder/Version5.0_CloudScreened/5day/FullRes/2009/2009361-2009365.pfv5sst_5day.hdf
    baseurl = 'http://data.nodc.noaa.gov/thredds/dodsC/pathfinder/Version5.0_CloudScreened/5day/FullRes';
  end;

  % First Pentad of 1985 started on jday 004 and would need special handling
  yrs = 1986:2009;

  % %%%% UNCOMMENT AND MODIFY THE FOLLOWING LINE TO CONTINUE A PARTIAL EXTRACTION
  % yrs = 2007:2009;

  alldts = [];
  for yr=yrs(:)'
    alldts = [ alldts datenum(yr,1,0)+[1:5:365] ];
  end;
  alldts = alldts';

  % Grid-point radii for time-series fields (e.g., seatemp_field)
  % Choose one axis slightly larger to guard against transposition errors
  % yrad = 10;
  % xrad = 11;
  yrad = 4;
  xrad = 5;


  needix = 1:numel(stns);

  for ix=1:numel(stns)
    matfname = fullfile(datapath,[lower(stns(ix).station_name) '_pathfinder_pentad.mat']);
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

    maxRetries = 5;

    for yr=yrs(:)'
      %DEBUG:
      disp(yr);      tic,

      for jd=1:5:365

        tryCount = 1;
        dayDone = false;

        while ( ~dayDone )
          
          dt = datenum(yr,1,0) + jd;
          [ig,dtix] = min(abs(alldts - dt));

          endjd = jd+4;
          % Stupid leap years!
          if ( jd==361 && mod(yr,4) == 0 )
            endjd = 366;
          end;

          url = sprintf('%s/%04d/%04d%03d-%04d%03d.pfv5sst_5day.hdf',...
                        baseurl,yr,yr,jd,yr,endjd);
          %DEBUG:          disp(url);
          dat = [];
          nc = mDataset(url);
          if ( isempty(nc) )
            warning('Ecoforecasts:Pathfinder:BadURL',...
                    'Bad URL "%s"',url);
          else
            try
              if ( isempty(lats) || isempty(lons) )
                lons = cast(nc{'Longitude'}(:),'double');
                lats = cast(nc{'Latitude'}(:),'double');

                for ix = needix(:)'
                  [yerr,yix(ix)] = min(abs(lats - stations(ix).lat));
                  [xerr,xix(ix)] = min(abs(lons - stations(ix).lon));
                  if ( yerr > max(diff(unique(lats))) || xerr > max(diff(unique(lons))) )
                    % Global product - this should never happen!
                    warning('Ecoforecasts:Pathfinder:StationOutsideDomain',...
                            'Station "%s" at %g,%g is outside Pathfinder domain?!',...
                            stations(ix).station_name,stations(ix).lat,stations(ix).lon);
                    needix(needix == ix) = [];
                  else
                    lat{ix} = lats(yix(ix)-2:yix(ix)+2);
                    lon{ix} = lons(xix(ix)-3:xix(ix)+3);

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
                    dat = cast(nc{var}(yix(ix)-yrad:yix(ix)+yrad,xix(ix)-xrad:xix(ix)+xrad),'double');
                    dat(dat<0) = nan;
                    stations(ix).(fld).field(dtix,:,:) = dat;
                  else
                    dat = cast(nc{var}(yix(ix)-2:yix(ix)+2,xix(ix)-3:xix(ix)+3),'double');
                    dat(dat<0) = nan;
                    stations(ix).(fld).data(dtix,1) = ...
                        interp_field(lat{ix},lon{ix},dat,stations(ix).lat,stations(ix).lon,interpMethod,[],true);
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
              warning('Ecoforecasts:Pathfinder:NCError',...
                      'Failure to extract data! Will retry "%s"',url);
              pause(3*tryCount);
            else
              error('Ecoforecasts:Pathfinder:NCError',...
                    'Max retries exceeded! Gave up at "%s"',url);
            end; %if ( tryCount <= maxRetries ) else
          end; %if ( isempty(dat) )

        end; %while ( ~dayDone )

      end; %for jd

      %DEBUG:
      toc,

      % Save to MAT files at end of every year as a fail-safe
      for ix=needix(:)'
        matfname = fullfile(datapath,[lower(stations(ix).station_name) '_pathfinder_pentad.mat']);
        disp(['Saving MAT file ' matfname]);
        station = stations(ix);
        save(matfname,'station');
        station=[]; clear station;
      end; %for ix=needix(:)'

      % % Give remote server a breather
      % pause(20);

    end; %for yr

    for ix=needix(:)'
      matfname = fullfile(datapath,[lower(stations(ix).station_name) '_pathfinder_pentad.mat']);
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
        warning('Ecoforecasts:Pathfinder:MissingField',...
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
                interp_field(f.lat,f.lon,f.field,stations(ix).lat,stations(ix).lon,interpMethod,'w',true);
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
