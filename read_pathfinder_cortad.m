function stns = read_pathfinder_cortad(stns_or_stnms,mindt,maxdt,vars,flds,interpMethod,baseurl)
%function stns = read_pathfinder_cortad(stns_or_stnms,mindt,maxdt,vars,flds,interpMethod,baseurl)
%
% Add fields for RSMAS/NODC Pathfinder v5.0 Coral Reef Temperature Anomaly
% Database ("CoRTAD") Sea Surface Temperature to each member of struct array
% (or scalar) STNS. If a station name string or cell array of strings STNMS
% (five characters, e.g., 'mlrf1', or {'smkf1','tnrf1'}) is given instead of
% structs, or if the structs have no valid .lon and .lat fields, but all have
% .station_name fields, tries to retrieve locations with GET_STATION_COORDS.
%
% Documentation for Pathfinder was accessed online as of 2011 Mar 24 from:
%  http://www.nodc.noaa.gov/SatelliteData/pathfinder4km/
% Documentation for CoRTAD was accessed online as of 2011 Apr 03 from:
%  http://www.nodc.noaa.gov/SatelliteData/Cortad/
%
% Optional VARS is a string or cellstr specifying which variables to extract.
% If VARS is specified, optional FLDS may also be specified as the corresp.
% list of field names to be added to structs in STNS, e.g., 'sst'.
%
% INTERPMETHOD is passed to INTERP_FIELD (DEFAULT 'linear') to do subsets.
%
% DEFAULT VARS and corresponding FLDS cell arrays are as follows:
%     vars = { 'FilledSST', 'MedfillSST', 'WeeklySST', 'SSTA', 'SSTA_DHW', 'TSA', 'TSA_DHW', ...
%              'FilledSST', };
%     flds = { 'sst',       'sst_med',    'sst_raw',   'ssta', 'ssta_dhw', 'tsa', 'tsa_dhw', ...
%              'sst_field',      };
%
% NOTE: String 'pfv5c_cortad_' is prepended to each of the names in FLDS.
%
% DEFAULT for BASEURL:
%  'http://data.nodc.noaa.gov/opendap/cortad/Version3'
%
% CALLS: MDATASET (netCDF-Java); INTERP_FIELD, GET_STATION_COORDS (Ecoforecasts).
%
% Last Saved Time-stamp: <Thu 2011-04-14 10:20:10  Lew.Gramer>

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
  PFX = 'pfv5c_cortad_';

  if ( ~exist('vars','var') || isempty(vars) )
    vars = { 'FilledSST', 'MedfillSST', 'WeeklySST', 'SSTA', 'SSTA_DHW', 'TSA', 'TSA_DHW', ...
             'FilledSST', };
  end;
  if ( ~exist('flds','var') || isempty(flds) )
    flds = { 'sst',       'sst_med',    'sst_raw',   'ssta', 'ssta_dhw', 'tsa', 'tsa_dhw', ...
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
    %http://data.nodc.noaa.gov/opendap/cortad/Version3/cortadv3_row04_col00.h5
    baseurl = 'http://data.nodc.noaa.gov/opendap/cortad/Version3';
  end;

  % % Full span of dataset - variable not used currently, just for reference
  % yrs = 1982:2009;

  % Grid-point radii for time-series fields (e.g., seatemp_field)
  % Choose one axis slightly larger to guard against transposition errors
  yrad = 4;
  xrad = 5;


  needix = 1:numel(stns);

  for ix=1:numel(stns)
    matfname = fullfile(datapath,[lower(stns(ix).station_name) '_pathfinder_cortad.mat']);
    % If we already did this before - do not need to load again.
    if ( exist(matfname,'file') )
      disp(['Reloading MAT file ' matfname]);
      load(matfname,'station');
      stations(ix) = station;
      station = []; clear station;
      if ( isempty(stns(ix).lon) || isempty(stns(ix).lat) )
        stns(ix).lon = stations(ix).lon; stns(ix).lat = stations(ix).lat;
      end;
      needix(needix == ix) = [];

    else
      if ( isempty(stns(ix).lon) || isempty(stns(ix).lat) )
        [stns(ix).lon,stns(ix).lat,stns(ix).depth] = get_station_coords(stns(ix).station_name);
      end;

      stations(ix).station_name = stns(ix).station_name;
      stations(ix).lon = stns(ix).lon; stations(ix).lat = stns(ix).lat;
      stations(ix).([PFX 'interp_method']) = lower(char(interpMethod));
      % for vix=1:numel(vars)
      %   fld = [PFX flds{vix}];
      %   % Preallocate for fast execution
      %   stations(ix).(fld) = [];
      %   if ( ~isempty(strfind(fld,'_field')) )
      %     stations(ix).(fld).date = repmat(nan,[length(alldts),1]);
      %     stations(ix).(fld).field = repmat(nan,[length(alldts),(2*yrad)+1,(2*xrad)+1]);
      %   else
      %     stations(ix).(fld).date = repmat(nan,[length(alldts),1]);
      %     stations(ix).(fld).data = repmat(nan,[length(alldts),1]);
      %   end;
      % end;
    end;
  end;


  if ( ~isempty(needix) )

    disp('Loading original data from online HD5 file...');

    %DEBUG:
    tic,

    %%%% ??? Later, the following can be replaced by logic to group all
    %stations in STNS(NEEDIX) by tile, accessing each tile HD5 separately!

    % Calculate the CoRTAD tile corresponding to stations,
    % also checking that all stations are in the same tile.

    tileno(1) = floor((90 - min([stns(needix).lat])) / 22.5);
    tileno(2) = floor((min([stns(needix).lon]) + 180) / 22.5);

    chktileno(1) = floor((90 - max([stns(needix).lat])) / 22.5);
    chktileno(2) = floor((max([stns(needix).lon]) + 180) / 22.5);
    if ( any(tileno ~= chktileno) )
      error('For now all stations must lie in same CoRTAD tile! %s vs. %s',...
            num2str(tileno),num2str(chktileno));
    end;

    % switch (region)
    %  case 'asam',           tileno = [4, 0];  % American Samoa region
    %  case 'ecarib',         tileno = [3, 5];  % eastern Caribbean
    %  case 'freef',          tileno = [2, 4];  % Florida and Bahamas reefs
    %  case 'gbr',            tileno = [4,14];  % Great Barrier Reef region
    %  case 'nind',           tileno = [3,11];  % north central Indian Ocean reefs
    %  case 'nred',           tileno = [2, 9];  % northern Red Sea and W Med
    %  case 'ntri',           tileno = [3,13];  % northern Coral Triangle
    %  case 'papa',           tileno = [2, 0];  % Papahanaumokuakea (NWHI)
    %  case 'sind',           tileno = [4,11];  % south central Indian Ocean reefs
    %  case 'sred',           tileno = [3, 9];  % southern Red Sea
    %  case 'swtri',          tileno = [4,13];  % SW Coral Triangle (incl Ningaloo)
    %  case 'wcarib',         tileno = [3, 4];  % western Caribbean
    % end;


    for ix = needix(:)'
      stations(ix).([PFX 'tileno']) = tileno;
    end;

    url = sprintf('%s/cortadv3_row%02d_col%02d.h5',baseurl,tileno(1),tileno(2));
    nc = mDataset(url);
    if ( isempty(nc) )
      error('Unable to open "%s"!',url);
    end;

    try
      dts = getTimes(nc{'FilledSST'});
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

      for vix=1:numel(vars)
        var = vars{vix};
        fld = [PFX flds{vix}];
        %DEBUG:
        disp(fld);        tic,
        for ix=needix(:)'
          stations(ix).(fld).date = dts(:);
          if ( ~isempty(strfind(fld,'_field')) )
            dat = cast(nc{var}(yix(ix)-yrad:yix(ix)+yrad,xix(ix)-xrad:xix(ix)+xrad,:),'double');
            stations(ix).(fld).field = permute(dat,[3 1 2]);
          else
            dat = cast(nc{var}(yix(ix)-2:yix(ix)+2,xix(ix)-3:xix(ix)+3,:),'double');
            dat = permute(dat,[3 1 2]);
            stations(ix).(fld).data = ...
                interp_field(lat{ix},lon{ix},dat,stations(ix).lat,stations(ix).lon,interpMethod);
          end; %if ( ~isempty(strfind(fld,'_field')) )
        end; %for ix=needix(:)'
        %DEBUG:
        toc,
      end; %for vix=1:numel(vars)
    catch
      if ( exist('nc','var') && ~isempty(nc) ); close(nc); end;
      clear nc
      rethrow(lasterror);
    end; %try catch
    if ( exist('nc','var') && ~isempty(nc) ); close(nc); end;
    clear nc

    %DEBUG:
    toc,

    for ix=needix(:)'
      matfname = fullfile(datapath,[lower(stations(ix).station_name) '_pathfinder_cortad.mat']);
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

        % Anomaly SST fields present - negative values possible
        ixes = find(-20 >= stations(ix).(fld).field | stations(ix).(fld).field >= 40);
        stations(ix).(fld).field(ixes) = nan;

        ixes = find(mindt <= stations(ix).(fld).date & stations(ix).(fld).date <= maxdt);
        stations(ix).(fld).date = stations(ix).(fld).date(ixes);
        stations(ix).(fld).field = stations(ix).(fld).field(ixes,:,:);

        % CoRTAD is supposed to be WEEKLY - but there was a glitch for 1990!
        gapix = find(diff(stations(ix).(fld).date) < (7-eps));
        stations(ix).(fld).date(gapix) = [];
        stations(ix).(fld).field(gapix,:,:) = [];
      else
        % Anomaly SST fields present - negative values possible
        ixes = find(-20 <= stations(ix).(fld).data & stations(ix).(fld).data <= 40);
        stations(ix).(fld).date = stations(ix).(fld).date(ixes);
        stations(ix).(fld).data = stations(ix).(fld).data(ixes);

        ixes = find(mindt <= stations(ix).(fld).date & stations(ix).(fld).date <= maxdt);
        stations(ix).(fld).date = stations(ix).(fld).date(ixes);
        stations(ix).(fld).data = stations(ix).(fld).data(ixes);

        % CoRTAD is supposed to be WEEKLY - but there was a glitch for 1990!
        gapix = find(diff(stations(ix).(fld).date) < (7-eps));
        stations(ix).(fld).date(gapix) = [];
        stations(ix).(fld).data(gapix) = [];
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
