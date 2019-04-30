function stns = get_nwps_stations(stns_or_locs,dataset,dts,fillSpatialGaps)
%function stns = get_nwps_stations(stns_or_locs,dataset,dts,fillSpatialGaps)
%
% Extract wave data at locations in STNS for dates DTS, from netCDF files for
% Nearshore Wave Prediction System (NWPS) DATASET (DEFAULT: 'mfl_nwps_CG1').
% Returns struct STNS with fields that were passed in, but with additional
% subfields raw_nwps_sigwavehgt, raw_nwps_sigswellhgt, raw_nwps_primwaveper,
% raw_nwps_primwavedir, raw_nwps_currdir, raw_nwps_currspeed, raw_nwps_winddir,
% raw_nwps_windspeed, raw_nwps_stokes_drift_u, raw_nwps_stokes_drift_v (with
% only values extracted from files) and their corresponding "non-raw" fields,
% with time series interpolated to hourly using INTERP_TS with MAXGAP 18 h.
%
% First arg may be a STRUCT of station STRUCTs, a matrix of station STRUCTS,
% or a CELL of numeric vectors of LON,LAT and optionally STNM (a CELLSTR).
% The returned value is always a STRUCT of STRUCTs with fieldnames STNM.
%
% If FILLSPATIALGAPS (DEFAULT: False), call INTERP_FIELD (v.) with @NANMEAN
% to attempt to fill NWPS values for stations near land.
%
% NOTE: If STNS.(stnm) already contains any of the above-named wave fields,
% their previous contents will be WIPED OUT in this call.
%
% Last Saved Time-stamp: <Fri 2016-12-02 09:55:47 Eastern Standard Time gramer>

  % Handle all four optional forms for first arg
  if ( ischar(stns_or_locs) )
    % 1st arg is a matrix of STRUCTs with .station_name,.lon,.lat
    stns_or_locs = get_station_from_station_name(stns_or_locs);
    if ( ~isfield(stns_or_locs,'lon') )
      error('Unknown station name %s',stns_or_locs);
    end;
  end;

  % 1st arg is a matrix of STRUCTs with .station_name,.lon,.lat
  if ( isstruct(stns_or_locs) )
    stnms = {}; stlons = []; stlats = [];
    if ( isfield(stns_or_locs,'station_name') && isfield(stns_or_locs,'lon') )
      % 1st arg is a matrix of STRUCTs with .station_name,.lon,.lat
      stnmat = stns_or_locs;
      stns = [];
      for stix = 1:numel(stnmat)
        stnms{end+1} = stnmat(stix).station_name;
        stlons(end+1,1) = stnmat(stix).lon;
        stlats(end+1,1) = stnmat(stix).lat;
        stns.(stnms{end}) = stnmat(stix);
      end;
      stnmat = []; clear stnmat;
    else
      % 1st arg is a STRUCT of STRUCTs with .lon,.lat
      stns = stns_or_locs;
      cstnms = fieldnames(stns);
      for stix = 1:numel(cstnms)
        stnm = cstnms{stix};
        if ( isfield(stns.(stnm),'lon') )
          stnms{end+1} = stnm;
          stlons(end+1,1) = stns.(stnm).lon;
          stlats(end+1,1) = stns.(stnm).lat;
        else
          disp(['SKIPPING ',stnm]);
        end;
      end; %for stix = 1:numel(cstnms)
      if ( isempty(stlons) )
        error('Is first arg not a STRUCT of STRUCTs??');
      end;
    end;
  elseif ( iscell(stns_or_locs) )
    % 1st arg is a cell of {LONS,LATS} or {LONS,LATS,STNMS{:}}
    locs = stns_or_locs;
    if ( 2 > numel(locs) || numel(locs) > 4 )
      error('If cell, first arg must have 2 or 3 elements');
    else
      stlons = locs{1};
      stlats = locs{2};
      if ( numel(locs) == 3 )
        stnms = locs{3};
      else
        stnms = matlab.lang.makeValidName(cellstr(strcat('stn_',num2str(stlons),'__',num2str(stlats))));
      end;
    end;
    stns = [];
    for stix = 1:numel(stnms)
      stnm = stnms{stix};
      stns.(stnm).station_name = stnm;
      stns.(stnm).lon = stlons(stix);
      stns.(stnm).lat = stlons(stix);
    end;
  else
    error('First arg must either be a struct of station structs or a CELL');
  end;
  clear stns_or_locs;

  if ( ~isnumeric(stlons) || ~ismatrix(stlons) || numel(stlons) ~= numel(stlats) )
    error('Must have equal-length numeric matrices of LON and LAT');
  end;

  if ( ~exist('dataset','var') || isempty(dataset) )
    %dataset = 'key_nwps_CG2';
    dataset = 'mfl_nwps_CG1';
  end;
  %% % Testing
  %% all_dts = [datenum(2016,07,31,6,0,0):0.5:datenum(2016,08,02,18,0,0)];
  all_dts = [datenum(2016,07,31,6,0,0):0.5:datenum(2016,10,1,6,0,0)];
  %all_dts = [datenum(2016,07,31,6,0,0):0.5:datenum(2016,10,31,18,0,0)];
  if ( ~exist('dts','var') || isempty(dts) )
    dts = all_dts;
  else
    dts = all_dts(min(dts) <= all_dts && all_dts <= max(dts));
  end;

  if ( ~exist('fillSpatialGaps','var') || isempty(fillSpatialGaps) )
    fillSpatialGaps = false;
  end;

  if ( ~exist('datapath','var') )
    datapath = get_ecoforecasts_path('data');
  end;

  flds = { 'currspeed', 'currdir', 'windspeed', 'winddir', ...
           'sigwavehgt', 'sigswellhgt', 'primwaveper', 'primwavedir', ...
           'ardhuin_surface_drift_u', 'ardhuin_surface_drift_v' };


  lon = [];
  lat = [];

  for dtix = 1:numel(dts)
    dt = dts(dtix);
    [y,m,d,H,M,S] = datevec(dt);
    basefile = sprintf('%s_%04d%02d%02d_%02d00',dataset,y,m,d,H);
    clear y m d H M S;

    ncfname = fullfile(datapath,[basefile,'.nc']);
    if ( ~exist(ncfname,'file') )
      warning('Ecoforecasts:NWPS:NoFile','FILE NOT FOUND: %s',ncfname);
      %%%% EARLY RETURN
      continue;
    end;
    %DEBUG:    disp(ncfname);

    nc = mDataset(ncfname);
    try,
      if ( isempty(lat) )
        lat = nc{'lat'}(:); 
        lon = nc{'lon'}(:); 
      end;
      hrs = nc{'time'}(:);
      nhrs = numel(hrs);
      t = datenum(1,1,1) + (hrs/24); 
      for fldix = 1:numel(flds)
        fld = flds{fldix};
        nval = nc{fld}(:,:,:);
        val = interp_field(lat,lon,nval,stlats,stlons,'*linear');
        if ( fillSpatialGaps )
          badix = find(all(isnan(val)));
          if ( ~isempty(badix) )
            val(:,badix) = interp_field(lat,lon,nval,stlats(badix),stlons(badix),@nanmean);
          end;
        end;
        % Careful with our memory
        nval = []; clear nval;
        for stix = 1:numel(stnms)
          stnm = stnms{stix};
          rawfld = ['raw_nwps_',fld];
          %if ( ~isfield(stns.(stnm),rawfld) )
          if ( dtix == 1 )
            % NOTE: This will wipe out any previous contents of RAWFLD!
            stns.(stnm).(rawfld).date = [];
            stns.(stnm).(rawfld).data = [];
          end;
          stns.(stnm).(rawfld).date(end+1:end+nhrs,1) = t;
          stns.(stnm).(rawfld).data(end+1:end+nhrs,1) = val(:,stix);
        end; %for stix = 1:numel(stnms)
        val=[]; clear val;
      end; %for fldix = 1:numel(flds)
    catch ME,
      % Make sure we close the netCDF file whenever possible
      try,
        close(nc); clear nc
      catch,
      end;
      rethrow(ME);
    end;
    close(nc); clear nc

  end; %for dtix = 1:numel(dts)


  % Interpolate station time series to hourly values
  for stix = 1:numel(stnms)
    stnm = stnms{stix};
    badflds = 0;
    for fldix = 1:numel(flds)
      fld = ['nwps_',flds{fldix}];
      rawfld = ['raw_',fld];
      if ( ~is_valid_ts(stns.(stnm).(rawfld)) )
        badflds = badflds + 1;
        % In case FLD existed previously, don't leave inconsistent data in it
        stns.(stnm).(fld) = empty_ts;
      else
        % If more than one file was missing (>12 h gap), do not gap-fill
        stns.(stnm).(fld) = interp_ts(stns.(stnm).(rawfld),[],[],'linear',(18/24));
      end;
    end; %for fldix = 1:numel(flds)
    if ( badflds > 0 )
      warning('Ecoforecasts:NWPS:BadFields','Invalid raw time series in STNS.%s',stnm);
    end;
  end; %for stix = 1:numel(stnms)

return;
