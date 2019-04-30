function stn = get_avhrr_weekly_field(stn_or_stnm,calcFieldTerms,interpMethod,npts,keepBadDates)
%function stn = get_avhrr_weekly_field(stn_or_stnm,calcFieldTerms,interpMethod,npts,keepBadDates)
%
% Extract a (default) 17x17 1km pixel field around the position of station
% contained in the struct, or named in the string STN_OR_STNM. Station struct
% is updated with a new field STN.avhrr_weekly_sst_field, with the subfields
% described in EXTRACT_AVHRR_WEEKLY_FIELD (v.), which this function calls. If
% fields .lat and .lon are present, interpolate SST field to time series. If
% optional CALCFIELDTERMS is True, calculate gradient and Laplacian, AND if
% .lon,.lat are also present, interpolate the field, gradients, and Laplacian
% to STN.avhrr_weekly_sst, _x, _y, and _l, resp. If optional arg INTERPMETHOD
% also given (DEFAULT: 'linear'), pass to INTERP_FIELD and CALC_FIELD_TERMS.
% If optional NPTS (DEFAULT: 5) is given, also pass to CALC_FIELD_TERMS. If
% optional KEEPBADDATES (DEFAULT: false), call STATION_FILTER_BAD_DATES (v.)
%
% Last Saved Time-stamp: <Wed 2012-07-04 17:47:29  lew.gramer>

  set_more off;

  datapath = get_ecoforecasts_path('data');

  stn = get_station_from_station_name(stn_or_stnm);
  stn_or_stnm=[]; clear stn_or_stnm;
  stnm = lower(stn.station_name);

  if ( ~exist('calcFieldTerms','var') || isempty(calcFieldTerms) )
    calcFieldTerms = false;
  end;
  if ( ~exist('interpMethod','var') || isempty(interpMethod) )
    interpMethod = 'linear';
  end;
  % If user wants CALC_FIELD_TERMS, how many points to use in TEMPLATE
  if ( ~exist('npts','var') || isempty(npts) )
    npts = 5;
  end;
  % If user does NOT want us to call STATION_FILTER_BAD_DATES
  if ( ~exist('keepBadDates','var') || isempty(keepBadDates) )
    keepBadDates = false;
  end;


  %% Extract weekly composite field (or re-load if already extracted)

  matfname = fullfile(datapath,[stnm '_avhrr_weekly_field.mat']);

  if ( exist(matfname,'file') )

    disp(['Loading from presaved ' matfname]);
    load(matfname,'station');

  else

    disp('Extracting from raw AVHRR Weekly PNG...');
    station = extract_avhrr_weekly_field(stn);
    disp(['Saving to ' matfname]);
    save(matfname,'station');

  end;


  %% Let user specify BAD weekly composites in each station's AVHRR data
  % (See STATION_FILTER_BAD_DATES.m)

  % stn.avhrr_weekly_sst_field = station.avhrr_weekly_sst;
  stn.raw_avhrr_weekly_sst_field = station.avhrr_weekly_sst;
  % First week of the year always begins on 01 Jan
  wkix = find(ismember(get_jday(station.avhrr_weekly_sst.date),[1:7:358]));
  % Offset raw data to time midway between start and end of week (+3.5)
  stn.raw_avhrr_weekly_sst_field.date = station.avhrr_weekly_sst.date(wkix);
  stn.raw_avhrr_weekly_sst_field.field = station.avhrr_weekly_sst.field(wkix,:,:);
  % Field .date should be a column vector (Nx1)
  stn.raw_avhrr_weekly_sst_field.date = stn.raw_avhrr_weekly_sst_field.date(:);
  station = []; clear station;


  if ( ~keepBadDates )
    stn = station_filter_bad_dates(stn);
  end;

  % Offset raw data to time midway between start and end of week (+3.5)
  stn.raw_avhrr_weekly_sst_field.date = stn.raw_avhrr_weekly_sst_field.date + 3.5;


  %% Spatially interpolate time series from weekly field

  if ( isfield(stn,'lat') && isfield(stn,'lon') )
    stn.raw_avhrr_weekly_sst.date = stn.raw_avhrr_weekly_sst_field.date;
    stn.raw_avhrr_weekly_sst.data = ...
        interp_field(stn.raw_avhrr_weekly_sst_field.lat,stn.raw_avhrr_weekly_sst_field.lon,...
                     stn.raw_avhrr_weekly_sst_field.field,stn.lat,stn.lon,interpMethod);
  end;


  %% Interpolate daily field from weekly field

  % Redo: EXTRACT_AVHRR_WEEKLY_FIELD did not call STATION_FILTER_BAD_DATES
  dts = stn.raw_avhrr_weekly_sst_field.date;
  sst = stn.raw_avhrr_weekly_sst_field.field;
  alldts = [dts(1):dts(end)];
  % Quick and dirty way of finding all non-land pixels
  meansst = squeeze(nanmean(sst,1));
  seaix = find(isfinite(meansst));
  nwks = size(sst,1);

  stn.avhrr_weekly_sst_field.date = alldts(:);
  stn.avhrr_weekly_sst_field.lon = stn.raw_avhrr_weekly_sst_field.lon;
  stn.avhrr_weekly_sst_field.lat = stn.raw_avhrr_weekly_sst_field.lat;
  stn.avhrr_weekly_sst_field.N = repmat(nan, [size(sst,2) size(sst,3)]);
  stn.avhrr_weekly_sst_field.pctused = repmat(nan, [size(sst,2) size(sst,3)]);
  stn.avhrr_weekly_sst_field.field = repmat(nan, [numel(alldts) size(sst,2) size(sst,3)]);
  for ix = seaix(:)'
    usedix = find(isfinite(sst(:,ix)));
    stn.avhrr_weekly_sst_field.N(ix) = numel(usedix);
    stn.avhrr_weekly_sst_field.pctused(ix) = (numel(usedix)/nwks);
    stn.avhrr_weekly_sst_field.field(:,ix) = ...
        interp1(dts(usedix),sst(usedix,ix),alldts,'pchip',nan);
  end;
  dts=[]; sst=[]; alldts=[]; clear dts sst alldts;


  %% Spatially interpolate time series from daily field

  if ( isfield(stn,'lat') && isfield(stn,'lon') )
    stn.avhrr_weekly_sst.date = stn.avhrr_weekly_sst_field.date;
    stn.avhrr_weekly_sst.data = ...
        interp_field(stn.avhrr_weekly_sst_field.lat,stn.avhrr_weekly_sst_field.lon,...
                     stn.avhrr_weekly_sst_field.field,stn.lat,stn.lon,interpMethod);
  end;


  %% Calc gradients, Laplacians, and interp hourly time series from daily fields

  if ( calcFieldTerms )
    if ( isfield(stn,'lat') && isfield(stn,'lon') )
      stn = calc_field_terms(stn,'raw_avhrr_weekly_sst_field',...
                             'raw_avhrr_weekly_sst',interpMethod,stn.lat,stn.lon,npts);
      stn = calc_field_terms(stn,'avhrr_weekly_sst_field',...
                             'avhrr_weekly_sst',interpMethod,stn.lat,stn.lon,npts);
      % Interpolate hourly time series from daily time series
      stn.hourly_avhrr_weekly_sst = interp_ts(stn.avhrr_weekly_sst);
    else
      stn = calc_field_terms(stn,'raw_avhrr_weekly_sst_field',[],[],[],[],npts);
      stn = calc_field_terms(stn,'avhrr_weekly_sst_field',[],[],[],[],npts);
      warning('No hourly time series');
    end;
  end;

  set_more;

return;
