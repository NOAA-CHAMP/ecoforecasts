function stn = get_avhrr_weekly_anomaly(stn_or_stnm,calcFieldTerms,interpMethod)
%function stn = get_avhrr_weekly_anomaly(stn_or_stnm,calcFieldTerms,interpMethod)
%
% Extract a (default) 17x17 1km pixel ANOMALY field around the position of
% station contained in struct, or named in string STN_OR_STNM. Station struct
% is updated with a new field STN.avhrr_weekly_anomaly_field, with the subfields
% described in EXTRACT_AVHRR_WEEKLY_FIELD (v.), which this function calls. If
% STN .lat and .lon present, interpolate SST anomaly field to time series. If
% optional CALCFIELDTERMS is True, calculate gradient and Laplacian, AND if
% .lon,.lat are also present, interpolate the field, gradients, and Laplacian
% to STN.avhrr_weekly_anomaly, _x, _y, _l, resp. If optional arg INTERPMETHOD
% also given (DEFAULT: 'linear'), pass to INTERP_FIELD and CALC_FIELD_TERMS.
%
% Last Saved Time-stamp: <Tue 2011-10-25 14:06:12  lew.gramer>

  set_more off;

  datapath = get_ecoforecasts_path('data');

  if ( ischar(stn_or_stnm) )
    stnm = lower(stn_or_stnm);
    stn.station_name = stnm;
  else
    stn = stn_or_stnm;
    stnm = lower(stn.station_name);
  end;

  % If user wants CALC_FIELD_TERMS, how many points to use in TEMPLATE
  npts = 5;

  if ( ~exist('calcFieldTerms','var') || isempty(calcFieldTerms) )
    calcFieldTerms = false;
  end;
  if ( ~exist('interpMethod','var') || isempty(interpMethod) )
    interpMethod = 'linear';
  end;


  matfname = fullfile(datapath,[stnm '_avhrr_weekly_anomaly.mat']);

  if ( exist(matfname,'file') )

    disp(['Loading from presaved ' matfname]);
    load(matfname,'station');

  else

    disp('Extracting from raw AVHRR Weekly PNG...');
    station = extract_avhrr_weekly_field(stnm,[],[],'anomaly');
    disp(['Saving to ' matfname]);
    save(matfname,'station');

  end;

  stn.avhrr_weekly_anomaly_field = station.avhrr_weekly_sst;
  % Field .date should be a column vector (Nx1)
  stn.avhrr_weekly_sst_field.date = stn.avhrr_weekly_sst_field.date(:);
  station = []; clear station;

  if ( isfield(stn,'lat') && isfield(stn,'lon') )
    stn.avhrr_weekly_sst.date = stn.avhrr_weekly_sst_field.date;
    stn.avhrr_weekly_sst.data = ...
        interp_field(stn.avhrr_weekly_sst_field.lat,stn.avhrr_weekly_sst_field.lon,...
                     stn.avhrr_weekly_sst_field.field,stn.lat,stn.lon,interpMethod);
  end;

  if ( calcFieldTerms )
    if ( isfield(stn,'lat') && isfield(stn,'lon') )
      stn = calc_field_terms(stn,'avhrr_weekly_sst_field',...
                             'avhrr_weekly_sst',interpMethod,stn.lat,stn.lon,npts);
    else
      stn = calc_field_terms(stn,'avhrr_weekly_sst_field',[],[],[],[],npts);
    end;
  end;

  set_more;

return;
