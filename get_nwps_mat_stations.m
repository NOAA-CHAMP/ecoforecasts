function stns = get_nwps_mat_stations(stns_or_coords,fldnms,mdl,dts,interpMethod)
%function stns = get_nwps_mat_stations(stns_or_coords,fldnms,mdl,dts,interpMethod)
%
% Read year-month .MAT files with analyzed Nearshore Wave Prediction System
% (NWPS) data, and subset into a continuous time series STRUCT for each point
% specified in the arg STNS_OR_COORDS (see GET_COORDS_FROM_ARGS). DEFAULTS:
% FLDNMS={'nwps_sigwavehgt','nwps_primwavedir','nwps_primwaveper','nwps_swellhgt'},
% MDL='mfl_nwps_CG1', DTS=[2016/6/1 - 2019/3/1], INTERPMETHOD='linear' (see INTERP_FIELD).
%
% Last Saved Time-stamp: <Sun 2019-03-03 17:36:53 Eastern Standard Time gramer>

  set_more off;
  
  %datapath = get_ecoforecasts_path('data');
  datapath = get_datapath('NWPS');
  
  if ( ~exist('fldnms','var') || isempty(fldnms) )
    fldnms = {'nwps_sigwavehgt','nwps_primwavedir','nwps_primwaveper','nwps_swellhgt'};
  end;
  if ( ischar(fldnms) )
    fldnms = {fldnms};
  end;
  if ( ~exist('mdl','var') || isempty(mdl) )
    %mdl = 'key_nwps_CG2';
    mdl = 'mfl_nwps_CG1';
  end;
  if ( ~exist('dts','var') || isempty(dts) )
    dts = datenum(2016,6,1):(1/24):datenum(2019,3,1);
  end;
  if ( ~exist('interpMethod','var') || isempty(interpMethod) )
    interpMethod = 'linear';
  end;
  
  stns = get_stations_from_args(stns_or_coords,{'station_name'});
  
  % Convert DTS to a 3-hourly DATENUM vector with whole-day boundaries
  begix = find(get_hour(dts)==0,1);
  endix = find(get_hour(dts)>=21,1,'last');
  dts = dts(begix:endix);
  dts = dts(ismember(get_hour(dts),[0:3:24]));
  
  for ix=1:numel(stns)
    lons(ix) = stns(ix).lon;
    lats(ix) = stns(ix).lat;
    for fldix=1:numel(fldnms)
      fldnm = fldnms{fldix};
      stns(ix).(fldnm).date = dts(:);
      stns(ix).(fldnm).data = repmat(nan,size(stns(ix).(fldnm).date));
    end;
  end;
  
  yrmos = unique(get_yearmonth(dts));
  for yrmo = yrmos(:)';
    yr = get_year(yrmo);
    mo = get_month(yrmo);
    
    nwps_matfname = fullfile(datapath,sprintf('NWPS_%s_%04d_%02d.mat',mdl,yr,mo));
    if ( ~exist(nwps_matfname,'file') )
      disp(['MISSING ',nwps_matfname]);
      %%%% EARLY CONTINUE
      continue;
    end;
    %disp(sprintf('Extracting NWPS for %s_%04d_%02d',mdl,yr,mo));
    fld = load(nwps_matfname);
    [dtix,fldix] = intersect_dates(dts,fld.date,(1.5/24));
    for fldix=1:numel(fldnms)
      fldnm = fldnms{fldix};
      dats = interp_field(fld.lat,fld.lon,fld.(fldnm).field,lats,lons,interpMethod);
      for ix=1:numel(stns)
        stns(ix).(fldnm).data(dtix) = dats(fldix,ix);
      end;
    end;
  end;
  
  set_more;

return;
