function tses = get_nwps_mat_tses(stns_or_coords,fldnm,mdl,dts,interpMethod)
%function tses = get_nwps_mat_tses(stns_or_coords,fldnm,mdl,dts,interpMethod)
%
% Read year-month .MAT files with analyzed Nearshore Wave Prediction System
% (NWPS) data, and subset into a continuous time series STRUCT for each point
% specified in the arg STNS_OR_COORDS (see GET_COORDS_FROM_ARGS). DEFAULTS:
% FLDNM='nwps_sigwavehgt', MDL='mfl_nwps_CG1', DTS=[2016/6/1 - 2019/3/1],
% INTERPMETHOD='linear' (see INTERP_FIELD).
%
% RETURNS: TSES, a vector of time series STRUCT (each with .date,.data).
%
% Last Saved Time-stamp: <Mon 2019-02-18 11:38:14 Eastern Standard Time gramer>

  set_more off;

  datapath = get_ecoforecasts_path('data');
  
  if ( ~exist('fldnm','var') || isempty(fldnm) )
    fldnm = 'nwps_sigwavehgt';
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
  
  [lons,lats,stnms] = get_coords_from_args(stns_or_coords,'station_name');
  
  % Convert DTS to a 3-hourly DATENUM vector with whole-day boundaries
  begix = find(get_hour(dts)==0,1);
  endix = find(get_hour(dts)>=21,1,'last');
  dts = dts(begix:endix);
  dts = dts(ismember(get_hour(dts),[0:3:24]));
  
  for ix=1:numel(lons)
    tses(ix).date = dts(:);
    tses(ix).data = repmat(nan,size(tses(ix).date));
    tses(ix).lon = lons(ix);
    tses(ix).lat = lats(ix);
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
    dats = interp_field(fld.lat,fld.lon,fld.(fldnm).field,lats,lons,interpMethod);
    for ix=1:numel(tses)
      tses(ix).data(dtix) = dats(fldix,ix);
    end;
  end;

  set_more;

return;
