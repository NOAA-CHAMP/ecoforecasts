function stns = get_nmc_winds(stns_or_stnms)
%function stns = get_nmc_winds(stns_or_stnms)
%
% Subset U and V wind components in M/S (and calculate speed in KTS and
% direction degT) from NMC Reanalysis, for all sites in struct(s) STNS.
% First arg may also be a string or cell array of strings of site names.
%
% Last Saved Time-stamp: <Fri 2011-03-04 07:10:12  Lew.Gramer>

  %DEBUG:  tic,

  datapath = get_ecoforecasts_path('data');
  nmcpath = fullfile(datapath,'nmc');

  if ( ischar(stns_or_stnms) )
    stns_or_stnms = {stns_or_stnms};
  end;
  if ( iscellstr(stns_or_stnms) )
    stns = [];
    for ix=1:length(stns_or_stnms)
      stns(ix).station_name = stns_or_stnms{ix};
    end;    
  elseif ( isstruct(stns_or_stnms) )
    stns = stns_or_stnms;
  else
    error('First arg must be (cell array of) string(s) or (vector of) struct(s)');
  end;

  needExtracted = [];

  vars = {'uwnd','vwnd'};
  flds = {'nmc_wind_u','nmc_wind_v'};

  for ix=1:length(stns(:))
    if ( isfield(stns(ix),'station_name') && ~isempty(stns(ix).station_name) )
      stnm = stns(ix).station_name;
    else
      stnm = stns(ix).name;
    end;

    matfname{ix} = fullfile(datapath,[stnm '_nmc_winds.mat']);
    if ( ~exist(matfname{ix},'file') )
      needExtracted(end+1) = ix;
      if ( ~isfield(stns(ix),'lon') || isempty(stns(ix).lon) )
        [stns(ix).lon,stns(ix).lat,stns(ix).depth] = get_station_coords(stnm);
      end;
      for varix = 1:length(vars)
        fld = ['native_' flds{varix}];
        stns(ix).(fld) = struct('date',[],'data',[]);
      end;
    else
      %DEBUG:
      disp(['Loading ' matfname{ix}]);
      load(matfname{ix},'station');
      sflds = fieldnames(station);
      for fldix = 1:length(sflds)
        fld = sflds{fldix};
        stns(ix).(fld) = station.(fld);
      end;
      station = []; clear station;
    end;
  end;

  if ( isempty(needExtracted) )
    %DEBUG:    disp('All requested stations were loaded');
    %%%%%%%%%%%%%%%
    % EARLY RETURN
    %%%%%%%%%%%%%%%
    return;
  end;


  for yr=2002:2009
    for varix = 1:length(vars)
      var = vars{varix};
      fld = ['native_' flds{varix}];
      %DEBUG:
      disp({yr,var});
      fname = fullfile(nmcpath,[var '.10m.gauss.' num2str(yr) '.nc']);
      nc = mDataset(fname);
      if ( isempty(nc) )
        warning('Skipping bad or missing file "%s"',fname);
        break;
      end;
      try
        lon=nc{'lon'}(:,:);
        lat=nc{'lat'}(:,:);
        hrs = nc{'time'}(:);
        dts = datenum(1,1,1,0,0,0) + (hrs/24);
        ndts = length(dts);
        for ix=needExtracted(:)'
          mylon = stns(ix).lon;
          mylon(mylon<0) = mylon(mylon<0) + 360;
          lonix=round(interp1(lon,1:length(lon),mylon));
          latix=round(interp1(lat,1:length(lat),stns(ix).lat));
          if ( 1<latix && latix<length(lat) && 1<lonix && lonix<length(lon) )
            mylatix = latix-1:latix+1; mylonix = lonix-1:lonix+1;
            rawdat = nc{var}(:,mylatix,mylonix);
            dat = interp_field(lat(mylatix),lon(mylonix),rawdat,stns(ix).lat,mylon);
          else
            dat = nc{var}(:,latix,lonix);
          end;
          stns(ix).(fld).date(end+1:end+ndts,1) = dts;
          stns(ix).(fld).data(end+1:end+ndts,1) = dat;
        end;
      catch
        close(nc); clear nc
        rethrow(lasterror);
      end;
      close(nc); clear nc
    end;
  end;

  for varix = 1:length(vars)
    rawfld = ['native_' flds{varix}];
    fld = flds{varix};
    for ix=needExtracted(:)'
      stns(ix).(fld).date = ...
          [stns(ix).(rawfld).date(1):(1/24):stns(ix).(rawfld).date(end)]';
      stns(ix).(fld).data = ...
          interp1(stns(ix).(rawfld).date,stns(ix).(rawfld).data,stns(ix).(fld).date,'spline');
    end;
  end;

  for ix=needExtracted(:)'
    stns(ix).native_nmc_wind_speed.date = stns(ix).native_nmc_wind_u.date;
    stns(ix).native_nmc_wind_speed.data = ...
        mps2kts(uv_to_spd(stns(ix).native_nmc_wind_u.data,stns(ix).native_nmc_wind_v.data));
    stns(ix).native_nmc_wind_dir.date = stns(ix).native_nmc_wind_u.date;
    stns(ix).native_nmc_wind_dir.data = ...
        uv_to_dir(stns(ix).native_nmc_wind_u.data,stns(ix).native_nmc_wind_v.data);

    stns(ix).nmc_wind_speed.date = stns(ix).nmc_wind_u.date;
    stns(ix).nmc_wind_speed.data = ...
        mps2kts(uv_to_spd(stns(ix).nmc_wind_u.data,stns(ix).nmc_wind_v.data));
    stns(ix).nmc_wind_dir.date = stns(ix).nmc_wind_u.date;
    stns(ix).nmc_wind_dir.data = ...
        uv_to_dir(stns(ix).nmc_wind_u.data,stns(ix).nmc_wind_v.data);
  end;

  for ix=needExtracted(:)'
    %DEBUG:
    disp(['Saving ' matfname{ix}]);
    station = stns(ix);
    save(matfname{ix},'station');
    station = []; clear station;
  end;

  %DEBUG:  toc,

return;
