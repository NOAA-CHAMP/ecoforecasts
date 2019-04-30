function flds = get_nwps_stations(fldnms,dataset,dts,nwpspath)
%function flds = get_nwps_stations(fldnms,dataset,dts,nwpspath)
%
% Extract wave fields FLDNMS for dates DTS, from netCDF files for Nearshore
% Wave Prediction System (NWPS) DATASET (DEFAULT: 'mfl_nwps_CG1'). Returns
% struct FLDS with fields .(fldnms{ix}) for each field named in the CHAR or
% CELLSTR FLDNMS that was successfully retrieved. Each field includes the
% subfields .lon,.lat,.date, and .field. Elements of FLDNMS may include
% 'currdir', 'currspeed', 'primwavedir', 'primwaveper', 'sigwavehgt',
% 'sigswellhgt', 'winddir', 'windspeed', and calculated fields of Stokes
% drift and combined Stokes drift and Ekman surface current:
% 'ardhuin_surface_drift_u', 'ardhuin_surface_drift_v',
% 'ardhuin_surface_drift', 'ardhuin_surface_drift_dir',
% 'monismith_surface_drift_u', 'monismith_surface_drift_v',
% 'monismith_surface_drift', or 'monismith_surface_drift_dir'.
%
% NWPSPATH defaults to GET_ECOFORECASTS_PATH('data/NWPS').
%
% Refer to help and comments in the Python script 'PROCESS_NWPS.PY' for more
% information on both the netCDF files interpolated by NOAA/AOML from NWS
% GRiB2 files, and the calculated surface drift fields referred to above.
%
% Last Saved Time-stamp: <Thu 2018-11-15 17:46:12 Eastern Standard Time gramer>

  all_fldnms = { ...
      'currdir', 'currspeed', 'primwavedir', 'primwaveper', 'sigwavehgt', ...
      'sigswellhgt', 'winddir', 'windspeed', 'ardhuin_surface_drift_u', ...
      'ardhuin_surface_drift_v', 'ardhuin_surface_drift', ...
      'ardhuin_surface_drift_dir', 'monismith_surface_drift_u', ...
      'monismith_surface_drift_v', 'monismith_surface_drift', ...
      'monismith_surface_drift_dir', ...
             };

  flds = [];
  if ( ~exist('fldnms','var') || isempty(fldnms) )
    fldnms = all_fldnms;
  end;
  fldnms = lower(cellstr(fldnms));

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
    dts = all_dts(min(dts) <= all_dts & all_dts <= max(dts));
  end;

  if ( ~exist('nwpspath','var') )
    nwpspath = get_ecoforecasts_path('data/NWPS');
  end;

  for dtix = 1:numel(dts)
    dt = dts(dtix);
    [y,m,d,H,M,S] = datevec(dt);
    basefile = sprintf('%s_%04d%02d%02d_%02d00',dataset,y,m,d,H);
    clear y m d H M S;

    ncfname = fullfile(nwpspath,[basefile,'.nc']);
    if ( ~exist(ncfname,'file') )
      warning('Ecoforecasts:NWPS:NoFile','FILE NOT FOUND: %s',ncfname);
      %%%% EARLY RETURN
      continue;
    end;
    %DEBUG:    disp(ncfname);

    nc = mDataset(ncfname);
    try,
      hrs = nc{'time'}(:);
      nhrs = numel(hrs);
      t = datenum(1,1,1) + (hrs/24); 
      %DEBUG:      disp(datestr(t));

      if ( ~exist('lat','var') )
        lat = nc{'lat'}(:); 
        lon = nc{'lon'}(:); 
      end;

      for fldix = 1:numel(fldnms)
        fldnm = fldnms{fldix};
        if ( ~isfield(flds,fldnm) )
          flds.(fldnm).lat = lat;
          flds.(fldnm).lon = lon;
          flds.(fldnm).date = repmat(nan,[numel(dts)*4,1]);
          flds.(fldnm).field = repmat(nan,[numel(dts)*4,numel(lat),numel(lon)]);
        end;

        dtixen = (((dtix-1)*4)+1):(((dtix-1)*4)+4);
        flds.(fldnm).date(dtixen) = t;
        try,
          flds.(fldnm).field(dtixen,:,:) = nc{fldnm}(:,:,:);
        catch miniME,
          fldnms(fldix) = [];
          flds = rmfield(flds,fldnm);
          disp(['REMOVING FLDS.(',fldnm,'): not found while processing ',ncfname]);
          catchwarn(miniME);
        end;
      end; %for fldix = 1:numel(fldnms)
    catch ME,
      % Make sure we close the netCDF file whenever possible
      try,
        disp(['While processing ',ncfname]);
        close(nc); clear nc
      catch,
      end;
      rethrow(ME);
    end;
    close(nc); clear nc

  end; %for dtix = 1:numel(dts)

  % Remove missing timestamps from fields
  for fldix = 1:numel(fldnms)
    fldnm = fldnms{fldix};
    missingix = find(all(isnan(flds.(fldnm).field(:,:)),2));
    if ( ~isempty(missingix) )
      disp([upper(fldnm),': removing ',num2str(numel(missingix)),' dates of ',num2str(numel(flds.(fldnm).date))]);
      flds.(fldnm).date(missingix) = [];
      flds.(fldnm).field(missingix,:,:) = [];
    end;
  end;

return;
