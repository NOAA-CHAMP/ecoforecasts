function stn = load_aims_ncs(stn_or_stnm,qcFlagsToKeep)
%function stn = load_aims_ncs(stn_or_stnm,qcFlagsToKeep)
%
% Load all data for station STN_OR_STNM from netCDF files. If QCFLAGSTOKEEP
% is a nonempty vector (DEFAULT: 0:2) of "quality_control" flags (see AIMS
% metadata), only RETURN data points associated with one of those flags.
%
% Last Saved Time-stamp: <Sun 2011-11-13 16:26:22  Lew.Gramer>

  stn=[];

  %DEBUG:
  tic,

  set_more off;

  datapath = get_ecoforecasts_path('data');
  aimspath = fullfile(datapath,'aims');

  if ( ~exist('qcFlagsToKeep','var') )
    % qcFlagsToKeep = [];
    qcFlagsToKeep = [0,1,2];
  end;

  stn = get_station_from_station_name(stn_or_stnm);
  stnm = lower(stn.station_name);

  matfname = fullfile(datapath,[stnm,'_aims.mat']);
  if ( exist(matfname,'file') )
    disp(['Loading ',matfname]);
    load(matfname,'station');

  else
    disp(['Extracting from raw netCDF files for ',stnm]);

    station = [];
    for cfld={'station_name','lat','lon','depth'};
      if ( isfield(stn,cfld{:}) )
        station.(cfld{:}) = stn.(cfld{:});
      end;
    end;

    vars = {...
        'ATMP',...
        'AIRT',...
        'Photosynthetically_Active_Radiation',...
        'TEMP',...
        'Wind_Direction',...
        'Wind_Speed',...
           };

    flds = {...
        'aims_barom',...
        'aims_air_t',...
        'aims_par',...
        'aims_sea_t',...
        'aims_wind1_dir',...
        'aims_wind1_speed_kph',...
           };

    nams = {...
        'Barometric Pressure hPa',...
        'Air Temperature ^oC',...
        'Photosynthetically Active Radiation umol^.s^-1.m^-2',...
        'Sea Water Temperature ^oC',...
        'Wind From Direction ^oTrue',...
        'Wind Speed km/hr',...
           };

    stnpath = fullfile(aimspath,lower(stnm));
    ds = dir(fullfile(stnpath,'*.nc'));
    if ( isempty(ds) )
      warning('No netCDF files found in "%s"!',stnpath);
    end;
    for ix=1:numel(ds);
      res = [];

      fname=fullfile(stnpath,ds(ix).name);
      %DEBUG:
      disp(fname);
      nc=mDataset(fname);
      if ( isempty(nc) )
        warning('Unable to open "%s"',fname);
        continue;
      end;
      try,
        % Seconds since the UNIX epoch
        secs = cast(nc{'time'}(:),'double');

        for varix=1:numel(vars)
          var=vars{varix};
          fld=flds{varix};
          flgfld = [fld '_flags'];
          nam=nams{varix};
          nv=[];
          try,
            [ig,nv]=evalc(['nc{','''',var,'''}']);
          catch;
          end;
          if ( ~isempty(nv) )
            res.(fld).long_name = nam;
            res.(fld).date = datenum(1970,1,1) + (secs./(24*3600));
            res.(fld).data = cast(nv(:),'double');
            nfv=[];
            try,
              [ig,nfv]=evalc(['nc{','''',var,'_quality_control''}']);
            catch;
            end;
            if ( ~isempty(nfv) )
              res.(flgfld).date = res.(fld).date;
              res.(flgfld).data = cast(nfv(:),'double');
            end; %if ( ~isempty(nfv) )
          end; %if ( ~isempty(nv) )
        end; %for varix=1:numel(vars)
      catch,
        close(nc);
        rethrow(lasterror);
      end;
      close(nc);

      if ( ~isempty(res) )
        station = merge_station_data(station,res);
      end;
      res = []; clear res;
    end; %for ix=1:numel(ds);

    disp(['Saving ',matfname]);
    save(matfname,'station');
  end;


  flds = grepstruct(station,'^aims_');

  % Filter on QC flags if requested to
  if ( ~isempty(qcFlagsToKeep) )
    for fldix=1:numel(flds)
      fld = flds{fldix};
      if ( isempty(strfind(fld,'_flags')) )
        flgfld = [fld '_flags'];
        if ( isfield(station,flgfld) )
          badix = find(~ismember(station.(flgfld).data,qcFlagsToKeep));
          if ( ~isempty(badix) )
            %DEBUG:
            disp(['Filtering ',num2str(numel(badix)),' from ',fld]);
            station.(fld).date(badix) = [];
            station.(fld).data(badix) = [];
          end;
        end;
      end;
    end;
  end;

  % Transfer AIMS data to return STN struct
  for fldix=1:numel(flds)
    fld = flds{fldix};
    stn.(fld) = station.(fld);
  end;
  station=[]; clear station;

  if ( isfield(stn,'aims_wind1_speed_kph') )
    NM_PER_KM=0.539956803455724;
    % Convert [km/hr] -> [kts]
    stn.aims_wind1_speed.long_name = ...
        [stn.aims_wind1_speed_kph.long_name,' (cvt KTS)'];
    stn.aims_wind1_speed.date = stn.aims_wind1_speed_kph.date;
    stn.aims_wind1_speed.data = stn.aims_wind1_speed_kph.data.*NM_PER_KM;
  end;

  set_more;

  %DEBUG:
  toc,


return;
