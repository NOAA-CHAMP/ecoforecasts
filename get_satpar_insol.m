function stn = get_satpar_insol(stnm_or_stn)
%function stn = get_satpar_insol(stnm_or_stn)
%
% Extract satellite insolation and other radiative flux (and total cloud
% cover) data from TOMS GSIP CSV files (which must have been previously
% extracted from TOMS GSIP and placed in DATAPATH), for station struct STN or
% station named STNM. Multiple new fields are added to returned STN struct.
%
% Last Saved Time-stamp: <Mon 2011-10-24 16:31:55  lew.gramer>

  datapath = get_ecoforecasts_path('data');

  stn = get_station_from_station_name(stnm_or_stn);
  stnm = lower(stn.station_name);

  matfname = fullfile(datapath,[stnm '_satpar_insol.mat']);
  if ( exist(matfname,'file') )
    disp(['Reloading from ' matfname]);
    load(matfname,'station');

  else
    disp([stnm ': Extracting insolation from raw CSV files...']);

    station.station_name = stnm;

    dat = read_all_insolation_csvs(stnm);

    if ( isempty(dat) )
      warning('No satellite insolation data could be loaded...');

    else
      station.sat_par.date = dat.date(:);
      station.sat_par.data = dat.par(:);
      station.sat_par_watt.date = dat.date(:);
      station.sat_par_watt.data = dat.parW(:);

      station.sat_insol_in.date = dat.date(:);
      station.sat_insol_in.data = dat.sd(:);
      station.sat_insol_out.date = dat.date(:);
      station.sat_insol_out.data = dat.su(:);
      station.sat_net_insol.date = dat.date(:);
      station.sat_net_insol.data = dat.sd(:) - dat.su(:);

      station.sat_longwave_in.date = dat.date(:);
      station.sat_longwave_in.data = dat.ld(:);
      station.sat_longwave_out.date = dat.date(:);
      station.sat_longwave_out.data = dat.lu(:);
      station.sat_net_longwave.date = dat.date(:);
      station.sat_net_longwave.data = dat.ld(:) - dat.lu(:);

      station.sat_cloud_cover.date = dat.date(:);
      station.sat_cloud_cover.data = dat.cld(:);
    end;

    dat = []; clear dat;

    disp(['Saving to ' matfname]);
    save(matfname,'station');

  end;

  flds = fieldnames(station);
  for fldix = 1:numel(flds)
    stn.(flds{fldix}) = station.(flds{fldix});
  end;
  station=[]; clear station;

return;
