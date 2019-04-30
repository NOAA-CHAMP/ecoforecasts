function stn = get_misst_station(stn_or_stnm)
%function stn = get_misst_station(stn_or_stnm)

  datapath = get_ecoforecasts_path('data');

  stn = get_station_from_station_name(stn_or_stnm);
  if ( ~isfield(stn,'station_name') )
    error('Need a station_name field to extract MISST data!');
  end;

  misstfname = fullfile(datapath,[lower(stn.station_name),'_misst_freef.mat']);
  if ( exist(misstfname,'file') )
    disp(['Loading ' misstfname]);
    load(misstfname,'station');
    flds = fieldnames(station);
    for fldix=1:length(flds)
      fld = flds{fldix};
      stn.(fld) = station.(fld);
    end;
    station=[]; clear station;

  else

    error('NO MISST DATA .MAT FILE STORED FOR %s YET!',stn.station_name);

  end;

return;
