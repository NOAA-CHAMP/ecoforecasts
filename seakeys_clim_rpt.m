function stn = seakeys_clim_rpt(stn_or_stnm)
%function stn = seakeys_clim_rpt(stn_or_stnm)

  datapath = get_ecoforecasts_path('data');

  stn = get_station_from_station_name(stn_or_stnm);
  stnm = lower(stn.station_name);
  if ( ~isfield(stn,'ndbc_wind1_speed') )
    stn = load_all_ndbc_data(stn);
  end;

  flds = grepstruct(stn,'^ndbc_');
  for fldix = 1:numel(flds)
    fld = flds{fldix};
    disp(['Calculating stats for ',fld]);


    fmg; boxplot_ts(stn.(fld)

    gotData = false;
    try
      [dly,dly_tid] = grp_ts(stn.(fld).data,stn.(fld).date,@floor,@nanmean,16);
      [mea,mea_tid] = grp_ts(dly,dly_tid,@get_month,@nanmean,20);
      [ sd, sd_tid] = grp_ts(dly,dly_tid,@get_month,@nanstd,20);
      [med,med_tid] = grp_ts(dly,dly_tid,@get_month,@nanmedian,20);
      [ iq, iq_tid] = grp_ts(dly,dly_tid,@get_month,@iqr,20);
      if ( numel(mea_tid)==12 && numel(sd_tid)==12 && numel(med_tid)==12 && numel(iq_tid)==12 )
        gotData = true;
      end;
    catch,
      % Ignore errors for now
    end;
    if ( ~gotData )
      warning('Missing month(s) for "%s"??',fld);
    else
      stn.ndbc_stats.(fld).date   = mea_tid;
      stn.ndbc_stats.(fld).mean   = mea;
      stn.ndbc_stats.(fld).sd     =  sd;
      stn.ndbc_stats.(fld).median = med;
      stn.ndbc_stats.(fld).iqr    =  iq;
    end;
  end;

  fname = fullfile(datapath,[lower(stn.station_name),'_ndbc_stats.csv']);
  fid = fopen(fname,'w');
  if ( fid < 0 )
    error('Unable to open for writing "%s"',fname);
  end;
  try
    flds = fieldnames(stn.ndbc_stats);
    for fldix = 1:numel(flds)
      fld = flds{fldix};
      fprintf(fid,'stn,fld,JAN_lower,JAN_center,JAN_upper,FEB_lower,FEB_center,FEB_upper,MAR_lower,MAR_center,MAR_upper,APR_lower,APR_center,APR_upper,MAY_lower,MAY_center,MAY_upper,JUN_lower,JUN_center,JUN_upper,JUL_lower,JUL_center,JUL_upper,AUG_lower,AUG_center,AUG_upper,SEP_lower,SEP_center,SEP_upper,OCT_lower,OCT_center,OCT_upper,NOV_lower,NOV_center,NOV_upper,DEC_lower,DEC_center,DEC_upper\n');
      fprintf(fid,'%s,%s,mean',stnm,fld);
      for mo=1:12
        fprintf(fid,',%g,%g,%g',stn.ndbc_stats.(fld).mean(mo)-(1*stn.ndbc_stats.(fld).sd(mo)),stn.ndbc_stats.(fld).mean(mo),stn.ndbc_stats.(fld).mean(mo)+(1*stn.ndbc_stats.(fld).sd(mo)));
      end;
      fprintf(fid,'\n');
      fprintf(fid,'%s,%s,median',stnm,fld);
      for mo=1:12
        fprintf(fid,',%g,%g,%g',stn.ndbc_stats.(fld).median(mo)-(stn.ndbc_stats.(fld).sd(mo)*2.5),stn.ndbc_stats.(fld).median(mo),stn.ndbc_stats.(fld).median(mo)+(stn.ndbc_stats.(fld).iqr(mo)*2.5));
      end;
      fprintf(fid,'\n');
    end;
  catch
    fclose(fid);
    rethrow lasterror;
  end;
  fclose(fid);
  disp(['Statistics written to ',fname]);


return;
