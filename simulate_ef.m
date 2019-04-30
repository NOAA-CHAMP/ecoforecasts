1;
% function simulate_ef(stn, ef)


  tfz = stns(ix).factories.seatemp_7d.fuzzies;
  pfz = stns(ix).factories.par_7d.fuzzies;
  wfz = stns(ix).factories.wind_7d.fuzzies;

  tdh = find(stns(ix).misst_sst_7_day_maximum.data >= tfz{9}{2}(1));
  tvh = find(stns(ix).misst_sst_7_day_maximum.data >= tfz{8}{2}(1));
  thi = find(stns(ix).misst_sst_7_day_maximum.data >= tfz{7}{2}(1));

  tdhj = unique(round(stns(ix).misst_sst_7_day_maximum.date(tdh)));
  tvhj = unique(round(stns(ix).misst_sst_7_day_maximum.date(tvh)));
  thij = unique(round(stns(ix).misst_sst_7_day_maximum.date(thi)));


  wdl = find(stns(ix).qscat_speed_7_day_maximum.data >= wfz{1}{2}(1));
  wvl = find(stns(ix).qscat_speed_7_day_maximum.data >= wfz{2}{2}(1));
  wlo = find(stns(ix).qscat_speed_7_day_maximum.data >= wfz{3}{2}(1));
  wsl = find(stns(ix).qscat_speed_7_day_maximum.data >= wfz{4}{2}(1));

  wdlj = unique(round(stns(ix).qscat_speed_7_day_maximum.date(wdl)));
  wvlj = unique(round(stns(ix).qscat_speed_7_day_maximum.date(wvl)));
  wloj = unique(round(stns(ix).qscat_speed_7_day_maximum.date(wlo)));
  wslj = unique(round(stns(ix).qscat_speed_7_day_maximum.date(wsl)));





  tl_tdt = stns(ix).misst_sst_7_day_maximum.date(stns(ix).misst_sst_7_day_maximum.data >= tfz{8}{2}(1));
  tl_pdt = stns(ix).global_par_7_day_maximum.date(stns(ix).global_par_7_day_maximum.data >= pfz{8}{2}(1));
  tl_tjd = unique(round(tl_tdt));
  tl_pjd = unique(round(tl_pdt));

  tlw_tdt = stns(ix).misst_sst_7_day_maximum.date(stns(ix).misst_sst_7_day_maximum.data >= tfz{7}{2}(1));
  tlw_pdt = stns(ix).global_par_7_day_maximum.date(stns(ix).global_par_7_day_maximum.data >= pfz{7}{2}(1));
  tlw_wdt = stns(ix).qscat_speed_7_day_maximum.date(stns(ix).qscat_speed_7_day_maximum.data >= wfz{7}{2}(1));
  tlw_tjd = unique(round(tlw_tdt));
  tlw_pjd = unique(round(tlw_pdt));
  tlw_wjd = unique(round(tlw_wdt));

  length(t_tjd),
  length(intersect(tl_tjd,tl_pjd)),
  length(intersect(tlw_tjd,intersect(tlw_pjd,tlw_wjd))),
