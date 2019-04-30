1;

load(fullfile(get_ecoforecasts_path('data'),'T_tses_1.mat'),'ts');

calc_upwelling_ekman_flux
ts = copy_T_ts(face,ts,'FACE_HWD_');
ts = copy_T_ts(sefcri,ts,'SEFCRI_');
ts = copy_T_ts(sfomc,ts,'SFOMC_');
face=[]; sefcri=[]; sfomc=[]; clear face sefcri sfomc

ncore = get_ncore_2000;
[ncore.klgf1,ncore.mrtf1] = an_ncore;
ts = copy_T_ts(ncore,ts,'NCORE_');
ncore=[]; clear ncore

save(fullfile(get_ecoforecasts_path('data'),'T_tses.mat'),'ts');
