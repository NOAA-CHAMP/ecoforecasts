function res = get_sfomc_data
%function res = get_sfomc_data
%
% Get data from the South Florida Oil Monitoring Center field project of Nova
% Southeastern University's Oceanography Center (Soloviev et al.)
%
% Last Saved Time-stamp: <Thu 2016-08-18 10:14:04 Eastern Daylight Time gramer>

  datapath = get_ecoforecasts_path('data');

  sfo_matfname = fullfile(datapath,'SFOMC.mat');

  if ( exist(sfo_matfname,'file') )
    disp(['Loading ',sfo_matfname]);
    load(sfo_matfname);

  else
    disp('Extracting SFOMC');
    doPlots = false;
    saveMat = false;
    forceExtract = false;

    cwd = pwd;
    cd(get_coral_path('CRCP/Upwelling'));
    extract_sfomc_etc;
    cd(cwd);

    res.mc_sites = mc_sites;
    res.upwix = upwix;
    res.upwdt = upwdt;

    disp(['Saving ',sfo_matfname]);
    save(sfo_matfname,'res','-v7.3');
  end;

return;
