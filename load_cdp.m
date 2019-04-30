function [dts,sst] = load_cdp(fbasename)
%function [dts,sst] = load_cdp(fbasename)

  datapath = get_ecoforecasts_path('data');
  cdppath = fullfile(datapath,'cdp');

  % datapath = 'Tutuila_2006-2008';
  % % fname = fullfile(datapath, 'SST268015_SBE39306591265_20080219.cdp');
  fname = fullfile(cdppath, fbasename);
  dat = load('-ascii',fname);

  dts = datenum(dat(:,1:6));
  sst = dat(:,7);

return;
