function mask = get_sst_landmask(region_or_fname)
%function mask = get_sst_landmask(region_or_fname)
%
% Return landmask matrix for University of South Florida SST dataset REGION.
%
% REGION values currently supported: 'florida','keys_sst'. This argument may
% also specify a path and filename to an actual (e.g., PNG) SST image.
%
% Last Saved Time-stamp: <Fri 2011-03-04 12:18:30  Lew.Gramer>

  mydocs = get_ecoforecasts_path('../..');

  switch (region_or_fname)
   case 'keys_sst',
    fpath = fullfile(mydocs,'/coral/Bleach/ColdEvent/keys_sst/1995_2010_avhrr_clim/Jan_1to31_mean_avhrr_clim.png');
   case 'florida',
    fpath = fullfile(mydocs,'/RSMAS/Coastal/thesis/data/avhrr/all.2000001.2000007.sst.day_night.florida.mean.png');
   otherwise,
    fpath = region_or_fname;
  end;

  %DEBUG:  dir(fpath),

  sstbytes = imread(fpath);

  mask = repmat(false,size(sstbytes));

  % Per Brian Barnes, USF, 23 Feb 2011:
  %  coastlines are 255
  %  land is 254
  mask(sstbytes >= 254) = true;

  sstbytes = []; clear sstbytes

return;

