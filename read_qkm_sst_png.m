function sst = read_qkm_sst_png(fpath,isAnom)
%function sst = read_qkm_sst_png(fpath,isAnom)
%
% Read PNG (Portable Network Graphics) image file FPATH containing University
% of South Florida 1km or 4km AVHRR Sea Surface Temperature, and extract 2-D
% matrix of double SST values from it (nominal precision ~ 0.2 [oC]). Cloud
% mask pixels or those with unrealistic values are returned as NaN. FPATH may
% be a synoptic (snapshot) or a composite (daily, weekly, monthly) image. If
% FPATH is an ANOMALY image, specify optional second arg ISANOM as TRUE, to
% prevent anomaly values below +3 from being incorrectly returned as NaNs.
%
% CALLS: IMREAD
%
% Last Saved Time-stamp: <Mon 2011-08-29 16:56:57  lew.gramer>

  % if ( ~exist(fpath,'file') )
  %   warning('Assuming URL "%s"',fpath);
  % end;

  %    case 'mean',
  minval = +3.0; maxval = +35.0;
  if ( exist('isAnom','var') && isAnom )
    %    case 'anomaly',
    minval = -10.0; maxval = +10.0;
  end;

  sstbytes = imread(fpath);

  % Calculate SST from 1-byte PNG (colormap) values
  sst = (22.0 * cast(sstbytes,'double') / 235.0) + 10.0;

  % Per Chuanmin - values above 235 (32oC) are mask values(?)
  sst(sstbytes >= 236) = nan;

  % % Just in case our assumptions about cloud mask are wrong!
  % sst(minval > sst | sst > maxval) = nan;

  sstbytes = []; clear sstbytes

return;
