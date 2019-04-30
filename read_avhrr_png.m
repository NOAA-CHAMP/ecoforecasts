function sst = read_avhrr_png(fpath,isAnom)
%function sst = read_avhrr_png(fpath,isAnom)
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
% Last Saved Time-stamp: <Fri 2011-08-19 11:31:22  Lew.Gramer>

  if ( ~exist(fpath,'file') )
    error('Missing "%s"',fpath);
  end;

  %    case 'mean',
  minval = +3.0; maxval = +35.0;
  if ( exist('isAnom','var') && isAnom )
    %    case 'anomaly',
    minval = -10.0; maxval = +10.0;
  end;

  sstbytes = imread(fpath);

  % Calculate SST from 1-byte PNG (colormap) values
  sst = (cast(sstbytes, 'double') * 0.1992) - 2.1;

  % Assume that BOTH 254 and 255 are cloud mask values(??)
  % sst(sstbytes >= 254) = nan;

  % Per Brian Barnes, USF, 23 Feb 2011:
  %  coastlines are 255
  %  land is 254
  %  no data (clouds) are 251
  sst(sstbytes >= 251) = nan;

  % Just in case our assumptions about cloud mask are wrong!
  sst(minval > sst | sst > maxval) = nan;

  sstbytes = []; clear sstbytes

return;
