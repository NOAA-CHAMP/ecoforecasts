function [sst,mask] = read_sst_image(fpath,isAnom)
%function [sst,mask] = read_sst_image(fpath,isAnom)
%
% Read PNG (Portable Network Graphics) or similar image file in FPATH
% containing University of South Florida 1km or 4km AVHRR or MODIS Sea
% Surface Temperature (SST) data, and extract 2-D matrix of double SST values
% from it (nominal precision ~ 0.1992 [oC]). Land, coast, and cloud mask
% pixels or those with unrealistic values are returned as NaN. FPATH may be a
% synoptic (snapshot) or a composite (daily, weekly, monthly) image. If
% FPATH is an ANOMALY image, specify optional second arg ISANOM as TRUE, to
% prevent anomaly values below +3 from being incorrectly returned as NaNs.
%
% If FPATH,
%
% CALLS: IMREAD (ImageSci); URLWRITE (IOFun)
%
% Last Saved Time-stamp: <Fri 2011-03-25 09:05:58  lew.gramer>

  datapath = get_ecoforecasts_path('data');

  %http://cyclops.marine.usf.edu/modis/level3/husf/florida/2005/186/1km/pass/final/MODIS.2005186.035203.florida.sst.png
  if ( strncmpi(fpath,'http',4) || strncmpi(fpath,'ftp:',4) )
    url = fpath;
    [ig,fname] = fileparts(fpath);
    fpath = fullfile(datapath,fname);
    [fpath,downloadResult] = urlwrite(url,fpath);
  end;

  if ( ~exist(fpath,'file') )
    error('Missing "%s"',fpath);
  end;

  %    case 'mean',
  minval = +2.0; maxval = +38.0;
  if ( exist('isAnom','var') && isAnom )
    %    case 'anomaly',
    minval = -10.0; maxval = +10.0;
  end;

  sstbytes = imread(fpath);

  % Calculate SST from 1-byte PNG (colormap) values
  % (Code MAY need tweaking for other image formats!)
  sst = (cast(sstbytes, 'double') * 0.1992) - 2.1;
  mask = repmat(false,size(sstbytes));

  % Per Brian Barnes, USF, 23 Feb 2011:
  %  coastlines are 255
  %  land is 254
  %  no data (clouds) are 251
  sst(sstbytes >= 251) = nan;
  mask(sstbytes >= 254) = true;

  % Also remove suspect data
  sst(minval > sst | sst > maxval) = nan;

  sstbytes = []; clear sstbytes

return;
