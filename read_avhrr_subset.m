function sst = read_avhrr_subset(sstbytes_or_fname)
%function sst = read_avhrr_subset(sstbytes_or_fname)
%
% Convert unsigned-char image SSTBYTES (or IMREAD image file named FNAME) in
% USF AVHRR format, into DOUBLE matrix SST of sea-surface temperatures [oC].
%
% Last Saved Time-stamp: <Fri 2011-08-19 11:31:30  Lew.Gramer>

  if ( ischar(sstbytes_or_fname) && exist(sstbytes_or_fname, 'file') )
    sstbytes = imread(sstbytes_or_fname);
  elseif ( strcmp(class(sstbytes_or_fname), 'uint8') )
    sstbytes = sstbytes_or_fname;
  else
    error('First arg must be an image filename, or an array of bytes ("uint8")');
  end;

  minval = +0.0; maxval = +40.0;

  % Calculate SST from 1-byte PNG (colormap) values
  sst = (cast(sstbytes, 'double') * 0.1992) - 2.1;
  % Assume that BOTH 254 and 255 are cloud mask values(??)
  sst(sstbytes >= 254) = nan;
  % Just in case our assumptions about cloud mask are wrong!
  sst(minval > sst | sst > maxval) = nan;

return;
