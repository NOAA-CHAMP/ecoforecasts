function fld = read_usf_modis_png(fldbytes_or_fname,typ)
%function fld = read_usf_modis_png(fldbytes_or_fname,typ)
%
% Convert unsigned-char image FLDBYTES (or IMREAD image file named FNAME) in
% USF Optics MODIS format, into DOUBLE matrix of FLD [oC], chlor_a, or NTU.
%
% Last Saved Time-stamp: <Wed 2013-10-23 13:00:19 Eastern Daylight Time gramer>

  fname = [];
  % Is first arg a filename or an array of byte (uint8) values?
  if ( ischar(fldbytes_or_fname) && exist(fldbytes_or_fname, 'file') )
    fname = fldbytes_or_fname;
    fldbytes = imread(fname);
  elseif ( strcmp(class(fldbytes_or_fname), 'uint8') )
    fldbytes = fldbytes_or_fname;
  else
    error('First arg must be an image filename, or an array of bytes ("uint8")');
  end;

  % Second arg: What type of image or dataset is this?
  if ( ~exist('typ','var') || isempty(typ) )
    % typ = 'sst';
    % Guess data type from filename
    if ( isempty(fname) )
      error('Cannot guess PNG file data type!');
    elseif ( ~isempty(regexp(fname,'[.]SST[.]png')) )
      typ = 'sst';
    elseif ( ~isempty(regexp(fname,'[.]CHL[.]png')) )
      typ = 'chlor_a';
    elseif ( ~isempty(regexp(fname,'[.]NTU[.]png')) )
      typ = 'ntu';
    else
      error('Unknown PNG file data type!');
    end;
  end;

  %% Calculate SST, chlorophyll, etc. from 1-byte PNG (colormap) values

  %On Oct 22, 2013 12:16 AM, "Chuanmin Hu" <huc@usf.edu> wrote:
  % If you do want to use the png files for whatever reason, here is how
  % to convert them to the physical values:
  %
  % 255: coastline
  % 254: landmask
  % 251: clouds
  %
  % For SST,  0 => 10 degree C; 235 => 32 degree C. Between 0 and 235,
  % use linear interpolation.
  %
  % For Chl and NTU, load the image in Google Earth you'll see the color
  % legend showing their range (say, x_min, x_max. For Chl, this may be
  % x_min = 0.05, x_max = 10 (you need to check); for NTU, this may be
  % x_min = 0.5, x_max = 10).
  %
  % Then, perform the following:
  % y_min = alog10(x_min); y_max = alog10(x_max)
  % Then, 0 => y_min; 235 => y_max. For any pixel value between 0 and
  % 235, use linear interpolation to find its corresponding value y.
  % Then, convert this y to the geophysical value (either chl or NTU)
  % as: x = 10^y

  fld = cast(fldbytes, 'double');
  switch ( lower(typ) )
   case 'sst',
    % oC
    fld = (fld * 0.0936) + 10;
   case 'chlor_a',
    % mg/m^3: Scale from colorbar of sample image A20132731815.QKM.SE_FL.PASS.L3D.CHL.png
    ymn=log10(0.10); ymx=log10(20.0);
    y = ymn + ((ymx-ymn).*(fld./235));
    fld = 10.^y;
   case 'ntu',
    % NTU ("relative"): Scale from colorbar of sample image A20132731815.QKM.SE_FL.PASS.L3D.NTU.png
    ymn=log10(0.50); ymx=log10(10.0);
    y = ymn + ((ymx-ymn).*(fld./235));
    fld = 10.^y;
   otherwise,
    error('Unknown PNG image type "%s"!',typ);
  end;

  % Anything >250 is some kind of mask value
  fld(fldbytes > 250) = nan;

return;
