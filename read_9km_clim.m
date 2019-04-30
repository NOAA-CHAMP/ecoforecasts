function [varargout] = read_9km_clim(fname)
%function [clim,climstd,climn] = read_9km_clim(fname)
%
% Load MISST 9km climatology into memory. All output arguments OPTIONAL:
%     clim    = Climatology. Last index indicates day/night: 1=day.
%     climstd = STD of SST at that point
%     climn   = Number of observations at that point.
%
% Adapted from 'RD2GZ', 'RD2GZUNIX', and sample code - all by Chelle Gentemann.
%
% Last Saved Time-stamp: <Thu 2009-09-24 16:46:50 Eastern Daylight Time lew.gramer>

  clim = [];
  climstd = [];
  climn = [];

  if ( ~exist('fname','var') || isempty(fname) )
    fname = '//cygnus/gramer/home/coral/misst/climatology/misst_climatology.dat';
  end;

  dchar='float32';
  idim = [4096 2048 2];
  ilen = prod(idim);

  % This should work for BOTH Windoze and Linux
  % fid = fopen(fname, 'r', 'ieee-be');
  % But for Windoze - you actually need this...
  fid = fopen(fname, 'r');
  if ( fid < 0 )
    error('UNABLE TO OPEN CLIMATOLOGY "%s"!', fname);
  end;

  for ix=1:nargout
    a = fread(fid, ilen, dchar);
    varargout{ix} = reshape(a, idim);
  end;

  fclose(fid);

return;
