function [flds,lon,lat] = read_usf_modis_hdf(fname_or_nc,varnms)
%function [flds,lon,lat] = read_usf_modis_hdf(fname_or_nc,varnms)
%
% Extracts data from HDF file or URL named FNAME or already opened (i.e., by
% MDATASET, v.) as NC: return vectors of LON and LAT if desired, and cell
% array of 2D data fields FLDS. Note: If VARNMS is a simple string (i.e., not
% a CELLSTR), then FLDS is a simple matrix of type DOUBLE.
%
% Last Saved Time-stamp: <Wed 2013-10-23 10:40:24 Eastern Daylight Time gramer>

  if ( ~exist('varnms','var') || isempty(varnms) )
    %chlor_a,epsilon,sst,sst4,qual_sst,flh,Kd_488_lee,Kd_547_lee,Zeu_lee
    varnms = 'sst';
  end;
  if ( ~iscell(varnms) )
    varnms = {varnms};
    doCells = false;
  else
    doCells = true;
  end;

  lon = [];
  lat = [];

  if ( isa(fname_or_nc,'mDataset') )
    fname = [];
    nc = fname_or_nc;
  else
    fname = fname_or_nc;
    nc = mDataset(fname);
  end;
  clear fname_or_nc;

  for varix = 1:numel(varnms)
    varnm = varnms{varix};
    v = nc{varnm};

    % fld = cast(v(:),'double');
    rawfld = cast(v(:),'double');
    flds{varix} = v.('Intercept') + (v.('Slope').*rawfld);

    if ( varix == 1 && nargout > 1 )
      lms = v.('Limit');
      [nlat,nlon] = size(rawfld);
      lon = linspace(lms(2),lms(4),nlon);
      lat = linspace(lms(3),lms(1),nlat);
    end;
  end;

  if ( ~isempty(fname) )
    close(nc);
    clear nc;
  end;

  if ( ~doCells )
    flds = flds{:};
  end;

return;
