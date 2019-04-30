function [lon,lat,sst,dx] = read_misst_region(region, yr_or_fname, jd)
%function [lon,lat,sst,dx] = read_misst_region(region, yr_or_fname, jd)
%
% Read a raw MISST binary file for one of the REGIONs 'asam', 'freef',
% 'ecarib', 'gbr', 'papa', or 'world' (global) for year YR, year-day JD.
% If 2nd arg is a char, instead read file [MISSTPATH FNAME] and ignore JD.
% Return 1xM vector of longitudes LON, 1xN vector of latitudes LAT, grid
% resolution DX (all in decimal degrees), and the NxM matrix of SST [oC].
%
% Alternatively, calling READ_MISST_REGION(REGION) with only one argument
% will just return the LON, LAT, and DX for that REGION; SST then is empty.
%
% Last Saved Time-stamp: <Sun 2011-04-10 19:33:02 Eastern Daylight Time gramer>

  datapath = get_ecoforecasts_path('data');
  misstpath = fullfile(datapath,'misst');

  lon = [];
  lat = [];
  sst = [];
  dx = 360/4096;

  %%%% Original code per Dr. Chelle Gentemann
  %         i1=2064;         i2=2319;
  %         j1=1153;         j2=1408;
  %         dx=360/4096;
  %         xlon=(i-1)*dx+dx/2;
  %         xlat=(j-1)*dx-(90-dx/2);

  % On 2/7/2011 5:00 PM, chelle gentemann wrote:
  % >         ireg1(1)=2018;ireg1(2)=2986;ireg1(3)=3190;ireg1(4)=1529;
  % >         ireg2(1)=2290;ireg2(2)=3258;ireg2(3)=3462;ireg2(4)=1801;
  % >         jreg1(1)=687;jreg1(2)=1109;jreg1(3)=1029;jreg1(4)=653;
  % >         jreg2(1)=959;jreg2(2)=1381;jreg2(3)=1301;jreg2(4)=925;
  % >         ireg1(5)=2055;
  % >         ireg2(5)=2327;
  % >         jreg1(5)=1144;
  % >         jreg2(5)=1416;
  % >
  % >         areg(1)='asam'; areg(2)='freef'; areg(3)='ecarib'; areg(4)='gbr'; areg(5)='Papa';
  % >
  % > land mask used is the modis landmask.

  switch ( lower(region) )
   case 'world',  lonoff =   -8;  latoff =   -8;  lonlen = 4096;  latlen = 2048;
   case 'asam',   lonoff = 2018;  latoff =  687;  lonlen =  256;  latlen =  256;
   case 'ecarib', lonoff = 3190;  latoff = 1029;  lonlen =  256;  latlen =  256;
   case 'freef',  lonoff = 2986;  latoff = 1109;  lonlen =  256;  latlen =  256;
   case 'gbr',    lonoff = 1529;  latoff =  653;  lonlen =  256;  latlen =  256;
   case 'papa',   lonoff = 2056;  latoff = 1145;  lonlen =  256;  latlen =  256;
   otherwise,     error('Unrecognized region "%s"!', region);
  end;

  % Why did Chelle do this -8 +8 trick? Dunno... just replicating her code.
  lonoff = lonoff + 8;
  latoff = latoff + 8;

  % Following logic matches GHRSST netCDF coords, e.g. (accessed 2011 Apr 10)
  %  http://dods.jpl.nasa.gov/opendap/sea_surface_temperature/GHRSST/data/L4/GLOB/REMSS/mw_ir_OI/2009/228/20090816-REMSS-L4HRfnd-GLOB-v01-fv02-mw_ir_OI.nc.gz
  lon = ( ([lonoff:(lonoff+lonlen-1)]) * dx ) + (dx/2);
  lat = ( ([latoff:(latoff+latlen-1)]) * dx ) + (dx/2) - 90;

  lon(lon > 180) = lon(lon > 180) - 360;

  if ( nargin < 2 )
    %%%%%%%%%%%%%%%
    % EARLY RETURN
    %%%%%%%%%%%%%%%
    return;
  end;


  if ( ischar(yr_or_fname) )
    fname = yr_or_fname;
  elseif ( isnumeric(yr_or_fname) && isnumeric(jd) )
    yr = yr_or_fname;
    if ( strcmpi(region,'world') ); region_dataset = 'fusion';
    else region_dataset = sprintf('%s.fusion', lower(region)); end;
    fname = fullfile(misstpath, ...
                     sprintf('mw_ir.%s.%04d.%03d.v03', region_dataset, yr, jd));
    if ( ~exist(fname,'file') )
      fname = fullfile(misstpath, ...
                       sprintf('mw_ir.%s.%04d.%03d.v02', region_dataset, yr, jd));
      if ( ~exist(fname,'file') )
        fname = fullfile(misstpath, ...
                         sprintf('mw_ir.%s.%04d.%03d.rt', region_dataset, yr, jd));
        if ( exist(fname,'file') )
          warning('MISST:UsingRTFile', 'Using RT file "%s"!', fname);
        end;
      end;
    end;
  else
    error('Second arg must either be a filename or year (1993-present)!');
  end;

  if ( exist(fname, 'file') )
    %DEBUG:    disp([mfilename ': File ' fname]);
    sst = read_misst(fname, lonlen, latlen);
  else
    warning('MISST:NoRawFile', 'No MISST file "%s" found!', fname);
  end;

return;
