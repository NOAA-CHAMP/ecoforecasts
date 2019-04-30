function stns = subset_misst_cfg(cfgfname,region)
%function stns = subset_misst_cfg(cfgfname,region)
%
% Load MISST configuration file named CFGFNAME, and subset all valid stations
% from this file that are contained within the subregion REGION. (If REGION
% is 'world', all valid stations in CFGFNAME will be returned.) Result STNS
% is a vector of structs, each having fields .station_name, .misst_latix,
% .misst_lonix, .misst_region (==REGION), and .misst_(REGION)_latix and
% .misst_(REGION)_lonix. Last two fields are indices for that station within
% the given REGION. A valid sample format for CFGFNAME would be as follows:
%
%   # The "#" if it is the first non-space character, makes that entire line a comment
%   # Config file for parse-misst executable
%   #  stnam,latidx1,latidx2,lonidx1,lonidx2
%   # Note below: Blank linkes are also OK
%  
%   # SEAKEYS
%   mlrf1,1309,1309,1134,1134
%
%
% Last Saved Time-stamp: <Wed 2010-11-10 21:24:23 Eastern Standard Time gramer>

  if ( ~exist('region','var') || isempty(region) )
    region = 'world';
  end;
  region = lower(region);

  % ( $station, $lat, $lon ) = split(/,/, $line);
  %     $lat += 90.0;
  %     $lat *= 11.37777;
  %     $lon += 180.0;
  %     $lon *= 11.37777;
  %     $lat++;
  %     $lon++;
  %     printf("%s,%.0f,%.0f,%.0f,%.0f\n", $station, $lat, $lat, $lon, $lon);

  %% Sample calls to READ_MISST that should produce identical figures
  % sst = read_misst('data/misst/mw_ir.freef.fusion.2009.365.v02',256,256);
  % figure; maxigraph; contourf(sst,[15:28]); colorbar; title('FReef');
  % ox=2986; oy=1109; ox=ox+8+1; oy=oy+8+1;
  % gsst = read_misst('data/misst/mw_ir.fusion.2009.365.v02',4096,2048);
  % figure; maxigraph; contourf(gsst(oy:oy+255,ox:ox+255),[15:28]); colorbar; title('Global');

  switch ( lower(region) )
   case 'world',  lonoff =    0;  latoff =    0;  lonlen = 4096;  latlen = 2048;
   case 'asam',   lonoff = 2018;  latoff =  687;  lonlen =  256;  latlen =  256;
   case 'freef',  lonoff = 2986;  latoff = 1109;  lonlen =  256;  latlen =  256;
   case 'ecarib', lonoff = 3190;  latoff = 1029;  lonlen =  256;  latlen =  256;
   case 'gbr',    lonoff = 1529;  latoff =  653;  lonlen =  256;  latlen =  256;
   otherwise,     error('Unrecognized region "%s"!', region);
  end;

  lonoff = lonoff+8+1;
  latoff = latoff+8+1;

  nstns = 0;
  stns = [];

  if ( ~exist(cfgfname,'file') )
    error('File not found! "%s"',cfgfname);
  end;

  fid = fopen(cfgfname,'r');
  while ( fid > 0 && isempty(ferror(fid)) )
    cfgline = fgetl(fid);
    if ( cfgline == -1 )
      break;
    elseif ( ~isempty(cfgline) && isempty(regexp(cfgline,'^ *[#]')) )
      flds = textscan(cfgline, '%[^,],%[^,],%[^,],%[^,],%[^,]');
      if ( numel(flds) ~= 5 )
        warning('Malformed config line: "%s"',cfgline);
      else
        latix = [str2num(flds{2}{:}):str2num(flds{3}{:})] + 1;
        lonix = [str2num(flds{4}{:}):str2num(flds{5}{:})] + 1;

        % My poor dyslexic brain is hurting... Freagin' grids, freagin' MATLAB
        lonix = lonix - 2048;
        lonix(lonix < 0) = lonix(lonix < 0) + 4096;

        % If station is inside this region, collect it
        rlatix = latix - latoff;
        rlonix = lonix - lonoff;
        if ( 0 < max(rlatix) && min(rlatix) <= latlen && 0 < max(rlonix) && min(rlonix) <= lonlen )
          rlatix(rlatix < 1) = [];          rlatix(rlatix > latlen) = [];
          rlonix(rlonix < 1) = [];          rlonix(rlonix > lonlen) = [];

          nstns = nstns + 1;
          stns(nstns).station_name = flds{1}{:};

          stns(nstns).misst_latix = latix;
          stns(nstns).misst_lonix = lonix;
% % ???
          stns(nstns).g2_misst_latix = latix;
          stns(nstns).g2_misst_lonix = lonix-1;
%           stns(nstns).g2_misst_latix = latix-1;
%           stns(nstns).g2_misst_lonix = lonix;
%           stns(nstns).g2_misst_latix = latix-1;
%           stns(nstns).g2_misst_lonix = lonix-1;

          stns(nstns).region = region;
          if ( ~strcmp(region,'world') )
            stns(nstns).(['misst_' region '_latix']) = rlatix;
            stns(nstns).(['misst_' region '_lonix']) = rlonix;
% % ???
% %             stns(nstns).(['g2_misst_' region '_latix']) = rlatix;
% %             stns(nstns).(['g2_misst_' region '_lonix']) = rlonix-1;
%             stns(nstns).(['g2_misst_' region '_latix']) = rlatix-1;
%             stns(nstns).(['g2_misst_' region '_lonix']) = rlonix;
            stns(nstns).(['g2_misst_' region '_latix']) = rlatix-1;
            stns(nstns).(['g2_misst_' region '_lonix']) = rlonix-1;
          end;
        end;
      end;
    end;
  end;
  fclose(fid);

return;
