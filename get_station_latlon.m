function stns = get_station_latlon(stns,cfgfname)
%function stns = get_station_latlon(stns,cfgfname)
%
% Load station lat/lon file CFGFNAME (DEFAULT: stations.txt in ecoforecasts
% PATH) and extract .lat and .lon fields for each station in vector of
% structs STNS (may also be a single struct). If STNS is not given or empty,
% load station coordinates from CFGFNAME and return them *all* in STNS.
%
% Last Saved Time-stamp: <Mon 2010-11-08 13:45:43 Eastern Standard Time gramer>

  rootpath = get_ecoforecasts_path();

  if ( ~exist('stns','var') || isempty(stns) )
    stns = [];
  end;
  if ( ~exist('cfgfname','var') || isempty(cfgfname) )
    cfgfname = fullfile(rootpath,'stations.txt');
  end;

  if ( ~exist(cfgfname,'file') )
    error('File not found! "%s"',cfgfname);
  end;

  nstns = 0;
  allstns = [];
  fid = fopen(cfgfname,'r');
  if ( fid < 0 )
    error('Unable to open stations file "%s"!',cfgfname);
    ferror(fid),
  end;
  while ( fid > 0 && isempty(ferror(fid)) )
    cfgline = fgetl(fid);
    if ( cfgline == -1 )
      break;
    elseif ( ~isempty(cfgline) && isempty(regexp(cfgline,'^ *[#]')) )
      flds = textscan(cfgline, '%[^,],%[^,],%[^,]');
      if ( numel(flds) ~= 3 )
        warning('Malformed config line: "%s"',cfgline);
      else
        nstns = nstns + 1;
        allstns(nstns).station_name = flds{1}{:};
        allstns(nstns).lat = str2num(flds{2}{:});
        allstns(nstns).lon = str2num(flds{3}{:});
      end;
    end;
  end;
  fclose(fid);

  if ( isempty(stns) )
    stns = allstns;
  else
    for stnix = 1:length(stns)
      stnm = stns(stnix).station_name;
      ix = find(strcmp(stnm,{allstns.station_name}));
      if ( ~isempty(ix) )
        stns(stnix).lat = allstns(ix).lat;
        stns(stnix).lon = allstns(ix).lon;
      else
        warning('Station "%s" not found in "%s"!',stnm,cfgfname);
      end;    
    end;
  end;

return;
