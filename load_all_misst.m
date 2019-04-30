function stns = load_all_misst(cfgfname,region)
%function stns = load_all_misst(cfgfname,region)
%
% Process all available MISST binary data files of a given REGION or type,
% and create a time series '.misst_sst' (and a time series '.g2_misst_sst',
% see comments below) for each station in config file CFGFNAME (which is
% processed by calling SUBSET_MISST_CFG, v.) The result is a vector of
% structs STNS, and a new MAT file for each station within REGION (or all
% stations, if REGION=='world') named ['misst_' REGION '_' STNM '.mat'].
%
% Last Saved Time-stamp: <Wed 2010-11-10 21:21:00 Eastern Standard Time gramer>

  set_more off
  %DEBUG:
  tic,

  datapath = get_ecoforecasts_path('data');
  misstpath = fullfile(datapath,'misst');

  if ( ~exist('region','var') || isempty(region) )
    region = 'world';
  end;

  stns = subset_misst_cfg(cfgfname,region);
  %DEBUG:
  disp(['Stations to process: ' num2str(numel(stns))]);
  % Get actual lat/lon coordinates
  stns = get_station_latlon(stns);
  for stnix = 1:length(stns)
    %DEBUG:  disp(stns(stnix).station_name);
    stns(stnix).misst_sst.date = [];
    stns(stnix).misst_sst.data = [];
    stns(stnix).g2_misst_sst.date = [];
    stns(stnix).g2_misst_sst.data = [];
  end;

  region = lower(region);

  switch ( region ),
   case 'world',
    % For world, accept best of whatever we have
    fnamepatt = 'mw_ir.fusion.*.(rt|v01|v02|v03)$';
   case 'asam',
    % For any region, accept only region or global 'v02/v03' files
    fnamepatt = ['mw_ir.(' region '.fusion.*.(v02|v03)|fusion.*(v02|v03))$'];
   case 'ecarib',
    fnamepatt = ['mw_ir.(' region '.fusion.*.(v02|v03)|fusion.*(v02|v03))$'];
   case 'freef',
    fnamepatt = ['mw_ir.(' region '.fusion.*.(v02|v03)|fusion.*(v02|v03))$'];
   case 'gbr',
    fnamepatt = ['mw_ir.(' region '.fusion.*.(v02|v03)|fusion.*(v02|v03))$'];
   otherwise,
    error('Unrecognized region "%s"!',region);
  end;

  ds = dir(fullfile(misstpath,'mw_ir.*'));  
  fnames = {ds.name};
  ds = []; clear ds;
  matchix = find(~cellfun(@isempty,regexp(fnames,fnamepatt)));
  fnames = fnames(matchix);

  %DEBUG:  disp(['Filenames to process: ' num2str(numel(fnames))]);
  %DEBUG:  disp(fnames);

  if ( isempty(fnames) )
    warning('No filenames found matching "%s"!',fnamepatt);
  else
    for yr = 2000:get_year(now)
      for jd = 1:366
        rgn = region;
        ix = strmatch(sprintf('mw_ir.%s.fusion.%04d.%03d.v0',region,yr,jd),fnames);
        if ( isempty(ix) )
          rgn = 'world';
          ix = strmatch(sprintf('mw_ir.fusion.%04d.%03d.v03',yr,jd),fnames);
        end;
        if ( isempty(ix) )
          ix = strmatch(sprintf('mw_ir.fusion.%04d.%03d.v02',yr,jd),fnames);
        end;
        if ( isempty(ix) )
          ix = strmatch(sprintf('mw_ir.fusion.%04d.%03d.v01',yr,jd),fnames);
        end;
        if ( isempty(ix) )
          ix = strmatch(sprintf('mw_ir.fusion.%04d.%03d.rt',yr,jd),fnames);
        end;
        if ( isempty(ix) )
          %DEBUG:          disp(sprintf('No file found for %04d.%03d',yr,jd));
          continue;
        end;
        % If by some foulup there are multiple matches, pick the first one
        ix = ix(1);

        if ( strcmp(rgn,'world') )
          lonfld = ['misst_lonix'];
          latfld = ['misst_latix'];
          g2lonfld = ['g2_misst_lonix'];
          g2latfld = ['g2_misst_latix'];
        else
          lonfld = ['misst_' region '_lonix'];
          latfld = ['misst_' region '_latix'];
          g2lonfld = ['g2_misst_' region '_lonix'];
          g2latfld = ['g2_misst_' region '_latix'];
        end;

        fname = fullfile(misstpath,fnames{ix});
        %DEBUG:        if (jd==1||jd>=365); disp(fname); end;
        %DEBUG:
        disp(fname);
        dt = datenum(yr,1,1) + jd - 1;
        [lon,lat,sst,dx] = read_misst_region(rgn, fname);
        if ( isempty(sst) )
          % Call to READ_MISST_REGION should have already produced a Warning
          warning('Skipping %04d.%03d...',yr,jd);
          continue;
        else
          for stnix = 1:length(stns)
            % Generate time series for correct pixel registration
            stnsst = sst(stns(stnix).(latfld),stns(stnix).(lonfld));
            stns(stnix).misst_sst.date(end+1) = dt;
            stns(stnix).misst_sst.data(end+1) = nanmean(stnsst(:));
            % Reproduce the time series being loaded into G2 also
            stnsst = sst(stns(stnix).(g2latfld),stns(stnix).(g2lonfld));
            stns(stnix).g2_misst_sst.date(end+1) = dt;
            stns(stnix).g2_misst_sst.data(end+1) = nanmean(stnsst(:));
          end;
        end;
      end;
    end;
  end;

  disp('Saving all MISST stations to respective MAT files...');
  for stnix = 1:length(stns)
    stnm = stns(stnix).station_name;
    %DEBUG:
    disp(stnm);
    fname = fullfile(datapath,sprintf('misst_%s_%s.mat',region,stnm));
    %DEBUG:    disp(fname);
    station = stns(stnix);
    save(fname,'station');
    station = []; clear station;
  end;

  %DEBUG:
  toc,
  set_more;

return;
