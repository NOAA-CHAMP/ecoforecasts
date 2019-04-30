function [stns,RGN] = get_misst_stns(region, radiuskm)
%function [stns,RGN] = get_misst_stns(region, radiuskm)
%
% Load positions of all in situ stations in REGION (e.g., 'asam', 'freef'),
% and construct multi-annual time series of medians over MISST SST for all
% pixels within RADIUSKM surrounding each station location. If RADIUSKM is
% not specified or is zero, indices from the file [REGION '.cfg'] are used:
% if a station name is not found in that file, the nearest neighbor is used.
%
% RETURNS vector STNS of station structs; and RGN struct with .MEANSST
% field, monthly and weekly climatologies .MCLIMSST and .WCLIMSST, and
% vectors of all .LONS and .LATS and grid size .DX (in degs). See HELP
% READ_MISST_REGION for further details on each of these RGN fields.
%
% Before return, STNS and RGNS are SAVEd (qv.) in a .MAT file in datapath.
%
% Last Saved Time-stamp: <Tue 2010-09-28 17:37:20 Eastern Daylight Time gramer>

  datapath = 'data';

  if ( ~ischar(region) )
    error('First arg REGION should be a name string!');
  end;

  if (~exist('radiuskm','var') || isempty(radiuskm))
    radiuskm = 0.0;
  end;


  stns = [];
  RGN = [];


  %asam_stns_9km.mat
  matfname = fullfile(datapath, sprintf('%s_stns_%gkm.mat', region, radiuskm));

  % If we already loaded data one MISST file at a time, don't do it again
  if ( exist(matfname, 'file') )

    disp(['Reloading existing ' matfname]);
    load(matfname, 'stns', 'RGN');

  else

    % Load basic metadata (name,lon,lat,depth) for each station of interest
    % Then load appropriate thermistor and/or station data
    switch ( lower(region) ),
     case {'asam','gbr'},
      stns = load_misst_region_metadata(region);
      % Load thermistor data from local CSV files - where present
      stns = load_clean_thermistor(stns);


     case {'freef','ecarib'},
      warning('off','MISST:NoRegionMetadata');
      stns = load_misst_region_metadata(region);
      warning('on','MISST:NoRegionMetadata');
      if ( isempty(stns) )
        all_stns = get_all_station_metadata;
        for ix = 1:length(all_stns.codes)
          stns(ix).name  = all_stns.codes{ix};
          stns(ix).fname = '';
          stns(ix).lon   = all_stns.lons{ix};
          stns(ix).lat   = all_stns.lats{ix};
          stns(ix).depth = all_stns.depths{ix};
        end;
        all_stns = []; clear all_stns;
      end;
      for ix = 1:length(stns)
        x = load_all_ndbc_data([],stns(ix).name);
        flds = grepstruct(x,'ndbc_');
        for fldix = 1:length(flds)
          fld = flds{fldix};
          stns(ix).(fld) = x.(fld);
        end;
        x = []; clear x;
      end;

     otherwise,
      error('Region "%s" not currently handled!', region);

    end;


    disp(['Loading raw station data for ' region]);

    % Load bleaching observations for each station (if any)
    stns = load_bleach_data(stns);


    disp(['Loading MISST data for ' region]);

    % Load preselected MISST gridpoint indices for each station (if any)
    stns = load_misst_cfg(stns, region);

    warning('OFF', 'MISST:NoRawFile');

    % Get region boundaries for later reference
    [lons,lats,ig,dx] = read_misst_region(region, 0, 0);

    meansst = repmat(0,[numel(lats) numel(lons)]);
    nsst = 0;

    mclimsst = repmat(0,[12 numel(lats) numel(lons)]);
    nmclim = repmat(0, [12 1]);

    wclimsst = repmat(0,[52 numel(lats) numel(lons)]);
    nwclim = repmat(0, [52 1]);

    ix = 0;
    for yr = 2002 : 2009
      % Ignore leap years until we get more data
      for jd = 1 : 365
        ix = ix + 1;
        dt = datenum(yr,1,1, 12,0,0) + jd - 1;
        [ig,mn,dy] = datevec(dt);
        wk = floor((jd-1)/7) + 1; wk(wk > 52) = 52;

        [lon,lat,sst] = read_misst_region(region, yr, jd);
        if ( ~isempty(sst) )
          meansst(~isnan(sst)) = meansst(~isnan(sst)) + sst(~isnan(sst));
          nsst = nsst + 1;
          mclimsst(mn, ~isnan(sst)) = mclimsst(mn, ~isnan(sst)) + sst(~isnan(sst))';
          nmclim(mn) = nmclim(mn) + 1;
          wclimsst(wk, ~isnan(sst)) = wclimsst(wk, ~isnan(sst)) + sst(~isnan(sst))';
          nwclim(wk) = nwclim(wk) + 1;
        end;
        for stnix = 1:length(stns)
          stns(stnix).misst_sst.date(ix,1) = dt;
          stns(stnix).misst_sst.jday(ix,1) = jd;
          if ( isempty(sst) )
            stns(stnix).misst_sst.data(ix,1) = NaN;
          else
            [lonix,latix] = get_misst_indices(lon,lat,stns(stnix),radiuskm);
            dat = sst(latix,lonix);
            stns(stnix).misst_sst.data(ix,1) = nanmedian(dat(:));
          end;
        end;
        sst = []; clear sst;
      end;
      %DEBUG:
      fprintf(2, 'Finished MISST year %d\n', yr);
    end;
    warning('ON', 'MISST:NoRawFile');

    if ( nsst > 0 )
      meansst = meansst ./ nsst;
      meansst(meansst == 0) = nan;
    else
      meansst = repmat(nan, size(meansst));
    end;

    for mn = 1:12
      if ( nmclim(mn) > 0 )
        mclimsst(mn, :) = mclimsst(mn, :) ./ nmclim(mn);
        mclimsst(mn, mclimsst(mn,:) == 0) = nan;
      else
        mclimsst(mn, :) = repmat(nan, size(mclimsst(mn, :)));
      end;
    end;

    for wk = 1:52
      if ( nwclim(wk) > 0 )
        wclimsst(wk, :) = wclimsst(wk, :) ./ nwclim(wk);
        wclimsst(wk, wclimsst(wk,:) == 0) = nan;
      else
        wclimsst(wk, :) = repmat(nan, size(wclimsst(wk, :)));
      end;
    end;


    %RGN: meansst,mclimsst,wclimsst,lons,lats,dx
    RGN.name = region;

    RGN.meansst = meansst;
    RGN.mclimsst = mclimsst;
    RGN.wclimsst = wclimsst;
    RGN.lons = lons;
    RGN.lats = lats;
    RGN.dx = dx;



    %
    % Calculate individual station temperature climatologies
    %

    disp(['Calculating various derived time series for sites in ' region]);

    stns = multi_verify_variable(stns, 'misst_sst_anom');
    stns = multi_verify_variable(stns, 'sea_t_anom');

    disp('Sea temperature climatologies');

    [stns.misst_sst_weekly_clim] = deal([]);
    [stns.misst_sst_weekly_anom] = deal([]);
    [stns.misst_sst_monthly_clim] = deal([]);
    [stns.misst_sst_monthly_anom] = deal([]);
    stns = calc_clims(stns, 'misst_sst');

    [stns.sea_t_weekly_clim] = deal([]);
    [stns.sea_t_weekly_anom] = deal([]);
    [stns.sea_t_monthly_clim] = deal([]);
    [stns.sea_t_monthly_anom] = deal([]);
    stns = calc_clims(stns, 'sea_t');


    % Calculate Degree Heating Weeks

    disp('MISST DHWs');
    [stns.misst_sst_dhw] = deal([]);
    stns = calc_cum_anom(stns, 'misst_sst');

    if ( isfield(stns(1),'sea_t') )
      disp('Thermistor DHWs');
      [stns.sea_t_dhw] = deal([]);
      stns = calc_cum_anom(stns, 'sea_t');
    end;


    disp(['Saving results to ' matfname]);
    save(matfname, 'stns', 'RGN');

  end;


  % Basic QC for the caller - just prints warnings right now
  for ix = 1:length(stns)
    if ( all(isnan(stns(ix).misst_sst.data(:))) )
      warning('NO SST DATA for station %d "%s"...', ix, stns(ix).name);
    end;
  end;


return;



%%%%%%%%%%
%%%%%%%%%% PRIVATE FUNCTIONS
%%%%%%%%%%

function [lonix,latix] = get_misst_indices(lon,lat,stn,radiuskm)

  lonix = stn.misst_lonix;
  latix = stn.misst_latix;
  if ( radiuskm > 0 || isempty(lonix) || isempty(latix) )
    [lonix, latix] = gridnbhd_km(lon,lat,stn.lon,stn.lat,radiuskm);
  end;

return;
