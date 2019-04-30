function [sm, pc] = ndbc2sst(ndbcsm, n, framedts, sm)
%function [sm, pc] = ndbc2sst(ndbcsm, n, framedts, sm)
%
% Produce a multiplot display showing Best Matching Unit(s) 'n' from a Self
% Organizing Map (SOM) analysis, stacked with the raw data upon which the SOM
% was based: allows easy cross-identification of raw-data features (both in
% the atmospheric forcing and in hydrographic "response" variables) that have
% elicited a given set of SOM patterns. 'n' is an integer scalar or vector.
% NOTE: Struct 'ndbcsm' is assumed to have fields 'station_name' and 'period'
% detailing exactly which station and which seasonal period produced the SOM.
%
% Last Saved Time-stamp: <Fri 2009-04-10 17:23:24 Eastern Daylight Time gramer>

  more_status = get(0, 'More'); more('off');

  % Store data and figures in this M-file's local directory
  [pathroot, ig, ig, ig] = fileparts(mfilename('fullpath'));
  if ( ~exist('datapath', 'var') || isempty(datapath) )
    datapath = fullfile(pathroot, 'data', '');
  end;
  if ( ~exist('figspath', 'var') || isempty(figspath) )
    figspath = fullfile(pathroot, 'figs', '');
  end;

  stanm = ndbcsm.station_name;
  period = ndbcsm.period;

  % Which BMUs from the prior run of ANNDBCS are we interested in?
  % (This vector 'n' is the set of Best Matching Units chosen by human review
  %  of a previous Self-Organizing Map analysis on a station's time series.)
  frameix = find(ismember(ndbcsm.bmus, n));
  framestr = sprintf('%d', n);

  [emptysst, LONS, LATS] = ansst(stanm, 'bounds');
  sstdims = size(emptysst);
  clear emptysst;

  dataset = 'mean';

  if ( ~exist('sst', 'var') )

    clear sst;
    sstfname = sprintf('%s-NDBC-%04d-%04d-%s-SST-%s-%s-%04d-%04d.mat', ...
                       stanm, ndbcsm.topol.msize([1 2]), framestr, ...
                       dataset, lower(period), sstdims([1 2]));
    sstfname = fullfile(datapath, sstfname);

    if ( exist(sstfname, 'file') )
      load(sstfname, 'sst');
      fprintf('Reloading SST weekly composites from "%s"\n', sstfname);

    else
      fprintf('Loading weekly SSTs');

      sstix = 1;
      for ix = 1:length(frameix)
        dt = framedts(frameix(ix),1);
        [yr mo dy] = datevec(dt);
        % Satellite data set starts in 1993!
        if ( yr > 1992 )
          jd = datenum(yr,mo,dy) - datenum(yr,1,1) + 1;
          wk = floor((jd-1)/7) + 1;
          x = ansst(stanm, yr, wk, dataset);
          % Turn matrix/map of pixels into a vector of variables
          sst(sstix, 1:numel(x)) = x(:);      sstix = sstix + 1;

          % ???HACK!!! Get second week of record also
          [yr mo dy] = datevec(dt+10);
          jd = datenum(yr,mo,dy) - datenum(yr,1,1) + 1;
          wk = floor((jd-1)/7) + 1;
          x = ansst(stanm, yr, wk, dataset);
          % Turn matrix/map of pixels into a vector of variables
          sst(sstix, 1:numel(x)) = x(:);      sstix = sstix + 1;
        end;
      end;
      fprintf('\n');
      save(sstfname, 'sst');

    end;

  end;

  % A week when every pixel was NaN is a missing or bad image - remove it
  rmrows = all(isnan(sst), 2);
  sst(rmrows, :) = [];

  % If a pixel is NaN for every week in our dataset, it's probably LAND
  landmask = all(isnan(sst), 1)';


  % 'goodsst' is 'sst' with all NaNs (masked values) set to 0
  goodsst = sst;
  goodsst(~isfinite(sst)) = 0;

%   meandim = 'none';
%   tsmean = repmat( 0, size(goodsst) );

  % Remove the time mean from SSTs (a la Mariano et al, 2006)
  meandim = 'space'; % i.e., the mean is a FUNCTION OF x,y
  tsmean = repmat( mean(goodsst,1), [size(goodsst,1) 1] );
  sst = sst - tsmean;
  goodsst = goodsst - tsmean;

%   % Remove the spatial mean from SSTs (a la Mariano et al, 2006)
%   meandim = 'time'; % i.e., the mean is a FUNCTION OF time
%   tsmean = repmat( mean(goodsst,2), [1 size(goodsst,2)] );
%   sst = sst - tsmean;
%   goodsst = goodsst - tsmean;

%   % Remove spatial AND temporal mean from SSTs (a la Weisberg and He, 2006)
%   meandim = 'space(time)';
%   tsmean = repmat( mean(goodsst,1), [size(goodsst,1) 1] );
%   sst = sst - tsmean;
%   goodsst = goodsst - tsmean;
%   tsmean = repmat( mean(goodsst,2), [1 size(goodsst,2)] );
%   sst = sst - tsmean;
%   goodsst = goodsst - tsmean;


  mapdims = [3 2];

  nmaps = mapdims(2) * mapdims(1);

  if ( ~exist('sm','var')  || isempty(sm) || ~isstruct(sm) )
    fprintf('Computing SOM, calculating Best-Match-Units\n');
    % sm = som_make(goodsst, 'msize',mapdims, 'lininit', 'shape','sheet', 'tracking',0, ...
    %               'mask',double(~landmask), 'training','long', 'neigh','ep');

    % Random initialization - slower, less robust, but LESS MEMORY (and NaNs OK)
    sm = som_make(sst, 'msize',mapdims, 'randinit', 'shape','sheet', 'tracking',0, ...
                  'mask',double(~landmask), 'training','long', 'neigh','gaussian');

    sm.station_name = stanm;
    sm.period = period;

    % Calculate the best-matched Unit or "mode" for each week of data
    [sm.bmus, sm.qerrs] = som_bmus(sm, sst);
    % [sm.bmus, sm.qerrs] = som_bmus(sm, goodsst);
  end;

  % Reset all landmasks to NaN in our 'results'
  results = sm.codebook;
  results(1:size(results,1), landmask) = nan;

  fprintf('SOM completed\n');


  fprintf('Plotting SOM "modes"\n');

  % figure;
  % % Do not try to show all 'variables' - there are thousands of them!
  % som_show(sm, 'umat', 'all');
  % som_trajectory(sm.bmus);


  minc = min(results(isfinite(results)));
  maxc = max(results(isfinite(results)));

  % Calculate number of matched weeks for each "mode"
  wksmatched = som_hits(sm, sst);
  % wksmatched = som_hits(sm, goodsst);
  pctmatched = wksmatched * (100 / size(sst,1));

  % Sort BMUs by number of weeks, from most-matched to least
  [ign, mode_order] = sort(wksmatched, 'descend');


  %
  % Plot SOM results in order from most "hits" (matching samples) to least
  %

  figure;
  set(gcf, 'units', 'normalized', 'outerposition', [0 0 1 1]);

  cmin = min(results(:));  cmax = max(results(:));
  for orderix = 1:length(mode_order)
    % Plot "modes" from most frequently matched to least
    ix = mode_order(orderix);
    subplot(mapdims(2), mapdims(1), orderix);

    pcolor(LONS, LATS, reshape(results(ix, :), sstdims));
    shading('flat');
    % For pcolor and similar, let user see SST values when in datacursor mode
    set_pcolor_cursor;

    set(gca, 'clim', [cmin cmax]);

    title(sprintf('#%d (%d wks - %0.1f%%)', ...
                  ix, wksmatched(ix), pctmatched(ix)));
  end;
  % suptitle or suplabel, from MATLAB online contrib
  ax = suplabel(sprintf('%s SOM "Modes" vs. mean(%s), NDBC %s "frames" [%s] (%d weeks)', ...
                        upper(stanm), meandim, period, framestr, size(sst,1)), 't');
  pfname = fullfile(figspath, ...
                    sprintf('%s-SOM-NDBC-%s-%dx%d-%s', stanm, ...
                            period, mapdims(1), mapdims(2), framestr));
  print('-dtiff', '-r300', [pfname '.tiff']);


  %
  % Compute Principal Components / EOFs
  %

  fprintf('Doing complementary Principal Component Analysis\n');

  % [v, d, z] = eof(goodsst');
  covm = cov(goodsst');
  % % Use correlation coefficients instead of covariance matrix
  % covm = corr(goodsst');
  [evecs, evals, explained] = pcacov(covm);
  pc = goodsst.' * evecs;
  pc = pc';
  % PCA scaling is arbitrary - this "fix" suggested by Arthur Mariano
  if ( max(pc(1,:)) < max(tsmean(:)) )
    pc = -pc;
  end;

  pc(~isfinite(sst)) = nan;

  fprintf('Plotting PCA modes\n');

  figure;
  set(gcf, 'units', 'normalized', 'outerposition', [0 0 1 1]);
  % Try to show the same number of PCA modes as we did SOM maps
  totpct = 0;
  ncomps = min(size(pc,1), (mapdims(2)*mapdims(1)));
  cmin = min(min(pc(1:ncomps,:)));  cmax = max(max(pc(1:ncomps,:)));
  for ix = 1:ncomps
    subplot(mapdims(2), mapdims(1), ix);

    pcolor(LONS, LATS, reshape(pc(ix, :), sstdims));
    shading('flat');
    set(gca, 'clim', [cmin cmax]);
    % For pcolor and similar, let user see SST values when in datacursor mode
    set_pcolor_cursor;

    title(sprintf('#%d (%0.1f%% %s)', ix, explained(ix), '\sigma^2'));
    totpct = totpct + explained(ix);
  end;
  ax = suplabel(sprintf('%s PCA Modes vs. mean(%s), NDBC %s "frames" [%s] (%3.1f%% of %s)', ...
                        upper(stanm), meandim, period, framestr, totpct, '\sigma^2'), 't');
  pfname = fullfile(figspath, ...
                    sprintf('%s-PCA-NDBC-%s-%dx%d-%s', stanm, ...
                            period, mapdims(1), mapdims(2), framestr));
  print('-dtiff', '-r300', [pfname '.tiff']);


  more(more_status);

return;
