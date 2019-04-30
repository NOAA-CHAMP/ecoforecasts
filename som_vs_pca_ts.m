function [station,sm,pc,sc,framedts,fldnms,fhs] = som_vs_pca_ts(stnam_or_station,periods,ndays,mapdims,normalrng,nodisplay,fldnms,fldabbr)
%function [station,sm,pc,sc,framedts,fldnms,fhs] = som_vs_pca_ts(stnam_or_station,periods,ndays,mapdims,normalrng,nodisplay,fldnms,fldabbr)
%
% Analyze station time series for correlations between atmospheric forcing
% and sea temperature variability (e.g., defined as 40HLP of 1-day s.dev.).
% If STNAM_OR_STATION is a struct, assume it contains all our time series.
% Otherwise, it must be a 5-char station code string (eg 'smkf1') to load.
% Only data within season named PERIODS (optionally a cellstr) are shown.
% E.g., if PERIODS='Summer', only frames of Summer data are examined; if
% PERIODS='Annual', all data are used. If PERIODS is instead a cell array of
% strings, then multiple periods are analyzed. ('Wet' and 'Dry' periods are
% also available for Florida "wet" and "dry" seasons of the year, resp.)
% Data are analyzed not as one long time series, but in NDAYS-long frames.
% A Self-Organizing Map of MAPDIMS (DEFAULT: [3 3]) nodes is constructed and
% displayed, along with a Principal Component Analysis showing the same
% number of EOF modes. If NORMALRNG is given and is a 2-value vector, then
% frames are only analyzed if our "response" variable lies OUTSIDE the range
% NORMALRNG during any day of that frame. Specify scalar NORMALRNG and that
% PERCENTILE will be used as an (upper) cutoff for NORMALRNG, i.e., NORMALRNG
% will be [0 PRCTILE(data,NORMALRNG]'. Otherwise, if NORMALRNG is a vector or
% a FUNCTION_HANDLE accepting a TS and returning indices, analyze only those
% indices.  If NODISPLAY is true (DEFAULT: false), then export graphs in TIFF
% file format, but do not display them to screen.
%
% RETURN VALUES: STATION is a struct with all the raw and processed time
% series loaded for the given station from NDBC files; SM is a cell array
% of structs returned by calls to SOM_MAKE (qv.); PC is a cell array of
% Principal Components (products of PCA eigenvectors and data); SC is a cell
% array of PC scores (PC vector components for each date frame); and finally,
% FRAMEDTS is a cell array of N x nhrs matrices, containing complete hourly
% datenums of all frames corresponding to each PERIOD (see above). The size
% of SM, PC, and FRAMEDTS will equal the size of the input PERIODS: if
% PERIODS is a simple string rather than a cell array of strings, SM will be
% a simple struct, and PC and FRAMEDTS simple numeric matrices.  Return value
% FHS is a vector of all figure handles created for displayed plots.
%
% NOTE: For SOM Unit (map-node, "mode") N, FRAMEDTS((SM.bmus==N),[1 end])
% are the start- and end-dates for all the frames best matched by that Unit.
%
% EXAMPLE: Time series and all dates of matching BMUs for modes [N1 ... NM]
% may be plotted by, e.g., ANNOTBMU(STATION,SM,[1 4 7 8],FRAMEDTS,FLDNMS).
%
% CALLS:
%   PRCTILE, COV *or* CORR, PCACOV *or* PRINCOMP (Statistics Toolbox)
%   SOM_BMUS,SOM_HITS,SOM_MAKE (http://www.cis.hut.fi/projects/somtoolbox)
%   PLOTS,MULTIPLOT,SUPLABEL (MATLAB Central - File Exchange)
%   LOAD_ALL_NDBC_DATA,VERIFY_VARIABLE,DATETICK3 (Ecoforecasts toolbox)
%
% Last Saved Time-stamp: <Wed 2013-07-31 16:54:16 Eastern Daylight Time Lew.Gramer>


%%%% DEBUG:
  doSave = false;
  doPrint = false;

set_more off

% Store data and figures in this M-file's local directory
datapath = get_ecoforecasts_path('data');
figspath = get_ecoforecasts_path('figs');


%
% Initial processing of input and output ARGUMENTS
%

station = [];
sm = [];
pc = [];
framedts = [];

fhs = [];
nfigs = 0;


if ( ~exist('stnam_or_station', 'var') || isempty(stnam_or_station) )
  error('Missing first arg: must be a station name or data struct!');
elseif ( isstruct(stnam_or_station) )
  station = stnam_or_station;
  stnam = station.station_name;
elseif ( ischar(stnam_or_station) )
  stnam = stnam_or_station;
  station = []; % We'll have to load data below...
else
  error('Bad first arg: must be a station name or data struct!');
end;
clear stnam_or_station;

if ( ~exist('periods', 'var') || isempty(periods) )
  periods = 'Annual';
end;
if ( ~iscell(periods) )
  periods = { periods };
end;
periods = lower(periods);

if ( ~exist('ndays', 'var') || isempty(ndays) )
  ndays = 14; % Days per "frame"
end;
nhrs = (ndays*24);

if ( ~exist('mapdims', 'var') || isempty(mapdims) )
  mapdims = [3 3];
end;
nmaps = mapdims(1) * mapdims(2);

%if ( ~exist('normalrng', 'var') || ~isnumeric(normalrng) )
if ( ~exist('normalrng', 'var') )
  normalrng = [];
end;

if ( ~exist('nodisplay', 'var') || isempty(nodisplay) )
  nodisplay = false;
end;
graphvis_status = get(0, 'DefaultFigureVisible');
if ( nodisplay )
  % Disable graph display to the screen - good for batch runs
  set(0, 'DefaultFigureVisible', 'off');
end;



%
% Station (NDBC) environmental variables on which to train SOM and do PCA.
% (NOTE **LAST** variable given is assumed to be the "response" variable.)
%


if ( ~exist('fldnms', 'var') || isempty(fldnms) )
  idx = 0;
  %simple_ndbc_erai_erai_30a_net_flux_term_3_d_sum
  idx=idx+1; fldnms{idx} = 'ndbc_wind1_speed_3_day_average';
  fldabbr{idx} = '\mu_3_d(W)';
  idx=idx+1; fldnms{idx} = 'ndbc_wind1_u_3_day_deviation_sum_ndbc_wind1_v';
  fldabbr{idx} = '\sigma_3_d(W_U)+\sigma_3_d(W_V)';
  idx=idx+1; fldnms{idx} = 'ndbc_air_t_1_day_deviation_3_day_average';
  fldabbr{idx} = '\mu_3_d(\sigma_1_d(T_a_i_r))';
  idx=idx+1; fldnms{idx} = 'ndbc_sea_t_1_day_deviation_3_day_average';
  fldabbr{idx} = '\mu_3_d(\sigma_1_d(T_s_e_a))';
elseif ( strcmpi(fldnms,'default') )
  idx = 0;
  idx=idx+1; fldnms{idx} = 'ndbc_wind1_speed_3_day_average';
  fldabbr{idx} = '\mu_3_d(W)';
  idx=idx+1; fldnms{idx} = 'ndbc_wind1_u_3_day_deviation_sum_ndbc_wind1_v';
  fldabbr{idx} = '\sigma_3_d(W_U)+\sigma_3_d(W_V)';
  idx=idx+1; fldnms{idx} = 'ndbc_air_t_1_day_deviation_3_day_average';
  fldabbr{idx} = '\mu_3_d(\sigma_1_d(T_a_i_r))';
  idx=idx+1; fldnms{idx} = 'ndbc_sea_t_1_day_deviation_3_day_average';
  fldabbr{idx} = '\mu_3_d(\sigma_1_d(T_s_e_a))';
end;

if ( ~exist('fldabbr', 'var') || isempty(fldabbr) )
  fldabbr = strrep(fldnms,'_','\_');
end;

nvars = numel(fldnms);

if ( nvars ~= numel(fldabbr) )
  error('Field abbreviations specified did not match fields');
end;

%
% Load station data if it does not already exist in memory
%

if ( ~exist('station', 'var') || isempty(station) )

  station = load_all_ndbc_data(station, stnam, 1981:2009);

end;


%
% Make sure our variables have been loaded or CALCULATED
%

for fidx = 1:nvars
  station = verify_variable(station, fldnms{fidx});
end;

% OPTIONAL CODE: PLOT INPUT TIME SERIES FOR REVIEW...
doMultiplot = false;
if ( doMultiplot && ~nodisplay )
  for fidx = 1:nvars
    multidts{fidx} = station.(fldnms{fidx}).date;
    multidat{fidx} = station.(fldnms{fidx}).data;
  end;
  multiplot(multidts, multidat, ...
            'YLabel', fldabbr, ...
            'LineSpec', {'b-','g-','r-','c-'}, ...
            'Title', [stnam ': UNprocessed data for analysis']);
  datetick3; drawnow;
  disp('Hit any key to continue...'); pause;
  clear multidts;
  clear multidat;
end;

% OPTIONAL CODE: DO STATISTICAL DISTRIBUTIONS (BOXPLOT ETC.)
doStatPlots = false;


%
% Restrict ourselves to date ranges where all variables overlap
%

dtmins = [];
dtmaxs = [];
for fidx = 1:nvars
  dtmins = [dtmins min(station.(fldnms{fidx}).date)];
  dtmaxs = [dtmaxs max(station.(fldnms{fidx}).date)];
end;

begdt = max(dtmins);
enddt = min(dtmaxs);

[begyr,ign,ign,ign,ign,ign] = datevec(begdt);
[endyr,ign,ign,ign,ign,ign] = datevec(enddt);

for fidx = 1:nvars
  begix = find( (station.(fldnms{fidx}).date >= begdt), 1, 'first' );
  endix = find( (station.(fldnms{fidx}).date <= enddt), 1, 'last' );

  dts{fidx} = station.(fldnms{fidx}).date(begix:endix);
  dat{fidx} = station.(fldnms{fidx}).data(begix:endix);
end;




% If user gave us a "non-critical range" for our LAST variable (e.g., [0 0.4]
% sea temperature variability range OUTSIDE WHICH "onshore flux" events are
% considered to occur), find dates when any hourly value lies OUTSIDE range.

% If user specified a single scalar in NORMALRNG, we assume they want
% a single (upper) cutoff based on the PERCENTILE of that scalar
outlierstr = '';

keepdts = [];

switch ( numel(normalrng) )
 case 0,
  normalrngstr = '';
  normalrngfname = '';

 case 1,
  if ( isa(normalrng,'function_handle') )
    % Keep only dates whose indices are returned by the FUNCTION_HANDLE NORMALRNG
    x.date = dts{end};
    x.data = dat{end};
    keepix = normalrng(x);
    x=[]; clear x
    if ( isempty(keepix) )
      error('Function NORMALRNG (%s) returned empty set!',char(normalrng));
    end;
    keepdts = unique(floor(dts{end}(keepix)));
    %DEBUG:
    disp(['Unique Dates ',num2str(numel(keepdts))]);
    normalrngstr = sprintf('(%s)', char(normalrng));
    normalrngfname = sprintf('-%s', char(normalrng));

  else
    if ( 51 > normalrng || normalrng > 100 )
      error('If arg NORMALRNG is scalar, it must be a PERCENTAGE > 50(%)');
    end;
    normalrngstr = sprintf('(>%g%%)', normalrng);
    outlier = prctile(dat{end}, normalrng);
    outlierstr = sprintf(' (>%0.3g)', outlier);
    fprintf('Considering only %s > %0.3g...\n', fldnms{end}, outlier);
    normalrng = [min(dat{end}) outlier];
    normalrngstr = sprintf('(%g,%g)', normalrng(1), normalrng(2));
    normalrngfname = sprintf('-%0.3f-%0.3f', 0.0, normalrng(2));
  end;

 case 2,
  normalrngstr = sprintf('(%g,%g)', normalrng(1), normalrng(2));
  normalrngfname = sprintf('-%0.3f-%0.3f', normalrng(1), normalrng(2));

 otherwise,
  if ( numel(normalrng) > numel(dat{end}) )
    error('Index NORMALRNG has too many elements!');
  else
    keepdts = unique(floor(normalrng));
  end;
end;

if ( numel(normalrng) == 2 )
  % Keep only dates where our LAST variable's data lie outside NORMALRNG
  keepix = find(normalrng(1) >= dat{end} | dat{end} >= normalrng(2));
  keepdts = unique(floor(dts{end}(keepix)));
end;



% Plot color-order - stupid, stupid MATLAB
color_order = 'bgrcmykw';

dependent_fldnmstr = sprintf('%s (%c) %s', fldnms{end}, color_order(nvars), normalrngstr);
independent_fldnmstr = ' vs.';
for fidx = (nvars-1):-1:1
  independent_fldnmstr = sprintf('%s %s (%c)', independent_fldnmstr, fldnms{fidx}, color_order(fidx));
  if ( fidx > 1 ); independent_fldnmstr = [independent_fldnmstr ',']; end;
end;

fprintf('Analyzing top %d modes, %s time series %s\n%s\n', ...
        nmaps, upper(stnam), dependent_fldnmstr, independent_fldnmstr);




%
% For each seasonal period user requested, subset data and do SOM and PCA
%

for perix = 1:length(periods)

 period = periods{perix};

 period(1) = upper(period(1));
 fprintf('Station %s: Period is %s, %dH (%dd) frames, [%d %d] map\n', ...
         stnam, period, nhrs, ndays, mapdims(1), mapdims(2));


 %
 % Subset each variable, and concatenate them to form our 'vectors'
 %

 fdat = [];
 fdts = [];
 for fidx = nvars:-1:1
   [fdt, fda] = som_vs_pca_ts_preprocess(dts{fidx}, dat{fidx}, period, ndays, keepdts, 'none');

   if ( doSave )
     fid = fopen([stnam '-stats.csv'], 'a');
     if ( fid > 0 )
       fprintf(fid, '"%s,"%s,"%s,%d,%g,%g,%g,%g,%g,%g\n', ...
               period, normalrngstr, fldnms{fidx}, ...
               size(fda,1),mean(fda(:)),std(fda(:)),median(fda(:)),min(fda(:)), ...
               max(fda(:)),iqr(fda(:)));
       fclose(fid);
     end;
   end;

   if ( fidx == nvars )
     framedts{perix} = fdt;
   end;

   fdts = [fdt , fdts];
   fdat = [fda , fdat];
   if ( doMultiplot )
     multidts{fidx} = fdt(:);
     multidat{fidx} = fda(:);
   end;

   if ( doStatPlots )
     boxdat(1:numel(fda), fidx) = fda(:);
   end;

   clear fdt;
   clear fda;
 end;
 if ( doMultiplot && ~nodisplay )
   multiplot(multidts, multidat, ...
             'YLabel', fldabbr, ...
             'LineSpec', {'b.','g.','r.','c.'}, ...
             'Title', [stnam ': Preprocessed data for analysis']);
   datetick3; drawnow;
   disp('Hit any key to continue...'); pause;
 end;
 clear multidts;
 clear multidat;

 % If requested, do statistical box plot of distributions for each variable
 if ( doStatPlots )
   disp('Doing statistical plots...');

   nfigs = nfigs + 1;
   fhs(nfigs) = figure;
   set(gcf, 'units', 'normalized', 'outerposition', [0 0 1 1]);
   set(gcf, 'DefaultAxesFontSize', [8]);
   boxplot(boxdat, 'notch','on', 'whisker',3, 'labels',fldabbr);
   title(sprintf( 'Box Plots %s %d-%d %s (N=%d)', ...
                  upper(stnam), begyr, endyr, period, size(boxdat,1) ));
   pfname = sprintf('%s/boxplots-%s-%s%s', figspath, ...
                    lower(stnam), lower(period), normalrngfname);
   if ( doPrint )
     print('-dtiff', '-r300', [pfname '.tiff']);
   end;
   if ( nodisplay )
     % No point in leaving figure around, if user can't see it
     close(gcf);
   end;

   boxdat = [];
   clear boxdat;
 end;



 %
 % Self-Organizing Map (Artificial Neural Network) training and analysis
 %

 fprintf('Training Self-Organizing Map...\n');

 sm{perix} = som_make(fdat, 'msize',mapdims, 'lininit', 'batch', 'tracking',0, ...
                      'shape','sheet', 'training','long', 'neigh','ep');

 sm{perix}.station_name = stnam;
 sm{perix}.period = period;
 sm{perix}.normalrng = normalrng;

 % Add Best-Matching Units vector as a member of each 'sm' struct.
 % NOTE: For SOM Unit (map-node, "mode") 'n', framedts((sm.bmus==n),[1 end])
 % are the start- and end-dates for all the frames best matched by that Unit.
 [sm{perix}.bmus, sm{perix}.qerrs] = som_bmus(sm{perix}, fdat);

 % Calculate number of matched sample frames for each "mode" (Unit).
 permatched = som_hits(sm{perix}, fdat);
 pctmatched = permatched * (100 / size(fdat,1));

 % Sort BMUs by number of frames, from most-matched to least.
 [ign, mode_order] = sort(permatched, 'descend');

 nfigs = nfigs + 1;
 fhs(nfigs) = figure;

 set(gcf, 'units', 'normalized', 'outerposition', [0 0 1 1]);
 set(gcf, 'DefaultAxesFontSize', [8]);
 cmins = min(sm{perix}.codebook); cmn = min(reshape(cmins, [nhrs nvars]));
 cmaxs = max(sm{perix}.codebook); cmx = max(reshape(cmaxs, [nhrs nvars]));
 cmin = min([-abs(cmn) ; -abs(cmx)]); cmax = max([abs(cmn) ; abs(cmx)]);
 cmin = cmn; cmax = cmx;
 %%%% ???
 cmin = repmat(0, size(cmn)); cmax = cmx;
 for orderix = 1:length(mode_order)
   % Plot "modes" from most frequently matched to least
   ix = mode_order(orderix);
   subplot(mapdims(2), mapdims(1), orderix);
   hold on;
   x = 1:nhrs;
   cdbk = reshape(sm{perix}.codebook(ix,:), [nhrs nvars]);
   do_plots(x, cdbk, cmin, cmax);
   hold off;
   title(sprintf('#%d (%d frames - %0.1f%%)', ...
                 ix, permatched(ix), pctmatched(ix)));
 end;
 suplabel(strrep(lower(independent_fldnmstr),'_','\_'), 'x');
 suplabel(sprintf( 'SOM "Modes" %s %d-%d %s, %dH frames (N=%d): %s%s', ...
                   upper(stnam), begyr, endyr, period, nhrs, size(fdat,1), ...
                   strrep(lower(dependent_fldnmstr),'_','\_'), outlierstr ), 't');
 pfname = sprintf('%s/SOM-%s-%s%s-%d-%d', figspath, ...
                  lower(stnam), lower(period), normalrngfname, mapdims(1), mapdims(2));
 if ( doPrint )
   print('-dtiff', '-r300', [pfname '.tiff']);
 end;
 if ( nodisplay )
   % No point in leaving figure around, if user can't see it
   close(gcf);
 end;


 fprintf('Doing complementary Principal Component Analysis...\n');

 % PCA should be done on anomaly values, rather than raw time series
 fanom = fdat;

 % Remove the TIME AVERAGE of each variable from each hour of our canonical
 % frame, e.g., for a 7-day (168-hour) "frame" with TWO variables, subtract
 % 168 x 2 = 336 mean values from the rows of each of our 336 columns.
 varmean = repmat(nanmean(fdat), [size(fdat,1) 1]);
 fanom = fdat - varmean;

 %approach = 'PCACOV(COV)';
 %approach = 'PCACOV(CORR)';
 %approach = 'PRINCOMP';
 approach = 'PRINCOMP(ZSCORE)';
 disp(approach);

 covm = [];
 if ( strcmpi(approach,'PCACOV(COV)') )
   % Calculate covariance matrix
   covm = cov(fanom');
   corresp = 'C';
 elseif ( strcmpi(approach,'PCACOV(CORR)') )
   % Use correlation coefficients instead of covariance matrix
   covm = corr(fanom');
   corresp = '\rho';
 end;

 if ( ~isempty(covm) && any(~isfinite(covm(:))) )

   warning('Non-finite correlation matrix! Skipping PCA...');

 else

   clear evecs score evals explained

   if ( strcmpi(approach,'EOF') )
     % % EOF imlementation of Rich Signell (RPSstuff)
     % [v, d, z] = eof(fanom');
     % EOF implementation of Martijn Hooimeijer (1998)
     [evals, evecs, score] = eof(fanom');
     corresp = 'C';

   elseif ( strcmpi(approach,'PCACOV(COV)') || strcmpi(approach,'PCACOV(CORR)') )
     [evecs, evals, explained] = pcacov(covm);
     score = [];

   elseif ( strcmpi(approach,'PRINCOMP') )
     [evecs, score, evals] = princomp(fanom');
     totalvar = nansum(evals);
     explained = 100.*evals./totalvar;
     corresp = 'C';

   elseif ( strcmpi(approach,'PRINCOMP(ZSCORE)') )
     % Use PRINCOMP to do PCA on normalized (ZSCORE) data vectors
     [evecs, score, evals] = princomp(zscore(fanom'));
     totalvar = nansum(evals);
     explained = 100.*evals./totalvar;
     corresp = '\rho';
   end;

   %vevecs = rotatefactors(evecs(:,1:nmaps));

   pc{perix} = fanom' * evecs;
   pc{perix} = pc{perix}';
   %pc{perix}(badix) = nan;

   %sc{perix} = score;
   sc{perix} = fanom * score;
   sc{perix} = sc{perix}';
   %sc{perix}(badix) = nan;

   ncomps = min(size(pc{perix},1), nmaps);
   totpct = 0;

   nfigs = nfigs + 1;
   fhs(nfigs) = figure;

   set(gcf, 'units', 'normalized', 'outerposition', [0 0 1 1]);
   set(gcf, 'DefaultAxesFontSize', [8]);
   cmins = min(pc{perix}(1:ncomps,:)); cmn = min(reshape(cmins, [nhrs nvars]));
   cmaxs = max(pc{perix}(1:ncomps,:)); cmx = max(reshape(cmaxs, [nhrs nvars]));
   cmin = min([-abs(cmn) ; -abs(cmx)]); cmax = max([abs(cmn) ; abs(cmx)]);
   cmin = cmn; cmax = cmx;
   for ix = 1:ncomps
     subplot(mapdims(2), mapdims(1), ix);
     hold on;
     x = 1:nhrs;
     pcadats = reshape(pc{perix}(ix,:), [nhrs nvars]);
     do_plots(x, pcadats, cmin, cmax);
     hold off;
     title(sprintf('#%d (%0.1f%% %s)', ix, explained(ix), '\sigma^2'));
     totpct = totpct + explained(ix);
   end;
   suplabel(strrep(lower(independent_fldnmstr),'_','\_'), 'x');
   suplabel(sprintf( 'PCA(%s) Modes %s %d-%d %s, %dH frames (%3.1f%% of %s): %s%s', ...
                     corresp, upper(stnam), begyr, endyr, period, ...
                     nhrs, totpct, '\sigma^2', ...
                     strrep(lower(dependent_fldnmstr),'_','\_'), outlierstr ), 't');
   pfname = sprintf('%s/PCA-%s-%s%s-%d-%d', figspath, ...
                    lower(stnam), lower(period), normalrngfname, mapdims(1), mapdims(2));
   if ( doPrint )
     print('-dtiff', '-r300', [pfname '.tiff']);
   end;
   if ( nodisplay )
     % No point in leaving figure around, if user can't see it
     close(gcf);
   end;

 end;

end;

% Fix up return values if necessary
if ( numel(sm) == 1 )
  sm = sm{1};
end;
if ( numel(pc) == 1 )
  pc = pc{1};
end;
if ( numel(sc) == 1 )
  sc = sc{1};
end;
if ( numel(framedts) == 1 )
  framedts = framedts{1};
end;


% Restore previous state of graph display to the screen
set(0, 'DefaultFigureVisible', graphvis_status);

% Done
set_more;

return;


%%
% PRIVATE function [hl, he] = do_plots(x, Y, cmin, cmax)
% CALLS: 'plots' from MATLAB Central Exchange
function [hl, ha] = do_plots(x, Y, cmin, cmax)
  % ylabs = strjust(strvcat('Hour','SeaT','AirT','WindU','WindV'),'center');
  % [hl, ha] = plots(x, Y, 'left', ylabs, 5); 
  [hl, ha] = plots(x, Y, 'left', 8); 
  % stylel = repmat(strvcat('--','-.','-',':'), [2 1]);
  % stylel = repmat(strvcat('--','-'), [4 1]);
  stylel = repmat(strvcat('-'), [8 1]);
  colorl = get(gcf, 'DefaultAxesColorOrder');
  set(ha(end), 'YAxisLocation', 'right');
  for ix = 1:numel(hl)
    set(hl(ix), 'LineStyle', stylel(ix,:), 'Color', colorl(ix,:));
    set(ha(ix), 'XGrid', 'on');
    set(ha(ix), 'YColor', colorl(ix,:))
    if ( exist('cmin', 'var') && exist('cmax', 'var') )
      % NOTE: If cmin and cmax are present but empty, NO YLim is set...
      if ( ~isempty(cmin) && ~isempty(cmax) )
        if ( cmin(ix) < cmax(ix) )
          ylim(ha(ix), [cmin(ix) cmax(ix)]);
        else
          % Needed when one component (e.g., latent heat flux) is always negative
          ylim(ha(ix), [min([cmin(ix) cmax(ix)]),max([cmin(ix) cmax(ix)])]);
        end;
      end;
    else
      ylm = max(abs(ylim(ha(ix))));
      ylim(ha(ix), [-ylm ylm]);
    end;
  end;
return;
