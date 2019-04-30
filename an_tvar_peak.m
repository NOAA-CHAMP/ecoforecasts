function [station,sm,pc,framedts,fhs] = an_tvar_peak(stnam_or_station,periods,ndays,mapdims,normalrng,nodisplay)
%function [station,sm,pc,framedts,fhs] = an_tvar_peak(stnam_or_station,periods,ndays,mapdims,normalrng,nodisplay)
%
% Analyze NDBC station data for correlations between atmospheric forcing and
% sea temperature variability (e.g., defined as 40HLP of 1-day std. dev.).
% If 'stnam_or_station' is a struct, assume it contains all our NDBC data.
% Otherwise, it must be a 5-char station code string (eg 'smkf1') to load.
% Only data within the seasonal 'periods' (optionally a cellstr) are shown.
% E.g., if periods='Summer', only frames of Summer data are examined; if
% periods='Annual', all data are used. If 'periods' is instead a cell array
% of strings, then multiple periods are analyzed. ('Wet' and 'Dry' periods
% are also available for Florida "wet" and "dry" seasons of the year, resp.)
% Data are analyzed not as one long time series, but in 'ndays'-long frames.
% A Self-Organizing Map with 'mapdims' (DEFAULT: [3 3]) nodes is constructed
% and displayed, along with a Principal Component Analysis showing the same
% number of EOF modes. If 'normalrng' is given and is a 2-value vector, then
% frames are only analyzed if our "response" variable lies OUTSIDE the range
% 'normalrng' during any day of that frame. Specify scalar 'normalrng' and
% that PERCENTILE will be used as an (upper) cutoff for 'normalrng', i.e., 
% 'normalrng = [0 prctile(data,normalrng]'. If 'nodisplay' is true (DEFAULT:
% false), then export graphs in TIFF file format, but do not display them.
%
% RETURN VALUES: 'station' is a struct with all the raw and processed time
% series loaded for the given station from NDBC files; 'sm' is a cell array
% of structs returned by calls to SOM_MAKE (qv.); 'pc' is a cell array of
% Principal Components (products of PCA eigenvectors and data); finally,
% 'framedts' is a cell array of N x nhrs matrices, containing complete hourly
% datenums of all frames corresponding to each 'period' (see above). The size
% of 'sm', 'pc', and 'framedts' will equal the size of the input 'periods':
% if 'periods' is a simple string rather than a cell array of strings, 'sm'
% will be a simple struct, and 'pc' and 'framedts' simple numeric matrices.
% Return value 'fhs' is a vector of all figure handles created by plots.
%
% NOTE: For SOM Unit (map-node, "mode") 'n', framedts((sm.bmus==n),[1 end])
% are the start- and end-dates for all the frames best matched by that Unit.
%
% Time series and all dates of matching BMUs for modes '[n1 ... nm]' may be
% displayed on screen by, e.g., 'annotbmu(station,sm,[1 4 7 8],framedts)'.
%
% CALLS:
%   PRCTILE - from MATLAB Statistics Toolbox
%   SOM_BMUS,SOM_HITS,SOM_MAKE - http://www.cis.hut.fi/projects/somtoolbox/
%   DATETICK2,PLOTS,MULTIPLOT,SUPLABEL - from MATLAB Central - File Exchange
%   LOAD_ALL_NDBC_DATA,VERIFY_VARIABLE - from ICON/MATLAB 'ecoforecasts' toolbox
%   SET_DATETICK_CURSOR - from Lew's personal stash
%
% Last Saved Time-stamp: <Sat 2009-08-01 14:28:15 Eastern Daylight Time gramer>


% Save MORE state, then turn MORE OFF
more_status = get(0, 'More'); more('off');

% Store data and figures in this M-file's local directory
[pathroot, ig, ig, ig] = fileparts(mfilename('fullpath'));
if ( ~exist('datapath', 'var') || isempty(datapath) )
  datapath = fullfile(pathroot, 'data', '');
end;
if ( ~exist('figspath', 'var') || isempty(figspath) )
  figspath = fullfile(pathroot, 'figs', '');
end;


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

if ( ~exist('normalrng', 'var') || ~isnumeric(normalrng) )
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
% NDBC Station environmental variables on which to train SOM and do PCA.
% (Note **LAST** variable given is assumed to be the "response" variable.)
%

idx = 0;

idx=idx+1; fldnms{idx} = 'wind1_speed_3_day_average';
  fldabbr{idx} = '\mu_3_d(W)';
idx=idx+1; fldnms{idx} = 'wind1_u_3_day_deviation_sum_wind1_v';
  fldabbr{idx} = '\sigma_3_d(W_U)+\sigma_3_d(W_V)';
idx=idx+1; fldnms{idx} = 'air_t_1_day_deviation_3_day_average';
  fldabbr{idx} = '\mu_3_d(\sigma_1_d(T_a_i_r))';
idx=idx+1; fldnms{idx} = 'sea_t_1_day_deviation_3_day_average';
  fldabbr{idx} = '\mu_3_d(\sigma_1_d(T_s_e_a))';

nvars = length(fldnms);

%
% Load station data if it does not already exist in memory
%

if ( ~exist('station', 'var') || isempty(station) )

  matfname = fullfile(datapath, [stnam '-ndbc.mat']);

  if ( exist(matfname, 'file') )
    disp(['Reloading station data from ' matfname]);
    load(matfname, 'station');
    if ( ~isfield(station, 'station_name') || ...
         ~strcmpi(station.station_name, stnam) )
      fprintf(2, '\n\n');
      warning('Station.station_name DID NOT MATCH "%s"', stnam);
      fprintf(2, '\n\n');
      station.station_name = lower(stnam);
    end;

  else
    disp('Reloading station data from raw NDBC files');
    station = load_all_ndbc_data(station, stnam, 1981:2009);
    save(matfname, 'station');

  end;

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
  datetick2; set_datetick_cursor; drawnow;
  disp('Hit any key to continue...'); pause;
  clear multidts;
  clear multidat;
end;

% OPTIONAL CODE: DO STATISTICAL DISTRIBUTIONS (BOXPLOT ETC.)
doStatPlots = true;


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

% If user specified a single scalar in 'normalrng', we assume they want
% a single (upper) cutoff based on the PERCENTILE of that scalar
outlierstr = '';

switch ( numel(normalrng) )
 case 1,
  if ( 51 > normalrng || normalrng > 100 )
    error('If arg "normalrng" is scalar, it must be a PERCENTAGE > 50(%)');
  end;

  normalrngstr = sprintf('(>%g%%)', normalrng);

  outlier = prctile(dat{end}, normalrng);
  outlierstr = sprintf(' (>%0.3g)', outlier);

  fprintf('Considering only %s > %0.3g...\n', fldnms{end}, outlier);
  normalrng = [min(dat{end}) outlier];

  normalrngfname = sprintf('-%0.3f-%0.3f', 0.0, normalrng(2));

 case 2,
  normalrngstr = sprintf('(%g,%g)', normalrng(1), normalrng(2));
  normalrngfname = sprintf('-%0.3f-%0.3f', normalrng(1), normalrng(2));

 otherwise,
  normalrngstr = '';
  normalrngfname = '';
end;


keepdts = [];

if ( numel(normalrng) == 2 )
  % Keep only dates where our LAST variable's data lie outside normalrng
  keepix = find(normalrng(1) >= dat{end} | dat{end} >= normalrng(2));
  keepdts = unique(fix(dts{end}(keepix)));
end;



% Plot color-order - stupid, stupid MATLAB
color_order = 'bgrcmykw';

dependent_fldnmstr = sprintf('%s (%c) %s', fldnms{end}, color_order(nvars), normalrngstr);
independent_fldnmstr = ' vs.';
for fidx = (nvars-1):-1:1
  independent_fldnmstr = sprintf('%s %s (%c)', independent_fldnmstr, fldnms{fidx}, color_order(fidx));
  if ( fidx > 1 ); independent_fldnmstr = [independent_fldnmstr ',']; end;
end;

fprintf('Analyzing top %d modes, NDBC %s time series %s\n%s\n', ...
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
   [fdt, fda] = anndbcs_preprocess(dts{fidx}, dat{fidx}, period, ndays, keepdts, 'none');

   fid = fopen([stnam '-stats.csv'], 'a');
   if ( fid > 0 )
     fprintf(fid, '"%s,"%s,"%s,%d,%g,%g,%g,%g,%g,%g\n', ...
             period, normalrngstr, fldnms{fidx}, ...
             size(fda,1),mean(fda(:)),std(fda(:)),median(fda(:)),min(fda(:)), ...
             max(fda(:)),iqr(fda(:)));
     fclose(fid);
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
   datetick2; set_datetick_cursor; drawnow;
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
   print('-dtiff', '-r300', [pfname '.tiff']);
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
 print('-dtiff', '-r300', [pfname '.tiff']);
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
 varmean = repmat(mean(fdat), [size(fdat,1) 1]);
 fanom = fdat - varmean;


 % [v, d, z] = eof(fanom');
 % covm = cov(fanom');
 % Use correlation coefficients instead of covariance matrix
 covm = corr(fanom');

 if ( any(~isfinite(covm(:))) )

   warning('Non-finite correlation matrix! Skipping PCA...');

 else

   [evecs, evals, explained] = pcacov(covm);
   pc{perix} = fanom' * evecs;
   pc{perix} = pc{perix}';
   %pc{perix}(badix) = nan;

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
                     '\rho', upper(stnam), begyr, endyr, period, ...
                     nhrs, totpct, '\sigma^2', ...
                     strrep(lower(dependent_fldnmstr),'_','\_'), outlierstr ), 't');
   pfname = sprintf('%s/PCA-%s-%s%s-%d-%d', figspath, ...
                    lower(stnam), lower(period), normalrngfname, mapdims(1), mapdims(2));
   print('-dtiff', '-r300', [pfname '.tiff']);
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
if ( numel(framedts) == 1 )
  framedts = framedts{1};
end;


% Restore previous state of graph display to the screen
set(0, 'DefaultFigureVisible', graphvis_status);

% Done - restore MORE state
more(more_status);

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
        ylim(ha(ix), [cmin(ix) cmax(ix)]);
      end;
    else
      ylm = max(abs(ylim(ha(ix))));
      ylim(ha(ix), [-ylm ylm]);
    end;
  end;
return;
