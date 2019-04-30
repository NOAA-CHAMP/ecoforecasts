1;

format short; format compact;

more_status = get(0, 'More'); more('off');


% stnam = 'fwyf1';
% stnam = 'mlrf1';
stnam = 'smkf1';

fldnm = 'sea_t';
% fldnm = 'wind1_speed';

filtOrder = 5;

mapdims = [3 3];
nmaps = mapdims(1) * mapdims(2);

fprintf('Analyzing top %d NDBC %s time series %s\n', nmaps, stnam, fldnm);


%
% Load station data if it does not already exist in memory
%

if ( ~exist('station', 'var') || isempty(station) )
  if ( ~exist('datapath', 'var') || isempty(datapath) )
    datapath = './data';
  end;

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

  clear datapath yrs ix yr fname matfname;
end;

% Make sure our variable has been loaded or CALCULATED
station = verify_variable(station, fldnm);

dmin = min(station.(fldnm).date(:));
dmax = max(station.(fldnm).date(:));
[ymin,m,d] = datevec(dmin);
[ymax,m,d] = datevec(dmax);


cmin = min(station.(fldnm).data(:));
cmax = max(station.(fldnm).data(:));


% Fill gaps in dates with interpolated hours, gaps in data with NaNs
[dts, rawdat] = gap_expand(station.(fldnm).date, station.(fldnm).data);

% Remove all leading and trailing NaNs
first_nonnan_idx = find(~isnan(rawdat), 1);
last_nonnan_idx = find(~isnan(rawdat), 1, 'last');
dts = dts(first_nonnan_idx:last_nonnan_idx);
rawdat = rawdat(first_nonnan_idx:last_nonnan_idx);

% Start at local midnight (GMT-4 or GMT-5) on next day
dvec = datevec(dts(1:25));
first_midnight_idx = find( (dvec(:,4) == 4), 1 );
dts = dts(first_midnight_idx:end);
rawdat = rawdat(first_midnight_idx:end);


% Linearly interpolate across any NaNs in the original
badix = (~isfinite(rawdat));
rawdat(badix) = interp1(dts(~badix), rawdat(~badix), dts(badix));


% % % Remove the daily temperature cycle from each day of time series
% % clim = [ ...
% %     0.030 0.036 0.058  0.106 0.169 0.246 ...
% %     0.314 0.367 0.410  0.426 0.427 0.398 ...
% %     0.355 0.295 0.263  0.234 0.207 0.177 ...
% %     0.142 0.112 0.096  0.070 0.051 0.036 ...
% %        ];
% % clim = repmat(clim, [??? ??? nrows ndays]);
% % dat = rawdat - clim;

% % Remove diurnal cycle from interpolated time series with a bandstop filter
% if ( strcmpi(fldnm, 'sea_t') || strcmpi(fldnm, 'tide') )
%   disp('Removing Fourier components of 10-14 and 18-30 hour periods');
%   [B,A] = butter(filtOrder, [(2/14) (2/10)], 'stop');
%   dat = filtfilt(B, A, rawdat);
%   [B,A] = butter(filtOrder, [(2/30) (2/18)], 'stop');
%   dat = filtfilt(B, A, dat);
% end;

dat = rawdat;


%
% Analyze our time series broken up into "frames" (e.g., weekly subsets)
%

ndays = 14; % Days per "frame"
ncols = (24*ndays);
nrows = floor(length(dat) / ncols);
nelts = (nrows * ncols);

fdts = reshape(dts(1:nelts)', [ncols nrows])';
fdat = reshape(dat(1:nelts)', [ncols nrows])';

% Remove the mean and trend OF EACH FRAME from all n of its hourly values
fdat = detrend(fdat', 'linear')';


%
% Subset frames further, to those that intersect a given seasonal period!
%
period = 'Annual';

switch ( lower(period) )
 case 'annual',
  %wks = [1:52];
  jdays = 1:366;
 case 'winter',
  %wks = [1:12 50:52];
  jdays = [1:84 (344-ndays):366];
 case 'spring',
  %wks = [11:25];
  jdays = [(71-ndays):175];
 case 'summer',
  %wks = [24:38];
  jdays = [(162-ndays):266];
 case 'autumn',
  %wks = [37:51];
  jdays = [(253-ndays):357];
 otherwise,
  error('Unrecognized seasonal period "%s"!', period);
end;
fprintf('Seasonal Period is %s\n', period);

[yy mm dd] = datevec(fdts(:, 1));
jd = datenum(yy,mm,dd) - datenum(yy,1,1) + 1;

goodfs = find(ismember(jd, jdays));
fdts = fdts(goodfs,:);
fdat = fdat(goodfs,:);



%
% Self-Organizing Map training and analysis
%

sm = som_make(fdat, 'msize',mapdims, 'lininit', 'batch', ...
              'shape','sheet', 'training','long', 'neigh','ep');

% Calculate number of matched sample frames for each "mode"
bmus = som_bmus(sm, fdat);
permatched = som_hits(sm, fdat);
pctmatched = permatched * (100 / length(bmus));

% Sort BMUs by number of frames, from most-matched to least
[ign, mode_order] = sort(permatched, 'descend');

figure;
set(gcf, 'units', 'normalized', 'outerposition', [0 0 1 1]);
for orderix = 1:length(mode_order)
  % Plot "modes" from most frequently matched to least
  ix = mode_order(orderix);
  subplot(mapdims(2), mapdims(1), orderix);
  plot(1:ncols, sm.codebook(ix,:));
  xlim([1 ncols]);
%   ylim([cmin cmax]);
  title(sprintf('#%d (%d pers - %0.1f%%)', ...
                ix, permatched(ix), pctmatched(ix)));
end;
ax = suplabel(sprintf('SOM "Modes" (%s %s), %04d-%04d %s, %d-day frames (%d frames)', ...
                      upper(stnam), fldabbr, ...
                      ymin, ymax, period, ...
                      ndays, size(fdat,1)), ...
              't', [0.075 0.075 0.850 0.870]);

print('-dpng', fullfile(figspath,['SOM-' lower(stnam) '-' lower(fldnm) '.png']));


fprintf('Doing complementary Principal Component Analysis...\n');

% [v, d, z] = eof(goodsst');
covm = cov(fdat');
[evecs, evals, explained] = pcacov(covm);
pcadat = fdat' * evecs;
pcadat = pcadat';

% (For now, use same scale as SOM outputs) NOT
minc = min(pcadat(:));
maxc = max(pcadat(:));

% pcadat(badix) = nan;


ncomps = min(size(pcadat,1), nmaps);
totpct = 0;

figure;
set(gcf, 'units', 'normalized', 'outerposition', [0 0 1 1]);
for ix = 1:ncomps
  subplot(mapdims(2), mapdims(1), ix);

  plot(1:ncols, pcadat(ix,:));
  xlim([1 ncols]);
%   ylim([cmin cmax]);
  title(sprintf('#%d (%0.1f%% %s)', ix, explained(ix), '\sigma^2'));
  totpct = totpct + explained(ix);
end;
ax = suplabel(sprintf('PCA Modes (%s %s), %04d-%04d %s, %d-day frames (%3.1f%% of %s)', ...
                      upper(stnam), fldabbr, ...
                      ymin, ymax, period, ...
                      ndays, totpct, '\sigma^2'), ...
              't', [0.075 0.075 0.850 0.870]);

print('-dpng', fullfile(figspath,['PCA-' lower(stnam) '-' lower(fldnm) '.png']));



fprintf('Power Spectral Density of signal...\n');

[Pxx, Fs] = pmtm(rawdat);
figure;
set(gcf, 'units', 'normalized', 'outerposition', [0 0 1 1]);
loglog(Fs, Pxx);
% Annotate with common periods
pers = {12.42, 24, 25.819, 27.5, (24*13.661), (24*183), (24*365)};
perstr = {'M_2', '24', '', 'T_i', 'M_f', '183d', '365d'};
for perix = 1:length(pers)
  annotline((2*pi)/pers{perix}, [], perstr{perix});
end;
title(sprintf('Multitaper Spectrum (%s %s), %04d-%04d %s', ...
              upper(stnam), fldabbr, ...
              ymin, ymax, period));

print('-dpng', fullfile(figspath,['PSD-' lower(stnam) '-' lower(fldnm) '.png']));


% Done
more(more_status);
clear more_status;
