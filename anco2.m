1;
% SCRIPT to analyze hourly pCO2 time series from Discovery Bay ICON station
%
% Last Saved Time-stamp: <Mon 2009-03-23 10:02:58 Eastern Daylight Time gramer>


datapath = './data';
figspath = './figs';

if ( ~exist('dts', 'var') || isempty(dts) )
  [dts, dat] = load_g2_data(fullfile(datapath, 'dbjm1-pco2.csv'), ...
                            [], datenum(2007,12,01));
end;


% [Pxx, Fn] = pwelch(dat);
% figure;
% loglog(Fn, Pxx);
% title('Spectrum - hourly pCO_2');


% Make an effort to remove all data corrupted by periodic "blanking"
badix = find(dat < 0);
bdix = [ badix ; badix+1 ; badix+2 ];
% Also remove any duplicate reports (i.e., nearly equal timestamps)
% (G2 occasionally has some trouble keeping track of time.)
dupes = find(diff(dts) < ((1/24)/2));
bdix = unique([bdix ; dupes]);

bdix(bdix <= 0 | bdix > numel(dat)) = [];

ddat = dat;
ddts = dts;
ddat(bdix) = [];
ddts(bdix) = [];


% Remove "trend" (drift) from cleaned-up time series
dddts = ddts;
dddat = detrend(ddat, 'linear');

% figure;
% plot(dddts, dddat);
% datetick; set_datetick_cursor;
% title('DBJM1 - pCO_2 anomaly (2007)');


% Expand hourly time series to fill every hour of record
% (Uses cubic spline fit to fill gaps in the data record.)
[fdts, fdat] = gap_expand(dddts, dddat);
fdat = gap_fill(fdts, fdat, 'cubic');

% Now we can apply a low-pass filter to isolate multi-day events
hlp = 96;

[B, A] = butter(10, (2/hlp));
fdat = filtfilt(B, A, fdat);


% Trim both ends of time series to remove edge effects (Gibbs)
cutoff = (hlp*2);
fdts = fdts((cutoff+1):end-cutoff);
fdat = fdat((cutoff+1):end-cutoff);


% Plot our result, along with +/- two standard deviation ranges
figure;
set(gcf, 'units', 'normalized', 'outerposition', [0 0 1 1]);
plot(fdts, fdat);
annotline([], [+20 +20], '+2\times\sigma_p_C_O_2')
annotline([], [-20 -20], '-2\times\sigma_p_C_O_2')
xlabel('2007');
ylabel('pCO_2 anomaly - \mu-atmospheres');
datetick; set_datetick_cursor;
xlim([min(fdts) max(fdts)]);
grid on;
title(sprintf('DBJM1 - %d hour mean (%dHLP) detrended pCO_2', ...
              hlp, hlp));
print( '-dbitmap', fullfile(figspath, ...
                            sprintf('Fig-DBJM1-SAMI-CO2-%dHLP.bmp', hlp)) );
print( '-dpng', fullfile(figspath, ...
                         sprintf('Fig-DBJM1-SAMI-CO2-%dHLP.png', hlp)) );


% dfdts = fdts(2:end);
% dfdat = diff(fdat);
% figure;
% plot(dfdts, dfdat);
% datetick; set_datetick_cursor;
% title(sprintf('DBJM1 - %dHLP detrend %s (2007)', hlp, '\deltapCO_2/\deltat'));
