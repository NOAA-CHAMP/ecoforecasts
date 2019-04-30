1;
% Perform wavelet analysis of station sea temperature in relation to other
% data (station wind, station tide height, satellite chlorophyll_a)

if ( ~exist('datpath', 'var') || isempty(datpath) )
  datpath = './data';
end;
if ( ~exist('figpath', 'var') || isempty(figpath) )
  figpath = './figs';
end;


% stanm = 'fwyf1';
stanm = 'mlrf1';
% stanm = 'smkf1';

if ( ~exist('station', 'var') )
  matfname = fullfile(datpath, [stanm '-ndbc.mat']);

  if ( exist(matfname, 'file') )
    disp(['Reloading station data from ' matfname]);
    load(matfname, 'station');

  else
    disp('Reloading station data from raw NDBC files');
    station = [];
    yrs = 1981:2009;
    for ix = 1:length(yrs)
      yr = yrs(ix);
      fname = sprintf('%s/%sh%d.txt', datpath, lower(stanm), yr);
      if ( exist(fname, 'file') )
        station = load_ndbc_data(station, fname);
      end;
    end;
    save(matfname, 'station');

  end;

  clear matfname yrs ix yr fname;
end;


%
% Process sea temperature, wind speed, and wind direction time series
%

[tdts, tdat] = gap_expand(station.sea_t.date, station.sea_t.data);

[sdts, sdat] = gap_expand(station.wind1_speed.date, station.wind1_speed.data);
[ddts, ddat] = gap_expand(station.wind1_dir.date, station.wind1_dir.data);

% Remove leading and trailing data where ANY variable had NaNs
first_nonnan_idx = find((~isnan(tdat) & ~isnan(sdat) & ~isnan(ddat)), 1);
last_nonnan_idx = find((~isnan(tdat) & ~isnan(sdat) & ~isnan(ddat)), 1, 'last');

tdts = tdts(first_nonnan_idx:end); tdat = tdat(first_nonnan_idx:last_nonnan_idx);
sdts = sdts(first_nonnan_idx:end); sdat = sdat(first_nonnan_idx:last_nonnan_idx);
ddts = ddts(first_nonnan_idx:end); ddat = ddat(first_nonnan_idx:last_nonnan_idx);


% Linearly interpolate across any NaNs within each original time series
badix = (~isfinite(tdat));
tdat(badix) = interp1(tdts(~badix), tdat(~badix), tdts(badix));

badix = (~isfinite(sdat));
sdat(badix) = interp1(sdts(~badix), sdat(~badix), sdts(badix));
badix = (~isfinite(ddat));
ddat(badix) = interp1(ddts(~badix), ddat(~badix), ddts(badix));

% Convert wind direction and speed [kts], to U/V components in [m/s]
wu = sdat .* sind(ddat) .* 0.51444;
wv = sdat .* cosd(ddat) .* 0.51444;

wdts = sdts;
wdat = wu + (i .* wv);

clear badix first_nonnan_idx station sdts sdat ddts ddat wu wv;


disp('Plotting wavelet spectra');

nyears = 2;
yearhours = (nyears*365*24) - 1;

wtw = [wdts(end-yearhours:end) ; wdat(end-yearhours:end)];
wtt = [tdts(end-yearhours:end) ; tdat(end-yearhours:end)];

figure;
set(gcf, 'units', 'normalized', 'outerposition', [0 0 1 1]);
xwt(wtw, wtt, 'BlackandWhite');
datetick('x', 'mmmyy', 'keeplimits','keepticks');
title(sprintf('WXT: wind vs. sea temp. (%s, last %d years)', upper(stanm), nyears));
print('-dpng', fullfile(figpath, sprintf('WXT_%s_wind_seatemp.png', stanm)));
saveas(gcf, fullfile(figpath, sprintf('WXT_%s_wind_seatemp.fig', stanm)));
close gcf

doWTs = true;
if ( doWTs )
  figure;
  set(gcf, 'units', 'normalized', 'outerposition', [0 0 1 1]);
  wt(wtw, 'BlackandWhite');
  datetick('x', 'mmmyy', 'keeplimits','keepticks');
  title(sprintf('WT: wind "u+vi" time series (%s, last %d years)', upper(stanm), nyears));
  print('-dpng', fullfile(figpath, sprintf('WT_%s_wind.png', stanm)));
  saveas(gcf, fullfile(figpath, sprintf('WT_%s_wind.fig', stanm)));
  close gcf

  figure;
  set(gcf, 'units', 'normalized', 'outerposition', [0 0 1 1]);
  wt(wtt, 'BlackandWhite');
  datetick('x', 'mmmyy', 'keeplimits','keepticks');
  title(sprintf('WT: sea temperature time series (%s, last %d years)', upper(stanm), nyears));
  print('-dpng', fullfile(figpath, sprintf('WT_%s_seatemp.png', stanm)));
  saveas(gcf, fullfile(figpath, sprintf('WT_%s_seatemp.fig', stanm)));
  close gcf
end;
