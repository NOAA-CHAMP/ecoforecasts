function [fdts, fdat] = som_vs_pca_ts_preprocess(dts, dat, period, ndays, keepdts, detrend_method)
%function [fdts, fdat] = som_vs_pca_ts_preprocess(dts, dat, period, ndays, keepdts, detrend_method)
%
% Preprocess data in time series 'dat' with timestamps 'dts', for seasonal
% 'period', breaking time series into 'ndays'-long frames. If 'keepdts' is
% a non-empty date vector, only return frames whose dates contain a member of
% 'keepdts'. If 'detrend_method' is a non-empty string, remove the trend from
% each FRAME of data using that method (DEFAULT: 'linear', q.v. DETREND).
% 
% Last Saved Time-stamp: <Wed 2010-11-24 17:18:32 Eastern Standard Time gramer>

  if ( ~exist('keepdts', 'var') )
    keepdts = [];
  end;
  if ( ~exist('detrend_method', 'var') || ~isstr(detrend_method) )
    detrend_method = 'linear';
  end;


  % Fill gaps in dates with interpolated hours, fill gaps in data with NaNs
  [dts, rawdat] = gap_expand(dts, dat);

  % Remove all leading and trailing NaNs
  first_nonnan_idx = find(~isnan(rawdat), 1);
  last_nonnan_idx = find(~isnan(rawdat), 1, 'last');
  dts = dts(first_nonnan_idx:last_nonnan_idx);
  rawdat = rawdat(first_nonnan_idx:last_nonnan_idx);

  % Start at local midnight (GMT-4 or GMT-5), on a 'jweek' boundary, i.e., a
  % jday that is evenly divisble by 7 (with days 365-366 lumped into week 52).
  dvec = datevec(dts);
  hrs = dvec(:,4);
  jdays = datenum(dvec(:,1),dvec(:,2),dvec(:,3)) - datenum(dvec(:,1),1,1) + 1;
  first_midnight_idx = find( ((mod(jdays, 7) == 1) & (jdays < 365) & (hrs == 4)), 1 );
  dts = dts(first_midnight_idx:end);
  rawdat = rawdat(first_midnight_idx:end);


  % Linearly interpolate across any NaNs in the original
  badix = (~isfinite(rawdat));
  rawdat(badix) = interp1(dts(~badix), rawdat(~badix), dts(badix));
  dat = rawdat;

  %
  % Analyze our time series broken up into "frames" (e.g., weekly subsets)
  %

  ncols = (24*ndays);
  nrows = floor(length(dat) / ncols);
  nelts = (nrows * ncols);

  fdts = reshape(dts(1:nelts)', [ncols nrows])';
  fdat = reshape(dat(1:nelts)', [ncols nrows])';

  % If the user requested it, remove the mean and/or trend OF EACH FRAME
  % from all 'ncols' of its hourly values
  switch ( detrend_method )
   case {'constant', 'linear'},
    fdat = detrend(fdat', detrend_method)';
  end;


  %
  % Subset frames further, to those that intersect a given seasonal period!
  %

  endperstr = min(3, length(period));

  switch ( lower(period(1:endperstr)) )
   case 'ann',    jdays_keep = 1:366;

   case 'dry',    jdays_keep = [1:84 (253-ndays):366];
   case 'wet',    jdays_keep = [(71-ndays):266];

   case 'win',    jdays_keep = [1:84 (344-ndays):366];
   case 'spr',    jdays_keep = [(71-ndays):175];
   case 'sum',    jdays_keep = [(162-ndays):266];
   case 'aut',    jdays_keep = [(253-ndays):357];

   case 'jan',    jdays_keep = [001:031];
   case 'feb',    jdays_keep = [032:059];
   case 'mar',    jdays_keep = [060:090];
   case 'apr',    jdays_keep = [091:120];
   case 'may',    jdays_keep = [121:151];
   case 'jun',    jdays_keep = [152:181];
   case 'jul',    jdays_keep = [182:212];
   case 'aug',    jdays_keep = [213:243];
   case 'sep',    jdays_keep = [244:273];
   case 'oct',    jdays_keep = [274:304];
   case 'nov',    jdays_keep = [305:334];
   case 'dec',    jdays_keep = [335:366];

   otherwise,
    error('Unrecognized seasonal period "%s"!', period);
  end;

  % If requested, again subset frames to those bracketing any element of
  % 'keepdts'. (Note: Only DATES are checked: hrs/min/sec are ignored.)
  if ( isempty(keepdts) )
    keepix = 1:size(fdts,1);
  else
    [keepix, ign] = find(ismember(fix(fdts),fix(keepdts)));
    keepix = unique(keepix);
  end;

  % In any case, only keep 'ndays'-long frames that overlap our period!
  % Get the date of the FIRST day of each frame
  [yy mm dd] = datevec(fdts(:, 1));
  jdays = datenum(yy,mm,dd) - datenum(yy,1,1) + 1;

  keepix = intersect(keepix, find(ismember(jdays, jdays_keep)));
  fdts = fdts(keepix,:);
  fdat = fdat(keepix,:);

return;
