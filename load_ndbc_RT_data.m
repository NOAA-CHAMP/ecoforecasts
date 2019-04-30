function stn = load_ndbc_RT_data(stn, fname)
%function stn = load_ndbc_RT_data(stn, fname)
%
% Load the file named FNAME, and append its data to the station data struct
% STN (optionally an empty matrix [], if we're just starting to load data).
% These archived NDBC C-MAN/SEAKEYS/ICON station data files can be obtained
% for a given the 5-character code 'stnam' from http://www.ndbc.noaa.gov.
%
% NOTE: 99 and 999 are used as both valid values (e.g., wind dir and barom,
% resp.) but also as "bad value" placeholders: this code is careful of that.
%
% N.B.: Unlike LOAD_NDBC_DATA, this function loads the "real-time" (45-day
% old) data files from NDBC: these are full of "MM" values for missing fields
% (those not yet judged to be QC'd by NDBC), and their headers differ from
% those of "regular" year and month QC'd data files from NDBC.
%
% Last Saved Time-stamp: <Mon 2018-10-01 12:51:07 Eastern Daylight Time gramer>


  if ( ~exist(fname, 'file') )
    warning('Skipping non-existent year file "%s"...', fname);
    return;
  end;

  x = importdata(fname);
  if ( ~iscellstr(x) || numel(x) < 3 || (~strncmp(x{1},'#YY',3) && ~strncmp(x{1},'YY',2)) )
    error('Unable to find header in "%s"!', fname);
  end;

  % NDBC headers have been remarkably consistent for >20 years! From...
  % YY MM DD hh WD   WSPD GST  WVHT  DPD   APD  MWD  BAR    ATMP  WTMP  DEWP  VIS
  %  to ...
  % #YY  MM DD hh mm WDIR WSPD GST  WVHT   DPD   APD MWD   PRES  ATMP  WTMP  DEWP  VIS  TIDE
  % #yr  mo dy hr mn degT m/s  m/s     m   sec   sec degT   hPa  degC  degC  degC   mi    ft

  chdr = textscan(x{1},'%s');
  hdr = chdr{1};

  if ( strncmp(x{2},'#yr',3) || strncmp(x{2},'yr',2) )
    % Stride -1 is because data appears in RT file "newest first"
    datrows = x(end:-1:3);
  else
    datrows = x(end:-1:2);
  end;
  
  datrows = strrep(datrows,...
                   '     MM', ...
                   ' 9999.0');
  datrows = strrep(datrows,...
                   '    MM',...
                   ' 999.0');
  datrows = strrep(datrows,...
                   '   MM',...
                   '  999');
  % CAREFUL: Some "direction" fields are only 4 chars wide!
  datrows = strrep(datrows,...
                   '  MM',...
                   ' 999');
  datrows = strrep(datrows,...
                   ' MM',...
                   ' 99');

  % This replaces the Pressure Tendency nonsense with NaN
  datrows = regexprep(datrows,' [+-][0-9][.][0-9][0-9]','  9999');
  datrows = regexprep(datrows,' [+-][0-9][.][0-9]',' 9999');

  % cdat = textscan(char(datrows(1:10)), '%f');
  % dat = cdat{:};
  dat = str2num(char(datrows));


  % Parse timestamp columns

  YYix = find((strcmpi(hdr, 'YEAR') | strcmpi(hdr, 'YR') | strncmpi(hdr, 'YY', 2) | strncmpi(hdr, '#YY', 3)), 1);
  MMix = find((strcmp(hdr, 'MM')), 1);
  DDix = find((strcmpi(hdr, 'DD')), 1);
  hhix = find((strcmpi(hdr, 'HH')), 1);
  mmix = find((strcmp(hdr, 'mm')), 1);
  if ( isempty(YYix) || isempty(MMix) || isempty(DDix) || isempty(hhix) )
    error('Needed timestamp header fields not found in "%s"!', fname);
  end;


  yrs = dat(:,YYix);
  % Oh those awful hairdos
  the90s = find(50 <= yrs & yrs < 100);
  yrs(the90s) = yrs(the90s) + 1900;
  % Mediocre mutual fund, but hopefully a great epoch
  newCentury = find(0 <= yrs & yrs < 50);
  yrs(newCentury) = yrs(newCentury) + 2000;

  if ( isempty(mmix) )
    mns = 0;
  else
    mns = dat(:,mmix);
  end;

  dts = datenum(yrs, dat(:,MMix), dat(:,DDix), dat(:,hhix), mns, 0);


  % Parse data columns (of interest)

  % (Fields generally found in NDBC "h" Standard meteorological data files)
  dirix = find((strcmpi(hdr, 'DIR') | strncmpi(hdr, 'WD', 2)), 1);
  % NOTE: Wind Speed is loaded in [m/s]; LOAD_ALL_NDBC_DATA (v.) will later
  % convert this value to [kts], in order to "be consistent with raw data".
  spdix = find(strcmpi(hdr, 'WSPD'), 1);
  gstix = find((strcmpi(hdr, 'GST') | strncmpi(hdr, 'WG', 2)), 1);
  baromix = find((strcmpi(hdr, 'PRES') | strncmpi(hdr, 'BAR', 3)), 1);
  airtix = find((strcmpi(hdr, 'ATMP') | strncmpi(hdr, 'AIR', 3)), 1);
  seatix = find((strcmpi(hdr, 'WTMP') | strncmpi(hdr, 'SEA', 3)), 1);
  dewix = find(strncmpi(hdr, 'DEW', 3), 1);
  tidix = find(strncmpi(hdr, 'TIDE', 4), 1);

  swhix = find(strncmpi(hdr, 'WVHT', 4), 1);
  swpix = find(strncmpi(hdr, 'DPD', 3), 1);
  apdix = find(strncmpi(hdr, 'APD', 3), 1);
  mwdix = find(strncmpi(hdr, 'MWD', 3), 1);

  % (Fields generally found in NDBC "o" Ocean data files)
  ocndix = find(strcmpi(hdr, 'DEPTH'), 1);
  ocntix = find(strcmpi(hdr, 'OTMP'), 1);
  ocncix = find(strcmpi(hdr, 'COND'), 1);
  ocnsix = find(strcmpi(hdr, 'SAL'), 1);

  % (Fields generally found in NDBC "r" Solar Radiation data files)
  sradrix = find(strncmpi(hdr, 'SRAD', 4), 1);
  sradsix = find(strncmpi(hdr, 'SWRAD', 5), 1);
  sradlix = find(strncmpi(hdr, 'LWRAD', 5), 1);


  stn = update_ndbc_RT_field(stn, 'ndbc_wind1_dir', dts, dat, dirix, true);
  stn = update_ndbc_RT_field(stn, 'ndbc_wind1_speed', dts, dat, spdix);
  stn = update_ndbc_RT_field(stn, 'ndbc_wind1_gust', dts, dat, gstix);
  stn = update_ndbc_RT_field(stn, 'ndbc_barom', dts, dat, baromix, false, true);
  stn = update_ndbc_RT_field(stn, 'ndbc_air_t', dts, dat, airtix);
  stn = update_ndbc_RT_field(stn, 'ndbc_sea_t', dts, dat, seatix);
  stn = update_ndbc_RT_field(stn, 'ndbc_dew_t', dts, dat, dewix);
  stn = update_ndbc_RT_field(stn, 'ndbc_tide', dts, dat, tidix);

  stn = update_ndbc_RT_field(stn, 'ndbc_sigwavehgt', dts, dat, swhix);
  stn = update_ndbc_RT_field(stn, 'ndbc_sigwaveper', dts, dat, swpix);
  stn = update_ndbc_RT_field(stn, 'ndbc_avgwaveper', dts, dat, apdix);
  stn = update_ndbc_RT_field(stn, 'ndbc_avgwavedir', dts, dat, mwdix, true);

  stn = update_ndbc_RT_field(stn, 'ndbc_ct_shallow_seatemp', dts, dat, ocntix);
  stn = update_ndbc_RT_field(stn, 'ndbc_ct_shallow_cond', dts, dat, ocncix);
  stn = update_ndbc_RT_field(stn, 'ndbc_ct_shallow_salin', dts, dat, ocnsix);
  stn = update_ndbc_RT_field(stn, 'ndbc_ct_shallow_i_depth', dts, dat, ocndix);

  stn = update_ndbc_RT_field(stn, 'ndbc_drf',  dts, dat, sradrix,true,true);
  stn = update_ndbc_RT_field(stn, 'ndbc_dsrf', dts, dat, sradsix,true,true);
  stn = update_ndbc_RT_field(stn, 'ndbc_dlrf', dts, dat, sradlix,true,true);

return;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% INTERNAL FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function stn = update_ndbc_RT_field(stn, fld, dts, dat, fldix, valid99, valid999)

  % Is 99 a valid value for this field?
  if ( ~exist('valid99','var') || isempty(valid99) )
    valid99 = false;
  end;
  % Is 999 a valid value for this field?
  if ( ~exist('valid999','var') || isempty(valid999) )
    valid999 = false;
  end;
  if ( ~isempty(fldix) )

    if ( size(dat,2) < fldix )
      if ( ~isempty(regexp(fld,'_adcp_')) )
      else
        warning('Possible partially populated column(s) for "%s"!', fld);
      end;
      return;
    end;

    dat = dat(:,fldix);

    % Idiotic IMPORTDATA, idiotic NDBC: If a field appears in the header, but
    % only begins to have non-blank data in it midway through the file, then
    % IMPORTDATA will parse that extra column as a line all by itself, with all
    % other columns NaN-filled. This tries to fix that, by ignoring such lines.
    % (And in the process, it also ignores any data from that column... Sorry!)
    dat(isnan(dts)) = [];
    dts(isnan(dts)) = [];

    % Some fields only available intermittently over the years, e.g., wave
    % height or dew-point temperature. Wasteful to reuse time series dates
    % from the more consistent variables, for these intermittent ones.
    if ( ~valid99 )
      dat(dat == 99.0) = nan;
    end;
    if ( ~valid999 )
      dat(dat == 999.0) = nan;
    end;
    dat(dat == 9999.0) = nan;
    dts(isnan(dat)) = [];
    dat(isnan(dat)) = [];

    % Ensure dates are monotonically increasing (and unique)
    [dts, dtsi] = unique(dts);
    dat = dat(dtsi);

    npts = numel(dts);

    % Only create/append field, if we have SOME data
    if ( npts > 1 )
      if ( ~isfield(stn, fld) )
        stn.(fld).date = dts(:);
        stn.(fld).data = dat(:);
      else
        % DO NOT overwrite values for duplicated dates, i.e., *keep* whole
        % DATE BLOCKS with any QC'd data that was passed in with STN
        dat(dts <= stn.(fld).date(end)) = [];
        dts(dts <= stn.(fld).date(end)) = [];
        npts = numel(dts);

        % Basic sanity check
        if ( npts == 0 )
          warning('No new valid data in "%s"!',fld);
        else
          dgap = dts(1) - stn.(fld).date(end);
          if ( dgap > 3 )
            warning('Apparent %g-day gap in "%s"!',dgap,fld);
          end;
          
          stn.(fld).date(end+1:end+npts,1) = dts;
          stn.(fld).data(end+1:end+npts,1) = dat;
        end; %if ( npts == 0 ) else
      end; %if ( ~isfield(stn, fld) ) else
    end; %if ( npts > 1 )

  end; %if ~isempty(fldix)

return;
