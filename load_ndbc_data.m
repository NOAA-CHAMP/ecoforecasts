function stn = load_ndbc_data(stn, fname)
%function stn = load_ndbc_data(stn, fname)
%
% Load the file named FNAME, and append its data to the station data struct
% STN (optionally an empty matrix [], if we're just starting to load data).
% These archived NDBC C-MAN/SEAKEYS/ICON station data files can be obtained
% for a given the 5-character code 'stnam' from http://www.ndbc.noaa.gov.
%
% NOTE: 99 and 999 are used as both valid values (e.g., wind dir and barom,
% resp.) but also as "bad value" placeholders: this code is careful of that.
%
% FUTURE IMPROVEMENT NEED: Automatically download "STNAM[ho]YYYY.txt" files!
%
% Last Saved Time-stamp: <Thu 2013-12-26 17:46:08 Eastern Standard Time gramer>


  if ( ~exist(fname, 'file') )
    warning('Skipping non-existent year file "%s"...', fname);
    return;
  end;

  x = importdata(fname);
  if ( ~isfield(x, 'textdata') )
    error('Unable to find header in "%s"!', fname);
  end;
  if ( ~isfield(x, 'data') )
    error('Unable to find data in "%s"!', fname);
  end;

  % NDBC headers have been remarkably consistent for >20 years! From...
  % YY MM DD hh WD   WSPD GST  WVHT  DPD   APD  MWD  BAR    ATMP  WTMP  DEWP  VIS
  %  to ...
  % #YY  MM DD hh mm WDIR WSPD GST  WVHT   DPD   APD MWD   PRES  ATMP  WTMP  DEWP  VIS  TIDE
  % #yr  mo dy hr mn degT m/s  m/s     m   sec   sec degT   hPa  degC  degC  degC   mi    ft

  % But sometimes, flaky IMPORTDATA() fails to parse header lines
  if ( size(x.textdata) == [1 1] )
    hdr = textscan(x.textdata{:}, '%s');
    hdr = hdr{:};
  elseif ( size(x.textdata,1) > 1 )
    hdr = textscan(x.textdata{1,1}, '%s');
    hdr = hdr{:};
  else
    hdr = x.textdata;
  end;

  dat = x.data;


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

  % (Fields generally found in NDBC "a" ADCP data files)
  adcphix = find(strncmpi(hdr, 'DEP', 3));
  adcpdix = find(strncmpi(hdr, 'DIR', 3));
  adcpsix = find(strncmpi(hdr, 'SPD', 3));

  % (Fields generally found in NDBC "r" Solar Radiation data files)
  sradrix = find(strncmpi(hdr, 'SRAD', 4), 1);
  sradsix = find(strncmpi(hdr, 'SWRAD', 5), 1);
  sradlix = find(strncmpi(hdr, 'LWRAD', 5), 1);


  if ( ~isempty(regexp(fname,'a[12][0-9][0-9][0-9][.]txt')) )
    if (isempty(adcphix) || numel(adcphix)~=numel(adcpsix) || numel(adcpdix)~=numel(adcpsix))
      warning('Skipping misformated ADCP file "%s"...', fname);
      return;
    else
      %NOTE: This code does NOT handle it well when the number of ADCP bins
      %      changes during the course of a one-year file! Work-around: break
      %      "a" files up so that each "year" has a consistent number of bins. 
      for ix=1:numel(adcphix)
        ixs = num2str(ix,'%02d');
        if ( size(dat,2) >= max([adcpdix(ix),adcpsix(ix),adcphix(ix)]) )
          % Bad speed or direction also means bad depth
          badix = find(dat(:,adcpdix(ix))==999 | dat(:,adcpsix(ix))==999 | dat(:,adcphix(ix))==999);
          dat(badix,adcpdix(ix)) = 999;
          dat(badix,adcpsix(ix)) = 999;
          dat(badix,adcphix(ix)) = 999;
        end;
        stn = update_ndbc_field(stn, ['ndbc_adcp_dir_',ixs],dts,dat,adcpdix(ix),true);
        stn = update_ndbc_field(stn, ['ndbc_adcp_spd_',ixs],dts,dat,adcpsix(ix),true);
        stn = update_ndbc_field(stn, ['ndbc_adcp_dep_',ixs],dts,dat,adcphix(ix),true);
      end;
    end;
  else
    stn = update_ndbc_field(stn, 'ndbc_wind1_dir', dts, dat, dirix, true);
    stn = update_ndbc_field(stn, 'ndbc_wind1_speed', dts, dat, spdix);
    stn = update_ndbc_field(stn, 'ndbc_wind1_gust', dts, dat, gstix);
    stn = update_ndbc_field(stn, 'ndbc_barom', dts, dat, baromix, false, true);
    stn = update_ndbc_field(stn, 'ndbc_air_t', dts, dat, airtix);
    stn = update_ndbc_field(stn, 'ndbc_sea_t', dts, dat, seatix);
    stn = update_ndbc_field(stn, 'ndbc_dew_t', dts, dat, dewix);
    stn = update_ndbc_field(stn, 'ndbc_tide', dts, dat, tidix);

    stn = update_ndbc_field(stn, 'ndbc_sigwavehgt', dts, dat, swhix);
    stn = update_ndbc_field(stn, 'ndbc_sigwaveper', dts, dat, swpix);
    stn = update_ndbc_field(stn, 'ndbc_avgwaveper', dts, dat, apdix);
    stn = update_ndbc_field(stn, 'ndbc_avgwavedir', dts, dat, mwdix, true);

    stn = update_ndbc_field(stn, 'ndbc_ct_shallow_seatemp', dts, dat, ocntix);
    stn = update_ndbc_field(stn, 'ndbc_ct_shallow_cond', dts, dat, ocncix);
    stn = update_ndbc_field(stn, 'ndbc_ct_shallow_salin', dts, dat, ocnsix);
    stn = update_ndbc_field(stn, 'ndbc_ct_shallow_i_depth', dts, dat, ocndix);

    stn = update_ndbc_field(stn, 'ndbc_drf',  dts, dat, sradrix,true,true);
    stn = update_ndbc_field(stn, 'ndbc_dsrf', dts, dat, sradsix,true,true);
    stn = update_ndbc_field(stn, 'ndbc_dlrf', dts, dat, sradlix,true,true);
  end;

return;


function stn = update_ndbc_field(stn, fld, dts, dat, fldix, valid99, valid999)

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
        stn.(fld).date = [];
        stn.(fld).data = [];
      end;

      stn.(fld).date(end+1:end+npts,1) = dts;
      % Overwrite values for duplicated dates
      stn.(fld).data(end+1:end+npts,1) = dat;
    end;

  end; %if ~isempty(fldix)

return;
