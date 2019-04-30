function stn = load_nodc_data(stn, fname)
%function stn = load_nodc_data(stn, fname)
%
% Load the file named FNAME, and append its data to the station data struct
% STN (optionally an empty matrix [], if we're just starting to load data).
% These archived CoRIS/NODC C-MAN/SEAKEYS/ICON station data files can be
% found at http://data.nodc.noaa.gov/coris/data/NOAA/oar/SEAKEYSstations
% and http://data.nodc.noaa.gov/coris/data/NOAA/oar/ICONstations.
%
% Last Saved Time-stamp: <Fri 2010-03-12 15:21:58 Eastern Standard Time gramer>


  if ( ~exist(fname, 'file') )
    warning('Skipping non-existent file "%s"...', fname);
    return;
  end;

  x = importdata(fname);
  if ( ~isfield(x, 'textdata') )
    error('Unable to find header in "%s"!', fname);
  end;
  if ( ~isfield(x, 'data') )
    error('Unable to find data in "%s"!', fname);
  end;

  % CHAMP-QA'd file headers are somewhat consistent across time and space

  % MLRF1: FROM 1992...
  % Station Name	Year	Julian Day	Hour	AirTemp	Dew Point	Wind Dir	Max Wind Speed	Wind Gust	Wind Speed	Barometer	Sea Temp NDBC
  % TO 2008...
  % Date	Name	Jyear	Jday	Hour	WDir1	WS1	WDir2	WS2	WGDir1	WG1	WGDir2	WG2	AirT	DewPt	Baro1	SeaTempNDBC	PARS	SEaTempFIO1m	conductivity	salinity

  % SMKF1: FROM 2003...
  % Date	Name	Jyear	Jday	Hour	WDir1	WS1	WDir2	WS2	WGDir1	WG1	WGDir2	WG2	AirT	DewPt	Baro1	SeaTempNDBC	PARS	SEaTempFIO1m	conductivity	salinity
  % TO 2005...
  % Date	Name	Jyear	Jday	Hour	WDir1	WS1	WDir2	WS2	WGDir1	WG1	WGDir2	WG2	AirT	DewPt	Baro	SeaTempNDBC	SeaTempFIO1m	conductivity	salinity


  % But sometimes, flaky IMPORTDATA() fails to parse header lines
  if ( size(x.textdata) == [1 1] )
    hdr = textscan(x.textdata{:}, '%s');
    hdr = hdr{:};
  elseif ( size(x.textdata,1) > 1 )
    hdr = x.textdata(1,:);
  else
    hdr = x.textdata;
  end;

  for ix = 1:length(x.textdata(2,:))
    if ( ~isempty(x.textdata{2,ix}) )
      skipToIx = ix;
    else
      break;
    end;
  end;

  dat = x.data;


  % Parse timestamp columns

  DTix = find((strcmpi(hdr, 'DATE')), 1);
  YRix = find((strcmpi(hdr, 'YEAR') | strcmpi(hdr, 'JYEAR')), 1);
  JDix = find((strcmpi(hdr, 'JDAY') | strcmpi(hdr, 'JULIAN DAY')), 1);
  MOix = find((strcmpi(hdr, 'MONTH') | strcmpi(hdr, 'MO')), 1);
  DYix = find((strcmpi(hdr, 'DAY')), 1);
  HRix = find((strcmpi(hdr, 'HOUR') | strcmpi(hdr, 'HR') | strcmpi(hdr, 'HH')), 1);
  MNix = find((strcmpi(hdr, 'MINUTE') | strcmpi(hdr, 'MIN') | strcmpi(hdr, 'MN')), 1);

  YRix = YRix - skipToIx;
  JDix = JDix - skipToIx;
  MOix = MOix - skipToIx;
  DYix = DYix - skipToIx;
  HRix = HRix - skipToIx;
  MNix = MNix - skipToIx;

  dts = [];
  if ( ~isempty(HRix) )
    hrs = dat(:,HRix);
    if ( ~isempty(MNix) )
      mns = dat(:,MNix);
    else
      mns = 0;
    end;
    if ( ~isempty(DTix) )
      dts = datenum(x.textdata(2:end,DTix)) + (hrs/24) + (mns/24/60);
    elseif ( ~isempty(YRix) && ~isempty(MOix) && ~isempty(DYix) )
      dts = datenum(dat(:,YRix),dat(:,MOix),dat(:,DYix),hrs,mns);
    elseif ( ~isempty(YRix) && ~isempty(JDix) )
      dts = datenum(dat(:,YRix),1,1,hrs,mns) + dat(:,JDix) - 1;
    end;
  end;
  if ( isempty(dts) )
    error('Needed timestamp header field(s) not found in "%s"!', fname);
  end;

  % Parse data columns (of interest)

  dirix = find((strncmpi(hdr,'WDIR',4) | strncmpi(hdr,'WIND DIR',8)), 1);
  spdix = find((strncmpi(hdr,'WS',2) | strncmpi(hdr,'WSPD',4) | strncmpi(hdr,'WIND SPEED',10)), 1);
  gstix = find((strncmpi(hdr,'WG',2) | strncmpi(hdr,'WGST',4) | strncmpi(hdr,'WIND GUST',9)), 1);
  baromix = find((strncmpi(hdr,'PRES',4) | strncmpi(hdr,'BAR',3)), 1);
  airtix = find((strncmpi(hdr,'ATMP',4) | strncmpi(hdr,'AIR',3)), 1);
  dewix = find((strncmpi(hdr,'DEW',3) | strncmpi(hdr,'DWPT',4)), 1);
  seatix = find((strncmpi(hdr, 'SEA',3) | strncmpi(hdr, 'WT', 2) | strncmpi(hdr, 'WATER', 5)), 1);
  salix = find(strncmpi(hdr, 'SAL', 3), 1);
  parix = find(strncmpi(hdr, 'PAR', 3), 1);

  stn = update_nodc_field(stn, 'nodc_wind1_dir', dts, dat, dirix, skipToIx);
  stn = update_nodc_field(stn, 'nodc_wind1_speed', dts, dat, spdix, skipToIx);
  stn = update_nodc_field(stn, 'nodc_wind1_gust', dts, dat, gstix, skipToIx);
  stn = update_nodc_field(stn, 'nodc_barom', dts, dat, baromix, skipToIx);
  stn = update_nodc_field(stn, 'nodc_air_t', dts, dat, airtix, skipToIx);
  stn = update_nodc_field(stn, 'nodc_dew_t', dts, dat, dewix, skipToIx);
  stn = update_nodc_field(stn, 'nodc_sea_t', dts, dat, seatix, skipToIx);
  stn = update_nodc_field(stn, 'nodc_salin', dts, dat, salix, skipToIx);
  stn = update_nodc_field(stn, 'nodc_surf_par', dts, dat, parix, skipToIx);

return;


function stn = update_nodc_field(stn, fld, dts, dat, fldix, skipToIx)

  if ( ~isempty(fldix) )

    fldix = fldix - skipToIx;

    if ( size(dat,2) < fldix )
      warning('Possible partially populated column(s) for "%s"!', fld);
      return;
    end;

    dat = dat(:,fldix);

    % Idiotic IMPORTDATA: If a particular field appears in the header, but
    % only begins to have non-blank data in it midway through the file, then
    % IMPORTDATA will parse that extra column as a line all by itself, with all
    % other columns NaN-filled. This tries to fix that, by ignoring such lines.
    % (And in the process, it also ignores any data from that column... Sorry!)
    dat(isnan(dts)) = [];
    dts(isnan(dts)) = [];

    % Some fields only available intermittently over the years, e.g., wave
    % height or dew-point temperature. Wasteful to reuse time series of dates
    % from the more consistent variables, for these intermittent ones.
    dat(dat == -9.0) = nan;
    % dat(dat == 99.0) = nan;
    dat(dat == 999.0) = nan;
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
