function stn = qa_ts(stn,varname,sds,doReview,workFromMAT,matfname)
%function stn = qa_ts(stn,varname,sds,doReview,workFromMAT,matfname)
% 
% Apply various quality-assurance procedures to time series STN.(VARNAME),
% returning a new pair of time series STN.([VARNAME '_qc']) containing only
% quality-controlled data, and STN.([VARNAME '_qc_flags']) containing a
% 4-byte bit-vector of QC flags. Flags==0 implies "best possible data":
% 
%     Bit 01 = outside known allowable range for this kind of data
%     Bit 02 = outside SDS standard deviations for this kind of data
%              (if SDS not given, DEFAULT is 5 standard deviations)
%        NOTE: PRCTILE is used to check ranges, with percentile range
%         set based on SDS st.devs. for a Gaussian distribution (DEFAULT
%         99.9%). Thus we do NOT assume a 'true' distribution of t.s.!
%     Bit 03 = outside seasonal range for this kind of data (if known)
%     Bit 04 = outside (SDS - 1) standard deviations from DESEASONALIZED data
%     Bit 05 = outside daily variance from previous 24 hours' data
%     Bit 06 = outside daily SEASONAL variance from previous 24 hours
% 
%     Bit 16 = marked as BAD by human reviewer
%     Bit 17 = marked as GOOD (despite above flags) by human reviewer
% 
%     Bit 32 = status unknown (not good)
% 
% Result STN.([VARNAME '_qc']) will only include data with flags 0 or 17
% 
% NOTE: Human review is performed (with a call to QA_TS_REVIEW) if DOREVIEW
% is given and true. If true and QA data already exists, user prompted to
% edit or replace it. If DOREVIEW, *and* WORKFROMMAT also present and True,
% and either MATFNAME specifies a useable pathname string or STN.station_name
% is present, then save the resulting QA'd time series in a .MAT file. If the
% named MAT file already exists, LOAD IT before calling QA_TS_REVIEW.
%
% Last Saved Time-stamp: <Fri 2011-12-30 21:13:02  lew.gramer>

  set_more off;

  datapath = get_ecoforecasts_path('data');

  if ( ~isstruct(stn) || ~isfield(stn, varname) )
    error('No field "%s" (2nd arg) in struct STN (1st arg)', varname);
  end;
  if ( ~isfield(stn.(varname), 'date') || ~isfield(stn.(varname), 'data') )
    error('"STN.%s" not a time series? (no date or no data)', varname);
  end;
  if ( ~exist('sds', 'var') || isempty(sds) )
    sds = 5;
  end;
  if ( ~exist('doReview', 'var') || isempty(doReview) )
    doReview = false;
  end;
  if ( ~exist('workFromMAT', 'var') || isempty(workFromMAT) )
    workFromMAT = false;
  end;
  if ( doReview && workFromMAT )
    if ( ~exist('matfname', 'var') || isempty(matfname) || ~ischar(matfname) )
      if ( ~isfield(stn,'station_name') )
        error('WORKFROMMAT when neither STN.station_name nor valid MATFNAME given');
      else
        matfname = fullfile(datapath,[lower(stn.station_name),'_',varname,'.mat']);
      end;
    end;
  end;


  flgname = [varname '_qc_flags'];
  qcname = [varname '_qc'];

  if ( doReview && workFromMAT )
    if ( ~exist(matfname,'file') )
      disp(['Results will be saved to ',matfname]);
    else
      disp(['Results will be loaded from/saved to ',matfname]);
      load(matfname,'station');
      stn.(flgname) = station.(flgname);
      stn.(qcname) = station.(qcname);
      station=[]; clear station
    end;
  end;

  dts = stn.(varname).date;
  dat = stn.(varname).data;

  %
  % Build bit-vector vector, with QC result flags for each value in 'varname'
  %

  % Roger, Roger. Over, Over. What's our vector, Victor?

  flg = zeros([1 1], 'uint32');
  if ( isfield(stn, flgname) )

    % Did user say DOREVIEW=true on a TS that has already been QA'd?
    if ( doReview )
      quest = sprintf('Edit QA data *already present* for "%s"?', varname);
      answr = questdlg(quest,[mfilename ' query']);
      if ( strcmp(answr,'Yes') )
        % Perform (optional) human graphical review of "good" vs "bad" data
        % (Also clears any variables derived from this QC variable...)
        [stn,changed] = qa_ts_review(stn, varname);
        if ( changed && workFromMAT )
          qa_ts_save_to_mat(stn,varname,flgname,qcname,matfname);
        end;
        set_more;
        return;
      elseif ( strcmp(answr,'Cancel') )
        warning('Exiting QA_TS("%s") at user request...', varname);
        set_more;
        return;
      end;
    end;

    stn = rmfield(stn, flgname);
  end;
  stn.(flgname).date = dts;
  stn.(flgname).data = repmat(flg, size(dat));


  nonnan = dat(~isnan(dat));

  mn = min(nonnan);
  av = mean(nonnan);
  mx = max(nonnan);

  % sd = std(nonnan) * sds;
  % lowersd = av - sd;
  % uppersd = av + sd;
  switch (sds),
   case 1, sdspct = 33;     climsdspct = 33;
   case 2, sdspct = 17;     climsdspct = 33;
   case 3, sdspct =  3;     climsdspct = 17;
   case 4, sdspct =  1;     climsdspct =  3;
   case 5, sdspct =  0.1;   climsdspct =  1;
   case 6, sdspct =  0.01;  climsdspct =  0.1;
   case 7, sdspct =  0.001; climsdspct =  0.01;
   otherwise,
    error('Variable %s, "sds" must be an INTEGER between 1 and 7 incl.!', varname);
  end;
  lowersd = prctile(nonnan, sdspct);
  uppersd = prctile(nonnan, (100-sdspct));

  rg = valid_var_range(varname);

  lo = max(max(lowersd, mn), rg(1));
  hi = min(min(uppersd, mx), rg(2));



  % Absolute range and variance checks
  outrng = find( rg(1) > dat | dat > rg(2) );
  %     Bit 01 = outside known allowable range for this kind of data
  stn.(flgname).data(outrng) = bitset(stn.(flgname).data(outrng), 01);

  extrem = find( lo > dat | dat > hi );
  %     Bit 02 = outside 'sds' standard deviations for this kind of data
  stn.(flgname).data(extrem) = bitset(stn.(flgname).data(extrem), 02);


  % Seasonal range and variance checks
  [clim, climsd, climlopct, climuppct] = ...
      climatologize_time_series(dts, dat, varname, climsdspct);

  [Y M D] = datevec(dts);
  jdays = fix(dts) - datenum(Y, 1, 1) + 1;
  % Darn leap years
  jdays(jdays == 366) = 365;


  %     Bit 03 = outside seasonal range for this kind of data (if known)
  % ...


  for jday = 1:365
    nowix = find(jdays == jday);
    % seasextrem = find( abs(dat(nowix) - clim(jday)) > (sds*climsd(jday)) );
    % Use percentile ranges in place of SD (which would assume Gaussian distrib.)
    seasextrem = find( (climlopct(jday) > dat(nowix)) | (dat(nowix) > climuppct(jday)) );
    %     Bit 04 = outside 'sds' standard deviations from DESEASONALIZED data
    stn.(flgname).data(nowix(seasextrem)) = bitset(stn.(flgname).data(nowix(seasextrem)), 04);
  end;


  % Daily variance checks
  [dailysddts, dailysd] = window_func(dts, dat, 'std', 24);
  dailysdmin = prctile(dailysd, 1);
  dailysdmax = prctile(dailysd, 99);
  baddts = dailysddts(dailysdmin > dailysd | dailysd > dailysdmax);
  flgdts = find(ismember(floor(dts), floor(baddts)));
  %     Bit 05 = outside daily variance from previous 24 hours' data
  stn.(flgname).data(flgdts) = bitset(stn.(flgname).data(flgdts), 05);


  % Absolute jumps (ds / dt > unphysical change bound)
  dtsdiff = abs(diff(dts));
  datdiff = abs(diff(dat));
  % Calculate a daily trend between points
  dsdt = datdiff ./ dtsdiff;
  dailysdextreme = prctile(dailysd, 99.9);
  jmpidx = find(dsdt > dailysdextreme) + 1;
  %     Bit 06 = outside daily SEASONAL variance from previous 24 hours
  stn.(flgname).data(jmpidx) = bitset(stn.(flgname).data(jmpidx), 06);


  %
  % Pass through QC'd data in a new time series stn.([varname '_qc'])
  %

  if ( isfield(stn, qcname) )
    stn = rmfield(stn, qcname);
  end;
  stn.(qcname).date = stn.(flgname).date;
  stn.(qcname).data = repmat(nan, size(dat));

  gooddata = find(stn.(flgname).data == 0);
  stn.(qcname).data(gooddata) = dat(gooddata);

  % No sense returning all the NaNs in our "good" vector
  stn.(qcname).date = stn.(qcname).date(~isnan(stn.(qcname).data));
  stn.(qcname).data = stn.(qcname).data(~isnan(stn.(qcname).data));


  % Perform (optional) human graphical review of "good" vs "bad" data
  if ( doReview )
    % (Also clears any variables derived from this QC variable...)
    [stn,changed] = qa_ts_review(stn, varname);
    if ( changed && workFromMAT )
      qa_ts_save_to_mat(stn,varname,flgname,qcname,matfname);
    end;

  else
    % Clear any variables DERIVED from this QC variable before we finish!
    flds = fieldnames(stn);
    for fidx = 1:length(flds)
      fld = flds{fidx};
      if ( ~strcmpi(fld, qcname) && ~strcmpi(fld, flgname) )
        if ( ~isempty(strfind(fld, qcname)) )
          saved_ts = stn.(fld);
          stn = rmfield(stn, fld);
          stn = verify_variable(stn, fld);
          if ( isfield(stn,fld) )
            fprintf('Recalculated %s!\n', fld);
          else
            stn.(fld) = saved_ts;
          end;
        end;
      end;
    end;

  end;

  set_more;

return;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INTERNAL FUNCTIONS 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function qa_ts_save_to_mat(stn,varname,flgname,qcname,matfname)
  station.(flgname) = stn.(flgname);
  station.(qcname) = stn.(qcname);
  save(matfname,'station');
  disp(['Saved to ',matfname]);
  station=[]; clear station
return;

