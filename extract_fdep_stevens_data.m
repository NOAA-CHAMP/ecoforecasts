function stn = extract_fdep_stevens_data(stcd,dts)
%function stn = extract_fdep_stevens_data(stcd,dts)
% 
% Run Web query to extract whole month(s) worth of environmental data for a
% Florida Dept. of Environmental Protection (FDEP) Stevens Water Monitoring
% Systems coastal monitoring station with code STCD (DEFAULT: 'k' for their
% "Station K" near St. Lucie Inlet). See: http://www.fldep-stevens.com
% 
% If return value STN is specified, also calls READ_FDEP_STEVENS_DATA (v.)
% 
% Last Saved Time-stamp: <Wed 2016-08-10 17:28:00 Eastern Daylight Time gramer>

  set_more off;

  datpath = get_ecoforecasts_path('data');

  if ( ~exist('stcd','var') || isempty(stcd) )
    stcd = 'k';
  end;
  if ( ~exist('dts','var') || isempty(dts) )
    [stn,dts] = get_fdep_stevens_metadata(stcd);
  else
    [stn,ig] = get_fdep_stevens_metadata(stcd);
    ig=[]; clear ig
  end;

  yrmos = unique(get_yearmonth(dts));

  for begdt = yrmos(:)';
    yr = get_year(begdt);
    mo = get_month(begdt);
    enddt = datenum(yr,mo+1,0,23,59,59);

    fname = fullfile(datpath,sprintf('%s-%04d%02d.csv',stn.station_name,yr,mo));

    % Was there a partial load of this month previously? If so, delete the
    % stub file "reminder" and the old data file: SEE COMMENTS BELOW.
    partialfname = [fname,'.PARTIAL'];
    if ( exist(fname,'file') && exist(partialfname,'file') )
      delete(fname);
      delete(partialfname);

      % A .MAT file is also created by READ_FDEP_STEVENS_DATA (v.)
      w = warning('OFF','MATLAB:DELETE:FileNotFound');
      matfname = fullfile(datpath,sprintf('%s-%04d%02d.mat',stn.station_name,yr,mo));
      delete(matfname);
      warning(w);
    end;

    % Do we need to extract this month from the Web site?
    if ( ~exist(fname,'file') )
      disp(['Extracting ',fname]);

      % If we are extracting the current month, then this file may be partial
      % later: leave a stub file as "reminder" for us to reload it next time.
      if ( begdt == get_yearmonth(now) )
        fid = fopen(partialfname,'w+');
        fclose(fid);
      end;

      % Meshuggahs because Stevens interface doesn't seem to support specifying
      % query dates in UTC: We want an exact month worth in UTC, so for EST/EDT
      % transition months, we will get slightly more or less than 24 hours...???
      begdati = datetime(begdt,'TimeZone','America/New_York','ConvertFrom','datenum');
      begdt = begdt + ( hour(tzoffset(begdati)) / 24);
      enddati = datetime(enddt,'TimeZone','America/New_York','ConvertFrom','datenum');
      enddt = enddt + ( hour(tzoffset(enddati)) / 24);

      [begyr,begmo,begdy,beghr,begmi,begse] = datevec(begdt);
      [endyr,endmo,enddy,endhr,endmi,endse] = datevec(enddt);

      % SAMPLE URL (for station K):
      %  http://www.fldep-stevens.com/export-8722375.php?t=c&txtFDate=02/01/2013&selFTime=00&txtTDate=02/28/2013&selTTime=23&selChannel=&rdoDateTime=1&rdoWTemp=1&rdoWSpeed=2&rdoATemp=1&rdoPressure=1&rdoRainfall=1&rdoWLevel=1

      url = sprintf('http://www.fldep-stevens.com/export-%s.php?t=c&txtFDate=%02d/%02d/%04d&selFTime=%02d:%02d&txtTDate=%02d/%02d/%04d&selTTime=%02d:%02d&selChannel=&rdoDateTime=1&rdoWTemp=1&rdoWSpeed=2&rdoATemp=1&rdoPressure=1&rdoRainfall=1&rdoWLevel=1',stn.station_id,begmo,begdy,begyr,beghr,begmi,endmo,enddy,endyr,endhr,endmi);
      %DEBUG:      disp(url);
      [fname,status] = urlwrite(url,fname,'Timeout',60);
      if ( status ~= 1 || numel(fname) < 3 )
        warning('No data from %s',url);
        continue;
      end; %if ( status ~= 1 || numel(fname) < 3 )

      clear begdati begyr begmo begdy beghr begmi begse enddati endyr endmo enddy endhr endmi endse status url
    end; %if ( ~exist(fname,'file') )

  end; %for begdt = yrmos(:)';

  if ( nargout > 0 )
    stn = read_fdep_stevens_data(stn,dts);
  end; %if ( nargout > 0 )

  set_more;

return;
