function download_ndbc_files(stnm,yrs,doRT,doOverwrite,datapath,thisyr)
%function download_ndbc_files(stnm,yrs,doRT,doOverwrite,datapath,thisyr)
%
% Download quality-controlled data files from http://www.ndbc.noaa.gov for
% station STNM and years YRS (DEFAULT: 1987 to now). If YRS includes this
% year, also download monthly QC'd files for the current year. NOTE: This
% will overwrite older versions of files, but does not combine files (e.g.,
% won't create a "mlrf1h2013.txt" from monthly files "mlrf1[12]2013.txt").
% If DORT (DEFAULT: False), also download the most recent data, specifically
% the real-time file for the past 45 days. If DOOVERWRITE (Default: True),
% download files whether or not they already exist on the data path. NOTE:
% Default DATAPATH is the result of calling GET_ECOFORECASTS_PATH('data').
%
% Last Saved Time-stamp: <Mon 2019-02-25 16:40:25 Eastern Standard Time gramer>

  set_more off;
  w = warning('OFF','MATLAB:DELETE:FileNotFound');

  if ( ~exist('doRT','var') || isempty(doRT) )
    doRT = false;
  end;
  if ( ~exist('doOverwrite','var') || isempty(doOverwrite) )
    doOverwrite = true;
  end;
  if ( ~exist('datapath','var') || strcmpi(datapath,'default') )
    datapath = get_ecoforecasts_path('data');
  end;

  % NDBC QC and annual archiving ultimately takes about 45 d
  if ( ~exist('thisyr','var') || isempty(thisyr) )
    thisyr = get_year(now-45);
  end;
  doThisYear = false;
  if ( ~exist('yrs','var') || isempty(yrs) )
    yrs = 1987:(thisyr-1);
    doThisYear = true;
  else
    ix = find(yrs >= thisyr);
    if ( ~isempty(ix) )
      %yrs(ix) = [];
      doThisYear = true;
    end;
    % NDBC cannot gaze into the future
    yrs(yrs > get_year(now)) = [];
  end;

  %% WAS:
  %http://www.ndbc.noaa.gov/view_text_file.php?filename=41009h1988.txt.gz&dir=data/historical/stdmet/
  %http://www.ndbc.noaa.gov/view_text_file.php?filename=41009o2008.txt.gz&dir=data/historical/ocean/

  % NOW:
  %<p>The document has moved <a href="https://www.ndbc.noaa.gov/view_text_file.php?filename=pgbp7h2008.txt.gz&amp;dir=data/historical/stdmet/">here</a>.</p>

  for yr=yrs(:)'
    %DEBUG:
    disp(yr);
    for ltr = 'hoar';
      fname = sprintf('%s%s%04d.txt',lower(stnm),ltr,yr);
      %DEBUG:    disp(fname);
      fpath = fullfile(datapath,fname);
      if ( exist(fpath,'file') && ~doOverwrite )
        disp(['Keeping existing ',fpath]);
      else
        %DEBUG:      disp(['Downloading ',fname]);
        switch ( ltr ),
         case 'h',
          % Historical meteorology (often incl. water temperature and tide)
          url = sprintf('https://www.ndbc.noaa.gov/view_text_file.php?filename=%s.gz&dir=data/historical/stdmet/',fname);
         case 'o',
          % Ocean Data - sea temperature, salinity, etc.
          url = sprintf('https://www.ndbc.noaa.gov/view_text_file.php?filename=%s.gz&dir=data/historical/ocean/',fname);
         case 'a',
          % Ocean Currents (ADCP) data - historical format
          url = sprintf('https://www.ndbc.noaa.gov/view_text_file.php?filename=%s.gz&dir=data/historical/adcp/',fname);
         case 'r',
          % Solar Radiation data (from, e.g., St. Augustine Buoy 41012)
          url = sprintf('https://www.ndbc.noaa.gov/view_text_file.php?filename=%s.gz&dir=data/historical/srad/',fname);
         otherwise,
          error('Do not yet know how to download files of type "%s"',ltr);
        end;
        if ( ltr == 'h' )
          disp(fname);
        end;
        try, F=websave(fpath,url); catch ME, try, delete(fpath); catch ME, ...
               end; try, delete([fpath,'.html']); catch ME, end; F=[]; end;
        if ( isempty(F) )
          if ( ltr == 'h' )
            warning('Failed download! %s',url);
          end;
        end;
      end;
    end;
  end;

  if ( doThisYear )
    %yr = get_year(now);
    yr = thisyr;
    disp(yr);

    mos = '123456789abc';
    % NDBC cannot gaze into the future
    if ( thisyr == get_year(now) )
      mos(get_month(now)+1:end) = [];
    end;

    need_mo = [];
    need_moix = [];
    for moix = 1:numel(mos);
      mo = mos(moix);
      monm = datestr(datenum(yr,moix,1),'mmm');
      disp(monm);

      % QC'd monthly meteorology (often incl. water temperature and tide)
      fname = sprintf('%s%s%04d.txt',lower(stnm),mo,yr);
      %DEBUG:    disp(fname);
      fpath = fullfile(datapath,fname);
      if ( exist(fpath,'file') && ~doOverwrite )
        disp(['Keeping existing ',fpath]);
      else
        url = sprintf('https://www.ndbc.noaa.gov/view_text_file.php?filename=%s.gz&dir=data/stdmet/%s/',fname,monm);
        disp(fname);
        try, F=websave(fpath,url); catch ME, try, delete(fpath); catch ME, ...
               end; try, delete([fpath,'.html']); catch ME, end; F=[]; end;
        if ( isempty(F) )
          warning('Failed download! %s',url);
          need_mo(end+1) = mo;
          need_moix(end+1) = moix;
        end;
      end;

      % QC'd monthly ocean data (often incl. salinity, tide, DO, etc.)
      o_fname = sprintf('%s%s%04d_o.txt',lower(stnm),mo,yr);
      %DEBUG:    disp(fname);
      o_fpath = fullfile(datapath,o_fname);
      if ( exist(o_fpath,'file') && ~doOverwrite )
        disp(['Keeping existing ',o_fpath]);
      else
        % Historical meteorology (often incl. water temperature and tide)
        url = sprintf('https://www.ndbc.noaa.gov/view_text_file.php?filename=%s.gz&dir=data/ocean/%s/',fname,monm);
        %disp(o_fname);
        try, F=websave(o_fpath,url); catch ME, try, delete(o_fpath); catch ME, ...
               end; try, delete([o_fpath,'.html']); catch ME, end; F=[]; end;
        % if ( isempty(F) )
        %   warning('Failed download! %s',url);
        % end;
      end;
    end; %for moix = 1:numel(mos);

    % There is usually one last file whose name is not in the format above!
    for ix = 1:numel(need_mo)
      mo = need_mo(ix);
      moix = need_moix(ix);
      % E.g., https://www.ndbc.noaa.gov/data/stdmet/Dec/lonf1.txt
      monm = datestr(datenum(yr,moix,1),'mmm');
      fname = sprintf('%s%s%04d.txt',lower(stnm),mo,yr);
      fpath = fullfile(datapath,fname);
      if ( exist(fpath,'file') && ~doOverwrite )
        disp(['Keeping existing ',fpath]);
      else
        disp(fname);
        url = sprintf('https://www.ndbc.noaa.gov/data/stdmet/%s/%s.txt',monm,lower(stnm));
        try, F=websave(fpath,url); catch ME, try, delete(fpath); catch ME, ...
               end; try, delete([fpath,'.html']); catch ME, end; F=[]; end;
        if ( isempty(F) )
          warning('Failed download! %s',url);
        end;
      end;
    end;

  end;

  if ( doRT )
    %https://www.ndbc.noaa.gov/data/realtime2/LONF1.txt
    %https://www.ndbc.noaa.gov/data/realtime2/LONF1.cwind
    %https://www.ndbc.noaa.gov/data/derived2/LONF1.dmv

    % Only do Meteorology files for now (often incl. water temperature and tide)
    fname = sprintf('%sRT.txt',lower(stnm));
    %DEBUG:    disp(fname);
    fpath = fullfile(datapath,fname);
    if ( exist(fpath,'file') && ~doOverwrite )
      disp(['Keeping existing ',fpath]);
    else
      disp(fname);
      url = sprintf('https://www.ndbc.noaa.gov/data/realtime2/%s.txt',upper(stnm));
      try, F=websave(fpath,url); catch ME, try, delete(fpath); catch ME, ...
             end; try, delete([fpath,'.html']); catch ME, end; F=[]; end;
      if ( isempty(F) )
        warning('Failed download! %s',url);
      end;
    end;
  end;

  warning(w);
  set_more;

return;
