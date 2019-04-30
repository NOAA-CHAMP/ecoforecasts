function stn = load_portal_data(stn,fnames)
%function stn = load_portal_data(stn,fnames)
%
% Load data query CSV result file(s) FNAMES downloaded from CHAMP Portal, and
% merge those time-series data with fields in STRUCT STN. If FNAMES is a dir.
% name or the string 'CWD', load all likely-named CSV files from that dir.,
% and save the result as an appropriately named .MAT file.
%
% Last Saved Time-stamp: <Thu 2017-04-27 12:56:54 Eastern Daylight Time lew.gramer>

  if ( ~exist('stn','var') || (~ischar(stn) && ~isfield(stn,'station_name')) )
    error('First arg must be a station name or STRUCT with field .station_name');
  end;
  if ( ischar(stn) )
    % stnm = stn; clear stn
    % stn.station_name = stnm;
    stn = get_station_from_station_name(stn);
  end;
  stnm_prefix = [stn.station_name,'_'];

  fdir = [];
  doALL = false;

  if ( ischar(fnames) )
    % Convenience argument - use Current Working Directory
    if ( strcmpi(fnames,'CWD') || strcmpi(fnames,'PWD') )
      fnames = pwd;
      % We will then fall into the IF immediately below...
    end;
    % If a directory, find all files likely to be loadable in that dir.
    if ( exist(fnames,'dir') )
      doALL = true;
      fdir = fnames;
      fdirpath = fullfile(fdir,[stnm_prefix,'*.csv']);
      x = dir(fdirpath);
      if ( isempty(x) )
        error('No files found! %s',fdirpath);
      end;
      fnames = strcat(fdir,filesep,{x.name});
    elseif ( exist(fnames,'file') )
      fnames = {fnames};
    end;
  end;
  if ( ~iscellstr(fnames) )
    error('Second arg must be ''CWD'', valid pathname, qualified filename, or CELLSTR of qualified filenames');
  end;

  % If caller passed a directory, load all its files and SAVE as a MAT file
  matfname = fullfile(fdir,[stnm_prefix,'portal_ALL.mat']);
  if ( doALL && exist(matfname,'file') )
    disp(['Loading ',matfname]);
    load(matfname);

  else

    res = [];
    for fnix = 1:numel(fnames)
      fname = fnames{fnix};
      disp(fname);
      if ( ~exist(fname,'file') )
        warning('Ecoforecasts:Portal:NoSuchFile',...
                'Skipped! File not found %s',fname);
        continue;
      end;

      x = importdata(fname);
      if ( isfield(x,'textdata') && ~isempty(x.textdata) && isfield(x,'data') && ~isempty(x.data) )
        % For most Portal CSV download files
        dts = datenum(x.textdata(2:end,1));
        hdrs = x.textdata(1,:);
        data = x.data;
      elseif ( iscellstr(x) )
        % Weirdo special case - usually when there is just one data column
        cols = split(x,',');
        dts = datenum(char(cols(2:end,1)));
        hdrs = cols(1,:);
        data = double(cols(2:end,2:end));
      else
        error('Invalid CHAMP Portal file download format?? %s',fname);
      end;
      x=[]; clear x

      if ( ~strcmpi(hdrs{1,1},'datetime') )
        error('Invalid CHAMP Portal file: No timestamp in column one?? %s',fname);
      else
        hdrs = hdrs(:,2:end);
      end;

      if ( numel(hdrs) ~= size(data,2) )
        error('Ecoforecasts:Portal:ColumnMismatch',...
              'Expected %d columns but found %d in %s',numel(hdrs),size(data,2),fname);
      end;
      if ( numel(dts) ~= size(data,1) )
        error('Ecoforecasts:Portal:RowMismatch',...
              'Expected %d rows but found %d in %s',numel(dts),size(data,1),fname);
      end;

      for hdrix = 1:numel(hdrs)
        hdr = hdrs{hdrix};
        if ( ~strncmpi(hdr,stnm_prefix,length(stnm_prefix)) )
          warning('Ecoforecasts:Portal:SuspectField',...
                  'Station name prefix %s not found in field %s',...
                  stnm_prefix,hdr);
        else
          hdr = lower(hdr(length(stnm_prefix)+1:end));
        end;

        fres.(hdr).date = dts;
        fres.(hdr).data = data(:,hdrix);
      end; %for hdrix = 1:numel(hdrs)

      res = merge_station_data(res,fres);
      fres=[]; clear fres;
      dts=[]; hdrs=[]; data=[]; clear dts hdrs data

    end; %for fnix = 1:numel(fnames)


    if ( doALL && ~isempty(res) )
      disp(['Saving ',matfname]);
      save(matfname,'res','-v7.3');
    end;

  end; %if ( doALL && exist(matfname,'file') ) else

  if ( isempty(res) )
    warning('Ecoforecasts:Portal:NoData','No data loaded!');
  else
    stn = merge_station_data(stn,res);
  end;

return;
