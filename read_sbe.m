function stn = read_sbe(fnames,stn)
%function stn = read_sbe(fnames,stn)
%
% Read one or more .ASC files, or a directory of .ASC files, containing SBE37
% or SBE39 ASC-formatted data. This supplants the seemingly FUBAR functions
% named SBE* (SBE37PARSE, etc.) contained in the CSIRO's IMOS-Toolbox.
%
% Last Saved Time-stamp: <Wed 2014-10-29 15:41:13 Eastern Daylight Time gramer>

  % Did user forget the first arg?
  if ( ~exist('fnames','var') )
    fnames = [];
  elseif ( ischar(fnames) )
    % Is it a directory?
    if ( exist(fnames,'dir') )
      fs = dir(fullfile(fnames,'*.asc'));
      %fnames = {fs.name};
      fnames = cellstr(strcat(fnames,filesep,char({fs.name})));
    % Is it a single filename?
    elseif ( exist(fnames,'file') )
      fnames = {fnames};
    % Oops
    else
      fnames = [];
    end;
  % Is it a cell array of filename strings?
  elseif ( ~iscellstr(fnames) )
    fnames = [];
  end;

  if ( ~exist('stn','var') )
    stn = [];
  end;

  % No matter what FNAMES was passed in as, it is now either empty, or a
  % cell array of (one or more) pathnames for data files.
  if ( isempty(fnames) )
      error('First arg must be CELLSTR, or  individual directory or filename string');
  end;


  for fnix = 1:numel(fnames)

    %DEBUG:    disp(fnames{fnix});

    fid = fopen(fnames{fnix},'r');
    if ( fid < 0 )
      warning('Unable to open %s',fnames{fnix});
      %%%% EARLY CONTINUE
      continue;
    end;

    ccstr = textscan(fid,'%[^\n]\n');
    fclose(fid);

    cstr=ccstr{:}; clear ccstr;

    begix = strmatch('start sample number',cstr); begix = begix + 1;
    while ( isempty(cstr{begix}) ); begix = begix + 1; end;

    ncoms = numel(strfind(cstr{begix},','));

    str = sprintf('%s\n',cstr{begix:end});

    tmp=[];
    prs=[];
    cnd=[];
    dts=[];

    % T,C,P,date,time
    if ( ncoms == 4 )
      dat = textscan(str,'%f,%f,%f,%f %[^ ] %f,%f:%f:%f\n');
      tmp = dat{1};
      % Convert S/m to mS/cm
      cnd = dat{2}*10.0;
      prs = dat{3};
      dts = datenum(dat{6},month_str2num(dat{5}),dat{4},dat{7},dat{8},dat{9});
    % T,C,date,time
    elseif ( ncoms == 3 )
      dat = textscan(str,'%f,%f,%f %[^ ] %f,%f:%f:%f\n');
      tmp = dat{1};
      % Convert S/m to mS/cm
      cnd = dat{2}*10.0;
      dts = datenum(dat{5},month_str2num(dat{4}),dat{3},dat{6},dat{7},dat{8});
    % T,date,time
    elseif ( ncoms == 2 )
      dat = textscan(str,'%f,%f %[^ ] %f,%f:%f:%f\n');
      tmp = dat{1};
      dts = datenum(dat{4},month_str2num(dat{3}),dat{2},dat{5},dat{6},dat{7});
    else
      error('Unknown format! First data line: "%s"',cstr{begix});
    end; %if ncoms else

    %DEBUG:    disp(datestr(dts([1,end])));

    %% Simple QA/QC

    %goodix = [1;find(diff(dts)>0)+1];
    % Cut off everything after FIRST bad date
    goodix = 1:length(dts);
    badix = find(diff(dts)<=0);
    if ~isempty(badix); goodix(badix(1)+1:end)=[]; end;

    % Remove any data point where temperature, pressure, OR conductivity is bad
    % NOTE: YES, this applies for coastal sub-temperate data ONLY... So sue me.
    badix = find(0 >= tmp(goodix) | tmp(goodix) >= 39);
    goodix(badix) = [];
    if ( ~isempty(cnd) )
      badix = find(10 >= cnd(goodix) | cnd(goodix) >= 100);
      goodix(badix) = [];
    end;
    if ( ~isempty(prs) )
      badix = find(0.2 >= prs(goodix) | prs(goodix) >= 2000);
      goodix(badix) = [];
    end;

    res=[]; clear res
    res.sbe_seatemp.date = dts(goodix);
    res.sbe_seatemp.data = tmp(goodix);
    if ( ~isempty(cnd) )
      res.sbe_seacond.date = dts(goodix);
      res.sbe_seacond.data = cnd(goodix);

      res.sbe_salin.date = dts(goodix);
      if ( isempty(prs) )
        % Approximate for 5db pressure (potential errors << 0.33psu)
        res.sbe_salin.data = gsw_SP_from_C(cnd(goodix),tmp(goodix),5);
      else
        res.sbe_salin.data = gsw_SP_from_C(cnd(goodix),tmp(goodix),prs(goodix));
      end;
    end;
    if ( ~isempty(prs) )
      res.sbe_seapres.date = dts(goodix);
      res.sbe_seapres.data = prs(goodix);
    end;

    stn = merge_station_data(stn,res);

  end; %for fnix

return;
