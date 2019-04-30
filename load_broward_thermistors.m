function brwd = load_broward_thermistors(facepath)
%function brwd = load_broward_thermistors(facepath)
%
% Load sea temperature data in an array of spreadsheets provided by Ken Banks
% of Broward County environmental protection. Saves data to .MAT file on
% first run, and reloads that .MAT file from then on.
%
% Last Saved Time-stamp: <Sat 2017-06-24 16:49:43 Eastern Daylight Time gramer>

  set_more off;

  if ( ~exist('facepath','var') )
    facepath = get_coral_path('FACE');
  end;

  browpath = fullfile(facepath,'Partners','Broward');

  matfname = fullfile(facepath,'broward_thermistors.mat');

  if ( exist(matfname,'file') )
    disp(['Loading ',matfname]);
    load(matfname);

  else
    disp(['Extracting Broward data from ',browpath]);

    brwd = [];
    stnms = {};

    [N,T,RAW] = xlsread(fullfile(browpath,'thermograph site locations.xlsx'));
    stnms = RAW(2:end-1,1);
    for stix=1:numel(stnms);
      stnm = stnms{stix};
      brwd.stns(stix).station_name = stnm;
      brwd.stns(stix).lon = dm2degrees(str2num(RAW{stix+1,5}));
      brwd.stns(stix).lat = dm2degrees(str2num(RAW{stix+1,4}));
      brwd.stns(stix).seatemp = empty_ts;
      brwd.(stnm) = brwd.stns(stix);
    end;


    ds = dir(fullfile(browpath,'*.xls'));

    monms = {'Jan','Feb','March','April','May','June','July','Aug','Sept','Oct','Nov','Dec'};

    for dix=1:numel(ds)

      C = textscan(ds(dix).name,'%[^-]-%f.xls');
      stnm = C{1}{:};
      yr = C{2};
      stix = find(strcmp(stnms,stnm));

      fname = fullfile(ds(dix).folder,ds(dix).name);

      if ( isempty(yr) || 1980 > yr || yr > 2050 || isempty(stix) )
        disp(['SKIPPING FILE: ',fname]);
        continue;
      end;

      disp(fname);
      % Read each monthly worksheet of each yearly file...
      for mo=1:12
        monm = monms{mo};


        % Format of Excel files is so hashed up, this is the best we can do
        % without lots of manual/special processing for each year and site :(
        try,
          [N,T,RAW]=xlsread(fname,monm);
        catch,
          RAW = {};
        end;
        if ( numel(RAW) == 1 && ~isempty(strfind(RAW{1},'NO DATA')) )
          disp(['Skipping NO DATA: ',num2str(yr),' ',monm]);
          continue;
        end;

        firstix = find(cellfun(@ischar,RAW(:,1)),1);

        res=[]; clear res
        if ( isnumeric(RAW{firstix,2}) )
          goodix = find(cellfun(@ischar,RAW(firstix:end,1))'&~isnan([RAW{firstix:end,2}]));
          if ( isempty(goodix) )
            % For some unaccountable reason, at certain stations beginning in
            % August 2003, the date format changes - and data gets much better!
            goodix = find(cellfun(@ischar,RAW(firstix:end,1))');
            goodix = goodix + firstix - 1;
            res.seatemp.date = datenum(RAW(goodix,1));
          else
            goodix = goodix + firstix - 1;
            res.seatemp.date = datenum(RAW(goodix,1)) + [RAW{goodix,2}]';
          end;
          tcolix = 3;
          if any([RAW{goodix,tcolix}]>50); tcolix = 4; end;
          res.seatemp.data = [RAW{goodix,tcolix}]';
        else
          res.seatemp.date = datenum(strcat(RAW(firstix:end,1),{' '},RAW(firstix:end,2)));
          res.seatemp.data = [RAW{firstix:end,3}]';
        end;

        %N=[]; T=[]; RAW=[]; clear N T RAW

        % Do some simple fix-ups
        badix = find(diff(res.seatemp.date)<(0.9/24))+1;
        badix = union(badix,find(isnan(res.seatemp.date)));
        if ( ~isempty(badix) )
          disp(['Skipping ',num2str(numel(badix)),' bad records: ',fname,': ',monm]);
          res.seatemp.date(badix) = [];
          res.seatemp.data(badix) = [];
        end;
        badix = find(diff(res.seatemp.date)<(0.9/24))+1;
        badix = union(badix,find(isnan(res.seatemp.date)));
        if ( ~isempty(badix) )
          disp(['Skipping ',num2str(numel(badix)),' MORE bad records: ',fname,': ',monm]);
          res.seatemp.date(badix) = [];
          res.seatemp.data(badix) = [];
        end;

        if ( ~is_valid_ts(res.seatemp) )
          disp(['INVALID TIME SERIES: ',fname,': ',monm]);
          %DEBUG:        keyboard;
          continue;
        elseif ( any(5>res.seatemp.data|res.seatemp.data>35) )
          disp(['INVALID TEMP. DATA: ',fname,': ',monm]);
          %DEBUG:        keyboard;
          continue;
        else
          w = warning('OFF','Ecoforecasts:mergedNonTS');
          brwd.stns(stix) = merge_station_data(brwd.stns(stix),res);
          brwd.(stnm) = merge_station_data(brwd.(stnm),res);
          warning(w);
        end; %if ( ~is_valid_ts(res.seatemp) ) else

        N=[]; T=[]; RAW=[]; clear N T RAW
        res=[]; clear res
      end; %for mo=1:12
    end; %for dix=1:numel(ds)

    disp(['Saving ',matfname]);
    save(matfname,'brwd');

  end; %if ( exist(matfname,'file') ) else

  set_more;

return;
