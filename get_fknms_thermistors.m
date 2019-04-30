function stns = get_fknms_thermistors(cfgfpath,matfname,doPlot)
%function stns = get_fknms_thermistors(cfgfpath,matfname,doPlot)
%
% Process and save metadata and sea temperature data from the Florida Keys
% National Marine Sanctuary thermistor network (deployed <1987 to present).
% Parses a configuration Excel file CFGFPATH for the number, location, and
% depth of individual thermistor sites, and directory where each site's data
% files are stored; attempts to process all data files found in directories.
%
% Returns a struct STNS with one field for each individual thermistor site;
% each thermistor site sub-struct has .lat,.lon,.depth (from CFGFPATH), and
% .fknms_seatemp,.hourly_fknms_seatemp,.hourly_fknms_seatemp_qc time series.
% DEFAULT CFGFPATH: FULLFILE(DATAPATH,'FKNMS','thermograph.xls')
%
% SIDE-EFFECT: After initial run, SAVEs resulting STNS struct in MATFNAME
% (DEFAULT: 'DATAPATH/FKNMS_thermograph.mat'), and loads it on later runs.
% If optional DOPLOT is True, plot thermistor sites map, and time series.
%
% Last Saved Time-stamp: <Fri 2011-10-14 09:03:05 Eastern Daylight Time gramer>

  set_more off;

  datapath = get_ecoforecasts_path('data');
  fknmspath = fullfile(datapath,'FKNMS');

  if ( ~exist('cfgfpath','var') || isempty(cfgfpath) )
    cfgfpath = fullfile(fknmspath,'thermograph.xls');
  end;
  if ( ~exist('matfname','var') || isempty(matfname) )
    [ig,matfbase,ig] = fileparts(cfgfpath);
    matfname = fullfile(datapath,['FKNMS_' matfbase '.mat']);
  end;
  if ( ~exist('doPlot','var') || isempty(doPlot) )
    doPlot = false;
  end;

  stns.lons = [];
  stns.lats = [];

  if ( exist(matfname,'file') )

    disp(['Loading ' matfname]);
    load(matfname,'stns');


  else

    disp(['Parsing configurations in ' cfgfpath]);
    cfgs = importdata(cfgfpath);

    % Parse each thermistor site from Excel file
    for ix = 1:size(cfgs.data)
      stnm = regexprep(cfgs.textdata{ix+1,1},'[^A-Z0-9_]','_');
      stnms{ix} = stnm;
      stns.(stnm).station_name = stnm;
      stns.(stnm).long_name = cfgs.textdata{ix+1,2};
      stns.(stnm).dir_name = cfgs.textdata{ix+1,3};

      stns.(stnm).lat = get_fknms_thermistors_dm2degree(cfgs.textdata{ix+1,9});
      stns.(stnm).lon = -get_fknms_thermistors_dm2degree(cfgs.textdata{ix+1,10});

      stns.lats(end+1,1) = stns.(stnm).lat;
      stns.lons(end+1,1) = stns.(stnm).lon;

      % Convert depth from feet to [m]
      stns.(stnm).depth = cfgs.data(ix) ./ 3.2808399;

      stns.(stnm).files = dir(fullfile(fknmspath,stns.(stnm).dir_name));
      stns.(stnm).processed_files = {};
    end; %for ix = 1:size(cfgs.data)

    stnms = grepstruct(stns,'^FKNMS_');
    for ix = 1:length(stnms);
      stnm = stnms{ix};

      % Find, read, and parse all thermistor files for this site
      for flix = 1:numel(stns.(stnm).files)
        fname = stns.(stnm).files(flix).name;

        x=[]; clear x;

        % TR is a combined header-data ("Record") file
        if ( regexp(fname,'^TR[0-9]') )

          fpath = fullfile(fknmspath,stns.(stnm).dir_name,fname);
          x = textread(fpath,'%s');
          if ( isempty(x) )
            warning('Empty record file "%s"',fpath);
            continue;
          end;
          %dtix = strmatch('Deployed:',x);
          dtix = strmatch('Began',x);
          if ( (numel(dtix)~=1) || (dtix+3>numel(x)) || ~strcmp(x{dtix+1},'Logging:') )
            warning('No unique valid "Began Logging" line found in "%s"',fpath);
            continue;
          end;
          begdt = datenum([x{dtix+2} ' ' x{dtix+3}]);
          openix = strmatch('Open',x);
          if ( isempty(openix) )
            warning('No "Open" tags found in "%s"',fpath);
            continue;
          end;
          x(openix) = {'-9.00'};
          tix = openix(1):length(x);
          valcelix = regexp(x(tix),'^[-]*[0-9]+[.][0-9]+$');
          valix = find(~cellfun(@isempty,valcelix));
          if ( isempty(valix) )
            warning('No sea temperature data found in "%s"',fpath);
            continue;
          end;
          dts = begdt + ([0:numel(valix)-1]*(2/24));
          dat = cellfun(@str2num,x(tix(valix)),'ErrorHandler',@get_fknms_thermistors_naner);
          ts.date = dts(:); ts.data = dat(:); dts=[]; dat=[]; clear dts dat;
          %DEBUG:
          oldts = ts;
          ts = despike_ts(ts,[],(2.1/24),[2,38]);
          ts.date(2>=ts.data | ts.data>=38 | ~isfinite(ts.data)) = [];
          ts.data(2>=ts.data | ts.data>=38 | ~isfinite(ts.data)) = [];
          if ( isempty(ts.date) )
            warning('No valid sea temperature data found in "%s"',fpath);
            ts=[]; clear ts;
            continue;
          end;
          %DEBUG:
          if ( length(ts.date) < 0.95*length(oldts.date) )
            fmg; plot_ts(oldts,ts,'r'); titlename(strrep([stnm ' ' fname],'_','\_'));
          end;
          result.fknms_seatemp.date = ts.date(:);
          result.fknms_seatemp.data = ts.data(:);
          warning('off','Ecoforecasts:mergedNonTS');
          stns.(stnm) = merge_station_data(stns.(stnm),result);
          warning('on','Ecoforecasts:mergedNonTS');
          result=[]; clear result;
          stns.(stnm).processed_files{end+1} = fname;
          %DEBUG:          disp(['Processed ' fpath]);


        % Just a header file - need to find corresponding "TA*" data file
        elseif ( regexp(fname,'^TH[0-9]') )

          fpath = fullfile(fknmspath,stns.(stnm).dir_name,fname);

          datfname = regexprep(fname,'^TH','TA');
          datfpath = fullfile(fknmspath,stns.(stnm).dir_name,datfname);
          if ( ~exist(datfpath,'file') )
            warning('Missing data file for header file "%s"',fpath);
            continue;
          end;
          x = textread(fpath,'%[^\n]\n');
          dtcelix = regexp(x,'^ *[0-9]+[/][0-9]+[/][0-9]+ +[0-9]+[:][0-9]+[:][0-9]+ ');
          dtix = find(~cellfun(@isempty,dtcelix));
          if ( numel(dtix) == 1 )
            begdt = datenum(x{dtix},'mm/dd/yy HH:MM:SS');
          else
            dtcelix = regexp(x,'^ *[0-9]+[/][0-9]+[/][0-9]+ +[0-9]+[ ][0-9]+[:][0-9]+ ');
            dtix = find(~cellfun(@isempty,dtcelix));
            if ( numel(dtix) == 1 )
              begdt = datenum(x{dtix},'mm/dd/yy HH MM');
            else
              dtcelix = regexp(x,'^ *[0-9]+[/][0-9]+[/][0-9]+ ');
              dtix = find(~cellfun(@isempty,dtcelix));
              if ( numel(dtix) == 1 )
                warning(['Found starting date "%s" but no time in "%s": ',...
                         'assuming starting time 16:00UT...'],...
                        x{dtix},fpath);
                begdt = datenum(x{dtix},'mm/dd/yy') + (16/24);
              else
                warning('Unique starting date-time not found in "%s"',fpath);
                continue;
              end;
            end;
          end;
          x=[]; clear x;

          x = textread(datfpath,'%s');
          valcelix = regexp(x,'^[-]*[0-9]+[.][0-9]+$');
          valix = find(~cellfun(@isempty,valcelix));
          if ( isempty(valix) )
            warning('No sea temperature data found in "%s"',datfpath);
            continue;
          end;
          dts = begdt + ([0:numel(valix)-1]*(2/24));
          %dat = str2num(x(valix});
          dat = cellfun(@str2num,x(valix),'ErrorHandler',@get_fknms_thermistors_naner);
          ts.date = dts(:); ts.data = dat(:); dts=[]; dat=[]; clear dts dat;
          %DEBUG:
          oldts = ts;
          ts = despike_ts(ts,[],(2.1/24),[2,38]);
          ts.date(2>=ts.data | ts.data>=38 | ~isfinite(ts.data)) = [];
          ts.data(2>=ts.data | ts.data>=38 | ~isfinite(ts.data)) = [];
          if ( isempty(ts.date) )
            warning('No valid sea temperature data found in "%s"',datfpath);
            ts=[]; clear ts;
            continue;
          end;
          %DEBUG:
          if ( length(ts.date) < 0.95*length(oldts.date) )
            fmg; plot_ts(oldts,ts,'r'); titlename(strrep([stnm ' ' fname],'_','\_'));
          end;
          result.fknms_seatemp.date = ts.date(:);
          result.fknms_seatemp.data = ts.data(:);
          warning('off','Ecoforecasts:mergedNonTS');
          stns.(stnm) = merge_station_data(stns.(stnm),result);
          warning('on','Ecoforecasts:mergedNonTS');
          result=[]; clear result;
          stns.(stnm).processed_files{end+1} = fname;
          %DEBUG:        disp(['Processed ' fpath]);

        end; %if ( regexp(fname,'^TR[0-9]') ) elseif

      end; %for flix = 1:numel(stns.(stnm).files)

      disp({stnm,numel(stns.(stnm).files),numel(stns.(stnm).processed_files)});

    end; %for ix = 1:length(stnms);


    %% Do basic sanity check, calculate derived time series
    stnms = grepstruct(stns,'^FKNMS_');
    for ix = 1:length(stnms);
      stnm = stnms{ix};
      if ( ~isfield(stns.(stnm),'fknms_seatemp') )
        warning('Missing STNS.%s.fknms_seatemp time series',stnm);
      elseif ( ~is_valid_ts(stns.(stnm).fknms_seatemp) )
        warning('Invalid STNS.%s.fknms_seatemp time series',stnm);
      else
        stns.(stnm).hourly_fknms_seatemp = ...
            interp_ts(stns.(stnm).fknms_seatemp,[],[],'linear',(4/24));
        stns.(stnm) = qa_ts(stns.(stnm),'hourly_fknms_seatemp',7,false);
      end;
    end; %for ix = 1:length(stnms);


    if ( ~isempty(matfname) )
      disp(['Saving to ' matfname]);
      save(matfname,'stns');
    end;

  end; %if ( exist(matfname,'file') ) else


  if ( doPlot )
    stnms = grepstruct(stns,'FKNMS_');
    stnms = strrep(stnms,'FKNMS_','');
    stnms = strrep(stnms,'_','\_');
    fmg;
    map_freef([min(stns.lons)-0.5,max(stns.lons)+0.5,min(stns.lats)-0.5,max(stns.lats)+0.5],[-5,-30,-200]);
    plot3(stns.lons,stns.lats,repmat(100,size(stns.lons)),'kp');
    text(stns.lons,stns.lats,repmat(100,size(stns.lons)),stnms,'Color','Red','FontSize',6);
    daspect([1,1,1]);
    titlename('FKNMS Thermistor Locations');

    tses = repmat(struct('date',[],'data',[]),0);
    legs = {};
    stnms = grepstruct(stns,'^FKNMS_');
    for ix = 1:length(stnms);
      stnm = stnms{ix};
      if ( isfield(stns.(stnm),'hourly_fknms_seatemp_qc') )
        tses(end+1) = stns.(stnm).hourly_fknms_seatemp_qc;
        stnm = strrep(stnm,'FKNMS_','');
        stnm = strrep(stnm,'_','\_');
        legs{end+1} = stnm;
      end;
    end;
    if ( ~isempty(tses) )
      fmg;
      plot_ts(tses);
      legend(legs, 'Location','Best');
      titlename('FKNMS Thermistor Time Series');
    end;
  end;

  set_more;

return;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PRIVATE FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function d = get_fknms_thermistors_dm2degree(str)
  [dstr,mstr] = strtok(str,'-');
  mstr(1) = [];
  d = dm2degrees([str2num(dstr),str2num(mstr)]);
return;

function val = get_fknms_thermistors_naner(varargin)
  val = NaN;
return;
