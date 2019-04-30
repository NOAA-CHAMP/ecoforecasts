function stns = get_lirman_thermistors(fname)
%function stns = get_lirman_thermistors(fname)
% Load thermistor location metadata and time series data from an Excel file
% named FNAME (DEFAULT: [DATAPATH 'Lirman Temperatures.Fla.2009-2010.xls']),
% populating a struct STNS whose fields are "station" structs, each with the
% 5-char site code, and each with its own '.seatemp' time series field.
%
% Last Saved Time-stamp: <Thu 2011-06-23 13:50:17  Lew.Gramer>

  datapath = get_ecoforecasts_path('data');

  global lirmanxls;

  if ( ~exist('fname','var') || isempty(fname) )
    fname = fullfile(datapath,'Lirman Temperatures.Fla.2009-2010.xls');
  end;

  matfname = strrep(fname,'.xls','.mat');
  if ( exist(matfname,'file') )
    disp(['Loading ' matfname]);
    load(matfname,'stns','lirmanxls');

  else

    disp(['Extracting from ' fname]);

    worksheets = { ...
        'TavRoxs',        'tavrk', 'Tav Rxs'                      ; ...
        'ShallowConch',   'consh', 'Shallow Conch'                ; ...
        'DeepConch',      'condp', 'Deep Conch'                   ; ...
        'NewNursery',     'bnpnn', 'BNP New Nursery'              ; ...
        'OldNursery',     'bnpon', 'BNP Old Nursery'              ; ...
        'BNPInshore',     'bnpin', 'BNP Inshore'                  ; ...
        'BNPMid',         'bnpmi', 'BNP Mid Channel'              ; ...
        'Pacific',        'bnppa', 'BNP Offshore (Pacific Reef)'  ; ...
                 };

    %DEBUG:
    if ( isempty(lirmanxls) || ~strcmpi(lirmanxls.filename,fname) )
      lirmanxls = [];
      clear lirmanxls;

      disp('Calling IMPORTDATA...');
      lirmanxls.filename = fname;
      lirmanxls = importdata(fname);
      if ( ~isstruct(lirmanxls) )
        error('Import FAILED on spreadsheet "%s"',fname);
      end;
    end;

    %DEBUG:
    disp('Extracting metadata...');
    ncoords = size(lirmanxls.textdata.Coords,1);
    nstns = size(worksheets,1);

    ngoodcoords = 0;
    for ix = 2:ncoords
      if ( ~isempty(lirmanxls.textdata.Coords{ix,2}) )
        latstr = lirmanxls.textdata.Coords{ix,2};
        latdm = sscanf(latstr,'%g %g');
        lonstr = lirmanxls.textdata.Coords{ix,3};
        londm = sscanf(lonstr,'%g %g');

        longnm = lirmanxls.textdata.Coords{ix,1};
        wkix = strmatch(longnm,worksheets(:,3));
        if ( numel(wkix) ~= 1 )
          warning('Found no unique site matching name "%s"!',longnm);
          continue;
        end;

        stnm = worksheets{wkix,2};

        stns.(stnm).station_name = stnm;
        stns.(stnm).long_name = longnm;
        stns.(stnm).lat = latdm(1) + (latdm(2) / 60.0);
        stns.(stnm).lon = -( londm(1) + (londm(2) / 60.0) );
        % Convert depth to [m]
        stns.(stnm).depth = lirmanxls.data.Coords((ix-1),1) * 0.3048;

        ngoodcoords = ngoodcoords + 1;
      end;
    end;

    if ( nstns ~= ngoodcoords )
      warning('Expected %g good coords, but found %g!',nstns,ngoodcoords);
    end;

    for ix = 1:nstns
      wknm = worksheets{ix,1};
      stnm = worksheets{ix,2};
      %DEBUG:
      disp(wknm);

      wks = lirmanxls.textdata.(wknm);
      wkn = lirmanxls.data.(wknm);

      stns.(stnm).seatemp.date = datenum(wks(2:end,2));
      stns.(stnm).seatemp.data = wkn(:,3);
    end;

    disp(['Save ' matfname]);
    save(matfname,'stns','lirmanxls');

  end; %if ( exist(matfname,'file') ) else

  % Make sure sea temperature is a COLUMN vector (Nx1)
  for cstnm=fieldnames(stns)';
    stnm=cstnm{:};
    stns.(stnm).seatemp.date = stns.(stnm).seatemp.date(:);
    stns.(stnm).seatemp.data = stns.(stnm).seatemp.data(:);
  end;

  % Manual quality-control (duplicate dates in record)
  stns.bnpin.seatemp.date([1:431]) = [];             stns.bnpin.seatemp.data([1:431]) = [];
  stns.bnpmi.seatemp.date([1:431]) = [];             stns.bnpmi.seatemp.data([1:431]) = [];
  stns.bnpnn.seatemp.date([1:431,2415]) = [];        stns.bnpnn.seatemp.data([1:431,2415]) = [];
  stns.bnpon.seatemp.date([1:431,22593:22644]) = []; stns.bnpon.seatemp.data([1:431,22593:22644]) = [];
  stns.bnppa.seatemp.date([1:431,10359:10361]) = []; stns.bnppa.seatemp.data([1:431,10359:10361]) = [];

return;
