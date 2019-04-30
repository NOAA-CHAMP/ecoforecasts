function stn = station_filter_bad_dates(stn,fname)
%function stn = station_filter_bad_dates(stn,fname)
%
% Remove "bad" dates from time series structs in struct STN. FNAME (DEFAULT:
% [ECOFORECASTS,'/',STNM,'-bad-dates.csv']) is the name of a CSV file with
% one line for each CONTIGUOUS range of bad dates for a given sensor. E.g.,
% the "bad dates file" might contain multiple lines like those shown below:
%   /* Sample "bad dates" file */
%   ndbc_sea_t,2001,5,24,2001,5,26
%   ndbc_sea_t,2003,5,24,2003,5,26
%   ndbc_air_t,2003,5,24,2003,5,26
% to specify a range of three dates (INCLUSIVE) in May of 2003 when both Air
% and Sea temperature from the NDBC were bad, and three dates in May of 2001
% when only the Sea temperature was bad. Parsing of FNAME ignores multi-line
% comments bracketed by "/* comment */" as in the C language.
%
% Last Saved Time-stamp: <Sat 2012-03-31 16:10:59  Lew.Gramer>


  if ( ~exist('fname','var') || isempty(fname) )
    fname = fullfile(get_ecoforecasts_path(),[stn.station_name,'-bad-dates.csv']);
  end;

  if ( ~exist(fname,'file') )
    warning('Ecoforecasts:BadDates:NoFile','No such file %s',fname);
    %%%%%%%%%%%%%%%
    % EARLY RETURN
    %%%%%%%%%%%%%%%
    return;
  end;

  fid = fopen(fname,'r');
  if ( fid <= 0 )
    error('Unable to read %s',fname);
    %%%%%%%%%%%%%%%
    % EARLY RETURN
    %%%%%%%%%%%%%%%
    return;
  end; %if ( fid <= 0 )

  try,
    C = textscan(fid,'%[^,],%f,%f,%f,%f,%f,%f', 'CommentStyle',{'/*','*/'});
    allflds = C{1};
    flds = unique(allflds);
    begyr = C{2};
    begmo = C{3};
    begdy = C{4};
    endyr = C{5};
    endmo = C{6};
    enddy = C{7};
    fclose(fid);
  catch,
    try, fclose(fid); catch, end;
    rethrow(lasterror);
  end; %try, catch

  if ( numel(allflds) ~= numel(begyr) || numel(allflds) ~= numel(endyr) )
    error('Invalid file format "%s"',fname);
  end;


  nRemoved = 0;

  for fldix = 1:numel(flds)
    fld = flds{fldix};
    if ( isfield(stn,fld) )
      for recix = find(strcmpi(fld,allflds))'
        begdt=datenum(begyr(recix),begmo(recix),begdy(recix));
        % Midnight of the day after ENDDY, minus 4 seconds
        enddt=datenum(endyr(recix),endmo(recix),enddy(recix)+1-(4/(24*3600)));
        if ( begdt > enddt )
          warning('Ecoforecasts:BadDates:BadDateRange',...
                  'Invalid dates for "%s": %s-%s%s', fld, datestr(begdt),...
                  datestr(enddt), sprintf('\n%s',fname));
        elseif ( ~isfield(stn.(fld),'date') )
          warning('Ecoforecasts:BadDates:NonTSField',...
                  'Field "%s" in STN has no "date" field',fld);
        else

          dtix = find(begdt<=stn.(fld).date & stn.(fld).date<=enddt);
          %DEBUG:          disp(['Removing ',num2str(length(dtix)),' records ',fld,': ',datestr(begdt),' to ',datestr(enddt)]);
          if ( ~isempty(dtix) )
            stn.(fld).date(dtix) = [];
            if ( isfield(stn.(fld),'data') )
              stn.(fld).data(dtix) = [];
              nRemoved = nRemoved + numel(dtix);
            end;
            % ADCP current or CTD profiles
            if ( isfield(stn.(fld),'prof') && ndims(stn.(fld).prof)==2 )
              stn.(fld).prof(dtix,:) = [];
              nRemoved = nRemoved + numel(dtix);
            end;
            % Two- and three-D fields, e.g., from satellite or model
            if ( isfield(stn.(fld),'field') )
              if ( ndims(stn.(fld).field)==3 )
                stn.(fld).field(dtix,:,:) = [];
                nRemoved = nRemoved + numel(dtix);
              elseif ( ndims(stn.(fld).field)==4 )
                stn.(fld).field(dtix,:,:,:) = [];
                nRemoved = nRemoved + numel(dtix);
              else
                warning('Ecoforecasts:BadDates:NonTimeField',...
                        'Field "%s" in STN has NDIMS<3 or >4',fld);
              end; %if ( ndims(stn.(fld).field)==3 ) elseif else
            end; %if ( isfield(stn.(fld),'field') )
          end; %if ( ~isempty(dtix) )

        end; %if ( begdt > enddt ) elseif else
      end; %for recix = find(strcmpi(fld,allflds))'

      %DEBUG:
      if (nRemoved); disp([upper(mfilename),' ',upper(stn.station_name),'.',fld,': Removed ',num2str(nRemoved),' bad points']); end;
      %DEBUG:
      nRemoved=0;

    end; %if ( isfield(stn,fld) )
  end; %for fldix = 1:numel(flds)

  %DEBUG:
  % if ( nRemoved>0 )
  %   disp([upper(mfilename),' ',upper(stn.station_name),...
  %         ': Removed ',num2str(nRemoved),' bad points']);
  % end;

return;
