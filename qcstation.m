function stn = qcstation(stn,varargin)
%function stn = qcstation(stn,[addfld,][varname|{varnames...}|jump_initial_frame|{'linked',fld1,fl2...}]...)
%
% Use BRUSH to allow caller to INTERACTIVELY  perform QC on all time-series
% (STRUCT fields with .date and .data subfields) in station STN. To limit
% work, User may specify a field name CHAR, or a CELLSTR of field names. Uses
% JUMPAX (v.) to jump through data graph of each variable; if caller passes
% in a numeric 2-vector JUMP_INITIAL_FRAME, initial XLIM of each data graph
% is set to that; if it is a numeric scalar < 1, initial XLIM is that pct. of
% the total record size (DEFAULT: 0.25); if >= 1, frame size is that many
% days; if +Inf, then JUMPAX is not called. Calls INPUT between fields.
% If a CELLSTR has first element 'linked', all fields named in that CELLSTR
% are linked together: only the first one is QC'd, and fields removed or
% NaN'd from that field are also removed from all the fields linked to it.
%
% NOTE: If optional first argument ADDFLD (Default: false) is Logical true or
% a cell array with Logical true first element and character second element,
% add new fields to STN rather than replacing existing ones. If ADDFLD is a
% cell array, the QC'd time series for STN.(FLD) will have the field name
% [ADDFLD{2},FLD]; otherwise, QC'd time series will have field ['qc_',FLD].
%
% Last Saved Time-stamp: <Tue 2017-05-30 16:13:36 Eastern Daylight Time lew.gramer>


  fieldsModified = 0;

  flds = {};

  % Default to one-quarter of total record size for JUMPAX (below)
  frame_size = 0.25;
  jump_initial_frame = [];

  args = varargin(:);
  nargs = numel(args);

  addfld = false;
  if ( nargs > 0 && ( islogical(args{1}) || ( iscell(args{1}) && islogical(args{1}{1}) ) ) )
    if ( islogical(args{1}) )
      addfld = 'qc_';
    elseif ( args{1} && ischar(args{1}{2}) )
      addfld = args{1}{2};
    else
      error('If specified, ADDFLD must be a Logical or a cell array with TRUE 1st element and CHAR 2nd element');
    end;
    args(1) = [];
    nargs = nargs - 1;
  end;

  if ( nargs == 0 )
    % DEFAULT: Quality-control all time-series fields in STN
    flds = fieldnames(stn);
  else
    while ( nargs > 0 )
      arg = args{1};
      if ( ischar(arg) )
        flds{end+1} = arg;
      elseif ( iscell(arg) )
        flds(end+1) = {arg};
      elseif ( isvector(arg) && isnumeric(arg) )
        if ( numel(arg) == 1 )
          frame_size = arg;
        elseif ( numel(arg) == 2 )
          jump_initial_frame = arg;
        end;
      else
        error('ecoforecasts:qcstation:UnknownArg',...
                'Optional arg must be a CHAR, CELLSTR, or numeric vector');
      end;
      args(1) = [];
      nargs = nargs - 1;
    end; %while ( nargs > 0 )
  end;


  % QC each selected field (with linked and display-buddy fields, if any)

  for fldix = 1:numel(flds)
    fld = flds{fldix};
    linkedflds = {};
    dispedflds = {};
    % Special handling for linked or display-buddy field arguments...
    if ( iscell(fld) )
      if ( numel(fld) > 2 )
        error('ecoforecasts:qcstation:UnknownArg',...
              'Each CELL arg must have exactly one or two elts.');
      end;

      fldcel = fld; fld=[];
      if ( ischar(fldcel{1}) )
        fld = fldcel{1};
      elseif ( ~iscellstr(fldcel{1}) )
        error('ecoforecasts:qcstation:UnknownArg',...
              'Each elt. of a CELL arg must be either CHAR or CELLSTR');
      elseif ( numel(fldcel{1}) == 1 )
        fld = fldcel{1}{1};
        disp([fld,' had no fields to link']);
      else
        fld = fldcel{1}{1};
        linkedflds = fldcel{1}(2:end);
        disp([fld,' HAS LINKED FIELDS: ']);
        disp(linkedflds);
      end;
      if ( numel(fldcel) == 2 )
        if ( ischar(fldcel{2}) )
          dispedflds = fldcel(2);
        elseif ( ~iscellstr(fldcel{2}) )
          error('ecoforecasts:qcstation:UnknownArg',...
                'Each elt. of a CELL arg must be either CHAR or CELLSTR');
        else
          dispedflds = fldcel{2};
        end;
        disp([fld,' will also display with: ']);
        disp(dispedflds);
      end;
    end;

    if ( ~ischar(fld) )
      error('ecoforecasts:qcstation:UnknownArg',...
            'Invalid field sub-argument (#%d)??',fldix);
    end;
    fldnm = strrep(fld,'_','\_');

    if ( ~isfield(stn,fld) )
      warning('ecoforecasts:qcstation:UnknownArg',...
              'Field not found! STN.%s',fld);
    elseif ( ~is_valid_ts(stn.(fld)) )
      warning('ecoforecasts:qcstation:UnknownArg',...
              'Skipping non-time-series STN.%s',fld);

    else

      dts = stn.(fld).date;
      dat = stn.(fld).data;

      fh = fmg;
      lh = plot_ts(stn.(fld),'k');
      %lh = plot(dts,dat); datetick3;

      ttl = ['BRUSH points and ''Del'' to remove from STN.',fldnm];
      titlename(ttl);
      if ( numel(linkedflds) > 0 )
        linkedfldnms = cellstr(strrep(linkedflds,'_','\_'));
        ttl = strvcat(ttl,['And field(s):',sprintf(' STN.%s',linkedfldnms{:})]);
      end;
      title(ttl);

      set(lh,'XDataSource','dts');
      set(lh,'YDataSource','dat');

      clrs = {'r','b','g','m','c',[.5,.5,.5]};
      for dispedix = 1:numel(dispedflds);
        dispedfld = dispedflds{dispedix};
        if ( ~isfield(stn,dispedfld) )
          disp(['DISPLAY field not found: STN.',dispedfld]);
        elseif ( ~is_valid_ts(stn.(dispedfld)) )
          disp(['DISPLAY field not a time-series: STN.',dispedfld]);
        else
          plot_ts(stn.(dispedfld),'Color',clrs{dispedix});
        end;
      end;
      if ( ~isempty(dispedflds) )
        dispedfldnms = cellstr(strrep(dispedflds,'_','\_'));
        legend(fldnm,dispedfldnms{:});
      end;

      % Whippy-dippy new MATLAB Brush Mode
      brush(fh,'on');
      h = brush(fh);

      if ( isinf(frame_size) )
        %disp('When done editing, hit Enter to continue to save prompt...'); pause;
        inp = 'r';
        while ( strncmpi(inp,'r',1) )
          inp = input('Enter when done, or ''r''edisplay X tick labels: ','s');
          datetick3;
        end;
      else
        if ( isempty(jump_initial_frame) )
          if ( frame_size >= 1 )
            jump_initial_frame = [dts(1),dts(1)+frame_size] - 1;
          elseif ( frame_size > 0 )
            jump_initial_frame = [dts(1),dts(floor(end*frame_size))] - 1;
          else
            error('Invalid FRAME SIZE %g',frame_size);
          end;
        end;

        xlim(jump_initial_frame); 
        qcstation_jumpaction;

        disp(['BRUSH (or shift-BRUSH) points and ''Del'' to remove from STN.',fld]);
        [ax,fh] = jumpax([],[],@qcstation_jumpaction);
      end; %if ( isinf(frame_size) ) else

      % Stupid DATALINK (set X/YDataSource above) does not seem to work in reverse??
      if ( ishandle(lh) )
        dts = get(lh,'XData');
        dat = get(lh,'YData');
        if ( isempty(dts) || numel(dts) ~= numel(dat) )
          error('Bad graph data! To remove field, use RMFIELD');
        end;
      end;

      if ( ishandle(fh) )
        close(fh);
      end;

      inp = input('ENTER to save and move to next field, "s"kip to next (ABANDONS edits), or "q"uit: ','s');

      inp = strip(inp);
      if ( strncmpi(inp,'q',1) )
        break;
      elseif ( ~strncmpi(inp,'s',1) )
        if ( numel(dts) ~= numel(stn.(fld).date) || ...
             ~all(dat(:) ~= stn.(fld).data | (isnan(dat(:)) & isnan(stn.(fld).data))) )

          origdts = stn.(fld).date;
          origdat = stn.(fld).data;

          if ( ischar(addfld) )
            qcfld = [addfld,fld];
            disp(['Adding field STN.',qcfld]);
            stn.(qcfld).date = dts(:);
            stn.(qcfld).data = dat(:);
          else
            disp(['Modifying STN.',fld]);
            stn.(fld).date = dts(:);
            stn.(fld).data = dat(:);
          end; %if ( ischar(addfld) ) else
          fieldsModified = fieldsModified + 1;

          % Which timestamps did the user delete or set data to NaN for?
          delix = find(~ismember(origdts,dts));
          if ( isempty(delix) )
            % If NONE MISSING, at which timestamps did user set data to NaN?
            alldat = origdat;
            alldat(ismember(origdts,dts)) = dat;
            delix = find(~isnan(origdat) & isnan(alldat));
          end;

          % Delete corresponding timestamps in each of the "linked" time
          % series (according to the time resolution of each time series)
          deldts = origdts(delix);
          for linkedix = 1:numel(linkedflds);
            linkedfld = linkedflds{linkedix};
            if ( ~isfield(stn,linkedfld) )
              disp(['LINKED field not found: STN.',linkedfld]);
            elseif ( ~is_valid_ts(stn.(linkedfld)) )
              disp(['LINKED field not a valid time series: STN.',linkedfld]);
            else
              % badix = find(ismember(stn.(linkedfld).date,deldts));
              [badix,ig] = intersect_dates(stn.(linkedfld).date,deldts,median(diff(stn.(linkedfld).date))/2);
              if( isempty(badix) )
                disp(['Nothing to modify in linked STN.',linkedfld]);
              else
                disp(['Modifying linked STN.',linkedfld]);
                stn.(linkedfld).date(badix) = [];
                stn.(linkedfld).data(badix) = [];
                fieldsModified = fieldsModified + 1;
              end; %if( isempty(badix) ) else
            end; %%if ( ~isfield(stn,linkedfld) ) elseif ( ~is_valid_ts(stn.(linkedfld)) ) else
          end; %for linkedix = 1:numel(linkedflds);

        end; %if ( numel(dts) ~= numel(stn.(fld).date) || ...
      end; %if ( strncmpi(inp,'q',1) ) elseif ( ~strncmpi(inp,'s',1) )

    end; %if ( ~is_valid_ts(stn.(fld)) ) else

  end; %for fldix = 1:numel(flds)

  if ( fieldsModified == 0 )
    disp('(No edits were saved...)');
  else
    if ( isfield(stn,'station_name') )
      STNM = upper(stn.station_name);
    else
      STNM = 'STNAM1'; % Just as an example
    end;
    if ( ~exist('fld','var') || ~ischar(fld) )
      fld = 'air_t_2_degc'; % Just as an example
    end;
    disp([num2str(fieldsModified),' fields were modified...']);
    disp(['You may now wish to save your work, e.g., ']);
    disp(['  qc.stn = qcstation(qc.stn,''',fld,''');']);
    disp(['  save(''',STNM,'_portal_ALL_qc.mat'',''-struct'',''qc'',''-v7.3'');']);
  end;

return;


function qcstation_jumpaction(ax)
  if ( exist('ax','var') && ishandle(ax) )
    %ylim(ax,'default');
    ylim(ax,'auto');
    datetick3(ax);
  else
    %ylim('default');
    ylim('auto');
    datetick3;
  end;
return;
