function [stn,changed] = qa_ts_review(stn,varname,forceNoQA)
%function [stn,changed] = qa_ts_review(stn,varname,forceNoQA)
% 
% Allow user to edit a graph by modifying STN.([VARNAME '_qc_flags']), a time
% series of 4-byte bit-vectors, QC flags for 'stn.(varname)'. Flags modified: 
% 
%     Bit 16 = marked as BAD by human reviewer
%     Bit 17 = marked as GOOD (despite above flags) by human reviewer
% 
% Expects QS_TS (v.) to have been previously run on STN.(VARNAME), or else
% QA_TS to be calling this function. If optional FORCENOQA is True however,
% allow user to do "pure human QA" (clear any QA flags prior to review).
% 
% Returns CHANGED==True if the user made any edits of the QA flags/data.
% 
% Last Saved Time-stamp: <Fri 2011-12-30 20:53:47  lew.gramer>

  changed = false;

  stname = 'Station';
  if ( isfield(stn,'station_name') )
    stname = stn.station_name;
  end;

  qcname = [varname '_qc'];
  flgname = [varname '_qc_flags'];

  if ( ~isstruct(stn) || ~isfield(stn, varname) || ...
       ~isfield(stn.(varname), 'date') || ~isfield(stn.(varname), 'data') )
    warning('No valid time-series field "%s" in "stn"!', varname);
    return;
  end;

  if ( ~exist('forceNoQA','var') || isempty(forceNoQA) )
    forceNoQA = false;
  end;

  if ( forceNoQA )

    YesNoCxl = questdlg('Removing all QA flags, starting review from scratch?');
    if ( ~strcmpi(YesNoCxl,'Yes') )
      warning('Please rerun with FORCENOQA argument False!');
      %%%%%%%%%%%%%%%
      % EARLY RETURN
      %%%%%%%%%%%%%%%
      return;
    end;
    if ( isfield(stn, qcname) )
      stn = rmfield(stn,qcname);
    end;
    if ( isfield(stn, flgname) )
      stn = rmfield(stn,flgname);
    end;
    stn.(qcname) = stn.(varname);
    stn.(flgname) = stn.(varname);
    stn.(flgname).data(:) = 0;

  elseif ( ~isfield(stn, qcname) || ~isfield(stn, flgname) || ...
           ~isfield(stn.(qcname), 'date') || ~isfield(stn.(qcname), 'data') || ...
           ~isfield(stn.(flgname), 'date') || ~isfield(stn.(flgname), 'data') )

    warning('Function QA_TS() has not been run on this data yet??');
    %%%%%%%%%%%%%%%
    % EARLY RETURN
    %%%%%%%%%%%%%%%
    return;

  end;

  fh = figure;
  set(fh, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);
  setappdata(fh, 'dates', stn.(flgname).date);
  setappdata(fh, 'flags', stn.(flgname).data);
  setappdata(fh, 'data', stn.(varname).data);

  th = title(sprintf('Human Review of QA: %s %s', stname, varname));
  % Variable names often have underscores in them, so...
  set(th, 'Interpreter', 'none');

  hold on;
  ax = gca;
  % set(ax, 'OuterPosition', [0.03 0.00 0.94 1.00]);
  % set(ax, 'OuterPosition', [0.03 0.00 0.98 1.00]);
  set(ax, 'OuterPosition', [0.04 0.00 0.98 1.00]);
  grid on; grid minor;
  rawph = plot(stn.(varname).date, stn.(varname).data, 'r.');
  set(rawph, 'MarkerSize', 5);
  qcph = plot(stn.(qcname).date, stn.(qcname).data, 'b.');
  set(qcph, 'MarkerSize', 5);
  set(qcph, 'HitTest', 'off');
  datetick;
  legend('Pre-QA', 'Post-QA', 'Location', 'Best');
  hold off;

  h = datacursormode(gcf);
  set(h, 'UpdateFcn', @qa_ts_datacursor);

  % Create the button group.
  bgh = uibuttongroup('Visible', 'off', 'Units', 'normalized', ...
                      'Position', [0.02 0.13 0.12 0.40]);

  % Create two pushbuttons for Confirming or Canceling edits
  okic = uicontrol('Style', 'pushbutton', 'String', 'Confirm', ...
                   'Units', 'normalized', 'Value', 1, 'Parent', bgh, ...
                   'Pos', [0.05 0.81 0.50 0.10]);
  set(okic, 'Callback', {@qa_confirm_cancel, fh, true});
  cxlic = uicontrol('Style', 'pushbutton', 'String', 'Cancel', ...
                    'Units', 'normalized', 'Value', 1, 'Parent', bgh, ...
                    'Pos', [0.05 0.71 0.50 0.10]);
  set(cxlic, 'Callback', {@qa_confirm_cancel, fh, false});

  % Create three radio buttons in the button group.
  % Actions for radio buttons controlled by buttongroup callback below.
  % (Or could have used a fancy uipushtool(...) instead somehow?)
  nop = uicontrol('Style', 'Radio', 'String', 'Review only', 'Parent', bgh, ...
                  'Units', 'normalized', 'Pos', [0.05 0.61 0.85 0.10]);
  mBd = uicontrol('Style', 'Radio', 'String', 'Flag as BAD', 'Parent', bgh, ...
                  'Units', 'normalized', 'Pos', [0.05 0.51 0.85 0.10]);
  mGd = uicontrol('Style', 'Radio', 'String', 'Flag as GOOD', 'Parent', bgh, ...
                  'Units', 'normalized', 'Pos', [0.05 0.41 0.85 0.10]);
  clr = uicontrol('Style', 'Radio', 'String', 'CLEAR User Flags', 'Parent', bgh, ...
                  'Units', 'normalized', 'Pos', [0.05 0.31 0.85 0.10]);

  % Create two checkboxes for toggling data display
  qcuic = uicontrol('Style', 'checkbox', 'String', 'Show QA', ...
                    'Units', 'normalized', 'Value', 1, 'Parent', bgh, ...
                    'Pos', [0.05 0.21 0.85 0.10]);
  set(qcuic, 'Callback', {@qa_ts_uicb, ax, qcph});
  rawuic = uicontrol('Style', 'checkbox', 'String', 'Show Raw', ...
                     'Units', 'normalized', 'Value', 1, 'Parent', bgh, ...
                     'Pos', [0.05 0.11 0.85 0.10]);
  set(rawuic, 'Callback', {@qa_ts_uicb, ax, rawph});

  % Make sure Figure toolbar is showing...
  set(gcf,'Toolbar','figure');

  % Initialize some button group properties. 
  set(bgh, 'SelectionChangeFcn', @qa_ts_goodbadcbk);
  % set(bgh, 'SelectedObject', mBd);
  set(bgh, 'Visible', 'on');


  txtcnt = ...
      uicontrol('Style', 'text', 'Units', 'normalized', ...
                'Position', [0.02 0.53 0.12 0.40], 'FontSize', 7, ...
                'HorizontalAlignment', 'left', 'String', ...
                sprintf('%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n', ...
                        '  Bit 01 = outside known allowable range for this kind of data',...
                        '  Bit 02 = outside SDS standard deviations for this kind of data',...
                        '  Bit 03 = outside seasonal range for this kind of data (if known)',...
                        '  Bit 04 = outside SDS standard deviations from DESEASONALIZED data',...
                        '  Bit 05 = outside daily variance from previous 24 hours data',...
                        '  Bit 06 = outside daily SEASONAL variance from previous 24 hours',...
                        '  Bit 16 = marked as BAD by human reviewer',...
                        '  Bit 17 = marked as GOOD (despite above flags) by human reviewer',...
                        '  Bit 32 = status unknown (not good)') ...
                );


  %%%% 
  %%%% UIWAIT - blocking call
  %%%% 
  uiwait(fh);

  OK = [];
  if ( ishandle(fh) )
    OK = getappdata(fh, 'confirmed');
  end;

  was_edited = [];
  if ( ishandle(fh) )
    was_edited = getappdata(fh, 'was_edited');
  end;

  if ( ~isempty(OK) && OK  && isappdata(fh, 'flags') )

    stn.(flgname).data = getappdata(fh, 'flags');
    ix = find(stn.(flgname).data == 0 | bitget(stn.(flgname).data, 17));
    stn.(qcname).date = stn.(varname).date(ix);
    stn.(qcname).data = stn.(varname).data(ix);
    if ( ~isempty(was_edited) && was_edited )
      disp('QC flags and human edits have been stored, QC time series updated...');
      changed = true;
    else
      disp('QC flags and QC time series stored: No human edits were attempted?');
    end;

    % Clear any variables DERIVED from this QC variable also!
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

  else
    disp('');
    if ( ~isempty(was_edited) && was_edited )
      disp('QC flags retained, but (new) human edits WERE NOT SAVED!');
    else
      disp('No changes to QC flags and QC time series');
    end;
    disp('');

  end;

  if ( ishandle(fh) )
    close(fh);
  end;

return;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PRIVATE FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function qa_ts_uicb(uich, eventdata, ax, datph)
%function qa_ts_uicb(uich, eventdata, ax, datph)
% Toggle display of either Raw or QA time series plot

  flg = get(uich, 'Value');
  xl = xlim(ax); yl = ylim(ax);
  if (flg)
    set(datph, 'Visible', 'on');
  else
    set(datph, 'Visible', 'off');
  end;
  xlim(ax, xl); ylim(ax, yl);

return;


function qa_confirm_cancel(uich, eventdata, fh, confirmed)
%function qa_confirm_cancel(uich, eventdata, fh, confirmed)
% Close figure window, either saving changes, or canceling

  if ( confirmed )
    setappdata(fh, 'confirmed', 1);
  end;
  %%%% 
  %%%% UIRESUME - unblocks uiwait call in primary function above
  %%%% 
  uiresume(fh);

return;


function output_txt = qa_ts_datacursor(obj, event_obj)
% Display the position of the data cursor (as a tooltip)
% obj          Currently not used (empty)
% event_obj    Handle to event object
% output_txt   Data cursor text string (string or cell array of strings).

  dates = getappdata(gcf, 'dates');
  flags = getappdata(gcf, 'flags');

  pos = get(event_obj, 'Position');
  xix = find(pos(1) == dates, 1);

  qabits = bitget(flags(xix), 1:32);
  qastr = '';
  for qix = 1:length(qabits)
    if ( qabits(qix) )
      qastr = [ qastr ' ' num2str(qix) ];
    end;
  end;

  dvec = datevec(pos(1));
  jday = num2str( fix(pos(1) - datenum(dvec(1), 1, 1) + 1) );
  output_txt = { ['X: ',datestr(pos(1)),' (',jday,')'],...
                 ['Y: ',sprintf('%f (pt# %d)', pos(2), xix)],...
                 ['QC Flags:' qastr] };

return;


function qa_ts_goodbadcbk(src, eventdata)
%function qa_ts_goodbadcbk(src, eventdata)
% User can select points as Good, or Bad, or just Explore the time series

  ah = gca;

  newval = get(eventdata.NewValue, 'String');

  if ( ~isempty(strfind(upper(newval), 'GOOD')) )
    set(ah, 'ButtonDownFcn', @qa_ts_good_btndn);
  elseif ( ~isempty(strfind(upper(newval), 'BAD')) )
    set(ah, 'ButtonDownFcn', @qa_ts_bad_btndn);
  elseif ( ~isempty(strfind(upper(newval), 'CLEAR')) )
    set(ah, 'ButtonDownFcn', @qa_ts_clr_btndn);
  else
    set(ah, 'ButtonDownFcn', '');
  end;

  drawnow;

return;


function qa_ts_goodbad_btndn(src, eventdata, good)
%function qa_ts_goodbad_btndn(src, eventdata, good)
% Let user rubber-band select range of points to mark as GOOD or BAD

  ah = gca;

  dates = getappdata(gcf, 'dates');
  flags = getappdata(gcf, 'flags');
  data = getappdata(gcf, 'data');
  setappdata(gcf, 'was_edited', true);

  point1 = get(ah, 'CurrentPoint');    % button down detected
  finalRect = rbbox;                   % return figure units
  point2 = get(ah, 'CurrentPoint');    % button up detected
  point1 = point1(1,1:2);              % extract x and y
  point2 = point2(1,1:2);
  p1 = min(point1,point2);             % calculate locations
  offset = abs(point1-point2);         % and dimensions
  x1 = p1(1);
  y1 = p1(2);
  x2 = p1(1) + offset(1);
  y2 = p1(2) + offset(2);

  edtix = find(x1 <= dates & dates <= x2 & y1 <= data & data <= y2);

  fprintf(1, 'Dates %s to %s (%d points)... ', ...
          datestr(x1), datestr(x2), numel(edtix));
  if ( good == 1 )
    flags(edtix) = bitset(flags(edtix), 16, 0);
    flags(edtix) = bitset(flags(edtix), 17);
    setappdata(gcf, 'flags', flags);
    disp(' are now GOOD');
  elseif ( good == 0 )
    flags(edtix) = bitset(flags(edtix), 17, 0);
    flags(edtix) = bitset(flags(edtix), 16);
    setappdata(gcf, 'flags', flags);
    disp(' are now BAD');
  elseif ( good == -1 )
    flags(edtix) = bitset(flags(edtix), 16, 0);
    flags(edtix) = bitset(flags(edtix), 17, 0);
    setappdata(gcf, 'flags', flags);
    disp(' are now CLEARED of human flags');
  end;

return;

function qa_ts_good_btndn(src, eventdata)
  qa_ts_goodbad_btndn(src, eventdata, 1);
return;
function qa_ts_bad_btndn(src, eventdata)
  qa_ts_goodbad_btndn(src, eventdata, 0);
return;
function qa_ts_clr_btndn(src, eventdata)
  qa_ts_goodbad_btndn(src, eventdata, -1);
return;
