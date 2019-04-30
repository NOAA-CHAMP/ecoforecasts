function remove_extra_axis_labels(varargin)
%function remove_extra_axis_labels(varargin)
%
% For a FIGURE with SUBPLOTs (v.), this function will remove all internal X-
% and Y-axis labels. I.e., only the left-most subplots will keep YLABELs and
% only the bottom-most subplots will keep XLABELs (v.) If first optional arg
% ISHANDLE (v.), it is the FIGURE handle (DEFAULT: GCF). If optional arg
% ISCHAR, it specifies which axes to remove internal labels for: may be any
% combination of the letters 'x', 'y', and 'z' (DEFAULT: 'xy'). Only AXES
% with Position coordinate matching the bottom-/left-/outer-most AXES will be
% changed: so for example, the XLABEL of SUBPLOT(3,2,3) will not be changed
% if the figure does not also contain a SUBPLOT(3,2,5) immediately below it.
%
% Last Saved Time-stamp: <Fri 2018-02-16 15:56:36 Eastern Standard Time gramer>

  args = varargin;
  while ( numel(args) > 0 )
    if ( ishandle(args{1}) )
      fh = args{1};
    elseif ( ischar(args{1}) )
      axis = args{1};
    end;
    args(1) = [];
  end;
  if ( ~exist('fh','var') )
    fh = gcf;
  end;
  if ( ~exist('axis','var') )
    axis = 'xy';
  end;
  axis = lower(axis);

  % Find all the AXES handles
  axs = [];
  cs = get(fh,'child');
  for cix = 1:numel(cs)
    if ( strcmpi(get(cs(cix),'Type'),'axes') )
      axs(end+1) = cs(cix);
    end;
  end;

  cpos = get(axs,'Position');
  if ( isnumeric(cpos) )
    pos = cpos;
  else
    pos = reshape([cpos{:}],[4,numel(axs)])';
  end;

  lefy = min(pos(:,1));
  botx = min(pos(:,2));

  if ( ~isempty(strfind(axis,'x')) )
    % MORE ADVANCED FUTURE FEATURE (still To Be Coded):
    % Keep the XTICKLABEL of the bottom-most AXES - all AXES above it whose
    % XLIM match its can share its XTICKLABEL; when we find an AXES with a
    % different XLIM, keep its XTICKLABEL and continue upward

    % botix = find( (pos(:,2) - botx) > 1e-6 );
    % [ig,aboveix] = sort(pos(botix,1));
    % aboveix = botix(aboveix);
    % for ix = aboveix(2:end)

    % [ig,lefix] = min(pos(botix,1));

    % for ix = 1:numel(axs)

    chgx = find( (pos(:,2) - botx) > 1e-6 );
    %% For linked AXES, this removes ALL X Tick Labels
    %set(axs(chgx),'XTickLabel',[]);
    set(axs(chgx),'XTickLabel','');
  end;

  if ( ~isempty(strfind(axis,'y')) )
    % Leave only the leftmost AXES unchanged
    chgy = find( (pos(:,1) - lefy) > 1e-6 );
    set(axs(chgy),'YTickLabel','');
  end;

  if ( ~isempty(strfind(axis,'z')) )
    % Leave only the outermost AXES unchanged
    innz = min(pos(:,3));
    chgz = find( (pos(:,3) - innz) > 1e-6 );
    set(axs(chgz),'ZTickLabel','');
  end;

return;
