function datetick3(varargin)
%function datetick3(varargin)
%
% By default, call DATETICK2_5('x',2,'keeplimits') if available, and then
% SET_DATETICK_CURSOR (qv. both). Otherwise, call DATETICK2 with only the
% arguments passed into DATETICK3, and then call SET_DATETICK_CURSOR. Sets
% DATETICK3 to be called *again* each time AXES is panned or zoomed.
%
% Formerly by default called DATETICK2_5('x',2,'keepticks','keeplimits').
%
% WARNING: This does *not* yet fix the annoying bug in DATETICK2 (Not My
% Code) that always puts focus to the first AXES in a MULTIPLOT or SUBPLOT.
%
% Last Saved Time-stamp: <Wed 2016-10-26 16:59:32 Eastern Daylight Time lew.gramer>

  args = varargin;

  if ( exist('datetick2_5','file') )
    %DEBUG:    disp('Calling DATETICK2_5');
    dtfn = @datetick2_5;
  elseif ( exist('datetick2','file') )
    %DEBUG:    disp('Calling DATETICK2');
    dtfn = @datetick2;
  elseif ( exist('datetickzoom','file') )
    %DEBUG:    disp('Calling DATETICKZOOM');
    dtfn = @datetickzoom;
  else
    %DEBUG:    disp('Calling DATETICK');
    dtfn = @datetick;
  end;

  %DEBUG:  disp(numel(args));

  if nargin > 0 & ishandle(args{1}) & isequal(get(args{1},'type'),'axes')
    % use the axes passed in
    axh = args{1};
    args(1) = [];
  else
    axh = gca;
  end;

  if ( isempty(args) )
    dtfn('x',2,'keeplimits');
  else
    dtfn(args{:});
  end;

  if ( exist('isOctave') && isOctave )
    warning('No ZOOM or PAN mode in Octave');
  else
    set(zoom(axh),'ActionPostCallback',@datetick3_panzoom)
    set(pan(get(axh,'parent')),'ActionPostCallback',@datetick3_panzoom)
  end;

  %addlistener(axh,'XLim','PostSet',@datetick3_panzoom);

  if ( exist('set_datetick_cursor','file') )
    set_datetick_cursor;
  end;

return;


%% PRIVATE FUNCTIONS

function datetick3_panzoom(varargin)
  %DEBUG:  disp('Calling DATETICK3_PANZOOM');
  %DEBUG:  disp(nargin);
  if nargin==2 && isstruct(varargin{2}) && isfield(varargin{2},'Axes') && isscalar(varargin{2}.Axes)
    datetick3(varargin{2}.Axes);
  end;
return;
