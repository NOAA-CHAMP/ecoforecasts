function varargout = xlim_datetick(varargin)
%function varargout = xlim_datetick([N_or_Y][AXS,][LM,...])
%
% Set XLIM of all axes AXS (DEFAULT: GCA) to LM. Call DATETICK with arguments
% specified in [...] (DEFAULT: 'x',2,'keeplimits'). If N_or_Y=='N' (DEFAULT)
% then remove *all* X-tick labels by calling SET(AXS,'XTICKLABEL','').
%
% Last Saved Time-stamp: <Mon 2018-02-19 18:26:33 Eastern Standard Time gramer>

  args = varargin;

  if ( ~isempty(args) && ischar(args{1}) )
    N_or_Y = args{1};
    args(1) = [];
  else
    N_or_Y = 'N';
  end;

  if ( ~isempty(args) && all(ishandle(args{1})) )
    axs = args{1};
    args(1) = [];
  else
    axs = gca;
  end;

  if ( ~isempty(args) && numel(args{1}) == 2 && all(isnumeric(args{1})) )
    xlim(axs,args{1});
    args(1) = [];
  else
    error('XLim limits must be numeric 2-vector');
  end;

  for ix=1:numel(axs)
    if ( ~isempty(args) )
      datetick(axs(ix),args{:});
    else
      datetick(axs(ix),'x',2,'keeplimits');
    end;
  end;
  if ( strncmpi(N_or_Y,'N',1) )
    set(axs,'XTickLabel','');
  end;

return;
