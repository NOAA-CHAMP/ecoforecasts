function ax = subplot_tight(n,m,ix,varargin)
%function ax = subplot_tight(n,m,ix,varargin)
%
% Identical to SUBPLOT (v., except for certain calling forms), but the gap
% or buffer space between subplots is significantly reduced. NOTE this will
% usually leave enough room on a screen-sized figure for XLABEL on each
% subplot, but not for titles on any subplots below row 1. The following
% calling forms for SUBPLOT are *NOT* supported by SUBPLOT_TIGHT:
%     SUBPLOT(m,n,p,'v6') % Does not work with SUBPLOT_TIGHT
%     SUBPLOT(m,n,p,'replace') % Nope
%
% Optional fourth arg (not passed on to SUBPLOT) may be any of 'sharex',
% 'sharey', or 'sharexy': this compresses buffer around all interior panes
% and removes the corresponding axis label(s): SET(AX,'X/YTICKLABEL',[]);
%
% Last Saved Time-stamp: <Tue 2018-10-02 08:24:02 Eastern Daylight Time gramer>

  nargs = nargin;
  args = varargin;

  ax = [];

  % Boundaries of whole figure (leaving room for, e.g., TITLE and lowest XLABEL)
  x = 0;
  w = 1;
  y = 0.97;
  h = 0.97; % Leave room for title
  %h = 0.96; % Leave room for LARGE title

  % Buffers around each subplot AXES
  %lbuf = 0.05; % Leave room for a SMALL left YLabel
  lbuf = 0.08; % Leave room for a left YLabel
  rbuf = 0.02;
  tbuf = 0.02;
  %bbuf = 0.03; % Leave room for a SMALL bottom XLabel
  bbuf = 0.05; % Leave room for a bottom XLabel
  %bbuf = 0.07; % Leave room for a LARGE bottom XLabel

  % Special (and ugly) single-argument handling - see SUBPLOT
  isHandleArg = false;
  if ( nargs == 1 )
    if ( isnumeric(n) )
      if ( ishandle(n) && strcmpi(get(n,'type'),'axes') )
        isHandleArg = true;
      elseif ( isscalar(n) && 121 <= n && n <= 999 )
        n = num2str(n);
      end;
    end;
    if ( ischar(n) && ~isempty(regexp(n,'^[1-9][1-9][1-9]$')) )
      ix = str2num(n(3));
      m = str2num(n(2));
      n = str2num(n(1));
      nargs = 3;
    end;
  end;


  % Standard argument handling
  if ( nargs == 0 )
    error('Ecoforecasts:subplot_tight:BadCall',...
          'Unsupported calling form.');

  elseif ( nargs == 1 )
    if ( isHandleArg )
      ax = subplot(n);
    elseif ( n == 111 )
      subplot(n);
    else
      error('Ecoforecasts:subplot_tight:BadCall',...
            'Unsupported calling form.');
    end;

  elseif ( nargs == 2 )
    if ( strcmpi(n,'position') && isnumeric(m) && numel(m) == 4 )
      ax = subplot(n,m);
    else
      error('Ecoforecasts:subplot_tight:BadCall',...
            'Unsupported calling form.');
    end;

  else
    if ( nargs > 3 && (strcmpi(args{1},'replace') || strcmpi(args{1},'v6')) )
      warning('Ecoforecasts:subplot_tight:CannotReplace',...
              'Argument ''%s'' ignored.',args{1});
      nargs = nargs - 1;
      args = args(2:end);
    end;

    sharex = false;
    sharey = false;
    if ( nargs > 3 && strncmpi(args{1},'share',length('share')) )
      if ( strcmpi(args{1},'sharex') )
        sharex = true;
      elseif ( strcmpi(args{1},'sharey') )
        sharey = true;
      else
        sharex = true;
        sharey = true;
      end;
      nargs = nargs - 1;
      args = args(2:end);
    end;

    % Tom and Dick
    [mix,nix] = ind2sub([m,n],ix(:));
    if ( max(nix) > n )
      error('Ecoforecasts:subplot_tight:SubplotIndexTooLarge',...
            'Index exceeds number of subplots.');
    end;

    % Pane total width (subplot AXES plus buffer)
    pw = numel(unique(mix)) * (w/m);
    % Pane X
    px = x + ( (min(mix)-1)*(w/m) );
    % Pane total height (subplot AXES plus buffer)
    ph = numel(unique(nix)) * (h/n);
    % Pane Y
    py = y - ph - ( (min(nix)-1)*(h/n) );

    if ( sharex && (min(mix)>1) )
      sw = pw - rbuf;
      sx = px;
    else
      % Suplot AXES width
      sw = pw - lbuf - rbuf;
      % Subplot AXES X
      sx = px + lbuf;
    end;

    if ( sharey && (max(nix)<n) )
      sh = ph - tbuf;
      sy = py;
    else
      % Subplot AXES height
      sh = ph - tbuf - bbuf;
      % Subplot AXES Y
      sy = py + bbuf;
    end;


    ax = subplot('position',[sx,sy,sw,sh]);

    if ( sharex && (min(mix)>1) ); set(ax,'YTickLabel',[]); end;
    if ( sharey && (max(nix)<n) ); set(ax,'XTickLabel',[]); end;

    if ( nargs > 3 )
      set(ax,args{:});
    end;
  end;

  if ( ishandle(ax) )
    set(ax,'FontSize',20);  % For publication-ready fonts when printing
  end;

return;
