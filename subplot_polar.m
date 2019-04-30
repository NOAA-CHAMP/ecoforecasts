function ax = subplot_polar(n,m,ix,varargin)
%function ax = subplot_polar(n,m,ix,varargin)
%
% Identical to SUBPLOT (v., except for certain calling forms), but the gap
% or buffer space between subplots is significantly reduced. NOTE this
% function differs from SUBPLOT_TIGHT in that room must be left above each
% subplot for titles, as there is no "XLABEL" property for POLARPLOTS (v.)
%
% Last Saved Time-stamp: <Mon 2017-06-05 13:47:55 Eastern Daylight Time gramer>

  nargs = nargin;
  args = varargin;

  ax = [];

  % Buffers around entire figure
  lbuf = 0.08;
  rbuf = 0.02;
  tbuf = 0.02;
  bbuf = 0.05;

  % Buffer around each subplot
  % buf=0.03;
  xbuf=0.010;
  ybuf=0.065;


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

    dy = ((1-(tbuf+bbuf))/n) - (2*ybuf);
    dx = ((1-(lbuf+rbuf))/m) - (2*xbuf);
    % Tom and Dick
    [mix,nix] = ind2sub([m,n],ix(:));
    if ( max(nix) > n )
      error('Ecoforecasts:subplot_tight:SubplotIndexTooLarge',...
            'Index exceeds number of subplots.');
    end;
    nm = max(mix)-min(mix)+1;
    nn = max(nix)-min(nix)+1;
    ax = subplot('position',...
                 [((min(mix)-1)*(dx+2*xbuf))+lbuf,((n-max(nix))*(dy+2*ybuf))+bbuf,(dx*nm),(dy*nn)]);
    if ( nargs > 3 )
      set(ax,args{:});
    end;
  end;

return;
