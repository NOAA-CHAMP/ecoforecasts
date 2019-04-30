function [cs,h] = plot_contour_cs(varargin)
%function [cs,h] = plot_contour_cs([curaxes,zlevel,]cs[,plotargs])
%
% PLOT lines for all contours in contour matrix CS (v. CONTOURC) on the
% CURAXES (DEFAULT: GCA). Optional PLOTARGS is passed onto PLOT.
% (Hard to believe MATLAB doesn't have a way to do this, but... ?)
% If ZLEVEL is specified, use PLOT3; if ZLEVEL is the logical true
% scalar (LOGICAL(1)), use PLOT3 to plot the value levels from CS.
% 
% Last Saved Time-stamp: <Mon 2015-05-25 13:51:08 Eastern Daylight Time gramer>

  cs=[];
  h=[];

  args = varargin;
  if ( ishandle(args{1}) )
    ax = args{1};
    args(1) = [];
  elseif ( isempty(get(0,'Children')) )
    ax = gca;
    hold on;
  else
    ax = gca;
  end;

  zlevel = [];
  if ( isscalar(args{1}) )
    zlevel = args{1};
    args(1) = [];
  end;

  cs = args{1};
  args(1) = [];


  % Handle additional PLOT/PLOT3 arguments
  plotargs = args;
  if ( isempty(zlevel) )
    if ( isempty(plotargs) )
      plotfun = @(x,y)(plot(ax,x,y));
    else
      plotfun = @(x,y)(plot(ax,x,y,plotargs{:}));
    end;
  elseif ( isscalar(zlevel) && (isnumeric(zlevel) || islogical(zlevel)) )
    if ( islogical(zlevel) )
      doScalarZlevel = zlevel;
    else
      doScalarZlevel = false;
    end;
    if ( isempty(plotargs) )
      plotfun = @(x,y,z)(plot3(ax,x,y,z));
    else
      plotfun = @(x,y,z)(plot3(ax,x,y,z,plotargs{:}));
    end;
  else
    error('If ZLEVEL is specified, it must be a double or logical scalar');
  end;


  % Step through segments of the CS matrix, plotting all LON,LAT pairs
  cursegix = 1;
  while ( cursegix < size(cs,2) )
    seglen = floor(cs(2,cursegix));
    if ( isempty(zlevel) )
      h(end+1) = plotfun(cs(1,cursegix+1:cursegix+seglen),cs(2,cursegix+1:cursegix+seglen));
    else
      if ( doScalarZlevel )
        zs=repmat(zlevel,[1,seglen]);
      else
        zs=repmat(cs(1,cursegix),[1,seglen]);
      end;
      h(end+1) = plotfun(cs(1,cursegix+1:cursegix+seglen),cs(2,cursegix+1:cursegix+seglen),zs);
      zs=[]; clear zs
    end;

    cursegix = cursegix+seglen+1;
  end;

return;
