function [ax,fh] = jumpax(ax,varargin)
%function [ax,fh] = jumpax(ax,[framesize,][jumpsize,][xyzc[,axfun[,exit_at_end]]])
%
% Keep jumping AX (DEFAULT: GCA) forward one screenful in response to user
% INPUT: "Enter" jump plot ahead length of XLIM; "b"ack; "r"eset; "q"uit.
% Optional arg XYZC (DEFAULT: 'x') specifies which axis to jump. Optional arg
% AXFUN (DEFAULT: DATETICK3) can be a function-handle that is called with AX
% as its sole arg after each jump; specify [] for no such call. EXIT_AT_END
% (DEFAULT: True) causes the function to return when the frame moves beyond
% the end of the data; otherwise, the frame is reset to the beginning again.
%
% Optional FRAMESIZE if a numerical vector is the initial frame to set XLIM
% (or YLIM or ZLIM) to: if FRAMESIZE is a numerical scalar between 0 and 1,
% set frame to that fraction of the total length of data. (If FRAMESIZE is
% NaN or is not set, initial X/Y/ZLIM is used as initial frame.)
%
% Optional JUMPSIZE if a numerical vector is the distance to jump between
% each FRAMESIZE frame; DEFAULT is set based on FRAMESIZE (see above). 
%
% Optionally, AX may also be an array of AXES handles (e.g., from different
% FIGURES): then JUMPAX advances the X/Y/ZLIM of each AXES in the array.
% Assumes that all AXES have the same range of ordinates as AX(1).
%
% RETURNS handles of (first) AXES AX and FIGURE FH, resp.: if FIGURE or AXES
% was deleted during the call, then both return values are empty ([]).
%
% Last Saved Time-stamp: <Fri 2017-06-09 16:07:37 Eastern Daylight Time gramer>

  if ( ~exist('ax','var') || isempty(ax) )
    ax = gca;
  end;
  if ( isscalar(ax) )
    firstax = ax;
  else
    firstax = ax(1);
  end;

  args = varargin;

  % Non-positional (presence-optional) arguments - may not be "[]"!
  framesize = [];
  if ( numel(args)>0 && ~isempty(args{1}) && isnumeric(args{1}) )
    framesize = args{1};
    args(1) = [];
  end;
  if ( isnan(framesize) ); framesize=[]; end;

  jumpsize = [];
  if ( numel(args)>0 && ~isempty(args{1}) && isnumeric(args{1}) )
    jumpsize = args{1};
    args(1) = [];
  end;

  % Positional arguments
  if ( numel(args)>=1 && ~isempty(args{1}) )
    xyzc = args{1};
  else
    xyzc = 'x';
  end;
  if ( numel(args)>=2 && ~isempty(args{2}) )
    axfun = args{2};
  else
    %axfun = [];
    if ( exist('datetick3') )
      axfun = @datetick3;
    elseif ( exist('datetick2_5') )
      axfun = @datetick2_5;
    elseif ( exist('datetick2') )
      axfun = @datetick2;
    elseif ( exist('datetick') )
      axfun = @datetick;
    else
      axfun = [];
    end;
  end;
  if ( numel(args)>=3 && ~isempty(args{3}) )
    exit_at_end = args{3};
  else
    exit_at_end = true;
  end;

  fh = ancestor(firstax,'figure');
  % Return focus to our figure before initial input
  figure(fh);

  if ( strncmpi(xyzc,'x',1) )
    lmfn = @xlim;
    ordattnm = 'XData';
  elseif ( strncmpi(xyzc,'y',1) )
    lmfn = @ylim;
    ordattnm = 'YData';
  elseif ( strncmpi(xyzc,'z',1) )
    lmfn = @zlim;
    ordattnm = 'ZData';
  elseif ( strncmpi(xyzc,'c',1) )
    lmfn = @caxis;
    ordattnm = 'CData';
  else
    error('Optional XYZC arg must be ''x'', ''y'', ''z'', or ''c''');
  end;

  % Go searching for a child of the AX which has a non-empty *Data property
  ordmin = -inf;
  ordmax = +inf;
  if ( exit_at_end || ~isempty(framesize) )
    ploth = get(firstax,'Children');
    for hix = 1:numel(ploth)
      h = ploth(hix);
      att = get(h);
      % Stop as soon as we find some data
      if ( isfield(att,ordattnm) && ~isempty(att.(ordattnm)) )
        orddat = att.(ordattnm);
        ordmin = nanmin(orddat(:));
        ordmax = nanmax(orddat(:));
        break;
      end; %if ( isfield(att,ordattnm) )
    end; %for hix = 1:numel(ploth)
  end; %if ( exit_at_end )

  if ( ~isempty(framesize) )
    if ( numel(framesize) == 2 )
      lmfn(ax,framesize);
    elseif ( numel(framesize) == 1 && 0<framesize && framesize<=1 )
      if ( isfinite(ordmin) )
        lmfn(ax,[ordmin,ordmin+(framesize*(ordmax-ordmin))]);
      end;
    else
      error('Optional FRAMESIZE arg must be a numeric 2-vector or scalar in (0,1]');
    end;
  end;

  lm = lmfn(firstax);
  orig_lm = lm;

  if ( isempty(jumpsize) )
    dlm = (lm(2)-lm(1))*0.9;
    disp('"Enter" jump plot ahead one screenful; "h"alf-jump; "b"ack; "r"eset; "q"uit');
  else
    dlm = jumpsize;
    disp('"Enter" jump plot ahead by JUMPSIZE; "h"alf-jump; "b"ack; "r"eset; "q"uit');
  end;

  done = false;
  while (~done)
    inp = input('','s');
    if ( isempty(inp) ); inp=' '; end;
    % Did user destroy figure while we weren't looking?
    if ( ~ishandle(fh) || any(~ishandle(ax)) )
      fh = [];
      ax = [];
      done = true;
    else
      switch ( lower(inp(1)) )
       case 'b',
        lmfn(ax,lm-dlm);
       case 'h',
        lmfn(ax,lm+(dlm/2));
       case 'r',
        lmfn(ax,orig_lm);
       case 'q',
        done = true;
       otherwise,
        lmfn(ax,lm+dlm);
      end;
    end; %if ( ~ishandle(fh) || any(~ishandle(ax)) ) else

    if ( done )
      disp('Done');

    else
      lm = lmfn(firstax);
      % Return focus to our figure
      figure(fh);
      if ( ~isempty(axfun) )
        for axix = 1:numel(ax)
          axfun(ax(axix));
        end;
      end;

      if ( ordmin > max(lm) || min(lm) > ordmax )
        if ( exit_at_end )
          disp('Returning - no more data');
          done = true;
        else
          lmfn(ax,orig_lm);
          lm = lmfn(firstax);
        end; %if ( exit_at_end ) else
      end; %if ( ordmin > max(lm) || min(lm) > ordmax )

    end; %if ( done ) else


  end; %while (~done)

return;
