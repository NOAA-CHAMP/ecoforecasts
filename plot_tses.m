function lhs = plot_tses(varargin)
%function lhs = plot_tses([ax,]tses,[nameN,{val1_1,...,val1_P},...],[nameN,valN,...])
%
% Call PLOT_TS (v.) for each element of matrix of STRUCTs TSES: This function
% differs from PLOT_TS itself, only in that just one TS matrix may be passed
% (right after the AXES handle, if any), and in that the user may also pass
% name-value pairs where the "value" is actually a cell array of values, one
% for each element of TSES. Thus, any "value" arg which is passed in that is
% a CELL must have the same number of elements as TSES does.
%
% Last Saved Time-stamp: <Wed 2018-05-02 17:13:49 Eastern Daylight Time gramer>

% % ALTERNATIVES TO THIS UGLY M-FUNCTION INCLUDE USE OF 'AXES.ColorOrder' OR...
% % From https://www.mathworks.com/matlabcentral/answers/49062-is-it-possible-to-increment-color-and-markers-automatically-for-a-plot-in-a-loop
% set(0,'DefaultAxesLineStyleOrder',{'+','o','*','.','x','s','d','^','v','>','<','p','h'});
% set(0,'defaultaxescolororder',[0 0 0]);


  lhs = [];

  args = varargin;
  if ( numel(args)>0 && isscalar(args{1}) && ishandle(args{1}) )
    ax = args{1};
    args(1) = [];
  else
    ax = gca;
  end;

  if ( numel(args)>0 && isfield(args{1},'date') && isfield(args{1},'data') )
    tses = args{1};
    args(1) = [];
  else
    error('First arg (after optional AXES) must be TS STRUCT or STRUCT matrix');
  end;

  ntses = numel(tses);

  for ix=1:ntses
    cellargs{ix} = {};
  end;
  if ( mod(numel(args),2) && ischar(args{1}) )
    % String PLOTSPEC
    for ix=1:ntses
      cellargs{ix}{end+1} = args{1};
    end;
    args(1) = [];
  end;

  while ( numel(args)>1 && ischar(args{1}) )
    for ix=1:ntses
      cellargs{ix}{end+1} = args{1};
    end;
    args(1) = [];

    % Second arg is the value or CELL of values
    if ( iscell(args{1}) )
      if ( numel(args{1}) ~= ntses )
        error('All cells should have %d values',ntses);
      end;
      for ix=1:ntses
        cellargs{ix}{end+1} = args{1}{ix};
      end;
    else
      for ix=1:ntses
        cellargs{ix}{end+1} = args{1};
      end;
    end;
    args(1) = [];
  end; %while ( numel(args)>1 && ischar(args{1}) )

  for ix=1:ntses
    lhs(ix) = plot_ts(ax,tses(ix),cellargs{ix}{:});
  end;

return;
