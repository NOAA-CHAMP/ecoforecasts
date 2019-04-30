function newargs = plotspec_check(varargin)
%function newargs = plotspec_check(args)
%
% Check argument list (cell array) ARGS for a plotspec, e.g., 'ro-.' (v.
% PLOT). Convert each plotspec into a string of additional Name, Value
% pairs, inserted into NEWARGS at the point where the plotspec appeared.  For
% example, 'ro-.' becomes 'Color','r','Marker','o','LineStyle','-.'. ARGS may
% also be a list of arguments, e.g., PLOTSPEC_CHECK('r','-.'). Note ARGS list
% may already include name-value pairs: these are simply passed through.
%
% Last Saved Time-stamp: <Fri 2016-08-26 11:56:34 Eastern Daylight Time gramer>

  newargs = {};

  args = {};
  if ( numel(varargin) > 0 )
    if ( iscell(varargin{1}) )
      args = varargin{1};
      if ( numel(varargin) > 1 )
        warning('First arg was a CELL: succeeding arguments will be ignored');
      end; %if ( numel(varargin) > 1 )
    elseif ( ischar(varargin{1}) )
      args = varargin;
    end; %if ( iscell(varargin{1}) ); elseif ( ischar(varargin{1}) )
  end; %if ( numel(varargin) > 0 )


  for argix = 1:numel(args)

    arg = args{argix};

    % Check if this seems to be a plot spec
    % This mess is due to PLOT's ridiculously flexible calling sequence! 
    if ( ischar(arg) && numel(arg) <= 4 && ~isempty(regexp(arg,'^[.ox+*sdv^<>ph:-]*[bgrcmykw]*[.ox+*sdv^<>ph:-]*$')) )

      % If caller gave a plot spec with a color in it
      clrix = regexp(arg,'[bgrcmykw]','once');
      if ( ~isempty(clrix) )
        newargs(end+1:end+2) = {'Color',arg(clrix) };
      end;
      % If caller gave a plot spec with a line-style in it
      lsix = regexp(arg,'[:-]','once');
      if ( ~isempty(lsix) )
        ls = regexprep(arg(lsix:end),'[^-]*(-[.]|--|-|:)[^-]*','$1','once');
        newargs(end+1:end+2) = {'LineStyle',ls};
        % Delete line style, so '.' in '-.' doesn't mistakenly also become a Marker 
        arg = strrep(arg,ls,'');
      end;
      % If caller gave a plot spec with a face marker in it
      fmix = regexp(arg,'[.ox+*sdv^<>ph]','once');
      if ( ~isempty(fmix) )
        newargs(end+1:end+2) = { 'Marker',arg(fmix) };
      end;

    else
      % Pass through all other arguments
      newargs(end+1) = args(argix);

    end; %if ( ischar(arg) && numel(arg) <= 4 && ... ); else

  end; %for argix = 1:numel(args)

return;
