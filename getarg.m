function [val,args,foundArg] = getarg(args,nam,varargin)
%function [val,args,foundArg] = getarg(args,nam[,'default',DFLT][,'matchFun',MATCHFUN])
%
% Get the value VAL of the name-value pair argument NAM, from the array ARGS
% (e.g., VARARGIN).  If ARGS is specified in the output arguments, then the
% name-value pair is also removed from ARGS before it is returned. If NAM is
% not found and a 'default',DFLT name-value pair is passed in the arguments,
% return DFLT instead. If NAM is not found and no DFLT is given, VAL is []:
% optional output FOUNDARG is then False. STRMATCH(LOWER(NAM),LOWER(ARGS))
% used to find NAM among the string elements of ARGS, unless name-value pair
% 'matchFun',MATCHFUN is among the inputs: MATCHFUN may be a function name or
% handle, in which case FIND(FEVAL(MATCHFUN,ARGS(charargs),NAM)) is used to
% find NAM among all string ARGS; or the string 'exact' (i.e., use STRCMP).
%
% SAMPLE CALL: [val,varargin] = getarg(varargin,'Color','black')
%
% Last Saved Time-stamp: <Tue 2018-10-16 12:01:12 EDT lew.gramer>

  val = [];
  foundArg = false;


  % Need a RECURSION-STOP, otherwise this block of code would be more elegant
  if ( ~isempty(varargin) )
    [dflt,ig,foundDflt] = getarg(varargin,'def');
    if ( foundDflt )
      val = dflt;
      foundArg = true;
    end;
  end;
  if ( isempty(args) )
    return;
  end;

  % Search all character string args for a match to NAM
  charix = find(cellfun(@ischar,args));

  if ( ~isempty(charix) )

    [matchfun,ig,foundMatchFun] = getarg(varargin,'matchFun');
    if ( ~foundMatchFun )
      argix = strmatch(lower(nam),lower(args(charix)));

    else
      if ( strcmpi(matchfun,'exact') )
        matchfun = @strcmp;
      elseif ( ischar(matchfun) )
        matchfun = str2func(matchfun);
      end;
      if ( ~isa(matchfun,'function_handle') )
        error('Optional arg MATCHFUN must be the string ''exact'' or a function name or handle');
      end;
      argix = feval(matchfun,args(charix),nam);
      if ( islogical(argix) && numel(argix)==numel(charix) )
        argix = find(argix);
      end;
    end;

    if ( ~isempty(argix) )
      ix = charix(argix);
      if ( numel(args) == ix )
        error('Arg "%s" had no value',nam);
      end;
      val = args{ix+1};
      foundArg = true;

      % Remove name-value pair from ARGS before returning it
      if ( nargout>1 )
        args = args([1:ix-1,ix+2:end]);
      end;
    end;
  end;

return;
