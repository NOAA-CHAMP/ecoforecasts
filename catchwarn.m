function catchwarn(varargin)
%function catchwarn([msg|ME[,msg]])
%
% Display the last error which occured (see LASTERROR) as a warning, with
% complete stack trace. Display of this "warning" may also be disabled by
% calling WARNING('OFF',IDENT) with the identifier of the error. Normally
% called, e.g., in the CATCH of a TRY,CATCH (v.) If first arg is an instance
% of class MException (v.), use that as the last error instead.
%
% Last Saved Time-stamp: <Tue 2016-12-27 22:59:57 Eastern Standard Time gramer>

  e = [];
  msg = '';
  if ( numel(varargin) > 0 )
    if ( isa(varargin{1},'MException') )
      e = varargin{1};
      if ( numel(varargin) > 1 )
        msg = varargin{2};
      end;
    elseif ( ischar(varargin{1}) )
      msg = varargin{1};
    end;
  end;
  if ( isempty(e) )
    e = lasterror;
  end;

  warnstr = sprintf('\nCAUGHT ERROR: %s\n%s\n%s\n',msg,e.identifier,e.message);
  for ix=1:numel(e.stack)
    %warnstr = sprintf('%s\n%s:%g\n',warnstr,strrep(e.stack(ix).file,'\','\\'),e.stack(ix).line);
    warnstr = sprintf('%s\n%s:%g',warnstr,e.stack(ix).file,e.stack(ix).line);
  end;
  warnstr = sprintf('%s\n',warnstr);

  if ( isempty(e.identifier) )
    warning(warnstr);
  else
    s = warning('OFF','MATLAB:printf:BadEscapeSequenceInFormat');
    warning(e.identifier,warnstr);
    warning(s.state,'MATLAB:printf:BadEscapeSequenceInFormat');
  end;

return;
