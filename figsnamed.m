function figs = figsnamed(patstr,varargin)
%function figs = figsnamed(patstr[,regexpfun][,figs])
%
% Return an array of all FIGUREs in FIGS (DEFAULT: all children of current
% display), whose Name property matches the regular expression pattern PATSTR
% (v. REGEXP). Useful, e.g., with REVIEWANIM (v.) If optional arg REGEXPFUN
% a function handle, use instead of REGEXP; if the string 'i', use REGEXPI.
%
% Last Saved Time-stamp: <Fri 2017-06-02 16:27:00 Eastern Daylight Time gramer>

  args = varargin;
  if ( numel(args) > 0 && isa(args{1},'function_handle') )
    regexpfun = args{1};
    args(1) = [];
  elseif ( numel(args) > 0 && strncmpi(args{1},'i',1) )
    regexpfun = @regexpi;
    args(1) = [];
  else
    regexpfun = @regexp;
  end;
  if ( numel(args) > 0 && isvector(args{1}) )
    figs = args{1};
    args(1) = [];
  else
    figs = get(0,'Children');
  end;
  if ( numel(args) > 0 )
    error('USAGE: figsnamed(patstr[,regexpfun][,figs])');
  end;

  nms = get(figs,'Name');
  ix = cellfun(@(x)(~isempty(x)),regexpfun(nms,patstr));
  figs = figs(ix);

return;
