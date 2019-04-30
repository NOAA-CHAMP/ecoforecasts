function th = appendtitlename(varargin)
%function th = appendtitlename(varargin)
%
% Append a string to an Axes title, and also to that Axes' Figure 'Name'
% string property. Otherwise, arguments are identical to TITLE (qv.).
% 
% Last Saved Time-stamp: <Mon 2011-04-11 13:20:48  Lew.Gramer>

  if ( ishandle(varargin{1}) )
    ah = varargin{1};
    txt = varargin{2};
  else
    ah = gca;
    txt = varargin{1};
  end;

  th = get(ah,'Title');
  ttl = get(th, 'String');
  txt = [ttl txt];

  th = title(ah,txt);
  txt = get(th, 'String');
  fh = getfig(th);
  % Just in case this TITLE was preprocessed to avoid "drooping letter
  % syndrome".  The underscore in MATLAB LaTex text means "subscript".
  txt = strrep(txt,'\_','_');
  set(fh, 'Name', txt);

return;
