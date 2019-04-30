function th = titlename(varargin)
%function th = titlename(varargin)
%
% Set Axes title, but also set Axes' Figure 'Name' property.
% Arguments are identical to TITLE (qv.).
% 
% Last Saved Time-stamp: <Tue 2016-12-13 16:05:20 Eastern Standard Time lew.gramer>

  th = title(varargin{:});
  txt = get(th, 'String');

  fh = getfig(th);
  % Just in case this TITLE was preprocessed to avoid "drooping letter
  % syndrome".  The underscore in MATLAB LaTex text means "subscript".
  try,
    txt = strrep(txt,'\_','_');
  catch,
    txt = strrep(string(txt),'\_','_');
  end;
  set(fh, 'Name', txt);

return;
