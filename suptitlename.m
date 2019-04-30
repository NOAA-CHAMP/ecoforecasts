function th = suptitlename(varargin)
%function th = suptitlename(varargin)
%
% Set Axes title, but also set Axes' Figure 'Name' property.
% Arguments are identical to SUPTITLE (qv., MATLAB Exchange).
% 
% Last Saved Time-stamp: <Mon 2018-02-26 14:37:48 EST lew.gramer>

  %th = suptitle(varargin{:});
  th = suptitle2(varargin{:});
  txt = get(th, 'String');
  fh = getfig(th);
  % Just in case this TITLE was preprocessed to avoid "drooping letter
  % syndrome".  The underscore in MATLAB LaTex text means "subscript".
  txt = strrep(txt,'\_','_');
  set(fh, 'Name', txt);

return;
