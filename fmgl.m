function fh = fmgl(fh)
%function fh = fmgl([fh])
%
% Create new full-screen Figure with HOLD ON, GRID ON, default FONTSIZE 20.
% If optional FH is passed in, pass focus to and do the above to it. Unlike
% FMG (v.), this version does NOT call TIGHTEN_AXES (v.) on the resulting
% figure! This is because e.g., for FIGUREs with right-labeled COLORBARs,
% this causes figure boundary to clip the colorbar label! In order to avoid
% this problem, TIGHTEN_AXES may still be called manually after calling
% COLORBAR, or try calling "COLORBAR_TIGHT".
%
% Last Saved Time-stamp: <Fri 2018-02-09 14:57:03 Eastern Standard Time gramer>

  if ( exist('fh','var') )
    fh = figure(fh);
  else
    fh = figure;
  end;
  maximize_graph;
  hold on;
  grid on;
  set(gca,'FontSize',20);  % For publication-ready fonts when printing

return;
