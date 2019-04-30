function fh = fmg(varargin)
%function fh = fmg([fh],[name1,param1,...])
%
% Create new full-screen Figure with HOLD ON, GRID ON, default FONTSIZE 20.
% If optional FH is passed in, pass focus to and do the above to it. From
% 2016 Nov 28 on, also calls TIGHTEN_AXES (v.) on the resulting figure. NOTE
% that after calling COLORBAR (v.), which also inherits 20 pt. fonts, this
% can cause the figure boundary to clip the colorbar label. To avoid this,
% call FMGL (v., "L" for "Loose") in place of FMG, then call COLORBAR_TIGHT
% in place of COLORBAR.
%
% Optional NAME,VALUE pairs are passed through to AXES of figure (GCA).
%
% Last Saved Time-stamp: <Sat 2018-11-24 15:00:36 Eastern Standard Time gramer>

  args = varargin;

  if ( numel(args)>0 && ~ischar(args{1}) )
    if ( ~ishandle(args{1}) && ~isnumeric(args{1}) )
      error('First arg must be a (potential) FIGURE handle or an attribute NAME');
    end;
    fh = figure(args{1});
    args(1) = [];
  else
    fh = figure;
  end;
  maximize_graph;
  hold on;
  grid on;
  grid minor;
  set(gca,'FontSize',20,args{:});  % For publication-ready fonts when printing

  % NOTE again: for right-labeled COLORBARs, figure may clip the label :(
  tighten_axes(fh);

return;
