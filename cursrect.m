function h = cursrect(fh)
%function h = cursrect(fh)
%
% Create a rectangle on the figure FH (DEFAULT: GCF) for use as a data
% cursor.  CALLS: ANNOTATION, GCF.
%
% Last Saved Time-stamp: <Fri 2013-01-18 15:44:36 Eastern Standard Time Lew.Gramer>

  if ( ~exist('fh','var') || isempty(fh) )
    %% If we preferred to do nothing (or error) when we had no Current Figure
    %fh = get(0,'CurrentFigure');
    fh = gcf;
  end;
  if ( ~isempty(fh) )
    annotation(fh,'rectangle',[0,0,0.5,1.0])
  end;

return;
