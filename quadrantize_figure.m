function quadrantize_figure(xcoords,lbls)
%function quadrantize_figure(xcoords,lbls)
%
% Place very light gray rectangles spaced horizontally across figure, between
% pairs of coordinates XCOORDS(:,1:2). Place optional label strings LBLS{:}
% at center top of each rectangle.

  for ix = 1:size(xcoords,1)
    coords = [xcoords(ix,1),miny,xcoords(ix,2)-xcoords(ix,1),maxy];
    rectpos = ds2nfu(coords);

    anh=annotation('rectangle',rectpos,'Color','m','LineWidth',0);
    set(anh,'FaceColor',[.8,.8,.8],'FaceAlpha',0.3,'Color','none');
  end;

return;
