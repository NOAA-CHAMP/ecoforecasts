function bbox = rect2bbox(rect)
%function bbox = rect2bbox(rect)
% Return a Bounding BOX [MINX,MAXX,MINY,MAXY] for RECTangle [X,Y,WIDTH,HEIGHT]

  if ( rect(3)<=0 || rect(4)<=0 )
    error('RECT should have form [X,Y,WIDTH,HEIGHT]');
  end;
  bbox = [rect(1),rect(1)+rect(3),rect(2),rect(2)+rect(4)];

return;
