function rect = bbox2rect(bboxes)
%function rect = bbox2rect(bboxes)
%
% Return a RECTangle [X,Y,WIDTH,HEIGHT] for bounding box [MINX,MAXX,MINY,MAXY]
% *OR* [MINX,MINY ; MAXX,MAXY]. If BBOXES is a 2x(N*2) matrix, assumes each
% BBOX(:,[1:2:end,2:2:end]) is a separate box, and returns N rectangles.
%
% Last Saved Time-stamp: <Thu 2018-10-11 17:30:52 Eastern Daylight Time gramer>

  if ( size(bboxes,1)==2 && mod(size(bboxes,2),2)==0 )
    boxix = 1:2:size(bboxes,2);
  elseif ( numel(bboxes)==4 )
    boxix = 1;
    if ( isvector(bboxes) )
      bboxes = [bboxes(1),bboxes(3) ; bboxes(2),bboxes(4)];
    else
      bboxes = [bboxes(1),bboxes(2) ; bboxes(3),bboxes(4)];
    end;
  else
    error('I do not understand BBOX shape (%d,%d)',size(bboxes,1),size(bboxes,2));
  end;

  % We could preallocate RECT, in case BBOXES was *really* big...
  for ixix=1:numel(boxix)
    ix = boxix(ixix);
    bbox = bboxes(:,[ix,ix+1]);
    if ( bbox(2,1)<=bbox(1,1) || bbox(2,2)<=bbox(1,2) )
      error('BBOX(%d) should have form [MINX,MAXX ; MINY,MAXY]',ixix);
    end;
    rect(ixix,1:4) = [bbox(1,1),bbox(1,2),bbox(2,1)-bbox(1,1),bbox(2,2)-bbox(1,2)];
  end;

return;
