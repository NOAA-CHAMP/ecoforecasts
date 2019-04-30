function align_subplot_to_subplot(ax,xy)
%function align_subplot_to_subplot(ax,xy)
%
% Align AX(2) to the XY-component(s) of the 'Position' of AX(1). XY can be
% 'x' (DEFAULT), 'y', or 'xy' (to overlay AX(2) on the position of AX(1)).
%
% Last Saved Time-stamp: <Sun 2018-02-18 14:51:59 Eastern Standard Time gramer>

  if ( ~exist('xy','var') || isempty(xy) )
    xy = 'x';
  end;
  if ( ~exist('posattr','var') || isempty(posattr) )
    posattr = 'Position';
  end;

  pos1 = get(ax(1),posattr);
  pos2 = get(ax(2),posattr);
  if ( findstr(xy,'x') )
    pos2(1) = pos1(1); pos2(3) = pos1(3);
  end;
  if ( findstr(xy,'y') )
    pos2(2) = pos1(2); pos2(4) = pos1(4);
  end;
  set(ax(2),posattr,pos2);

return;
