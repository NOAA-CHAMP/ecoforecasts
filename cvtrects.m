1;

if ( ~exist('LON','var') )
  [LON,LAT] = meshgrid(lon,lat);
end;

load('rota-rectangles.mat');

rs(1,1:4) = r;
rs(2,1:4) = r1;
rs(3,1:4) = r2;
rs(4,1:4) = r3;
rs(5,1:4) = r4;

clear Tops Rights Bottoms Lefts
for rix = 1:size(rs,1);
  b = rect2bbox(rs(rix,:));
  ixes = bbox2ind(b,LON,LAT);

  disp(sprintf('      %03d,%03d , %03d,%03d , %03d,%03d , %03d,%03d ; ...',...
               ixes(1,1),ixes(1,2), ixes(2,1),ixes(2,2), ixes(3,1),ixes(3,2), ixes(4,1),ixes(4,2) ));

  [ig,ix] = max(LAT(ixes(:,2)));
  Tops(rix) = ixes(ix,2);

  [ig,ix] = max(LON(ixes(:,1)));
  Rights(rix) = ixes(ix,1);

  [ig,ix] = min(LAT(ixes(:,2)));
  Bottoms(rix) = ixes(ix,2);

  [ig,ix] = max(LON(ixes(:,1)));
  Lefts(rix) = ixes(ix,1);
end;

clear b ig ixes r rix r1 r2 r3 r4 r5 rs

