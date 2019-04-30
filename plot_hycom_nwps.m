1;

if ( ~exist('lonix','var') )
  [ULON,ULAT] = meshgrid(ulon,ulat);
  lons = ULON(reef_line_ix(:));
  ULON=[]; ULAT=[]; clear ULON ULAT
  [ig,lonix] = ismember(lons,ulon);
  lons=[]; clear lons
end;

if ( ~exist('H','var') )
  H = repmat(ubat(reef_line_ix(:))',[numel(wdts)-2,1]);
end;

fmg;
surf(wdts(1:end-2),ulon(lonix),zwx(:,reef_line_ix(:))',H');
shading interp;
datetick3;
view(3);
